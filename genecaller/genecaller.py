import os
import sys
import argparse
import yaml
import shutil
from genecaller.logging import get_logger
from genecaller.util import get_cmd, IntervalList, load_bam, get_data_file
from genecaller.bam_process import Bam, CopyNumberModel, Phased_vcf
from genecaller.gene import Gene
import tempfile
import copy
from packaging.version import Version

logger = get_logger(__name__, "INFO")

DEFAULT_MIN_MAP_QUAL = 30

CALLING_MIN_VERSIONS = {
    "sentieon driver": Version("202503"),
    "whatshap": Version("2.3"),
    "bcftools": Version("1.10"),
    "samtools": Version("1.16"),
}


class GeneCaller:
    tools = []

    def __init__(self):
        self.tools = self.setup_paths()
        self.genes = []
        self.input = []
        self.read_data = {}
        self.params = {}
        self.ref = ""

    def setup_paths(self):
        cmds = {}
        for cmd, version in CALLING_MIN_VERSIONS.items():
            cmd_path = get_cmd(cmd, version)
            if not cmd_path:
                sys.exit(1)
            cmds[cmd.split()[0]] = cmd_path
        return cmds

    def init_read_data(self):
        short, long = self.input
        sr_params = copy.copy(self.params)
        sr_params["long_read"] = False
        sr_params["prefix"] = sr_params["sr_prefix"]
        short_bam = Bam(short, self.ref, sr_params)
        read_data = {"short_read": {"bam": short_bam, "params": sr_params}}
        if long:
            lr_params = copy.copy(self.params)
            lr_params["long_read"] = True
            lr_params["prefix"] = lr_params["lr_prefix"]
            long_bam = Bam(long, self.ref, lr_params)
            read_data["long_read"] = {"bam": long_bam, "params": lr_params}
        self.read_data = read_data

    def reset_read_data(self):
        for r in self.read_data.values():
            r["bam"].reset()

    def liftover(self, gene):
        for name, rd in self.read_data.items():
            logger.info(f"Performing liftover on {name}...")
            bam = rd["bam"]
            params = copy.deepcopy(rd["params"])
            params["prefix"] += ".liftover"
            liftover_fname = f"{params['prefix']}.bam"
            liftover_results = bam.liftover_combine(gene, liftover_fname)
            liftover_bam = Bam(liftover_fname, self.ref, params)
            rd["liftover"] = liftover_bam
            rd["liftover_valid"] = [
                IntervalList(region=gene.gene_regions[i])
                .subtract(liftover_results[i]["failed_region"])
                .to_region_str()
                for i in range(len(gene.gene_names))
            ]

    def save_data(self, data, outfile, fmt="auto"):
        if fmt == "auto":
            if outfile.endswith(".json"):
                fmt = "json"
            elif outfile.endswith(".yaml"):
                fmt = "yaml"
            else:
                fmt = "json"
        with open(outfile, "w") as fout:
            if fmt == "json":
                import json

                serialized = json.dumps(data, indent=4)
            else:
                serialized = yaml.dump(data, default_flow_style=False)
            fout.write(serialized)

    def process_gene(self, gene):
        logger.info(f"Processing gene {', '.join(gene.gene_names)}")
        self.read_data["short_read"]["params"]["prefix"] = (
            f"{self.params['tmpdir']}/{self.params['sample_name']}.{gene.gene_names[0]}.{self.params['sr_prefix']}"
        )
        if "long_read" in self.read_data:
            self.read_data["long_read"]["params"]["prefix"] = (
                f"{self.params['tmpdir']}/{self.params['sample_name']}.{gene.gene_names[0]}.{self.params['lr_prefix']}"
            )

        self.liftover(gene)

        cn_model = CopyNumberModel(self.read_data, gene, self.params)
        if gene.known_longdel:
            logger.info("Detecting likely long deletion...")
            long_dels = cn_model.detect_longdels(
                gene.liftover_nodel[0], gene.known_longdel
            )
            logger.info("Result:")
            for i, d in enumerate(gene.known_longdel):
                logger.info(f" * {d}: delta_CN = {-long_dels[i]}")
            gene.set_longdels(long_dels)

        data = gene.init_result_data(self.read_data)

        # logger.info(f'Performing given variant calling on liftover BAM in {gene.gene_names[0]}...')
        # for name, rd in self.read_data.items():
        #     logger.info(f' * Running {name}...')
        #     liftover_bam, liftover_valid = rd['liftover'], rd['liftover_valid']
        #     liftover_vcf = f'{rd["params"]["prefix"]}.liftover.vcf.gz'
        #     liftover_bam.call_dnascope(liftover_vcf, {"region": liftover_valid[0], 'given': gene.diff_vcf[0]})
        #     for i, region in enumerate(gene.gene_regions):
        #         given_vcf = f'{rd["params"]["prefix"]}.given{i}.vcf.gz'
        #         rd['bam'].call_dnascope(given_vcf, {"region": region, 'given': gene.diff_vcf[i]})
        #         for _, r in data[i]['orig'].items():
        #             r[name]['given_vcf'] = given_vcf
        #     for _, r in data[0]['liftover'].items():
        #         r[name]['given_vcf'] = liftover_vcf

        # process each region. Get high mapq regions
        logger.info("Analyzing read stats of each gene")
        for name, rd in self.read_data.items():
            bam, params = rd["bam"], rd["params"]
            logger.info(f" * Running {name}...")
            for i, gene_region in enumerate(gene.gene_regions):
                # get high-map-qual region
                segments = bam.get_segments(gene_region)
                high_qual_regions = []
                last_start, last_end = -1, -1
                for s in segments:
                    if s.mq >= params["min_map_qual"]:
                        if last_end != s.start:
                            high_qual_regions.append(
                                IntervalList.interval_str(s.chrom, s.start, s.end)
                            )
                            last_start, last_end = s.start, s.end
                        else:
                            high_qual_regions[-1] = IntervalList.interval_str(
                                s.chrom, last_start, s.end
                            )
                            last_end = s.end
                high_mq_interval = (
                    IntervalList(region=",".join(high_qual_regions))
                    if high_qual_regions
                    else None
                )
                for region_id, d in data[i]["orig"].items():
                    if name not in d:
                        d[name] = {}
                    if region_id == "nodel":
                        d[name]["region"] = (
                            high_mq_interval.intersect(
                                gene.nodel_region[i]
                            ).to_region_str()
                            if high_mq_interval
                            else ""
                        )
                    else:
                        d[name]["region"] = (
                            high_mq_interval.intersect(d["region"]).to_region_str()
                            if high_mq_interval
                            else ""
                        )

        # call Copy Numbers
        logger.info("Call copy number...")
        data = cn_model.call_cn(data, ["short_read"], params)
        for i, d in data[0]["orig"].items():
            if i == "nodel":
                msg = " * Overall: "
            else:
                msg = f" * Long-deletion {d['region']}: "
            msg += f"Total CN = {d['total_cn']}. "
            for j, dd in enumerate(data):
                cn = dd["orig"][i]["cn"]
                if cn == -1:
                    logger.error("Call copy number failed.")
                msg += f"CN_{gene.gene_names[j]} = {cn}. "
            logger.info(msg)
        # update result data based on CN
        for i, d in enumerate(data):
            for kk, dd in d.items():
                keys_remove = [
                    r
                    for r, ddd in dd.items()
                    if r != "nodel" and ddd["cn"] == dd["nodel"]["cn"]
                ]
                if keys_remove:
                    d[kk] = dict(
                        [(k, v) for k, v in dd.items() if k not in keys_remove]
                    )
                    for seq_key in self.read_data.keys():
                        nodel_region = IntervalList(
                            region=dd["nodel"][seq_key]["region"]
                        )
                        for k in keys_remove:
                            nodel_region = nodel_region.union(dd[k][seq_key]["region"])
                        dd["nodel"][seq_key]["region"] = nodel_region.to_region_str()

        # variant calling with short reads for each region
        logger.info("Call small variants...")
        for name, rd in self.read_data.items():
            logger.info(f" * Running {name}...")
            bam, params, liftover_bam = rd["bam"], rd["params"], rd["liftover"]
            bam.call_variant_all_regions(gene, data, "orig", params)
            liftover_bam.call_variant_all_regions(gene, data, "liftover", params)

        # assign phased liftover variants to genes
        logger.info("Create VCF for each gene")
        vcf_proc = Phased_vcf(data, gene, self.read_data)
        data = vcf_proc.process()

        # save result
        logger.info("Prepare final output")
        gene_data = gene.prepare_output(data)
        self.reset_read_data()
        return gene_data

    def print_param(self):
        logger.info("Input loaded:")
        logger.info(f" * Reference file: {self.ref}")
        logger.info(f" * Input short reads bam: {self.input[0]}")
        logger.info(
            f" * Input long reads bam: {self.input[1] if self.input[1] else 'None'}"
        )
        logger.info(f" * Genes: {' '.join([g.gene_names[0] for g in self.genes if g != 'main'])}")

    def load_params(self):
        default_cfg_fname = get_data_file("genes.yaml")
        with open(default_cfg_fname) as f:
            default_config = yaml.safe_load(f)
            default_genes = [g for g in default_config.keys() if g != 'main']
        parser = argparse.ArgumentParser(
            description="Targeted variant caller for genes with highly similar paralogs"
        )
        parser.add_argument(
            "--short",
            "-s",
            required=True,
            help="Input short-read BAM or CRAM (required)",
        )
        parser.add_argument(
            "--long", "-l", help="Input long-read BAM or CRAM (optional)"
        )
        parser.add_argument(
            "--sr_model", 
            required=True,
            help="Short read model bundle (required)"
        )
        parser.add_argument(
            "--lr_model", help="Long read model bundle (required if --long is available)"
        )
        parser.add_argument("--reference", "-r", required=True, help="Reference file")
        parser.add_argument(
            "--genes",
            "-g",
            required=True,
            help=f"List of genes to be called (comma seperated). \nSupported genes: {', '.join(default_genes)}",
        )
        parser.add_argument(
            "--sample_name",
            help="Sample name (default: SM tag in the input short read bam file will be used.)",
        )
        parser.add_argument(
            "--sr_prefix", help="short read result prefix", default="sr"
        )
        parser.add_argument("--lr_prefix", help="long read result prefix", default="lr")
        parser.add_argument("--config", help=argparse.SUPPRESS)
        parser.add_argument("--outdir", "-o", required=True, help="Output directory")
        parser.add_argument(
            "--keep_temp", help=argparse.SUPPRESS, action="store_true"
        )
        args = parser.parse_args()
        for input in (args.short, args.long):
            if not input:
                continue
            if input.endswith("bam") and not os.path.exists(input + ".bai"):
                raise Exception(f"Index file {input}.bai does not exist.")
            if input.endswith("cram") and not os.path.exists(input + ".crai"):
                raise Exception(f"Index file {input}.crai does not exist.")
        if args.long and not args.lr_model:
            raise Exception('Long read model is not specified. Please download the long read model from Sentieon')
        if not os.path.exists(args.reference + ".fai"):
            raise Exception(f"Index file {args.reference}.fai does not exist.")
        if args.config:
            with open(args.config) as f:
                config = yaml.safe_load(f)
        else:
            config = default_config
        config['main']['sr_model'] = args.sr_model + '/dnascope.model'
        if args.lr_model:
            config['main']['lr_model'] = args.sr_model + '/diploid_model'
        genes = [s.strip().upper() for s in args.genes.split(",")]
        missing = [g for g in genes if g not in config]
        if missing:
            raise Exception(
                f"Genes {', '.join(missing)} are not supported at this moment. Below is a list of supported regions: \n{', '.join(list(config.keys()))}"
            )
        self.genes = [Gene(config[g]) for g in genes]
        self.input = args.short, args.long
        self.ref = args.reference
        if not args.sample_name:
            bam = load_bam(args.short, args.reference)
            rg = bam.header.get("RG")
            sample_name = rg[0].get("SM") if rg else None
            if not sample_name:
                raise Exception(
                    f"SM tag not found in {args.short}. Please explicit set --sample_name."
                )
        else:
            sample_name = args.sample_name
        self.params = config["main"]
        os.makedirs(args.outdir, exist_ok=True)
        sent_tmpdir = os.environ.get("SENTIEON_TMPDIR")
        if sent_tmpdir:
            tmpdir = tempfile.mkdtemp(dir=sent_tmpdir, prefix="job")
        else:
            tmpdir = tempfile.mkdtemp(dir=args.outdir, prefix=f"{sample_name}_tmp")
        self.params.update(
            {
                "outdir": args.outdir,
                "tmpdir": tmpdir,
                "sample_name": sample_name,
                "tools": self.tools,
                "sr_prefix": args.sr_prefix,
                "lr_prefix": args.lr_prefix,
                "logger": logger,
                "keep_tmp": args.keep_temp,
            }
        )
        self.print_param()
        self.init_read_data()

    def call(self):
        self.load_params()
        result = {}
        for g in self.genes:
            gene_data = self.process_gene(g)
            result.update(gene_data)
        self.save_data(
            result, f"{self.params['outdir']}/{self.params['sample_name']}.yaml"
        )
        if not self.params["keep_tmp"]:
            shutil.rmtree(self.params["tmpdir"])
