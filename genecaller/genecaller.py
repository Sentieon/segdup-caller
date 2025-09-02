import os
import sys
import argparse
import yaml
import shutil
import logging
from typing import Dict, Any
from genecaller.logging import get_logger, set_global_log_level
from genecaller.util import get_cmd, IntervalList, load_bam, get_data_file
from genecaller.bam_process import Bam
from genecaller.gene import Gene
from genecaller import __version__
import tempfile
import copy
from packaging.version import Version

logger = get_logger(__name__)

CALLING_MIN_VERSIONS = {
    "sentieon driver": Version("202503"),
    "whatshap": Version("2.3"),
    "bcftools": Version("1.10"),
    "samtools": Version("1.16"),
}


class GeneCaller:
    tools = []

    def __init__(self) -> None:
        self.tools = {}
        self.genes = []
        self.input = []
        self.read_data = {}
        self.params = {}
        self.ref = ""

    def setup_paths(self) -> dict[str, str]:
        cmds = {}
        for cmd, version in CALLING_MIN_VERSIONS.items():
            cmd_path = get_cmd(cmd, version)
            if not cmd_path:
                sys.exit(1)
            cmds[cmd.split()[0]] = cmd_path
        return cmds

    def init_read_data(self) -> None:
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

    def reset_read_data(self, params: Dict[str, Any]) -> None:
        for r in self.read_data.values():
            r["bam"].reset()
        self.params = params

    def liftover(self, gene: "Gene") -> None:
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
                IntervalList(region=gene.liftover_target_regions[i])
                .subtract(liftover_results[i]["failed_region"])
                .to_region_str()
                for i in range(len(gene.liftover_region_names))
            ]

    def save_data(
        self, data: Dict[str, Dict[str, Any]], outfile: str, fmt: str = "auto"
    ) -> None:
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

    def process_gene(self, gene: "Gene") -> Dict[str, Any]:
        logger.info(f"Processing gene {', '.join(gene.gene_names)}")
        gene_config = gene.config
        orig_params = self.params
        if gene_config:
            self.params.update(gene_config)
            if "dp_norm" in gene_config:
                self.read_data["short_read"]["bam"].dp_norm = gene_config["dp_norm"]
        self.read_data["short_read"]["params"]["prefix"] = (
            f"{self.params['tmpdir']}/{self.params['sample_name']}.{gene.gene_names[0]}.{self.params['sr_prefix']}"
        )
        if "long_read" in self.read_data:
            self.read_data["long_read"]["params"]["prefix"] = (
                f"{self.params['tmpdir']}/{self.params['sample_name']}.{gene.gene_names[0]}.{self.params['lr_prefix']}"
            )

        self.liftover(gene)

        logger.info("Call copy number...")
        gene.call_cn(self.read_data, self.params)

        # variant calling with short reads for each region
        logger.info("Call small variants...")
        gene.call_small_var()
        gene.resolve_phased(self.params)

        # save result
        logger.info("Prepare final output")
        gene_data = gene.prepare_output()
        self.reset_read_data(orig_params)
        return gene_data

    def print_param(self) -> None:
        logger.info("Input loaded:")
        logger.info(f" * Reference file: {self.ref}")
        logger.info(f" * Input short reads BAM: {self.input[0]}")
        logger.info(
            f" * Input long reads BAM: {self.input[1] if self.input[1] else 'None'}"
        )
        logger.info(
            f" * Genes: {' '.join([g.gene_names[0] for g in self.genes if g != 'main'])}"
        )

    def load_params(self) -> None:
        default_cfg_fname = get_data_file("genes.yaml")
        if default_cfg_fname is None:
            logger.error("Failed to load genes.yaml file.")
            sys.exit(1)
        with open(default_cfg_fname) as f:
            default_config = yaml.safe_load(f)
            default_genes = [g for g in default_config.keys() if g != "main"]
        parser = argparse.ArgumentParser(
            description="Targeted variant caller for genes with highly similar paralogs",
            prog="segdup-caller",
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
            "--sr_model", required=True, help="Short read model bundle (required)"
        )
        parser.add_argument(
            "--lr_model",
            help="Long read model bundle (required if --long is available)",
        )
        parser.add_argument("--reference", "-r", required=True, help="Reference file")
        parser.add_argument(
            "--genes",
            "-g",
            required=True,
            help=f"List of genes to be called (comma separated). \nSupported genes: {', '.join(default_genes)}",
        )
        parser.add_argument(
            "--sample_name",
            help="Sample name (default: SM tag in the input short-read BAM file will be used)",
        )
        parser.add_argument(
            "--sr_prefix", help="Short read result prefix", default="sr"
        )
        parser.add_argument("--lr_prefix", help="Long read result prefix", default="lr")
        parser.add_argument(
            "--config", help="Custom gene configuration file (advanced users only)"
        )
        parser.add_argument("--outdir", "-o", required=True, help="Output directory")
        parser.add_argument(
            "--version", action="version", version=f"segdup-caller {__version__}"
        )
        parser.add_argument(
            "--keep_temp",
            action="store_true",
            help="Keep temporary files for debugging",
        )
        parser.add_argument(
            "--log_level",
            choices=["DEBUG", "INFO", "WARNING", "ERROR"],
            default="INFO",
            help=argparse.SUPPRESS,
        )
        args = parser.parse_args()

        # Setup tools after argument parsing (so --version and --help work without dependencies)
        self.tools = self.setup_paths()
        for input in (args.short, args.long):
            if not input:
                continue
            if input.endswith("bam") and not os.path.exists(input + ".bai"):
                raise Exception(f"Index file {input}.bai does not exist.")
            if input.endswith("cram") and not os.path.exists(input + ".crai"):
                raise Exception(f"Index file {input}.crai does not exist.")
        if args.long and not args.lr_model:
            raise Exception(
                "Long read model is not specified. Please download the long read model from Sentieon"
            )
        if not os.path.exists(args.reference + ".fai"):
            raise Exception(f"Index file {args.reference}.fai does not exist.")
        if args.config:
            with open(args.config) as f:
                config = yaml.safe_load(f)
        else:
            config = default_config
        config["main"]["sr_model"] = args.sr_model + "/dnascope.model"
        if args.lr_model:
            config["main"]["lr_model"] = args.lr_model + "/diploid_model"
        genes = [s.strip().upper() for s in args.genes.split(",")]
        missing = [g for g in genes if g not in config]
        if missing:
            raise Exception(
                f"Genes {', '.join(missing)} are not supported at this moment. Below is a list of supported regions: \n{', '.join(list(config.keys()))}"
            )
        Gene.cn_prior_std = config["main"].get("cn_prior_std", 2)
        self.genes = [Gene(config[g]) for g in genes]
        self.input = args.short, args.long
        self.ref = args.reference
        if not args.sample_name:
            bam = load_bam(args.short, args.reference)
            rg = bam.header.to_dict().get("RG", [])
            sample_name = rg[0]["SM"] if rg and "SM" in rg[0] else None
            if not sample_name:
                raise Exception(
                    f"SM tag not found in {args.short}. Please explicitly set --sample_name."
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
        # Set global logging level
        set_global_log_level(args.log_level)
        # Update the logger level since it was created before global level was set
        logger.setLevel(getattr(logging, args.log_level.upper()))

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
                "log_level": args.log_level,
            }
        )
        self.print_param()
        self.init_read_data()

    def _validate_input_files(self) -> None:
        """Validate that all required input files and their indices exist."""
        short_bam, long_bam = self.input

        # Check short read BAM file
        if not os.path.exists(short_bam):
            logger.error(f"Short read BAM file does not exist: {short_bam}")
            sys.exit(1)

        # Check short read BAM index
        if short_bam.endswith(".bam"):
            bam_index = short_bam + ".bai"
            if not os.path.exists(bam_index):
                logger.error(f"Short read BAM index file does not exist: {bam_index}")
                sys.exit(1)
        elif short_bam.endswith(".cram"):
            cram_index = short_bam + ".crai"
            if not os.path.exists(cram_index):
                logger.error(f"Short read CRAM index file does not exist: {cram_index}")
                sys.exit(1)

        # Check long read BAM file if provided
        if long_bam:
            if not os.path.exists(long_bam):
                logger.error(f"Long read BAM file does not exist: {long_bam}")
                sys.exit(1)

            # Check long read BAM index
            if long_bam.endswith(".bam"):
                bam_index = long_bam + ".bai"
                if not os.path.exists(bam_index):
                    logger.error(
                        f"Long read BAM index file does not exist: {bam_index}"
                    )
                    sys.exit(1)
            elif long_bam.endswith(".cram"):
                cram_index = long_bam + ".crai"
                if not os.path.exists(cram_index):
                    logger.error(
                        f"Long read CRAM index file does not exist: {cram_index}"
                    )
                    sys.exit(1)

        # Check reference file
        if not os.path.exists(self.ref):
            logger.error(f"Reference file does not exist: {self.ref}")
            sys.exit(1)

        # Check reference index
        ref_index = self.ref + ".fai"
        if not os.path.exists(ref_index):
            logger.error(f"Reference index file does not exist: {ref_index}")
            sys.exit(1)

        # Check model bundles
        sr_model_path = self.params["sr_model"]
        actual_sr_model_path = get_data_file(sr_model_path)
        if not actual_sr_model_path:
            logger.error(f"Short read model file does not exist: {sr_model_path}")
            sys.exit(1)

        if long_bam and "lr_model" in self.params:
            lr_model_path = self.params["lr_model"]
            actual_lr_model_path = get_data_file(lr_model_path)
            if not actual_lr_model_path:
                logger.error(f"Long read model file does not exist: {lr_model_path}")
                sys.exit(1)
        elif long_bam:
            logger.error("Long read model is required when long read BAM is provided")
            sys.exit(1)

        logger.info("All required input files and indices validated successfully")

    def call(self) -> None:
        self.load_params()

        # Check if required files exist before processing
        self._validate_input_files()

        result = {}
        for g in self.genes:
            result[g.gene_names[0]] = self.process_gene(g)
        self.save_data(
            result, f"{self.params['outdir']}/{self.params['sample_name']}.yaml"
        )
        if not self.params["keep_tmp"]:
            shutil.rmtree(self.params["tmpdir"])
