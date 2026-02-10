import os
import sys
import argparse
import yaml
import shutil
import logging
import time
import resource
from typing import Dict, Any
from genecaller.logging import get_logger, set_global_log_level
from genecaller.util import (
    get_cmd,
    IntervalList,
    load_bam,
    get_data_file,
    get_software_versions,
)
from genecaller.bam_process import Bam
from genecaller.gene import Gene
from genecaller.genes import get_gene_class
from genecaller import __version__
import tempfile
import copy
from packaging.version import Version
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from datetime import datetime

logger = get_logger(__name__)


CALLING_MIN_VERSIONS = {
    "sentieon driver": Version("202503"),
    "whatshap": Version("2.3"),
    "bcftools": Version("1.10"),
    "samtools": Version("1.16"),
}


def process_single_gene(gene_data: tuple) -> Dict:
    gene, input_files, ref, base_params, gene_index = gene_data

    read_data = _init_read_data(input_files, ref, base_params, gene.gene_names[0])

    gene_config = gene.config
    params = copy.deepcopy(base_params)
    if gene_config:
        params.update(gene_config)

    # Create gene-specific prefixes to avoid file conflicts
    gene_name = gene.gene_names[0]

    # Process the gene
    gene_logger = get_logger(f"Gene.{gene_name}")
    gene_logger.info(f"Processing gene {', '.join(gene.gene_names)}")

    try:
        # Apply gene-specific dp_norm if configured
        if gene_config and "dp_norm" in gene_config:
            read_data["short_read"]["bam"].dp_norm = gene_config["dp_norm"]

        # Liftover
        gene_logger.info("Performing liftover...")
        _perform_liftover(gene, read_data, ref)

        # Configure variant calling parallelism for this gene
        gene.variant_call_threads = min(
            params.get("threads", 4), 8
        )  # Cap at 8 threads for variant calling

        # Copy number calling
        gene_logger.info("Call copy number...")
        if gene.config.get("cnv_conversion_iters", 1) > 1:
            # Use iterative CNV-conversion detection
            gene.call_cn_with_conversion_iteration(read_data, params)
        else:
            # Original behavior
            gene.call_cn(read_data, params)
            gene.merge_regions_by_cn()
            gene.detect_gene_conversions()

        # Variant calling
        if not params.get("skip_variant_call", False):
            gene_logger.info("Call small variants...")
            gene.call_small_var()
            gene.resolve_phased(params)
        else:
            gene_logger.info("Skipping small variant calling (--skip-variant-call)")

        # Prepare output
        gene_logger.info("Prepare final output")
        gene_data_result = gene.prepare_output()

        gene_logger.info(f"Completed processing gene {gene_name}")
        return gene_data_result

    except Exception as e:
        gene_logger.error(f"Error processing gene {gene_name}: {e}")
        raise


def _perform_liftover(gene: "Gene", read_data: dict, ref: str) -> None:
    """Helper function to perform liftover for a gene."""
    for rd in read_data.values():
        bam = rd["bam"]
        liftover_params = copy.deepcopy(rd["params"])
        liftover_params["prefix"] += ".liftover"
        liftover_fname = f"{liftover_params['prefix']}.bam"
        liftover_results = bam.liftover_combine(gene, liftover_fname)
        liftover_bam = Bam(liftover_fname, ref, liftover_params)
        rd["liftover"] = liftover_bam
        rd["liftover_valid"] = [
            IntervalList(region=gene.liftover_target_regions[i])
            .subtract(liftover_results[i]["failed_region"])
            .to_region_str()
            for i in range(len(gene.liftover_region_names))
        ]


def _init_read_data(input_files, ref, params, gene_name) -> dict:
    short, long = input_files
    sr_params = copy.deepcopy(params)
    sr_params["long_read"] = False
    sr_params["prefix"] = (
        f"{params['tmpdir']}/{params['sample_name']}.{gene_name}.{params['sr_prefix']}"
    )
    sr_params["model"] = params["sr_model"]
    short_bam = Bam(short, ref, sr_params)
    read_data = {"short_read": {"bam": short_bam, "params": sr_params}}
    if long:
        lr_params = copy.copy(params)
        lr_params["long_read"] = True
        lr_params["prefix"] = (
            f"{params['tmpdir']}/{params['sample_name']}.{gene_name}.{params['lr_prefix']}"
        )
        if "lr_model" in params:
            lr_params["model"] = params["lr_model"]
        long_bam = Bam(long, ref, lr_params)
        read_data["long_read"] = {"bam": long_bam, "params": lr_params}
    return read_data


class GeneCaller:
    tools = []

    def __init__(self) -> None:
        self.tools = {}
        self.genes = []
        self.input = []
        self.read_data = {}
        self.params = {}
        self.ref = ""
        self.reset_tmpdir = False

    def setup_paths(self) -> dict[str, str]:
        cmds = {}
        for cmd, version in CALLING_MIN_VERSIONS.items():
            cmd_path = get_cmd(cmd, version)
            if not cmd_path:
                sys.exit(1)
            cmds[cmd.split()[0]] = cmd_path
        return cmds

    def reset_read_data(self, params: Dict[str, Any]) -> None:
        for r in self.read_data.values():
            r["bam"].reset()
        self.params = params

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
                # Configure YAML to use literal block style (|) for multiline strings
                def str_representer(dumper, data):
                    if "\n" in data:
                        return dumper.represent_scalar(
                            "tag:yaml.org,2002:str", data, style="|"
                        )
                    return dumper.represent_scalar("tag:yaml.org,2002:str", data)

                yaml.add_representer(str, str_representer)
                serialized = yaml.dump(
                    data, default_flow_style=False, allow_unicode=True
                )
            fout.write(serialized)

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
        logger.info(f" * short read model: {self.params['sr_model']}")
        if self.params.get("lr_model", None):
            logger.info(f" * long read model: {self.params['lr_model']}")
        logger.info(f" * Parallel threads: {self.params['threads']}")
        logger.info(f" * Parallel workers: {self.params['workers']}")
        logger.info(f" * Output directory: {self.params['outdir']}")
        logger.info(f" * Sample name: {self.params['sample_name']}")

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
            help=f"List of genes to be called (comma separated). If not specified, all supported genes will be called. \nSupported genes: {', '.join(default_genes)}",
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
            "--threads",
            "-t",
            type=int,
            default=cpu_count(),
            help="Number of parallel threads for gene processing (default: number of CPU count)",
        )
        parser.add_argument(
            "--workers",
            "-w",
            type=int,
            default=4,
            help="Number of genes processed concurrently (default: 4)",
        )
        parser.add_argument(
            "--log_level",
            choices=["DEBUG", "INFO", "WARNING", "ERROR"],
            default="INFO",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--skip-variant-call",
            action="store_true",
            default=False,
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

        # If no genes specified, use all available genes
        if args.genes:
            genes = [s.strip().upper() for s in args.genes.split(",")]
        else:
            genes = default_genes
            logger.info(
                f"No genes specified, will process all available genes: {', '.join(genes)}"
            )

        missing = [g for g in genes if g not in config]
        if missing:
            raise Exception(
                f"Genes {', '.join(missing)} are not supported at this moment. Below is a list of supported regions: \n{', '.join(default_genes)}"
            )
        Gene.cn_prior_std = config["main"].get("cn_prior_std", 2)

        # Set conversion database path if specified in config
        conversion_db_path = config["main"].get("conversion_db", None)
        if conversion_db_path:
            Gene.conversion_db = get_data_file(conversion_db_path)  # type: ignore

        # Instantiate genes using gene-specific classes if available
        self.input = args.short, args.long
        self.ref = args.reference
        self.genes = [get_gene_class(g)(config[g], self.ref) for g in genes]

        if args.threads > cpu_count():
            logger.warning(
                f"Requested {args.threads} threads, but only {cpu_count()} CPU cores available. Using {cpu_count()} threads."
            )
            args.threads = cpu_count()
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
            os.environ["SENTIEON_TMPDIR"] = (
                tmpdir  # this is needed for parallel jobs of Sentieon not conflicting each other
            )
            self.reset_tmpdir = True
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
                "threads": args.threads,
                "workers": args.workers,
                "skip_variant_call": args.skip_variant_call,
            }
        )
        self.args = args
        self.print_param()

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

    def process_genes(self) -> Dict[str, Any]:
        # Prepare task arguments for each gene
        gene_tasks = []
        for i, gene in enumerate(self.genes):
            gene_task = (
                gene,  # Gene object
                self.input,  # (short_bam, long_bam)
                self.ref,  # Reference file
                self.params,  # Base parameters
                i + 1,  # Gene index for logging
            )
            gene_tasks.append(gene_task)

        results = {}
        max_workers = min(self.params["workers"], len(self.genes))

        if max_workers > 1:
            # Use a process pool for parallel processing
            try:
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    # submit task
                    future_results = {
                        executor.submit(process_single_gene, task): task[0].gene_names[
                            0
                        ]
                        for task in gene_tasks
                    }
                # Collect results
                for returned, gene_name in future_results.items():
                    try:
                        gene_result = returned.result()
                        results[gene_name] = gene_result
                        logger.info(f"Completed {gene_name}")
                    except Exception as e:
                        logger.error(f"Failed to process {gene_name}: {e}")
                        raise

            except Exception as e:
                logger.error(f"Parallel processing failed: {e}")
                raise
        else:
            for task in gene_tasks:
                gene_name = task[0].gene_names[0]
                try:
                    gene_result = process_single_gene(task)
                    results[gene_name] = gene_result
                    logger.info(f"Completed {gene_name}")
                except Exception as e:
                    logger.error(f"Failed to process {gene_name}: {e}")
                    raise

        return results

    def call(self) -> None:
        self.load_params()
        self._validate_input_files()

        t0 = time.time()
        start_datetime = datetime.now()

        if len(self.genes) > 1 and self.params["threads"] * self.params["workers"] > 1:
            logger.info(
                f"Processing {len(self.genes)} genes in parallel using {self.params['workers']} workers and {self.params['threads']} threads"
            )
        result = self.process_genes()

        t1 = time.time()
        mm, ut, st = 0, 0, 0
        for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN):
            ru = resource.getrusage(who)
            mm += ru.ru_maxrss * 1024
            ut += ru.ru_utime
            st += ru.ru_stime

        logger.info(
            f"Resource usage: {mm} mem {ut:.3f} user {st:.3f} sys {t1 - t0:.3f} real"
        )

        config_section = {
            "software_versions": {
                "segdup-caller": __version__,
                **get_software_versions(CALLING_MIN_VERSIONS),
            },
            "run_date": start_datetime.strftime("%Y-%m-%d"),
            "run_time": start_datetime.strftime("%H:%M:%S"),
            "run_datetime": start_datetime.isoformat(),
            "resource_usage": {
                "memory_bytes": mm,
                "user_time_seconds": round(ut, 3),
                "system_time_seconds": round(st, 3),
                "wall_time_seconds": round(t1 - t0, 3),
            },
            "job_options": {
                "short_read_bam": self.args.short,
                "reference": self.args.reference,
                "sr_model": self.args.sr_model,
                "genes": self.args.genes if self.args.genes else "all",
                "outdir": self.args.outdir,
                "sample_name": self.params["sample_name"],
                "sr_prefix": self.args.sr_prefix,
                "lr_prefix": self.args.lr_prefix,
                "threads": self.args.threads,
                "keep_temp": self.args.keep_temp,
                "log_level": self.args.log_level,
            },
        }

        if self.args.long:
            config_section["job_options"]["long_read_bam"] = self.args.long
            config_section["job_options"]["lr_model"] = self.args.lr_model
        if self.args.config:
            config_section["job_options"]["config"] = self.args.config

        result["config"] = config_section

        self.save_data(
            result, f"{self.params['outdir']}/{self.params['sample_name']}.yaml"
        )
        if not self.params["keep_tmp"]:
            shutil.rmtree(self.params["tmpdir"])

        if self.reset_tmpdir:
            os.environ.pop("SENTIEON_TMPDIR", None)
