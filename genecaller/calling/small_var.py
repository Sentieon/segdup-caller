"""Small-variant calling: DNAscope, Haplotyper, Genotyper, and helpers.

Functions take a Bam instance as a state-holder for tool paths, reference, and
BAM file paths. They do not mutate the Bam's variant-call-irrelevant state.

When DNAscope is invoked with a population model bundle (one whose
bundle_info.json declares a SentieonVcfID), a transfer step is run between
DNAscope and DNAModelApply that annotates the raw VCF with INFO fields from
the bundled population VCF. The bundled pop VCF is looked up by
SentieonVcfID via data/pop_vcfs.yaml, then by the detected reference build
(hg38/hg19/b37). Non-population bundles run the existing two-step pipeline
unchanged.
"""

import os
import re
import subprocess
import sys
from importlib import resources
from typing import Any, Dict, List, Optional

import yaml

from .. import throttle
from ..logging import get_logger
from ..util import IntervalList, Reference, get_data_file, read_bundle_info

logger = get_logger(__name__)


# Caches (per process). The ProcessPoolExecutor in genecaller.process_genes
# means each worker pays the parse cost at most once; threads within a worker
# share the dict.
_POP_VCF_CACHE: Dict[tuple, Optional[str]] = {}  # (bundle_dir, build) -> pop_vcf path (None for non-pop)
_MERGE_RULES_CACHE: Dict[str, str] = {}         # pop_vcf path -> "ID:sum,ID:sum,..."
_MANIFEST_CACHE: Optional[Dict[str, dict]] = None
_BUILD_CACHE: Dict[str, str] = {}               # reference path -> genome build (hg38/hg19/b37)


def subset_input_vcf(
    bcftools: str,
    sentieon: str,
    input_vcf: str,
    region: str,
    output: str,
) -> None:
    cmd = (
        f"{bcftools} view -r {region} -f PASS,. {input_vcf} "
        f"| {sentieon} util vcfconvert - {output}"
    )
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    if result.returncode:
        raise Exception(result.stderr)


def _load_pop_vcf_manifest() -> Dict[str, dict]:
    """Load data/pop_vcfs.yaml. Returns {} when the file is absent or empty."""
    global _MANIFEST_CACHE
    if _MANIFEST_CACHE is not None:
        return _MANIFEST_CACHE
    path = get_data_file("pop_vcfs.yaml")
    if not path or not os.path.exists(path):
        _MANIFEST_CACHE = {}
        return _MANIFEST_CACHE
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    _MANIFEST_CACHE = data
    return _MANIFEST_CACHE


def _build_for_ref(ref: str) -> str:
    """Detect (and cache) the genome build for a reference FASTA."""
    b = _BUILD_CACHE.get(ref)
    if b is None:
        b = Reference(ref).build()
        _BUILD_CACHE[ref] = b
    return b


def _entry_pop_vcf(entry: dict, build: str) -> Optional[str]:
    """Relative pop_vcf path for ``build`` from a manifest entry, or None.

    Supports the per-build mapping ``pop_vcf: {hg38: ..., hg19: ...}`` and the
    legacy flat ``pop_vcf: <path>`` (which applies to its declared ``assembly``,
    defaulting to hg38). hg19 and b37 share GRCh37 coordinates but differ in
    contig naming, so they require separate entries.
    """
    pv = entry.get("pop_vcf")
    if isinstance(pv, dict):
        return pv.get(build)
    if pv:
        return pv if entry.get("assembly", "hg38") == build else None
    return None


def _resolve_pop_vcf_for_model(model_path: str, build: str) -> Optional[str]:
    """Return the bundled pop_vcf path for this model bundle and build, or None.

    Non-population bundles (no bundle_info.json or no SentieonVcfID) return
    None and the caller should run the existing flow unchanged.
    """
    bundle_dir = os.path.dirname(model_path)
    cache_key = (bundle_dir, build)
    if cache_key in _POP_VCF_CACHE:
        return _POP_VCF_CACHE[cache_key]
    info = read_bundle_info(bundle_dir)
    vcf_id = info.get("SentieonVcfID")
    if not vcf_id:
        _POP_VCF_CACHE[cache_key] = None
        return None
    entry = _load_pop_vcf_manifest().get(vcf_id)
    if entry is None:
        # Startup validation should have caught this; defensive only.
        raise RuntimeError(
            f"Pop VCF manifest is missing an entry for SentieonVcfID "
            f"'{vcf_id}' (required by model bundle '{bundle_dir}'). "
            f"Update segdup-caller or check data/pop_vcfs.yaml."
        )
    rel = _entry_pop_vcf(entry, build)
    if rel is None:
        raise RuntimeError(
            f"Pop VCF manifest entry for SentieonVcfID '{vcf_id}' has no "
            f"pop_vcf for reference build '{build}'. Add a '{build}' entry "
            f"under its 'pop_vcf' in data/pop_vcfs.yaml."
        )
    path = get_data_file(rel)
    if not path or not os.path.exists(path):
        raise RuntimeError(
            f"Pop VCF file '{rel}' for SentieonVcfID '{vcf_id}' (build "
            f"'{build}') was not found in the segdup-caller package data."
        )
    logger.info(
        f"Population model detected (SentieonVcfID={vcf_id}, build={build}); "
        f"transfer step will use {path}"
    )
    _POP_VCF_CACHE[cache_key] = path
    return path


def validate_pop_vcf_coverage(bundle_paths: List[str], build: str) -> None:
    """Confirm every population-flavor bundle has a manifest entry + file for ``build``.

    Bundles without bundle_info.json or without SentieonVcfID are
    non-population models and skipped. On any failure, exit(2) with a clear
    message before any per-gene processing begins.
    """
    manifest = _load_pop_vcf_manifest()
    for bundle in bundle_paths:
        info = read_bundle_info(bundle)
        vcf_id = info.get("SentieonVcfID")
        if not vcf_id:
            continue
        entry = manifest.get(vcf_id)
        if entry is None:
            logger.error(
                f"Model bundle '{bundle}' requires population VCF "
                f"'{vcf_id}' but no matching pop_vcf is bundled with "
                f"segdup-caller. Please update segdup-caller."
            )
            sys.exit(2)
        rel = _entry_pop_vcf(entry, build)
        if rel is None:
            logger.error(
                f"Model bundle '{bundle}' requires population VCF '{vcf_id}' "
                f"for reference build '{build}', but its manifest entry has no "
                f"'{build}' pop_vcf. Add one under 'pop_vcf' in data/pop_vcfs.yaml."
            )
            sys.exit(2)
        path = get_data_file(rel)
        if not path or not os.path.exists(path):
            logger.error(
                f"Manifest entry for '{vcf_id}' (build '{build}') points at "
                f"'{rel}' which is missing from the "
                f"segdup-caller package data."
            )
            sys.exit(2)
        if not os.path.exists(path + ".tbi"):
            logger.error(
                f"Pop VCF '{path}' is missing its tabix index "
                f"({path}.tbi)."
            )
            sys.exit(2)


def _merge_rules_for(pop_vcf: str, bcftools: str) -> str:
    """Derive bcftools merge-rules from the pop VCF header.

    Mirrors sentieon-cli/transfer.py: every ##INFO line with Number=A gets a
    "<ID>:sum" rule. Cached per pop_vcf path.
    """
    if pop_vcf in _MERGE_RULES_CACHE:
        return _MERGE_RULES_CACHE[pop_vcf]
    cmd = [bcftools, "view", "-h", pop_vcf]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        raise Exception(p.stderr)
    kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)')
    id_fields: List[str] = []
    for line in p.stdout.split("\n"):
        if not line.startswith("##INFO"):
            continue
        if ",Number=A" not in line:
            continue
        s = line.index("<")
        e = line.index(">")
        d = dict(kvpat.findall(line[s + 1 : e]))  # noqa: E203
        id_fields.append(d["ID"])
    if not id_fields:
        # No mergeable Number=A INFO fields; merge will be a no-op of sorts.
        rules = ""
    else:
        rules = ",".join(f"{x}:sum" for x in id_fields)
    _MERGE_RULES_CACHE[pop_vcf] = rules
    return rules


def _write_region_bed(region: str, bed_path: str) -> None:
    """Write BED rows (0-based half-open) for a comma-separated region string."""
    itv = IntervalList(region=region)
    with open(bed_path, "w") as fh:
        for chrom, coords in itv.regions.items():
            i = 0
            while i < len(coords):
                start, end = coords[i], coords[i + 1]
                fh.write(f"{chrom}\t{start}\t{end}\n")
                i += 2


def _trimalt_path() -> str:
    """Resolve the bundled trimalt.py path. Assumes unzipped install."""
    return str(resources.files("genecaller.scripts").joinpath("trimalt.py"))


def _transfer_pop_annotations(
    bam: Any,
    raw_vcf: str,
    pop_vcf: str,
    region: str,
    out_vcf: str,
) -> str:
    """Annotate raw DNAscope VCF with pop_vcf INFO fields. Returns out_vcf.

    The merge + trimalt steps match sentieon-cli's transfer pipeline
    (sentieon_cli/transfer.py + cmd_bcftools_merge_trim) — flag order and
    dynamic merge rules — so DNAModelApply sees the same annotations the
    model was trained against. The final write uses `sentieon util
    vcfconvert` (which also creates the .tbi) instead of
    `bcftools view -W=tbi`, matching the rest of segdup-caller and avoiding
    a bcftools 1.18+ requirement just for the `-W=tbi` option. No sharding:
    each call_dnascope already targets a small gene region.

    Runs the shell pipeline with `pipefail` so a failure in `bcftools
    merge` or `trimalt` surfaces an error here rather than yielding an
    empty transfer VCF that would only fail downstream in `DNAModelApply`.
    `--regions-overlap pos` requires bcftools 1.16+ (enforced at startup
    via `CALLING_MIN_VERSIONS`).
    """
    bed_path = out_vcf + ".region.bed"
    _write_region_bed(region, bed_path)
    merge_rules = _merge_rules_for(pop_vcf, bam.bcftools)
    merge_rules_opt = f"-i '{merge_rules}'" if merge_rules else ""
    trimalt = _trimalt_path()
    cmd = (
        f"set -o pipefail; "
        f"{bam.bcftools} merge --regions-file {bed_path} --no-version "
        f"--regions-overlap pos -m all {merge_rules_opt} "
        f"{raw_vcf} {pop_vcf} "
        f"| {sys.executable} {trimalt} "
        f"| {bam.sentieon} util vcfconvert - {out_vcf}"
    )
    bam.logger.debug(f"Running command: {cmd}")
    result = subprocess.run(
        cmd, capture_output=True, text=True, shell=True, executable="/bin/bash"
    )
    if result.returncode:
        raise Exception(result.stderr)
    return out_vcf


def call_variant(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    if (
        param.get("input_vcf")
        and param.get("ploidy") == 2
        and param.get("bam_type") == "orig"
        and param.get("seq_key") == "short_read"
        and not param.get("given")
    ):
        bam.logger.info(
            f"Using --input_vcf for {param['region']} (skipping DNAscope diploid call)"
        )
        subset_input_vcf(
            bam.bcftools,
            bam.sentieon,
            param["input_vcf"],
            param["region"],
            output,
        )
        return
    if "ploidy" in param and param["ploidy"] != 2:
        param["algo"] = "Haplotyper"
    call_dnascope(bam, output, param)


def call_dnascope(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    driver_opt = param.get("driver_opt", "") + f" -i {bam_path} -r {bam.ref}"
    algo_opt = param.get("algo_opt", "")
    apply_model = ""
    algo = "DNAscope" if "algo" not in param else param["algo"]
    ploidy = param.get("ploidy", 2)
    algo_opt += f" --ploidy {ploidy}"
    given = param.get("given", "")
    if given and os.path.exists(given):
        algo_opt += f" --given {given} -d {given}"
    if "model" in param and algo == "DNAscope" and ploidy == 2 and not given:
        algo_opt += f" --model {param['model']}"
        apply_model = param["model"]
    if "region" in param:
        driver_opt += f" --interval {param['region']}"
    if "dbsnp" in param:
        algo_opt += f" -d {param['dbsnp']}"
    if apply_model:
        tmp_output = output.replace("vcf", "tmp.vcf")
    else:
        tmp_output = output
    cmd = f"{bam.sentieon} driver {driver_opt} --algo {algo} {algo_opt} {tmp_output}"
    bam.logger.debug(f"Running command: {cmd}")
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    if result.returncode:
        raise Exception(result.stderr)
    if apply_model:
        pop_vcf = _resolve_pop_vcf_for_model(apply_model, _build_for_ref(bam.ref))
        apply_input = tmp_output
        if pop_vcf:
            transfer_output = output.replace(".vcf", ".transfer.vcf")
            apply_input = _transfer_pop_annotations(
                bam, tmp_output, pop_vcf, param["region"], transfer_output
            )
        cmd = (
            f"{bam.sentieon} driver -r {bam.ref} --algo DNAModelApply "
            f"-v {apply_input} --model {apply_model} {output}"
        )
        bam.logger.debug(f"Running command: {cmd}")
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if result.returncode:
            raise Exception(result.stderr)


def call_genotyper(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    driver_opt = param.get("driver_opt", "") + f" -i {bam_path} -r {bam.ref}"
    algo_opt = param.get("algo_opt", "")
    if "region" in param:
        driver_opt += f" --interval {param['region']}"
    if "ploidy" in param:
        algo_opt += f" --ploidy {param['ploidy']}"
    if "given" in param:
        algo_opt += f" --given {param['given']} -d {param['given']}"
    if "dbsnp" in param:
        algo_opt += f" -d {param['dbsnp']}"
    cmd = f"{bam.sentieon} driver {driver_opt} --algo Genotyper {algo_opt} {output}"
    bam.logger.debug(f"Running command: {cmd}")
    with throttle.slot():
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    if result.returncode:
        raise Exception(result.stderr)
