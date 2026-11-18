"""CNV/depth calling: CNVscope wrapper and depth/segment caching.

`call_cnvscope`, `call_depth`, and `trim_to_region_vcf` are pure helpers
taking a Bam instance as a state-holder. `call_depth` mutates Bam's depth/
segment caches as part of its caching contract. `source_cnvscope` runs the
source CNVscope of every chromosome once, before the genes start; `call_depth`
then reads those outputs by path (`staged_run`).
"""

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import vcflib

from .. import throttle
from ..util import IntervalList

# output paths blanked: two runs with the same template are the same command
_OUT, _NORM = "@CNV_OUT@", "@CNV_NORM@"


def cnvscope_cmd(bam: Any, output: str, param: Dict[str, Any]) -> str:
    # a lifted Bam may read several BAMs at once (cnv_inputs); an input Bam reads bam.bam
    inputs = getattr(bam, "cnv_inputs", None) or [bam.bam]
    driver_opt = param.get("driver_opt", "") + "".join(f" -i {p}" for p in inputs)
    driver_opt += f" -r {bam.ref}"
    algo_opt = param.get("algo_opt", "")
    if "region" in param:
        driver_opt += f" --interval {param['region']}"
    if "output_norm" in param:
        algo_opt += f" --output_normalized {param['output_norm']}"
    if "model" in param:
        algo_opt += f" --model {param['model']}"
    else:
        raise Exception("Model is not specified for CNVscope.")
    return f"{bam.sentieon} driver {driver_opt} --algo CNVscope {algo_opt} {output}"


def cnvscope_template(bam: Any, param: Dict[str, Any]) -> str:
    return cnvscope_cmd(bam, _OUT, dict(param, output_norm=_NORM))


def _run(bam: Any, cmd: str) -> None:
    bam.logger.debug(f"Running command: {cmd}")
    with throttle.slot():
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()


def call_cnvscope(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    _run(bam, cnvscope_cmd(bam, output, param))


def _segment_param(bam: Any, chrom: str, output_norm: str = _NORM) -> Dict[str, Any]:
    # what a gene's segmentation passes call_cnvscope for one chromosome
    param = {"region": chrom, "output_norm": output_norm}
    if bam.cnv_model:
        param["model"] = bam.cnv_model
    return param


def source_cnvscope(
    bams: Dict[str, Any], chroms: List[str], out_dir: str, max_parallel: int
) -> Dict[Tuple[str, str], Tuple[str, str]]:
    """Each input's source CNVscope per chromosome: {(input tag, chromosome): (VCF, TSV)}.

    The one place that decides how the stage runs CNVscope; callers see only these outputs. It runs
    the genes' own command per input and chromosome, max_parallel at a time: one call over several
    chromosomes normalises them together and changes the depths. A failed run is left out, and
    its genes run their own.
    """

    def run(job):
        tag, bam, chrom = job
        vcf, norm = (os.path.join(out_dir, f"{tag}.{chrom}{x}") for x in (".vcf.gz", ".norm.tsv"))
        try:
            call_cnvscope(bam, vcf, _segment_param(bam, chrom, norm))
        except Exception as e:
            bam.logger.error(f"Source CNVscope on {chrom} failed; its genes run their own: {e}")
            return None
        return (tag, chrom), (vcf, norm)

    jobs = [(tag, bam, chrom) for tag, bam in bams.items() for chrom in chroms]
    with ThreadPoolExecutor(max_workers=max(1, min(max_parallel, len(jobs)))) as pool:
        return dict(r for r in pool.map(run, jobs) if r)


def staged_table(
    bams: Dict[str, Any], outputs: Dict[Tuple[str, str], Tuple[str, str]]
) -> Dict[str, Tuple[str, str]]:
    """{command template: (VCF, TSV)}: each staged output under the gene command it stands for."""
    return {
        cnvscope_template(bams[tag], _segment_param(bams[tag], chrom)): paths
        for (tag, chrom), paths in outputs.items()
    }


def staged_run(bam: Any, param: Dict[str, Any]) -> Optional[Tuple[str, str]]:
    """(VCF, TSV) the up-front stage wrote for this exact command, or None to run it here."""
    table = getattr(bam, "source_cnv", None)
    if table is None:  # no stage for this Bam (a lifted Bam, or a run without one)
        return None
    staged = table.get(cnvscope_template(bam, param))
    bam.logger.info(
        f"Source CNVscope for {param.get('region')}: "
        + (f"staged output {staged[0]}" if staged else "no staged output, running its own")
    )
    return staged


def call_depth(bam: Any, output: str, param: Dict[str, Any]) -> None:
    # Import locally to avoid a circular import with bam_process at module load.
    from ..bam_process import Segment, DepthSegment

    orig_region = None
    if "region" in param:
        orig_region = param["region"]
        chrs = set()
        for parts in param["region"].split(","):
            chrs.add(parts.split(":")[0])
        param["region"] = ",".join(chrs)
    prefix = output[: output.index(".vcf")] if ".vcf" in output else output
    output_norm = param.get("output_norm", f"{prefix}.norm.tsv")
    if orig_region:
        tmp_output = prefix + ".tmp" + output[output.index(".vcf") :]
        tmp_output_norm = prefix + ".tmp.norm.tsv"
    else:
        tmp_output = output
        tmp_output_norm = output_norm
    param["output_norm"] = tmp_output_norm
    staged = staged_run(bam, param)
    if staged:  # the up-front stage ran this exact command
        tmp_output, tmp_output_norm = staged
    else:
        call_cnvscope(bam, tmp_output, param)
    if orig_region:
        param["region"] = orig_region
    vcf = vcflib.VCF(tmp_output)
    depth_fields = ["contig", "start", "end", "MQ", "DP", "DP0"]
    for v in vcf:
        dp = v.samples[0]["DPS"]
        dp0 = v.samples[0]["DP0S"]
        bam.segments.setdefault(v.chrom, []).append(
            Segment(
                v.chrom,
                v.pos,
                v.end,
                v.info["CNT"],
                v.samples[0]["MQS"],
                dp[0] / bam.dp_norm,
                dp[1] / bam.dp_norm,
                dp0[0] / bam.dp_norm,
                dp0[0] / bam.dp_norm,
            )
        )
    depth_df = pd.read_csv(tmp_output_norm, sep="\t")[depth_fields].copy()
    depth_df["contig"] = depth_df["contig"].astype(str)
    depth_df["DP"] /= bam.dp_norm
    depth_df["DP0"] /= bam.dp_norm
    high_mq_depth_df = depth_df[depth_df["MQ"] > 50]
    bam.overall_depth_std = high_mq_depth_df["DP"].std()
    # Optimize DepthSegment creation using vectorized operations
    bam.depths = {}
    bam._depth_coords_cache = {}
    for chr, group in depth_df.groupby("contig"):
        # Use vectorized numpy operations instead of itertuples
        values = group.values
        chr_depths = [DepthSegment(*row) for row in values]
        bam.depths[chr] = chr_depths
        # Cache coordinates directly from numpy arrays for better performance
        bam._depth_coords_cache[chr] = (
            values[:, 1].tolist(),  # start column
            values[:, 2].tolist(),  # end column
        )
    if orig_region:
        orig_region_itv = IntervalList(region=orig_region)

        # Vectorized filtering using the same logic as the earlier optimization
        mask = pd.Series(False, index=depth_df.index)

        for chr in orig_region_itv.regions:
            chr_mask = depth_df["contig"] == chr
            if not chr_mask.any():
                continue

            chr_regions = orig_region_itv.regions[chr]
            chr_rows = depth_df[chr_mask]

            i = 0
            while i < len(chr_regions):
                start, end = chr_regions[i], chr_regions[i + 1]

                # Vectorized overlap check
                overlap_mask = (chr_rows["start"] < end) & (chr_rows["end"] > start)
                mask[chr_rows[overlap_mask].index] = True

                i += 2

        trim_depth_df = depth_df[mask]
        trim_depth_df.to_csv(output_norm, sep="\t", index=None)
        trim_to_region_vcf(bam, orig_region, tmp_output, output)


def trim_to_region_vcf(bam: Any, region: str, in_vcf: str, out_vcf: str) -> None:
    if region.endswith(".bed"):
        region_opt = f"-T {region}"
    else:
        region_opt = f"-r {region}"
    cmd = (
        f"{bam.bcftools} view {region_opt} {in_vcf} "
        f"| {bam.sentieon} util vcfconvert - {out_vcf}"
    )
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()
