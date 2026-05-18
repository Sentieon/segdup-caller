"""CNV/depth calling: CNVscope wrapper and depth/segment caching.

`call_cnvscope`, `call_depth`, and `trim_to_region_vcf` are pure helpers
taking a Bam instance as a state-holder. `call_depth` mutates Bam's depth/
segment caches as part of its caching contract.
"""

import os
import subprocess
from typing import Any, Dict

import pandas as pd
import vcflib

from ..util import IntervalList


def call_cnvscope(bam: Any, output: str, param: Dict[str, Any]) -> None:
    if bam.use_existing and os.path.exists(output) and os.path.getsize(output):
        return
    driver_opt = param.get("driver_opt", "") + f" -i {bam.bam} -r {bam.ref}"
    algo_opt = param.get("algo_opt", "")
    if "region" in param:
        driver_opt += f" --interval {param['region']}"
    if "output_norm" in param:
        algo_opt += f" --output_normalized {param['output_norm']}"
    if "model" in param:
        algo_opt += f" --model {param['model']}"
    else:
        raise Exception("Model is not specified for CNVscope.")
    cmd = f"{bam.sentieon} driver {driver_opt} --algo CNVscope {algo_opt} {output}"
    bam.logger.debug(f"Running command: {cmd}")
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()


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
