"""VCF phasing dispatch.

Auto (default): ploidy==1 → symlink; short reads → k-modes; long reads → anchor.
The in-process phasers preserve the input genotype (only add phase, never
re-genotype) at every ploidy. `method` also accepts "kmodes", "anchor", and
"whatshap"; whatshap is an opt-in fallback resolved from PATH — not a dependency,
never auto-invoked.
"""

import os
import subprocess
from typing import Any, Optional


def phase_vcf(
    bam: Any,
    in_vcf: str,
    out_vcf: str,
    ploidy: int = 2,
    method: Optional[str] = None,
) -> None:
    if bam.use_existing and os.path.exists(out_vcf) and os.path.getsize(out_vcf):
        return

    if ploidy == 1:
        if os.path.exists(out_vcf):
            os.remove(out_vcf)
        os.symlink(os.path.abspath(in_vcf), out_vcf)
        in_vcf_idx = in_vcf + ".tbi"
        out_vcf_idx = out_vcf + ".tbi"
        if os.path.exists(in_vcf_idx):
            if os.path.exists(out_vcf_idx):
                os.remove(out_vcf_idx)
            os.symlink(os.path.abspath(in_vcf_idx), out_vcf_idx)
        else:
            cmd = f"{bam.sentieon} util vcfindex {out_vcf}"
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            result.check_returncode()
        return

    method = (method or "auto").lower()
    auto = method == "auto"
    # auto: short reads → k-modes, long reads → anchor (any ploidy >= 2)
    use_kmodes = method == "kmodes" or (auto and not bam.long_read)
    use_anchor = method == "anchor" or (auto and bam.long_read)

    if use_kmodes:
        from .polyphase_short import polyphase_short_reads
        polyphase_short_reads(bam, in_vcf, out_vcf, ploidy)
        return
    if use_anchor:
        from .polyphase_long import polyphase_long_reads
        polyphase_long_reads(bam, in_vcf, out_vcf, ploidy)
        return

    # method == "whatshap": opt-in fallback; whatshap resolved from PATH, not a dependency.
    whatshap = getattr(bam, "whatshap", None) or "whatshap"
    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    if ploidy == 2:
        cmd = (
            f"{whatshap} phase -o {out_vcf} {in_vcf} {bam_path} "
            f"--reference {bam.ref}"
        )
    else:
        cmd = (
            f"{whatshap} polyphase -o {out_vcf} {in_vcf} {bam_path} "
            f"-p {ploidy} --reference {bam.ref} -B 5 -t {bam.threads}"
        )
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()
    cmd = f"{bam.sentieon} util vcfindex {out_vcf}"
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()
