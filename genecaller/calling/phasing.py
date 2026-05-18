"""VCF phasing via whatshap."""

import os
import subprocess
from typing import Any


def phase_vcf(bam: Any, in_vcf: str, out_vcf: str, ploidy: int = 2) -> None:
    if bam.use_existing and os.path.exists(out_vcf) and os.path.getsize(out_vcf):
        return

    # Haploid (ploidy=1) does not require phasing - create symlink instead
    if ploidy == 1:
        # Create symlink from input VCF to output VCF (phasing not needed)
        if os.path.exists(out_vcf):
            os.remove(out_vcf)
        os.symlink(os.path.abspath(in_vcf), out_vcf)
        # Also symlink the index if it exists
        in_vcf_idx = in_vcf + ".tbi"
        out_vcf_idx = out_vcf + ".tbi"
        if os.path.exists(in_vcf_idx):
            if os.path.exists(out_vcf_idx):
                os.remove(out_vcf_idx)
            os.symlink(os.path.abspath(in_vcf_idx), out_vcf_idx)
        else:
            # Index doesn't exist, create it
            cmd = f"{bam.sentieon} util vcfindex {out_vcf}"
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            result.check_returncode()
        return

    bam_path = bam.clipped_bam if bam.clipped_bam else bam.bam
    if ploidy == 2:
        cmd = (
            f"{bam.whatshap} phase -o {out_vcf} {in_vcf} {bam_path} "
            f"--reference {bam.ref}"
        )
    else:
        # For ploidy > 2, use polyphase
        cmd = (
            f"{bam.whatshap} polyphase -o {out_vcf} {in_vcf} {bam_path} "
            f"-p {ploidy} --reference {bam.ref} -B 5 -t {bam.threads}"
        )
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()
    cmd = f"{bam.sentieon} util vcfindex {out_vcf}"
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    result.check_returncode()
