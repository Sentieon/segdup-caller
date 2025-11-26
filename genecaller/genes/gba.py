"""GBA gene with enhanced interpretation for GBAP1 deletions and gene conversions."""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
import vcflib


class GBA(Gene):
    """GBA gene with enhanced interpretation for GBAP1 deletions and gene conversions."""

    CONVERSION_VARIANTS = {
        "L483P": {"chrom": "chr1", "pos": 155235252, "ref": "A", "alt": "G"},
        "A495P": {"chrom": "chr1", "pos": 155235217, "ref": "C", "alt": "G"},
        "V499V": {"chrom": "chr1", "pos": 155235203, "ref": "C", "alt": "G"},
        "c.1263del": {
            "chrom": "chr1",
            "pos": 155235749,
            "ref": "GGGACTGTCGACAAAGTTACGCACCCAATTGGGTCCTCCTTCGGGGTTCAGGGCAA",
            "alt": "G",
        },
        "D409H": {"chrom": "chr1", "pos": 155235727, "ref": "C", "alt": "G"},
        "rs708606": {"chrom": "chr1", "pos": 155234903, "ref": "C", "alt": "T"},
    }

    def _check_conversion_variants(self, vcf_path: str) -> List[Dict[str, str]]:
        found_conversions = []

        if not os.path.exists(vcf_path):
            return found_conversions

        try:
            vcf = vcflib.VCF(vcf_path)
            for v in vcf:
                if not v.alt:
                    continue

                chrom = v.chrom if v.chrom.startswith("chr") else f"chr{v.chrom}"
                pos = v.pos + 1

                for var_name, var_info in self.CONVERSION_VARIANTS.items():
                    if (
                        chrom == var_info["chrom"]
                        and pos == var_info["pos"]
                        and v.ref == var_info["ref"]
                    ):
                        alts = v.alt if isinstance(v.alt, list) else [v.alt]
                        if var_info["alt"] not in alts:
                            continue

                        gt_str = v.samples[0].get("GT", "./.")
                        if not gt_str or gt_str in ["./.", "."]:
                            continue

                        gt_parts = gt_str.replace("|", "/").split("/")
                        try:
                            gt_indices = [int(x) for x in gt_parts if x != "."]
                        except (ValueError, TypeError):
                            continue

                        if not gt_indices or all(x == 0 for x in gt_indices):
                            continue

                        alt_count = sum(1 for x in gt_indices if x > 0)

                        if alt_count == 1:
                            allele_type = "mono-allelic"
                        elif alt_count >= 2:
                            allele_type = "bi-allelic"
                        else:
                            continue

                        found_conversions.append(
                            {
                                "variant": var_name,
                                "type": allele_type,
                                "pos": pos,
                                "ref": v.ref,
                                "alt": var_info["alt"],
                            }
                        )

            vcf.close()
        except Exception:
            pass

        return found_conversions

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()

        interpretation_lines = []

        gba1_cn = result_data["Copy numbers"].get("GBA1", 2)
        gbap1_cn = result_data["Copy numbers"].get("GBAP1", 2)

        interpretation_lines.append(f"GBA1: CN={gba1_cn}")
        interpretation_lines.append(f"GBAP1: CN={gbap1_cn}")

        vcf_path = result_data.get("Variants", "")
        conversion_variants = self._check_conversion_variants(vcf_path)

        if conversion_variants:
            interpretation_lines.append(
                f"Potential GBAP1→GBA1 conversion variants detected:"
            )
            for conv in conversion_variants:
                interpretation_lines.append(
                    f"  {conv['type']} {conv['variant']} "
                    f"(possible GBAP1→GBA1 gene-conversion)"
                )

        # 2. Check for gene conversions and fusion events
        # For GBA, only positive allele changes matter (GBAP1→GBA1 direction)
        relevant_conversions = [
            c for c in result_data["gene_conversions"] if c["converted_alleles"] > 0
        ]

        if relevant_conversions:
            interpretation_lines.append(
                f"Gene conversion/fusion events detected: {len(relevant_conversions)}"
            )

            for i, conv in enumerate(relevant_conversions, 1):
                region = conv["region"]
                converted_alleles = conv["converted_alleles"]
                conversion_sites = conv.get("conversion_sites", [])
                is_fusion = conv.get("is_fusion_candidate", False)
                event_type = "FUSION" if is_fusion else "CONVERSION"

                interpretation_lines.append(f"  Event {i} ({event_type}): {region}")
                interpretation_lines.append(
                    f"    Converted alleles: {converted_alleles}"
                )
                interpretation_lines.append(
                    f"    Conversion sites: {', '.join(conversion_sites) if conversion_sites else 'None'}"
                )

                # Add interpretation if available
                interp = conv.get("interpretation")
                if isinstance(interp, list):
                    # Known conversion events matched
                    for match in interp:
                        event_name = match["event_name"]
                        # Handle position-based fusion matches (no signature variants)
                        if "matched_variants" in match:
                            interpretation_lines.append(
                                f"    Interpretation: {event_name} "
                                f"(match score: {match['match_score']:.1%}, "
                                f"variants: {', '.join(match['matched_variants'])})"
                            )
                        else:
                            interpretation_lines.append(
                                f"    Interpretation: {event_name} "
                                f"(position-based match)"
                            )
                elif isinstance(interp, str):
                    # Novel or no database
                    interpretation_lines.append(f"    Interpretation: {interp}")

                # Add clinical significance for fusion events
                if is_fusion and event_type == "FUSION":
                    interpretation_lines.append(
                        f"    Clinical note: Gene fusion event (GBA1::GBAP1) detected. "
                        f"This may result in altered protein structure and function."
                    )
        elif not conversion_variants:
            interpretation_lines.append("No gene conversion or fusion events detected")

        # Add cn_interpretation field
        result_data["cn_interpretation"] = "\n".join(interpretation_lines)

        return result_data
