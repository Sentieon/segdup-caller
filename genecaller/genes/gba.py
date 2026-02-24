"""GBA gene with enhanced interpretation for GBAP1 deletions and gene conversions."""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.conversion_detector import GeneConversionDetector
import vcflib


class GBA(Gene):
    """GBA gene with enhanced interpretation for GBAP1 deletions and gene conversions."""

    def __init__(self, cfg: dict, ref_file: str, **kwargs) -> None:
        super().__init__(cfg, ref_file, **kwargs)
        # Load conversion variants from database instead of hardcoding
        self.conversion_variants = self._load_conversion_variants_from_db()

    def _load_conversion_variants_from_db(self) -> Dict[str, Dict[str, Any]]:
        """
        Load GBA-specific conversion variants from the conversion database.

        Returns a dictionary mapping variant names to their genomic coordinates and alleles.
        Only includes single-position variants (not complex events like RecTL or fusions).
        """
        if not self.conversion_db:
            self.logger.warning(
                "No conversion database configured; variant checking will be disabled"
            )
            return {}

        try:
            known_conversions = GeneConversionDetector.load_known_conversions(
                self.conversion_db
            )
        except Exception as e:
            self.logger.error(f"Failed to load conversion database: {e}")
            return {}

        variants = {}
        for kc in known_conversions:
            # Only load GBA1 variants with valid ref/alt alleles
            if kc.recipient_gene != "GBA1":
                continue
            if kc.ref == "NA" or kc.alt == "NA":
                continue  # Skip complex events without single ref/alt

            # Parse recipient_locus to extract chrom and position
            # Format: "chr1:155235252" or "chr1:155235749-155235805"
            try:
                locus_parts = kc.recipient_locus.split(":")
                if len(locus_parts) != 2:
                    continue

                chrom = locus_parts[0]
                pos_str = locus_parts[1]

                # For range loci (e.g., deletions), use the start position
                if "-" in pos_str:
                    pos = int(pos_str.split("-")[0])
                else:
                    pos = int(pos_str)

                variants[kc.event_name] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": kc.ref,
                    "alt": kc.alt,
                }

            except (ValueError, IndexError) as e:
                self.logger.warning(
                    f"Failed to parse locus for {kc.event_name}: {kc.recipient_locus} ({e})"
                )
                continue

        self.logger.info(
            f"Loaded {len(variants)} conversion variants from database for GBA1"
        )
        return variants

    def _check_conversion_variants(self, vcf_path: str) -> List[Dict[str, str]]:
        found_conversions = []

        if not os.path.exists(vcf_path):
            return found_conversions

        if not self.conversion_variants:
            return found_conversions

        try:
            vcf = vcflib.VCF(vcf_path)
            for v in vcf:
                if not v.alt:
                    continue

                chrom = v.chrom if v.chrom.startswith("chr") else f"chr{v.chrom}"
                pos = v.pos + 1

                for var_name, var_info in self.conversion_variants.items():
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
