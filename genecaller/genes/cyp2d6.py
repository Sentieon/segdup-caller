"""CYP2D6 gene with enhanced interpretation for gene conversions and fusion events."""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.util import get_data_file
import vcflib


class StarAlleleCaller:
    """Call CYP2D6 star alleles using phasing-aware consumption method."""

    def __init__(self, config: Dict[str, Any]):
        self.allele_defs = self._load_allele_definitions(config)
        self.backbone_vars = self._load_backbone_variants(config)

    def _load_allele_definitions(self, config: Dict[str, Any]) -> List[Dict[str, Any]]:
        definitions = []
        tsv_path = get_data_file(config["star_allele_file"])
        assert tsv_path is not None, (
            f"Star allele file not found: {config['star_allele_file']}"
        )
        with open(tsv_path, "r") as f:
            header = f.readline().strip().split("\t")
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < len(header):
                    continue
                allele_data = dict(zip(header, fields))
                allele_data["tier"] = int(allele_data["tier"])
                try:
                    allele_data["pos"] = int(allele_data["pos"])
                except ValueError:
                    pass
                allele_data["backbone"] = bool(int(allele_data["backbone"]))
                definitions.append(allele_data)
        return sorted(definitions, key=lambda x: x["tier"])

    def _load_backbone_variants(self, config: Dict[str, Any]) -> List[Dict[str, Any]]:
        backbone = []
        tsv_path = get_data_file(config["backbone_file"])
        assert tsv_path is not None, (
            f"Backbone file not found: {config['backbone_file']}"
        )
        with open(tsv_path, "r") as f:
            header = f.readline().strip().split("\t")
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < len(header):
                    continue
                var_data = dict(zip(header, fields))
                var_data["pos"] = int(var_data["pos"])
                backbone.append(var_data)
        return backbone

    def call_alleles(
        self, vcf_path: str, copy_number: int, conversions: List[Dict]
    ) -> Dict[str, Any]:
        if copy_number == 0:
            return {"alleles": ["*5", "*5"], "method": "CNV"}

        variants = self._extract_variants_from_vcf(vcf_path)

        structural_alleles = self._get_structural_alleles(conversions)
        for allele_name, count in structural_alleles.items():
            gt = "1/1" if count >= 2 else "0/1"
            key = f"STRUCTURAL:{allele_name}:.:."
            variants[key] = {"gt": gt, "count": count, "structural": True}

        variant_counts = {}
        for k, v in variants.items():
            variant_counts[k] = v.get("count", 1)

        called_alleles = []
        max_alleles = copy_number if copy_number > 0 else 2

        for allele_def in self.allele_defs:
            if len(called_alleles) >= max_alleles:
                break

            var_key = self._make_variant_key(allele_def)

            if variant_counts.get(var_key, 0) == 0:
                continue

            if allele_def["backbone"]:
                backbone_keys = [
                    self._make_variant_key(bv) for bv in self.backbone_vars
                ]
                backbone_present = all(
                    variant_counts.get(k, 0) > 0 for k in backbone_keys
                )

                if backbone_present and self._check_cis_phasing(
                    variants, var_key, backbone_keys
                ):
                    called_alleles.append(allele_def["allele"])
                    variant_counts[var_key] -= 1
                    for bb_key in backbone_keys:
                        variant_counts[bb_key] -= 1

                    if variant_counts[var_key] > 0:
                        backbone_enough = all(
                            variant_counts.get(k, 0) > 0 for k in backbone_keys
                        )
                        if backbone_enough and len(called_alleles) < max_alleles:
                            called_alleles.append(allele_def["allele"])
                            variant_counts[var_key] -= 1
                            for bb_key in backbone_keys:
                                variant_counts[bb_key] -= 1
            else:
                called_alleles.append(allele_def["allele"])
                variant_counts[var_key] -= 1

                if variant_counts[var_key] > 0 and len(called_alleles) < max_alleles:
                    called_alleles.append(allele_def["allele"])
                    variant_counts[var_key] -= 1

        while len(called_alleles) < max_alleles:
            called_alleles.append("*1")

        if copy_number == 1:
            called_alleles.append("*5")

        result = {
            "alleles": called_alleles,
            "method": "phased_consumption",
            "copy_number": copy_number,
        }

        if copy_number >= 3:
            result["note"] = "Duplication present - allele assignment may be ambiguous"

        return result

    def _extract_variants_from_vcf(self, vcf_path: str) -> Dict[str, Dict]:
        variants = {}
        if not os.path.exists(vcf_path):
            return variants

        try:
            vcf = vcflib.VCF(vcf_path)
            for v in vcf:
                if not v.alt:
                    continue

                gt_str = v.samples[0].get("GT", "./.")
                if not gt_str or gt_str in ["./.", "."]:
                    continue

                ps = v.samples[0].get("PS", None)
                phased = "|" in gt_str

                if "/" in gt_str or "|" in gt_str:
                    gt_parts = gt_str.replace("|", "/").split("/")
                else:
                    gt_parts = [gt_str]

                try:
                    gt_indices = [int(x) for x in gt_parts if x != "."]
                except (ValueError, TypeError):
                    continue

                if not gt_indices or all(x == 0 for x in gt_indices):
                    continue

                alts = v.alt if isinstance(v.alt, list) else [v.alt]

                for i, alt in enumerate(alts):
                    alt_index = i + 1
                    count = gt_indices.count(alt_index)

                    if count > 0:
                        chrom = v.chrom
                        if not chrom.startswith("chr"):
                            chrom = f"chr{chrom}"

                        key = f"{chrom}:{v.pos + 1}:{v.ref}:{alt}"

                        count_str = "1/1" if count >= 2 else "0/1"

                        haplotype = None
                        if count >= 2:
                            haplotype = "both"
                        elif phased and len(gt_indices) == 2 and count == 1:
                            for idx, allele_idx in enumerate(gt_indices):
                                if allele_idx == alt_index:
                                    haplotype = idx
                                    break

                        variants[key] = {
                            "gt": count_str,
                            "count": count,
                            "ps": ps,
                            "phased": phased,
                            "haplotype": haplotype,
                            "chrom": v.chrom,
                            "pos": v.pos + 1,
                            "ref": v.ref,
                            "alt": alt,
                        }
            vcf.close()
        except Exception:
            pass

        return variants

    def _make_variant_key(self, var_def: Dict) -> str:
        return f"{var_def['chrom']}:{var_def['pos']}:{var_def['ref']}:{var_def['alt']}"

    def _check_cis_phasing(
        self, variants: Dict, defining_key: str, backbone_keys: List[str]
    ) -> bool:
        def_var = variants.get(defining_key)
        if not def_var:
            return False

        def_ps = def_var.get("ps")
        def_phased = def_var.get("phased", False)
        def_hap = def_var.get("haplotype")

        if not def_phased or def_ps is None or def_hap is None or def_hap == "both":
            return True

        for bb_key in backbone_keys:
            bb_var = variants.get(bb_key)
            if not bb_var:
                continue

            bb_ps = bb_var.get("ps")
            bb_phased = bb_var.get("phased", False)
            bb_hap = bb_var.get("haplotype")

            if not bb_phased or bb_ps is None or bb_hap is None or bb_hap == "both":
                continue

            if def_ps == bb_ps and def_hap != bb_hap:
                return False

        return True

    def _get_structural_alleles(self, conversions: List[Dict]) -> Dict[str, int]:
        structural = {}
        for conv in conversions:
            star_allele = conv.get("star_allele", "")
            converted_alleles = conv.get("converted_alleles", 0)
            count = abs(converted_alleles) if converted_alleles != 0 else 1
            if "*13" in star_allele:
                structural["*13"] = structural.get("*13", 0) + count
            elif "*68" in star_allele:
                structural["*68"] = structural.get("*68", 0) + count
            elif "*36" in star_allele:
                structural["*36"] = structural.get("*36", 0) + count
        return structural


class CYP2D6(Gene):
    """CYP2D6 gene with interpretation for structural variants (*13, *36, *68, etc.)."""

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()

        interpretation_lines = []

        cyp2d6_cn = result_data["Copy numbers"].get("CYP2D6", 2)
        cyp2d7_cn = result_data["Copy numbers"].get("CYP2D7", 2)

        interpretation_lines.append(f"CYP2D6: CN={cyp2d6_cn}")
        interpretation_lines.append(f"CYP2D7: CN={cyp2d7_cn}")

        if result_data["gene_conversions"]:
            interpretation_lines.append(
                f"Structural variants detected: {len(result_data['gene_conversions'])}"
            )

            for i, conv in enumerate(result_data["gene_conversions"], 1):
                region = conv["region"]
                converted_alleles = conv["converted_alleles"]
                conversion_sites = conv.get("conversion_sites", [])
                is_fusion = conv.get("is_fusion_candidate", False)
                fusion_type = conv.get("fusion_type", None)
                interp = conv.get("interpretation")

                # Get star allele interpretation
                star_allele = self._interpret_star_allele(
                    converted_alleles, fusion_type, is_fusion, interp
                )

                # Update the interpretation field with the star allele determination
                # This ensures downstream tools see the refined interpretation, not just "Novel fusion"
                conv["star_allele"] = star_allele

                interpretation_lines.append(f"  Event {i}: {star_allele}")
                interpretation_lines.append(f"    Region: {region}")
                interpretation_lines.append(
                    f"    Allele change: {converted_alleles:+d}"
                )

                if conversion_sites:
                    interpretation_lines.append(
                        f"    Signature sites: {', '.join(conversion_sites)}"
                    )

                # Add database match info if available
                interp = conv.get("interpretation")
                if isinstance(interp, list):
                    for match in interp:
                        event_name = match["event_name"]
                        match_score = match["match_score"]
                        interpretation_lines.append(
                            f"    Database match: {event_name} ({match_score:.0%})"
                        )

                # Clinical note
                clinical_note = self._get_clinical_note(star_allele)
                if clinical_note:
                    interpretation_lines.append(f"    Clinical: {clinical_note}")

        else:
            interpretation_lines.append("No structural variants detected")

        star_caller = StarAlleleCaller(self.config)
        vcf_path = result_data["Variants"]
        for conv in result_data["gene_conversions"]:
            star_allele = conv.get("star_allele", "")
            converted_alleles = conv.get("converted_alleles", 0)
            if converted_alleles < 0 and any(x in star_allele for x in ["*13", "*68", "*36"]):
                cyp2d6_cn += abs(converted_alleles)
        star_result = star_caller.call_alleles(
            vcf_path, cyp2d6_cn, result_data.get("gene_conversions", [])
        )
        result_data["star_alleles"] = star_result

        if star_result["alleles"]:
            allele_str = " / ".join(star_result["alleles"])
            interpretation_lines.append(f"Star Alleles: {allele_str}")
            if "note" in star_result:
                interpretation_lines.append(f"  Note: {star_result['note']}")

        # 3. Phenotype prediction
        phenotype = self._predict_phenotype(
            cyp2d6_cn, result_data.get("gene_conversions", [])
        )
        if phenotype:
            interpretation_lines.append(f"\nPhenotype: {phenotype}")

        result_data["cn_interpretation"] = "\n".join(interpretation_lines)
        return result_data

    def _interpret_star_allele(
        self, converted_alleles: int, fusion_type: str, is_fusion: bool, interpretation
    ) -> str:
        """Interpret structural variant as CYP2D6 star allele.

        Key alleles:
        - *13: CYP2D6-CYP2D7 hybrid (D6 5' + D7 3'), non-functional
        - *68: CYP2D7-CYP2D6 hybrid (D7 5' + D6 3'), non-functional
        - *36: Exon 9 conversion from CYP2D7, non-functional
        """
        if not is_fusion:
            # Internal conversion (not at terminus)
            # Check if database matched *36
            if isinstance(interpretation, list):
                for match in interpretation:
                    if "*36" in match.get("event_name", ""):
                        return "*36 (Exon 9 conversion)"

            # No database match - use generic labels
            if converted_alleles > 0:
                return "Conversion (CYP2D7→CYP2D6)"
            else:
                return "Conversion (CYP2D6→CYP2D7)"

        # Fusion events - determine orientation from sign and position
        if fusion_type == "3_PRIME":
            if converted_alleles > 0:
                # Excess D7 at 3' end = D7-D6 hybrid
                return "*68 (CYP2D7-CYP2D6 hybrid)"
            else:
                # Deficit D6 at 3' end = D6-D7 hybrid
                return "*13 (CYP2D6-CYP2D7 hybrid)"

        elif fusion_type == "5_PRIME":
            if converted_alleles < 0:
                return "*68-like (5' D7-D6 fusion)"
            else:
                return "*13-like (5' D6-D7 fusion)"

        return "Fusion (unclassified)"

    def _get_clinical_note(self, star_allele: str) -> str:
        """Get clinical significance for star allele."""
        if "*13" in star_allele:
            return "Non-functional allele, poor metabolizer"
        elif "*68" in star_allele:
            return "Non-functional allele, poor metabolizer"
        elif "*36" in star_allele:
            return "Non-functional allele (frameshift), poor metabolizer"
        elif "Conversion" in star_allele:
            return "May affect enzyme function"
        return ""

    def _predict_phenotype(self, cyp2d6_cn: int, conversions: list) -> str:
        """Predict metabolizer phenotype based on CN and structural variants."""
        if cyp2d6_cn == 0:
            return "Poor metabolizer (no functional copies)"

        has_fusion = any(c.get("is_fusion_candidate", False) for c in conversions)
        has_conversion = len(conversions) > 0

        if cyp2d6_cn == 1:
            if has_fusion or has_conversion:
                return "Poor metabolizer (single copy with non-functional allele)"
            return "Intermediate metabolizer (single functional copy)"

        elif cyp2d6_cn == 2:
            if has_fusion or has_conversion:
                return "Variable (depends on zygosity of structural variant)"
            return "Normal metabolizer (likely *1/*1 or equivalent)"

        elif cyp2d6_cn >= 3:
            if has_fusion or has_conversion:
                return "Variable (CN≥3 with structural variant)"
            return "Potential ultra-rapid metabolizer (CN≥3)"

        return "Unknown"
