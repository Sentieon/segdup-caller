"""CYP2D6 gene with enhanced interpretation for gene conversions and fusion events."""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.util import get_data_file
import vcflib


class StarAlleleCaller:
    """Call CYP2D6 star alleles using phasing-aware consumption method."""

    def __init__(self, config: Dict[str, Any], ref: Any = None):
        # `ref` is the genome Reference (build-specific). It lets star calling be
        # reference-aware: on GRCh37/hg19 the CYP2D6 reference already encodes the
        # *2 backbone (rs16947/rs1135840), so those markers appear as REFERENCE
        # (no variant call) and must be reconstructed. On GRCh38 no marker is
        # reference-encoded, so this is a no-op and behavior is unchanged.
        self.ref = ref
        self.allele_defs = self._load_allele_definitions(config)
        self.backbone_vars = self._load_backbone_variants(config)

    @staticmethod
    def _canon_chrom(chrom: str) -> str:
        """Canonical contig form for variant keys: strip a leading 'chr'.
        Makes hg38 (chr22), hg19 (chr22) and b37 (22) tables all key-compatible
        with variants extracted from the result VCF."""
        return chrom[3:] if chrom.startswith("chr") else chrom

    @classmethod
    def _canon_av(cls, av: str) -> str:
        """Canonicalize a 'chrom:pos:ref:alt' associated-variant string."""
        parts = av.split(":")
        if parts:
            parts[0] = cls._canon_chrom(parts[0])
        return ":".join(parts)

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
                required_cols = len(header) - 1  # associated_variants is optional
                if len(fields) < required_cols:
                    continue
                allele_data = dict(zip(header, fields))
                allele_data["tier"] = int(allele_data["tier"])
                try:
                    allele_data["pos"] = int(allele_data["pos"])
                except ValueError:
                    pass
                allele_data["backbone"] = bool(int(allele_data["backbone"]))
                assoc = allele_data.get("associated_variants", "")
                allele_data["associated_variants"] = (
                    [self._canon_av(v.strip()) for v in assoc.split(";") if v.strip()]
                    if assoc
                    else []
                )
                definitions.append(allele_data)
        return sorted(definitions, key=lambda x: (x["tier"], not x["backbone"]))

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

        # Reference-aware reconstruction of markers the build reference encodes
        # (e.g. the *2 backbone on GRCh37). No-op when no marker is ref-encoded
        # (always the case on GRCh38), so hg38 calling is unchanged.
        self._apply_reference_encoding(variants, variant_counts, copy_number)

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
                    for av in allele_def.get("associated_variants", []):
                        if variant_counts.get(av, 0) > 0:
                            variant_counts[av] -= 1

                    while (
                        variant_counts[var_key] > 0
                        and len(called_alleles) < max_alleles
                    ):
                        backbone_enough = all(
                            variant_counts.get(k, 0) > 0 for k in backbone_keys
                        )
                        if not backbone_enough:
                            break
                        called_alleles.append(allele_def["allele"])
                        variant_counts[var_key] -= 1
                        for bb_key in backbone_keys:
                            variant_counts[bb_key] -= 1
                        for av in allele_def.get("associated_variants", []):
                            if variant_counts.get(av, 0) > 0:
                                variant_counts[av] -= 1
            else:
                called_alleles.append(allele_def["allele"])
                variant_counts[var_key] -= 1
                for av in allele_def.get("associated_variants", []):
                    if variant_counts.get(av, 0) > 0:
                        variant_counts[av] -= 1

                while variant_counts.get(var_key, 0) > 0 and len(called_alleles) < max_alleles:
                    called_alleles.append(allele_def["allele"])
                    variant_counts[var_key] -= 1
                    for av in allele_def.get("associated_variants", []):
                        if variant_counts.get(av, 0) > 0:
                            variant_counts[av] -= 1

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
                        chrom = self._canon_chrom(v.chrom)

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
        return f"{self._canon_chrom(var_def['chrom'])}:{var_def['pos']}:{var_def['ref']}:{var_def['alt']}"

    def _ref_base(self, chrom: str, pos: int):
        """Single reference base (1-based pos) on the build genome, or None.
        Tries the contig with and without a 'chr' prefix to match either naming."""
        if self.ref is None:
            return None
        for c in (chrom, f"chr{self._canon_chrom(chrom)}", self._canon_chrom(chrom)):
            try:
                b = self.ref.get(c, pos - 1, pos)
                if b:
                    return b.upper()
            except Exception:
                continue
        return None

    def _iter_marker_defs(self):
        """Yield (chrom, pos, ref, alt) for every SNV marker the consumption
        looks up: each allele's defining variant, the backbone variants, and all
        associated variants. Indels are skipped (single-base ref comparison is
        only valid for SNVs; the known reference-encoded markers are all SNVs)."""
        seen = set()

        def emit(chrom, pos, ref, alt):
            if (
                ref and alt and len(ref) == 1 and len(alt) == 1
                and ref in "ACGT" and alt in "ACGT"
            ):
                key = (self._canon_chrom(chrom), int(pos), ref, alt)
                if key not in seen:
                    seen.add(key)
                    yield key

        for d in self.allele_defs:
            if isinstance(d.get("pos"), int) and str(d.get("chrom", "")).startswith(("chr", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "X", "Y")):
                yield from emit(d["chrom"], d["pos"], d.get("ref", ""), d.get("alt", ""))
        for bv in self.backbone_vars:
            yield from emit(bv["chrom"], bv["pos"], bv.get("ref", ""), bv.get("alt", ""))
        for d in self.allele_defs:
            for av in d.get("associated_variants", []):
                p = av.split(":")
                if len(p) == 4:
                    try:
                        yield from emit(p[0], int(p[1]), p[2], p[3])
                    except ValueError:
                        continue

    def _apply_reference_encoding(self, variants, variant_counts, copy_number):
        """For each SNV marker whose build reference base equals the marker's ALT
        allele, the marker is reference-encoded: it is present on every haplotype
        that does NOT carry a back-mutation (ALT->REF) call. Reconstruct its
        presence so the variant-presence consumption logic still works.

        On GRCh38 no marker is reference-encoded, so nothing changes."""
        if self.ref is None:
            return
        cn = copy_number if copy_number and copy_number > 0 else 2
        for chrom, pos, ref, alt in self._iter_marker_defs():
            rb = self._ref_base(chrom, pos)
            if rb is None or rb != alt:
                continue  # normal (ref==ref) or third-allele -> leave as-is
            marker_key = f"{chrom}:{pos}:{ref}:{alt}"
            back_key = f"{chrom}:{pos}:{alt}:{ref}"  # ALT->REF back-mutation call
            back_var = variants.get(back_key)
            back_count = variant_counts.get(back_key, 0)
            marker_count = max(0, cn - back_count)
            if marker_count == 0:
                continue
            variant_counts[marker_key] = max(
                variant_counts.get(marker_key, 0), marker_count
            )
            # Best-effort phasing: marker sits opposite the back-mutation.
            if marker_count >= cn:
                hap = "both"
            elif back_var and back_var.get("haplotype") in (0, 1):
                hap = 1 - back_var["haplotype"]
            else:
                hap = None
            variants.setdefault(
                marker_key,
                {
                    "gt": "1/1" if marker_count >= 2 else "0/1",
                    "count": marker_count,
                    "ps": back_var.get("ps") if back_var else None,
                    "phased": back_var.get("phased", False) if back_var else False,
                    "haplotype": hap,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "reference_encoded": True,
                },
            )

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

        star_caller = StarAlleleCaller(self.config, ref=self.ref)
        vcf_path = result_data["Variants"]
        for conv in result_data["gene_conversions"]:
            star_allele = conv.get("star_allele", "")
            converted_alleles = conv.get("converted_alleles", 0)
            if converted_alleles < 0 and any(
                x in star_allele for x in ["*13", "*68", "*36"]
            ):
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

        phenotype = self._predict_phenotype(star_result.get("alleles", []))
        if phenotype:
            interpretation_lines.append(f"Phenotype: {phenotype}")

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

    def _predict_phenotype(self, star_alleles: List[str]) -> str:
        """Predict metabolizer phenotype based on star allele calls."""
        if not star_alleles:
            return "Unknown"

        allele_to_function = {}
        for allele_def in StarAlleleCaller(self.config).allele_defs:
            allele_to_function[allele_def["allele"]] = allele_def["function"]

        function_counts = {
            "No Function": 0,
            "Reduced": 0,
            "Uncertain Function": 0,
            "Normal": 0,
        }
        for allele in star_alleles:
            allele_clean = allele.split("/")[0].strip()
            if allele_clean == "*5":
                function_counts["No Function"] += 1
            elif allele_clean == "*1":
                function_counts["Normal"] += 1
            else:
                func = allele_to_function.get(allele_clean, "Uncertain Function")
                function_counts[func] += 1

        no_func = function_counts["No Function"]
        reduced = function_counts["Reduced"]
        uncertain = function_counts["Uncertain Function"]
        normal = function_counts["Normal"]
        total = len(star_alleles)

        if no_func == total:
            return "Poor metabolizer (no functional alleles)"
        elif no_func >= 1 and no_func + reduced == total:
            return "Poor metabolizer (all alleles non-functional or reduced)"
        elif reduced == total:
            return "Intermediate metabolizer (all reduced function)"
        elif no_func >= 1 or reduced >= 1:
            return "Intermediate metabolizer (mix of functional and impaired alleles)"
        elif total >= 3 and normal >= 2:
            return "Ultrarapid metabolizer (gene duplication with functional alleles)"
        elif normal == total:
            return "Normal metabolizer"
        elif uncertain > 0:
            return "Uncertain (contains alleles with unknown function)"

        return "Unknown"
