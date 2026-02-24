"""RCCX module interpretation for C4, CYP21, and TNX genes.

The RCCX (Repeat C4-CYP21-TNX) module is a complex tandem repeat on chromosome 6p21.3
containing Complement Component 4 (C4), Steroid 21-Hydroxylase (CYP21), and Tenascin X (TNX) genes.

This class interprets already-computed copy number states to provide clinical information about:
- Total module count (based on CYP21 copies)
- C4 isotypes (C4A vs C4B) and sizes (Long vs Short based on HERV-K)
- CYP21 functional status (CYP21A2 functional vs CYP21A1P pseudogene) and chimeras
- TNX status (TNXB functional vs TNXA pseudogene) and CAH-X syndrome detection

Note: This class does NOT re-compute copy numbers or re-identify genes.
The upstream CN solver has already distinguished C4A/C4B and CYP21A2/CYP21A1P
using differentiating variants. This class only aggregates and formats results.
"""

from typing import Dict, Any, List, Optional, Tuple
from genecaller.gene import Gene
from genecaller.util import IntervalList


class RCCX(Gene):
    """RCCX module interpretation for C4, CYP21, and TNX genes."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # CAH-X buffer: 2500bp covers ~half the NAHR homology track (~4kb Ten-X-Box)
        self.fusion_boundary_buffer = self.config.get("fusion_boundary_buffer", 2500)
        # CYP21 edge buffer: tight threshold for small genes (~2.9kb CYP21A2)
        self.edge_anchor_buffer = self.config.get("edge_anchor_buffer", 200)

    def cn_constraint_penalty(self) -> float:
        """Penalize biologically implausible CN states in RCCX modules.

        Each RCCX module is a structural unit: C4A pairs with CYP21A1P,
        C4B pairs with CYP21A2. Mismatched copy numbers between paired
        regions indicate a state that is unlikely without strong evidence.
        """
        penalty_weight = self.config.get("constraint_penalty_weight", 1000.0)
        penalty = 0.0

        # Build name -> cn_diff lookup
        cn_by_name = {v["name"]: v["cn_diff"] for v in self.all_vars["cns"].values()}

        c4a = cn_by_name.get("C4A", 0)
        cyp21a1p = cn_by_name.get("CYP21A1P", 0)
        tnxa = cn_by_name.get("TNXA", 0)
        c4b = cn_by_name.get("C4B", 0)
        cyp21a2 = cn_by_name.get("CYP21A2", 0)
        tnxb = cn_by_name.get("TNXB", 0)

        # Each RCCX module is a unit: so total C4 = total CYP21A
        penalty -= 2 * penalty_weight * abs(c4a + c4b - cyp21a1p - cyp21a2)
        penalty -= 2 * penalty_weight * abs(c4a + c4b - tnxa - tnxb)
        return penalty

    def prepare_output(self) -> Dict[str, Any]:
        """Generate comprehensive RCCX module interpretation."""
        result_data = super().prepare_output()
        copy_numbers = result_data["Copy numbers"]

        # Consolidate dual HERV-K variables into a single user-facing value.
        # Internally, two longdel variables (C4A_HERV-K, C4B_HERV-K) are used
        # to lift the CN constraint, but their individual values are degenerate.
        # Only the total Short count is meaningful.
        hervk_a = copy_numbers.pop("C4A_HERV-K", 0)
        hervk_b = copy_numbers.pop("C4B_HERV-K", 0)
        total_hervk = hervk_a + hervk_b
        if total_hervk != 0:
            copy_numbers["C4_HERV-K"] = total_hervk

        module_interpretation = self._interpret_module_structure(copy_numbers)
        c4_interpretation = self._interpret_c4(copy_numbers)
        cyp21_interpretation = self._interpret_cyp21(
            copy_numbers, result_data["gene_conversions"]
        )
        tnx_interpretation = self._interpret_tnx(
            copy_numbers, result_data["gene_conversions"]
        )

        result_data["module_count"] = module_interpretation["total_modules"]
        result_data["c4_analysis"] = c4_interpretation
        result_data["cyp21_analysis"] = cyp21_interpretation
        result_data["tnx_analysis"] = tnx_interpretation
        result_data["cn_interpretation"] = self._build_interpretation_text(
            module_interpretation,
            c4_interpretation,
            cyp21_interpretation,
            tnx_interpretation,
        )

        return result_data

    def _interpret_module_structure(self, copy_numbers: Dict) -> Dict[str, Any]:
        """Calculate total RCCX module count and describe structural haplotype.

        Module count is a structural property independent of gene fusions.
        Common haplotypes: bimodular (2 modules per chromosome) is the reference;
        monomodular (1 module) is very common and benign.
        """
        cyp21a2_cn = copy_numbers.get("CYP21A2", 0)
        cyp21a1p_cn = copy_numbers.get("CYP21A1P", 0)
        total_modules = cyp21a2_cn + cyp21a1p_cn

        if total_modules >= 5:
            haplotype = f"Expanded ({total_modules} modules)"
        elif total_modules == 4:
            haplotype = "Bimodular/Bimodular (Standard)"
        elif total_modules == 3:
            haplotype = "Bimodular/Monomodular"
        elif total_modules == 2:
            haplotype = "Monomodular/Monomodular"
        elif total_modules == 1:
            haplotype = "Monomodular/Null"
        else:
            haplotype = "Null/Null"

        return {"total_modules": total_modules, "haplotype": haplotype}

    def _interpret_c4(self, copy_numbers: Dict) -> Dict[str, Any]:
        """Interpret C4 isotypes, Long/Short status, and isotype switching.

        Three independent dimensions:
        1. Isotype (C4A/C4B) - from differentiating variants
        2. Size (Long/Short) - from HERV-K longdel cn_diff
        3. Isotype switching - inferred from module pairing mismatch

        Isotype switching detection: In standard RCCX modules, C4A pairs with
        CYP21A1P (module 1) and C4B pairs with CYP21A2 (module 2). CYP21
        identity is structurally stable (differences distributed across gene),
        so CYP21 CN serves as the structural ground truth for module counts.
        A mismatch between C4 isotype and CYP21 pairing indicates C4 isotype
        switching (gene conversion at exon 26).
        """
        c4a_cn = copy_numbers.get("C4A", 0)
        c4b_cn = copy_numbers.get("C4B", 0)
        total_c4 = c4a_cn + c4b_cn

        # HERV-K cn_diff: negative = that many copies have HERV-K deleted (Short)
        # Already consolidated from dual internal variables in prepare_output().
        short_count = abs(copy_numbers.get("C4_HERV-K", 0))
        long_count = total_c4 - short_count

        detail_parts = []
        if c4a_cn > 0:
            detail_parts.append(f"C4A (x{c4a_cn})")
        if c4b_cn > 0:
            detail_parts.append(f"C4B (x{c4b_cn})")
        if long_count > 0:
            detail_parts.append(f"Long (x{long_count})")
        if short_count > 0:
            detail_parts.append(f"Short (x{short_count})")

        # Isotype switching: C4A should pair with CYP21A1P, C4B with CYP21A2
        cyp21a1p_cn = copy_numbers.get("CYP21A1P", 0)
        cyp21a2_cn = copy_numbers.get("CYP21A2", 0)
        total_cyp21 = cyp21a1p_cn + cyp21a2_cn

        isotype_mismatch = 0
        isotype_switching = "None"

        if total_c4 == total_cyp21 and total_c4 > 0:
            # Module pairing intact — mismatch indicates isotype switching
            isotype_mismatch = c4a_cn - cyp21a1p_cn
            if isotype_mismatch > 0:
                if self._has_deletion_evidence(copy_numbers):
                    isotype_switching = (
                        f"{isotype_mismatch} C4B->C4A (may be structural - "
                        f"CYP21 chimera detected)"
                    )
                else:
                    isotype_switching = f"{isotype_mismatch} C4B->C4A isotype switching"
            elif isotype_mismatch < 0:
                n = abs(isotype_mismatch)
                if self._has_deletion_evidence(copy_numbers):
                    isotype_switching = (
                        f"{n} C4A->C4B (may be structural - CYP21 chimera detected)"
                    )
                else:
                    isotype_switching = f"{n} C4A->C4B isotype switching"

        return {
            "total_c4": total_c4,
            "c4a_count": c4a_cn,
            "c4b_count": c4b_cn,
            "c4_long_count": long_count,
            "c4_short_count": short_count,
            "isotype_mismatch": isotype_mismatch,
            "isotype_switching": isotype_switching,
            "detail": detail_parts,
        }

    def _interpret_cyp21(
        self, copy_numbers: Dict, conversions: List[Dict]
    ) -> Dict[str, Any]:
        """Interpret CYP21 functional/pseudo status with 3-tier classification.

        Tier 1 - Microconversion (Internal): patch entirely within gene body.
        Tier 2 - Gene Conversion (Boundary): touches gene edge but CN is intact.
        Tier 3 - Unequal Crossover (30kb Deletion): touches edge AND CN shows
                 missing C4/TNXA, confirming physical chromosome breakage.
        """
        functional_cn = copy_numbers.get("CYP21A2", 0)
        pseudo_cn = copy_numbers.get("CYP21A1P", 0)
        total_cyp21 = functional_cn + pseudo_cn

        boundary_events = []
        microconversion_events = []

        for conv in conversions:
            region = conv.get("region", "")
            converted_alleles = conv.get("converted_alleles", 0)
            if converted_alleles <= 0:
                continue
            if self._overlaps_gene(region, "CYP21A2") or self._overlaps_gene(
                region, "CYP21A1P"
            ):
                event = {
                    "region": region,
                    "converted_alleles": converted_alleles,
                }
                if self._is_cyp21_boundary_anchored(
                    region, "CYP21A2"
                ) or self._is_cyp21_boundary_anchored(region, "CYP21A1P"):
                    boundary_events.append(event)
                else:
                    event["type"] = "CYP21 Microconversion (Internal)"
                    microconversion_events.append(event)

        # Classify boundary events: Class 2 vs Class 3 based on CN evidence
        deletion_evidence = self._has_deletion_evidence(copy_numbers)
        chimera_class = "None"

        if boundary_events and deletion_evidence:
            chimera_class = "Class 3: Unequal Crossover (30kb Deletion)"
            for e in boundary_events:
                e["type"] = "CYP21 Chimera (30kb Deletion)"
        elif boundary_events:
            chimera_class = "Class 2: Gene Conversion (Boundary)"
            for e in boundary_events:
                e["type"] = "CYP21 Gene Conversion (Boundary)"

        return {
            "total_cyp21": total_cyp21,
            "functional_count": functional_cn,
            "pseudo_count": pseudo_cn,
            "chimera_detected": chimera_class != "None",
            "chimera_class": chimera_class,
            "conversion_events": boundary_events + microconversion_events,
            "detail": f"{functional_cn} Functional, {pseudo_cn} Pseudo",
        }

    def _has_deletion_evidence(self, copy_numbers: Dict) -> bool:
        """CN evidence for NAHR chimera: the physical array lost a module.

        A boundary gene conversion (copy-paste) leaves the chromosome intact
        — total modules remain at the expected diploid baseline (4 for
        bimodular/bimodular). A true NAHR chimera (cut-paste deletion)
        physically removes ~30kb, shrinking total modules below baseline.

        Total module count (CYP21A2 + CYP21A1P) is the structural anchor.
        Individual gene subtypes (C4B, TNXA) are too noisy for this purpose
        — they are subject to isotype switching and assembly artifacts.
        The HMM boundary location determines *which* chimera (CYP21 vs TNX);
        total module count determines *whether* it is a true deletion.
        """
        total_modules = copy_numbers.get("CYP21A2", 0) + copy_numbers.get("CYP21A1P", 0)
        expected_baseline = 2 * self.baseline_cn  # 4 for diploid
        return total_modules < expected_baseline

    def _interpret_tnx(
        self, copy_numbers: Dict, conversions: List[Dict]
    ) -> Dict[str, Any]:
        """Interpret TNX status with 3-tier classification (parallel to CYP21).

        Tier 1 - Microconversion (Internal): patch inside TNXB, not touching boundary.
        Tier 2 - Boundary Conversion: touches 5' TNXB boundary but total module
                 count == baseline, chromosome structurally intact.
        Tier 3 - CAH-X Chimeric Event: touches 5' boundary AND total module count
                 < baseline, confirming NAHR deletion (unequal crossover).
        """
        tnxb_cn = copy_numbers.get("TNXB", 0)
        tnxa_cn = copy_numbers.get("TNXA", 0)

        boundary_events = []
        microconversion_events = []

        for conv in conversions:
            region = conv.get("region", "")
            if self._overlaps_gene(region, "TNXB"):
                converted_alleles = conv.get("converted_alleles", 0)
                if converted_alleles > 0:
                    event = {
                        "region": region,
                        "converted_alleles": converted_alleles,
                    }
                    if self._is_tnxb_boundary_anchored(region):
                        boundary_events.append(event)
                    else:
                        event["type"] = "TNXB Microconversion (Internal)"
                        microconversion_events.append(event)

        # Classify boundary events: boundary conversion vs CAH-X chimera
        deletion_evidence = self._has_deletion_evidence(copy_numbers)
        cah_x_class = "None"

        if boundary_events and deletion_evidence:
            cah_x_class = "Class 3: CAH-X Chimeric Event (Unequal Crossover)"
            for e in boundary_events:
                e["type"] = "CAH-X Chimera (30kb Deletion)"
        elif boundary_events:
            cah_x_class = "Class 2: TNXB Boundary Conversion"
            for e in boundary_events:
                e["type"] = "TNXB Boundary Conversion"

        status = "Normal"
        clinical_note = ""

        if cah_x_class.startswith("Class 3"):
            status = "CAH-X Chimeric Event"
            clinical_note = (
                "Unequal crossover with C4/CYP21A2 CN loss confirmed. "
                "High risk for EDS; evaluate CYP21A2 for concurrent CAH."
            )
        elif cah_x_class.startswith("Class 2"):
            status = "TNXB Boundary Conversion"
            clinical_note = (
                "TNXA sequence overwrote TNXB 5' boundary without structural "
                "deletion. Chromosome intact."
            )
        elif tnxb_cn < 2:
            status = "TNXB Haploinsufficiency (Carrier)"
            clinical_note = "Reduced copy number without gene conversion signal."

        return {
            "functional_tnxb": tnxb_cn,
            "tnxa_count": tnxa_cn,
            "cah_x_class": cah_x_class,
            "cah_x_detected": cah_x_class.startswith("Class 3"),
            "haploinsufficiency": tnxb_cn < 2 and cah_x_class == "None",
            "status": status,
            "clinical_note": clinical_note,
            "events": boundary_events + microconversion_events,
        }

    def _build_interpretation_text(
        self,
        module_interpretation: Dict,
        c4_interpretation: Dict,
        cyp21_interpretation: Dict,
        tnx_interpretation: Dict,
    ) -> str:
        """Build human-readable interpretation text."""
        lines = []

        lines.append("RCCX Module Analysis:")
        lines.append(f"  Total Modules: {module_interpretation['total_modules']}")
        lines.append(f"  Haplotype: {module_interpretation['haplotype']}")
        lines.append("")

        lines.append("C4 (Complement Component 4):")
        lines.append(f"  Total C4 Copies: {c4_interpretation['total_c4']}")
        lines.append(
            f"  Isotypes: {c4_interpretation['c4a_count']} C4A, {c4_interpretation['c4b_count']} C4B"
        )
        lines.append(
            f"  Sizes: {c4_interpretation['c4_long_count']} Long (with HERV-K), "
            f"{c4_interpretation['c4_short_count']} Short (HERV-K deleted)"
        )
        isotype_switching = c4_interpretation.get("isotype_switching", "None")
        if isotype_switching != "None":
            lines.append(f"  Isotype Switching: {isotype_switching}")
        lines.append("")

        lines.append("CYP21 (Steroid 21-Hydroxylase):")
        lines.append(f"  Total CYP21 Copies: {cyp21_interpretation['total_cyp21']}")
        lines.append(
            f"  Functional (CYP21A2): {cyp21_interpretation['functional_count']}"
        )
        lines.append(f"  Pseudogene (CYP21A1P): {cyp21_interpretation['pseudo_count']}")

        chimera_class = cyp21_interpretation.get("chimera_class", "None")
        lines.append(f"  Chimera Classification: {chimera_class}")

        boundary = [
            e
            for e in cyp21_interpretation.get("conversion_events", [])
            if "Microconversion" not in e.get("type", "")
        ]
        microconversions = [
            e
            for e in cyp21_interpretation.get("conversion_events", [])
            if "Microconversion" in e.get("type", "")
        ]

        if boundary:
            lines.append("  Boundary Conversions:")
            for event in boundary:
                lines.append(
                    f"    {event['region']} (alleles: {event['converted_alleles']:+d})"
                )
        if microconversions:
            lines.append("  Microconversions:")
            for event in microconversions:
                lines.append(
                    f"    {event['region']} (alleles: {event['converted_alleles']:+d})"
                )
        if not boundary and not microconversions:
            lines.append("  Conversion Events: None")
        lines.append("")

        lines.append("TNX (Tenascin X):")
        lines.append(f"  Functional TNXB: {tnx_interpretation['functional_tnxb']}")
        lines.append(
            f"  CAH-X Classification: {tnx_interpretation.get('cah_x_class', 'None')}"
        )

        if tnx_interpretation.get("clinical_note"):
            lines.append(f"  Clinical Note: {tnx_interpretation['clinical_note']}")

        tnx_boundary = [
            e
            for e in tnx_interpretation.get("events", [])
            if "Microconversion" not in e.get("type", "")
        ]
        tnx_micro = [
            e
            for e in tnx_interpretation.get("events", [])
            if "Microconversion" in e.get("type", "")
        ]

        if tnx_boundary:
            lines.append("  Boundary Conversions:")
            for event in tnx_boundary:
                lines.append(
                    f"    {event['region']} (alleles: {event['converted_alleles']:+d})"
                )
        if tnx_micro:
            lines.append("  Microconversions:")
            for event in tnx_micro:
                lines.append(
                    f"    {event['region']} (alleles: {event['converted_alleles']:+d})"
                )
        if not tnx_boundary and not tnx_micro:
            lines.append("  Conversion Events: None")

        return "\n".join(lines)

    def _get_gene_coordinates(self, gene_name: str) -> Optional[Tuple[str, int, int]]:
        """Get coordinates from cn_regions configuration."""
        if isinstance(self.cn_regions, dict):
            coords = self.cn_regions.get(gene_name)
            if coords:
                chrom, start_end = coords.split(":")
                start, end = map(int, start_end.split("-"))
                return (chrom, start, end)
        else:
            for region_name, coords in self.cn_regions:
                if region_name == gene_name:
                    chrom, start_end = coords.split(":")
                    start, end = map(int, start_end.split("-"))
                    return (chrom, start, end)
        return None

    def _overlaps_gene(self, region_str: str, gene_name: str) -> bool:
        """Check if a region overlaps with a gene's coordinates."""
        gene_coords = self._get_gene_coordinates(gene_name)
        if not gene_coords:
            return False
        region_itv = IntervalList(region=region_str)
        gene_itv = IntervalList(
            region=f"{gene_coords[0]}:{gene_coords[1]}-{gene_coords[2]}"
        )
        return region_itv.intersect(gene_itv).size() > 0

    def _is_cyp21_boundary_anchored(self, conv_region: str, gene_name: str) -> bool:
        """Check if conversion anchors to the 5' (promoter) edge of a CYP21 gene.

        Only the 5' edge matters for chimera detection. The 3' edge faces
        TNXB/TNXA and conversions there are handled by TNX interpretation.
        The 5' edge is determined as the edge farthest from TNXB.
        """
        gene_coords = self._get_gene_coordinates(gene_name)
        tnxb_coords = self._get_gene_coordinates("TNXB")
        if not gene_coords or not tnxb_coords:
            return True

        try:
            chrom, coords = conv_region.split(":")
            conv_start, conv_end = map(int, coords.split("-"))
        except (ValueError, AttributeError):
            return True

        # 5' (promoter) edge = farthest from TNXB
        gene_start, gene_end = gene_coords[1], gene_coords[2]
        tnxb_mid = (tnxb_coords[1] + tnxb_coords[2]) // 2

        if abs(gene_start - tnxb_mid) > abs(gene_end - tnxb_mid):
            promoter_edge = gene_start
        else:
            promoter_edge = gene_end

        buf = self.edge_anchor_buffer
        return (
            abs(conv_start - promoter_edge) < buf or abs(conv_end - promoter_edge) < buf
        )

    def _is_tnxb_boundary_anchored(self, conv_region: str) -> bool:
        """Check if conversion anchors to the 5' TNXB boundary (facing CYP21A2).

        Boundary-anchored conversions are candidates for Class 2 (boundary
        conversion) or Class 3 (CAH-X chimera) — distinguished by CN evidence.
        Internal patches that don't touch this edge are Class 1 microconversions.
        Uses configurable buffer (default 2500bp) for NAHR homology block uncertainty.
        """
        tnxb_coords = self._get_gene_coordinates("TNXB")
        cyp21a2_coords = self._get_gene_coordinates("CYP21A2")

        if not tnxb_coords or not cyp21a2_coords:
            return True

        # Find the TNXB edge closest to CYP21A2
        tnxb_start, tnxb_end = tnxb_coords[1], tnxb_coords[2]
        cyp21_mid = (cyp21a2_coords[1] + cyp21a2_coords[2]) // 2

        if abs(tnxb_start - cyp21_mid) < abs(tnxb_end - cyp21_mid):
            fusion_boundary = tnxb_start
        else:
            fusion_boundary = tnxb_end

        try:
            chrom, coords = conv_region.split(":")
            conv_start, conv_end = map(int, coords.split("-"))
        except (ValueError, AttributeError):
            return True

        return (
            abs(conv_start - fusion_boundary) < self.fusion_boundary_buffer
            or abs(conv_end - fusion_boundary) < self.fusion_boundary_buffer
        )
