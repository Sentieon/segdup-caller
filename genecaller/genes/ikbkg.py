"""IKBKG gene-specific class for X-linked NEMO gene analysis.

IKBKG (inhibitor of kappa light polypeptide gene enhancer in B-cells, kinase gamma)
is an X-linked gene where loss-of-function is lethal in males. Key features:

- Sex-specific ploidy: baseline CN=1 (male) or CN=2 (female)
- Bidirectional 11.7kb CNV: the segmental duplication can be deleted OR duplicated
- Mosaic detection: living males with pathogenic events must be mosaic
- 6 cn_regions modeling: each gene split into flank/11.7kb/flank sub-regions
  with soft constraints tying flanks together (same allele)

The cn_diff_priors in config are sex-agnostic (keys = deviation from baseline).
This class converts them to absolute CN priors at runtime based on sample sex.
"""

from typing import Dict, Any
from genecaller.gene import Gene
from genecaller.gene_model import CategoricalPriors


# Internal sub-region names used for constraint and output collapsing
_FLANK_PAIRS = [
    ("IKBKG_5prime", "IKBKG_3prime"),
    ("IKBKGP1_5prime", "IKBKGP1_3prime"),
]

# Mapping: internal sub-region -> reported gene name
_FLANK_TO_GENE = {
    "IKBKG_5prime": "IKBKG",
    "IKBKG_3prime": "IKBKG",
    "IKBKGP1_5prime": "IKBKGP1",
    "IKBKGP1_3prime": "IKBKGP1",
}

# The 11.7kb segment regions — consolidated into a single reported value
# because no differentiating variants exist within the segment
_SEGMENT_REGIONS = ("IKBKGdel", "IKBKGP1del")


class IKBKG(Gene):
    """IKBKG gene-specific class with bidirectional CNV and mosaic support."""

    def _initialize_priors(self) -> None:
        """Convert cn_diff_priors to absolute CN priors based on baseline_cn."""
        if "cn_diff_priors" not in self.config:
            super()._initialize_priors()
            return

        cn_diff_priors = self.config["cn_diff_priors"]
        absolute_priors = {}

        for region_name, diff_probs in cn_diff_priors.items():
            absolute_priors[region_name] = {}
            for cn_diff, prob in diff_probs.items():
                cn_diff = int(cn_diff)
                abs_cn = self.baseline_cn + cn_diff
                if abs_cn >= 0:  # Skip negative absolute CN
                    absolute_priors[region_name][abs_cn] = prob

        # Inject as cn_priors so CategoricalPriors can consume it
        self.config["cn_priors"] = absolute_priors
        self._priors = CategoricalPriors(
            self.config, self.all_vars["cns"], self.baseline_cn
        )

    def cn_constraint_penalty(self) -> float:
        """Penalize states where flanking regions of the same gene differ.

        Each gene's 5' and 3' flanks represent the same contiguous allele
        outside the 11.7kb segment. Their CN must match — a mismatch is
        biologically implausible (would require two independent breakpoints).
        """
        penalty_weight = self.config.get("constraint_penalty_weight", 1000.0)
        penalty = 0.0

        cn_by_name = {v["name"]: v["cn_diff"] for v in self.all_vars["cns"].values()}

        for prime5, prime3 in _FLANK_PAIRS:
            diff = abs(cn_by_name.get(prime5, 0) - cn_by_name.get(prime3, 0))
            penalty -= penalty_weight * diff

        return penalty

    def _collapse_copy_numbers(self, copy_numbers: Dict[str, int]) -> Dict[str, int]:
        """Collapse 6 internal cn_regions into 3 reported values.

        Reported copy numbers:
        - IKBKG: CN of flanking regions (= baseline + cn_diff of IKBKG_5prime)
        - IKBKGP1: CN of flanking regions (= baseline + cn_diff of IKBKGP1_5prime)
        - IKBKGdel_region: sum of IKBKGdel + IKBKGP1del (total 11.7kb
          segment CN across both loci); the two segment regions contain no
          differentiating variants so gene-vs-pseudogene attribution is not
          possible.
        """
        collapsed = {}
        seen_genes = set()
        segment_total = 0

        for region_name, cn in list(copy_numbers.items()):
            if region_name in _SEGMENT_REGIONS:
                segment_total += cn
            elif region_name in _FLANK_TO_GENE:
                gene_name = _FLANK_TO_GENE[region_name]
                if gene_name not in seen_genes:
                    collapsed[gene_name] = cn
                    seen_genes.add(gene_name)
            else:
                collapsed[region_name] = cn

        collapsed["IKBKGdel_region (IKBKG+IKBKGP1 total)"] = segment_total
        return collapsed

    def prepare_output(self) -> Dict[str, Any]:
        """Collapse CN regions and add clinical interpretation.

        Interpretation covers:
        - Overall structural status (whole-gene vs segment-only events)
        - 11.7kb segment deletion/duplication with NAHR evidence
        - Gene conversion events (IKBKGP1→IKBKG pathogenic vs reverse benign)
        - Mosaic events with fraction estimates
        - Sex-specific clinical significance and summary
        """
        result_data = super().prepare_output()
        result_data["Copy numbers"] = self._collapse_copy_numbers(
            result_data["Copy numbers"]
        )

        lines = []
        copy_numbers = result_data["Copy numbers"]

        # Extract CN values — flanks are distinguishable, segment is consolidated
        ikbkg_cn = copy_numbers.get("IKBKG", self.baseline_cn)
        ikbkgp1_cn = copy_numbers.get("IKBKGP1", self.baseline_cn)
        segment_total = copy_numbers.get(
            "IKBKGdel_region (IKBKG+IKBKGP1 total)", 2 * self.baseline_cn
        )
        segment_expected = 2 * self.baseline_cn

        # Determine sex (requires_sex=true guarantees sample_sex is set)
        normalized = self._normalize_sex(self.sample_sex)
        is_male = normalized == "male"
        is_female = normalized == "female"

        # --- 1. Header ---
        lines.append(f"IKBKG (NEMO) Analysis ({normalized}):")

        # --- 2. CN summary ---
        lines.append(f"  IKBKG flanks: CN={ikbkg_cn}")
        lines.append(f"  IKBKGP1 flanks: CN={ikbkgp1_cn}")
        lines.append(
            f"  IKBKGdel_region — 11.7kb shared segment (exons 4-10): "
            f"total CN={segment_total} (expected {segment_expected})"
        )

        # --- 3. Structural event analysis ---
        has_pathogenic_event = False
        has_any_event = False
        # Track event types for female bi-allelic check
        has_pathogenic_deletion = False
        has_pathogenic_conversion = False

        segment_diff = segment_total - segment_expected

        # Check for whole-gene deletion (flanks reduced AND segment reduced)
        # Flank CN *is* distinguishable, so whole-gene attribution is valid
        ikbkg_whole_del = (
            ikbkg_cn < self.baseline_cn and segment_diff < 0
        )

        # --- 3a. Whole-gene deletion (flank-based, attributable) ---
        if ikbkg_whole_del:
            has_any_event = True
            has_pathogenic_event = True
            has_pathogenic_deletion = True
            lines.append("")
            lines.append(
                f"IKBKG whole-gene deletion detected "
                f"(flank CN={ikbkg_cn}, segment total CN={segment_total})"
            )
            lines.append(
                "  Flank CN reduction confirms IKBKG involvement."
            )
            if is_male:
                lines.append(
                    "  Clinical: Complete IKBKG loss in male — "
                    "non-viable unless mosaic."
                )
            elif is_female:
                if ikbkg_cn == 0:
                    lines.append(
                        "  Clinical: Homozygous whole-gene IKBKG deletion — "
                        "expected to be lethal or severely affected."
                    )
                else:
                    lines.append(
                        "  Clinical: Heterozygous whole-gene IKBKG deletion — "
                        "Incontinentia Pigmenti carrier."
                    )

        # --- 3b. Segment-only CNV (gene attribution unresolved) ---
        if segment_diff < 0 and not ikbkg_whole_del:
            has_any_event = True
            lines.append("")
            lines.append(
                f"[UNRESOLVED] 11.7kb segment deletion detected "
                f"(total CN={segment_total}, expected={segment_expected})"
            )
            lines.append(
                "  The 11.7kb segment is identical between IKBKG and "
                "IKBKGP1."
            )
            lines.append(
                "  Short reads cannot resolve which gene carries this "
                "event."
            )
            lines.append(
                "  Mechanism: NAHR between segmental duplications."
            )
            lines.append(
                "  If IKBKG: pathogenic (exon 4-10 loss). "
                "If IKBKGP1: benign."
            )
            lines.append(
                "  Recommend long-read confirmation for definitive "
                "assignment."
            )
            # Conservatively flag as potentially pathogenic
            has_pathogenic_event = True
            has_pathogenic_deletion = True
        elif segment_diff > 0:
            has_any_event = True
            lines.append("")
            lines.append(
                f"[UNRESOLVED] 11.7kb segment duplication detected "
                f"(total CN={segment_total}, expected={segment_expected})"
            )
            lines.append(
                "  The 11.7kb segment is identical between IKBKG and "
                "IKBKGP1."
            )
            lines.append(
                "  Short reads cannot resolve which gene carries this "
                "event."
            )
            lines.append(
                "  Mechanism: NAHR between segmental duplications."
            )
            lines.append(
                "  Duplication of exons 4-10 is usually benign "
                "regardless of gene."
            )

        # --- 4. Gene conversion events ---
        gene_conversions = result_data.get("gene_conversions", [])
        if gene_conversions:
            has_any_event = True
            lines.append("")
            lines.append(
                f"Gene Conversion Events: {len(gene_conversions)}"
            )
            for i, conv in enumerate(gene_conversions, 1):
                alleles = conv["converted_alleles"]
                if alleles > 0:
                    direction = "IKBKGP1\u2192IKBKG (pathogenic direction)"
                    has_pathogenic_event = True
                    has_pathogenic_conversion = True
                    if is_male:
                        clinical = (
                            "Pseudogene sequence in functional IKBKG — "
                            "male must be mosaic to survive."
                        )
                    else:  # is_female
                        clinical = (
                            "Pseudogene sequence in functional IKBKG — "
                            "Incontinentia Pigmenti carrier."
                        )
                else:
                    direction = "IKBKG\u2192IKBKGP1 (benign direction)"
                    clinical = (
                        "Functional sequence replacing pseudogene — "
                        "likely benign."
                    )
                lines.append(
                    f"  Event {i}: {direction}, alleles={alleles}"
                )
                lines.append(f"    Region: {conv['region']}")
                lines.append(f"    Clinical: {clinical}")

                # Include known conversion match if available
                interp = conv.get("interpretation")
                if isinstance(interp, list):
                    for match in interp:
                        event_name = match["event_name"]
                        if "matched_variants" in match:
                            lines.append(
                                f"    Known event: {event_name} "
                                f"(score: {match['match_score']:.1%})"
                            )
                        else:
                            lines.append(
                                f"    Known event: {event_name} "
                                f"(position-based)"
                            )

                # Include conversion sites
                sites = conv.get("conversion_sites", [])
                if sites:
                    lines.append(
                        f"    Variant sites: {', '.join(sites)}"
                    )

        # --- 5. Mosaic events ---
        mosaic = result_data.get("mosaic")
        if mosaic:
            has_any_event = True
            has_pathogenic_event = True
            lines.append("")
            mosaic_type = mosaic["type"]
            mosaic_frac = mosaic["fraction"]
            lines.append(
                f"Mosaic Event Detected: {mosaic['event']} ({mosaic_type}), "
                f"fraction={mosaic_frac:.1%}"
            )
            if is_male and mosaic_type in ("longdel", "conversion"):
                lines.append(
                    f"  Mosaicism explains viability in affected male — "
                    f"{mosaic_frac:.0%} of cells carry the pathogenic event."
                )

        # --- 6. Female bi-allelic warning ---
        if is_female and has_pathogenic_deletion and has_pathogenic_conversion:
            lines.append("")
            lines.append(
                "NOTE: Both deletion and pathogenic gene conversion detected. "
                "If on different alleles (compound heterozygous), this "
                "individual may be affected rather than a carrier."
            )

        # --- 7. Summary ---
        if not has_any_event:
            lines.append("")
            lines.append("No structural events detected.")

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data
