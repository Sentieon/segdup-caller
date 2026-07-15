"""CYP2A6 pharmacogene: star-allele calling, CYP2A7 hybrid/deletion interpretation,
and a nicotine-metabolism (NMR) metabolizer grouping.

CYP2A6 combines the two structural motifs the caller already handles:

1. **CYP2D6-shaped structural family.** Whole-gene events dominate: ``*4`` (whole-gene
   deletion, common in East Asians), ``*1x2`` (duplication), and the CYP2A7::CYP2A6 5'
   hybrids ``*12``/``*34``. So structural calling is copy-number driven (like CYP2D6
   ``*5`` / ``xN``), NOT the single-copy intergenic-probe model CYP2B6 uses for its
   *partial* ``*29``/``*30`` events.

2. **CYP2B6-shaped combination star matcher.** The exon-9 gene-conversion cluster
   (``*7``=I471T, ``*8``=R485L, ``*35``=N438Y, sitting on the shared ``*46`` 3'-conversion
   background) is combinatorial: ``*10`` = I471T+R485L in cis. So SNV matching REQUIRES
   all of an allele's variants (and cis when phased), most-specific-first -- reused from
   :class:`CYP2B6StarAlleleCaller`. Each cluster allele is distilled to its exon-9
   signature; the shared ``*46`` background is intentionally dropped because the alleles
   it distinguishes (``*46``/``*1B``) are increased/normal function, so reading them as
   ``*1`` never changes the metabolizer group.

CYP2A6 has NO CPIC guideline -- PharmVar/CPIC clinical allele function is N/A. Function
labels here are literature / enzyme-activity based and the reported phenotype is a
research nicotine-metabolism grouping (ultrarapid / rapid / normal / intermediate /
slow / poor), NOT a CPIC dosing recommendation.
"""

from typing import Dict, Any, List, Optional
from genecaller.gene import Gene
from genecaller.genes.cyp2b6 import (
    CYP2B6StarAlleleCaller,
    _FUNC_NORMAL,
    _FUNC_INCREASED,
    _FUNC_DECREASED,
    _FUNC_NONE,
    _FUNC_UNCERTAIN,
)

# Approximate max population allele frequency (PharmVar / literature / gnomAD, order-of-
# magnitude only). Used SOLELY to rank equally-consistent unphased diplotypes so the
# most likely is reported as primary -- NOT a clinical frequency estimate. CYP2A6*4/*9
# are common (esp. East Asian / global), the exon-9 combination alleles are rare.
_ALLELE_FREQ_DEFAULT = 0.001
_ALLELE_FREQ = {
    "*1": 0.55, "*9": 0.08, "*4": 0.05, "*2": 0.03, "*17": 0.03, "*35": 0.03,
    "*7": 0.02, "*8": 0.02, "*1x2": 0.006, "*13": 0.005, "*10": 0.005,
    "*12": 0.003, "*5": 0.001, "*20": 0.001, "*34": 0.001,
}

# Per-allele activity contribution for the nicotine-metabolism activity score (AS).
# Research grouping (not CPIC): none=0, decreased~0.5, normal=1.0, increased~1.5.
_FUNC_ACTIVITY = {
    _FUNC_NONE: 0.0,
    _FUNC_DECREASED: 0.5,
    _FUNC_NORMAL: 1.0,
    _FUNC_INCREASED: 1.5,
}


class CYP2A6StarAlleleCaller(CYP2B6StarAlleleCaller):
    """Require-all-cis SNV matcher (inherited from CYP2B6) wired to CYP2D6-style
    whole-gene copy-number structural calling: ``*4`` deletion, ``*1x2`` duplication,
    and CYP2A7::CYP2A6 ``*12``/``*34`` 5' hybrids from the fusion junction."""

    def _get_structural_alleles(self, conversions: List[Dict]) -> Dict[str, int]:
        """Map detected CYP2A7::CYP2A6 **5' hybrids** to ``*12``/``*34`` from the fusion
        junction. A named database match (``*34``/``*12`` in the event name) wins;
        otherwise a generic fusion with allele loss falls back to the common ``*12`` —
        but ONLY when the junction is at the **5' terminus** (``fusion_type == "5_PRIME"``).
        A ``3_PRIME`` event is the common **benign ``*1B``/``*46`` 3'-UTR conversion**, NOT a
        hybrid: calling ``*12`` there was a systematic false positive (cohort: ~29% of
        samples carry the benign 3' conversion), so 3' events add no structural allele and
        the SNV-based core call (``*1``) stands."""
        structural: Dict[str, int] = {}
        for conv in conversions:
            star = str(conv.get("star_allele", ""))
            interp = conv.get("interpretation")
            names = []
            if isinstance(interp, list):
                names = [str(m.get("event_name", "")) for m in interp]
            blob = star + " " + " ".join(names)
            converted = conv.get("converted_alleles", 0)
            count = abs(converted) if converted else 1
            is_5p = conv.get("fusion_type") == "5_PRIME"
            if "*34" in blob:
                structural["*34"] = structural.get("*34", 0) + count
            elif "*12" in blob:
                structural["*12"] = structural.get("*12", 0) + count
            elif conv.get("is_fusion_candidate") and converted < 0 and is_5p:
                structural["*12"] = structural.get("*12", 0) + count
        return structural

    def call_alleles(
        self, vcf_path: str, copy_number: int, conversions: List[Dict],
        intergenic_cn: int = 2,
    ) -> Dict[str, Any]:
        # Homozygous whole-gene deletion -> *4/*4 (no CYP2A6 sequence to genotype).
        if copy_number == 0:
            return {"alleles": ["*4", "*4"], "method": "CNV", "copy_number": 0}

        variants = self._extract_variants_from_vcf(vcf_path)
        variant_counts = {k: v.get("count", 1) for k, v in variants.items()}
        # No-op on GRCh38 (no CYP2A6 marker is reference-encoded); wired for a future
        # GRCh37 port where CYP2A7-shared bases may be reference-encoded.
        self._apply_reference_encoding(variants, variant_counts, copy_number)
        self._rescue_clustered_alleles(variants, variant_counts)

        cn = copy_number if copy_number and copy_number > 0 else 2
        called: List[str] = []

        # Tier 1: CYP2A7::CYP2A6 5' hybrids (*12/*34) from the fusion junction. A
        # hybrid occupies one of the present gene copies.
        hybrids = self._get_structural_alleles(conversions)
        for name in sorted(hybrids):
            for _ in range(hybrids[name]):
                if len(called) < cn:
                    called.append(name)

        # Tier 2+: SNV alleles, most-specific-first (require-all-cis consumption). The
        # remaining present gene copies are filled with SNV calls, then *1.
        defs = sorted(
            (d for d in self.allele_defs if d.get("chrom") != "STRUCTURAL"),
            key=lambda d: (d["tier"], -self._num_required(d)),
        )
        flagged: List[str] = []  # combination alleles called from unphased hets
        for d in defs:
            required = self._required_keys(d)
            if not required:
                continue
            while (
                len(called) < cn
                and all(variant_counts.get(k, 0) > 0 for k in required)
                and self._all_cis(variants, required)
            ):
                if len(required) > 1 and not self._phased_evidence(variants, required):
                    flagged.append(d["allele"])
                called.append(d["allele"])
                for k in required:
                    variant_counts[k] -= 1

        while len(called) < cn:
            called.append("*1")

        # Whole-gene deletion: each copy below diploid is a *4 (null) allele.
        n_del = max(0, 2 - copy_number)
        called.extend(["*4"] * n_del)

        # Phase-unresolved trans alternative: a combination allele called from unphased
        # hets could instead be the trans (compound-het) configuration.
        candidates = [list(called)]
        trans = self._trans_alternative(called, flagged)
        if trans and trans not in candidates:
            candidates.append(trans)
        candidates.sort(key=self._diplotype_freq_score, reverse=True)
        primary = candidates[0]
        alternatives = candidates[1:]

        result: Dict[str, Any] = {
            "alleles": primary,
            "method": "phased_consumption_require_all",
            "copy_number": copy_number,
        }
        notes: List[str] = []
        if hybrids:
            result["structural"] = dict(hybrids)
            for name in sorted(hybrids):
                notes.append(f"{name} CYP2A7::CYP2A6 hybrid (fusion junction detected)")
        if n_del:
            notes.append(f"*4 whole-gene deletion x{n_del} (CYP2A6 CN={copy_number})")
        if copy_number is not None and copy_number >= 3 and not hybrids:
            notes.append(
                "Duplication present (CN>=3) - if the extra copy is a functional allele "
                "this is a *1xN-type increased-function duplication; allele assignment "
                "may be ambiguous"
            )
        if alternatives:
            result["alt_alleles"] = alternatives
            result["phase_unresolved"] = True
            others = "; ".join(" / ".join(a) for a in alternatives)
            notes.append(
                "Phase unresolved from short reads - primary is the most likely "
                f"diplotype; also consistent with: {others} "
                "(compound-heterozygote configuration cannot be excluded)"
            )
        elif flagged:
            notes.append(
                "Combination allele called from unphased heterozygous variants - "
                "compound-heterozygote configuration cannot be excluded"
            )
        if notes:
            result["note"] = "; ".join(notes)
        return result

    @staticmethod
    def _diplotype_freq_score(alleles: List[str]) -> float:
        """Relative likelihood of a diplotype = product of its allele frequencies.
        Used only to order equally-consistent unphased calls."""
        score = 1.0
        for a in alleles:
            score *= _ALLELE_FREQ.get(a, _ALLELE_FREQ_DEFAULT)
        return score


class CYP2A6(Gene):
    """CYP2A6 pharmacogene: star-allele calling + nicotine-metabolism phenotype."""

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()
        lines: List[str] = []

        # The *4 whole-gene deletion RETAINS the 3'-flank (crossover), so the liftover
        # whole-gene CN tracks the flank and averages away the deleted exon body. Subtract
        # the deleted exon-body copies (CYP2A6_4del longdel, <=0) to recover the FUNCTIONAL
        # CYP2A6 copy number used for star calling and reporting (e.g. *4/*4 -> 0, not 1).
        whole_gene_cn = result_data["Copy numbers"].get("CYP2A6", 2)
        exon_del = result_data["Copy numbers"].get("CYP2A6_4del", 0)
        cyp2a6_cn = max(0, whole_gene_cn + exon_del)

        # PSV-consistency gate on a called deletion. The *4 exon-body deletion is inferred
        # from liftover DEPTH, but the exon body is the zone of highest CYP2A6/CYP2A7 identity
        # so its uniquely-mappable depth is fragile — it can dip on some platforms/libraries
        # (e.g. an Illumina mappability/coverage artifact) with no copy loss. The CYP2A6-vs-
        # CYP2A7 PSV allele ratio is the reliable gene-differentiation signal: a real het *4
        # depletes the CYP2A6 fraction to ~1/3 (HPRC cohort: true het *4 <=0.37, true hom ->0)
        # while a CN2 locus stays diploid ~0.5 (cohort CN2 >=0.48). If a deletion was called
        # but the observed exon-body CYP2A6 fraction is still diploid-level, the depth drop is
        # an artifact -> veto and revert to CN2. (HG002 Illumina: fraction 0.55 -> vetoed;
        # confirmed CN2 by its diploid assembly + AJ-trio.)
        gate_note: Optional[str] = None
        gate_cfg = self.config.get("psv_del_gate") if isinstance(self.config, dict) else None
        if cyp2a6_cn < 2 and gate_cfg:
            frac, n_sites = self._psv_del_gate_fraction(gate_cfg.get("region", ""))
            thresh = gate_cfg.get("ratio", 0.42)
            min_sites = gate_cfg.get("min_sites", 5)
            if frac is not None and n_sites >= min_sites and frac >= thresh:
                gate_note = (
                    f"*4 deletion call VETOED by CYP2A6/CYP2A7 PSV gate: observed CYP2A6 "
                    f"allele fraction {frac:.2f} over {n_sites} exon-body PSV sites is "
                    f"diploid-level (>= {thresh}), inconsistent with a deletion; the depth "
                    f"drop is a mappability/coverage artifact. CYP2A6 CN {cyp2a6_cn} -> 2."
                )
                cyp2a6_cn = 2
                exon_del = 0

        result_data["Copy numbers"]["CYP2A6"] = cyp2a6_cn
        result_data["Copy numbers"]["CYP2A6_4del"] = exon_del
        cyp2a7_cn = result_data["Copy numbers"].get("CYP2A7", 2)
        lines.append(f"CYP2A6: CN={cyp2a6_cn}")
        if exon_del:
            lines.append(
                f"  (*4 exon-body deletion {exon_del:+d}; 3'-flank retained, "
                f"whole-locus CN={whole_gene_cn})"
            )
        if gate_note:
            lines.append(f"  NOTE: {gate_note}")
        lines.append(f"CYP2A7: CN={cyp2a7_cn}")

        conversions = result_data.get("gene_conversions", []) or []
        if conversions:
            lines.append(f"Structural/conversion events detected: {len(conversions)}")
            for i, conv in enumerate(conversions, 1):
                event = self._interpret_event(conv)
                conv["star_allele"] = event
                lines.append(f"  Event {i}: {event}")
                lines.append(f"    Region: {conv.get('region', '')}")
                lines.append(f"    Allele change: {conv.get('converted_alleles', 0):+d}")
                sites = conv.get("conversion_sites", [])
                if sites:
                    lines.append(f"    Signature sites: {', '.join(sites)}")
        else:
            lines.append("No structural/conversion events detected")

        star_caller = CYP2A6StarAlleleCaller(self.config, ref=self.ref)
        star_result = star_caller.call_alleles(
            result_data.get("Variants", ""), cyp2a6_cn, conversions,
        )
        result_data["star_alleles"] = star_result

        alt_diplotypes = star_result.get("alt_alleles") or []
        if star_result["alleles"]:
            line = "Star Alleles: " + " / ".join(star_result["alleles"])
            if alt_diplotypes:
                line += " (most likely; also possible: " + \
                    "; ".join(" / ".join(a) for a in alt_diplotypes) + ")"
            lines.append(line)
        if star_result.get("note"):
            lines.append(f"  Note: {star_result['note']}")

        phenotype = self._predict_phenotype(star_result.get("alleles", []), star_caller)
        if phenotype:
            alt_phenos = sorted(
                {self._predict_phenotype(a, star_caller) for a in alt_diplotypes}
                - {phenotype}
            )
            if alt_phenos:
                lines.append(
                    f"Nicotine-metabolism phenotype: {phenotype} (primary); "
                    f"alternative configuration(s): {', '.join(alt_phenos)}"
                )
            else:
                suffix = " (all configurations agree)" if alt_diplotypes else ""
                lines.append(f"Nicotine-metabolism phenotype: {phenotype}{suffix}")
            lines.append(
                "  Note: research nicotine-metabolism grouping (activity score); "
                "CYP2A6 has no CPIC dosing guideline"
            )

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data

    def _psv_del_gate_fraction(self, region_str: str):
        """Observed CYP2A6 (gene1) allele fraction over the deletion region's PSV sites,
        from the per-region liftover allele depths (``r.ads`` = ``{pos: [CYP2A6, CYP2A7]}``)
        the CN solve already computed. Depth-weighted, matching the model's ``mean_ads_ratio``.
        Returns ``(fraction, n_sites)``; ``(None, 0)`` when no PSV evidence is available."""
        import re
        m = re.match(r"[^:]+:(\d+)-(\d+)", region_str or "")
        if not m:
            return None, 0
        lo, hi = int(m.group(1)), int(m.group(2))
        seen: Dict[int, Any] = {}
        for regions in getattr(self, "liftover_segdup_regions", []) or []:
            for r in regions:
                for pos, ad in (getattr(r, "ads", None) or {}).items():
                    if lo <= pos <= hi and isinstance(ad, (list, tuple)) and len(ad) >= 2:
                        seen[pos] = ad
        if not seen:
            return None, 0
        first = sum(ad[0] for ad in seen.values())
        both = sum(ad[0] + ad[1] for ad in seen.values())
        if both <= 0:
            return None, 0
        return first / both, len(seen)

    def _interpret_event(self, conv: Dict[str, Any]) -> str:
        """Label a detected CYP2A7<->CYP2A6 structural/conversion event."""
        interp = conv.get("interpretation")
        if isinstance(interp, list):
            for m in interp:
                name = str(m.get("event_name", ""))
                if any(t in name for t in ("*12", "*34", "*46", "*1B", "hybrid", "conversion")):
                    return name
        converted = conv.get("converted_alleles", 0)
        is_fusion = conv.get("is_fusion_candidate", False)
        is_5p = conv.get("fusion_type") == "5_PRIME"
        # Only a 5'-terminus fusion is a *12/*34 hybrid; a 3' event is the benign
        # *1B/*46 3'-UTR conversion (does not alter the core star call).
        if is_fusion and is_5p and converted < 0:
            return "*12 (CYP2A7::CYP2A6 5' hybrid / partial deletion)"
        if is_fusion and is_5p and converted > 0:
            return "CYP2A7::CYP2A6 hybrid duplication"
        return "3'-flank conversion (CYP2A7<->CYP2A6; *1B/*46-type, benign/normal function)"

    def _allele_function_map(self, star_caller: CYP2A6StarAlleleCaller) -> Dict[str, str]:
        m = {"*1": _FUNC_NORMAL}
        for d in star_caller.allele_defs:
            m[d["allele"]] = d.get("function", _FUNC_UNCERTAIN)
        return m

    def _predict_phenotype(
        self, alleles: List[str], star_caller: CYP2A6StarAlleleCaller
    ) -> str:
        """Nicotine-metabolism metabolizer grouping from an activity score (AS).

        AS = sum over ALL called alleles of the per-allele activity contribution
        (none 0.0, decreased 0.5, normal 1.0, increased 1.5). *4 deletion alleles
        contribute 0. Any uncertain-function allele makes the group indeterminate.
        This is a research grouping (nicotine metabolite ratio), NOT a CPIC phenotype;
        CYP2A6 has no CPIC guideline."""
        if not alleles:
            return "Indeterminate"
        fmap = self._allele_function_map(star_caller)
        funcs = [fmap.get(a.split("/")[0].strip(), _FUNC_UNCERTAIN) for a in alleles]
        if any(f == _FUNC_UNCERTAIN for f in funcs):
            return "Indeterminate (uncertain-function allele present)"
        activity = sum(_FUNC_ACTIVITY.get(f, 0.0) for f in funcs)

        if activity >= 3.0:
            return "Ultrarapid metabolizer"
        if activity >= 2.25:
            return "Rapid metabolizer"
        if activity >= 1.75:
            return "Normal metabolizer"
        if activity >= 1.0:
            return "Intermediate metabolizer"
        if activity >= 0.25:
            return "Slow metabolizer"
        return "Poor metabolizer"
