"""CYP2B6 pharmacogene: star-allele calling, CYP2B7P hybrid interpretation, and
CPIC efavirenz metabolizer-phenotype prediction.

CYP2B6 differs from CYP2D6 in two ways that shape this module:

1. **Combination alleles.** The clinically dominant allele *6 is defined by TWO
   variants in cis (c.516G>T rs3745274 + c.785A>G rs2279343); *9 = c.516 alone
   and *4 = c.785 alone. Many other alleles (*7, *13, *34, *36) are 3-4 variant
   combinations. So matching REQUIRES all of an allele's variants to be present
   (and in cis when phased) — the opposite of CYP2D6's opportunistic single-
   defining-variant model, which would over-call *6 on a lone *9 marker.

2. **Low structural burden.** Whole-gene CYP2B6 CNV is rare (typical CN=2). The
   only structural alleles are the reciprocal NAHR products with CYP2B7P: *29
   (partial deletion / 5' hybrid) and *30 (duplication / hybrid). There is no
   *5-style clean deletion or *13/*36/*68-style hybrid family as in CYP2D6.

Phenotype uses the CPIC CYP2B6/efavirenz diplotype->phenotype lookup (UM/RM/NM/
IM/PM), NOT a CYP2D6-style activity score.
"""

import os
from typing import Dict, Any, List, Optional, Tuple
import vcflib
from genecaller.gene import Gene
from genecaller.genes.cyp2d6 import StarAlleleCaller

# Approximate max population allele frequency (PharmVar/CPIC/gnomAD, order-of-
# magnitude only). Used SOLELY to rank equally-consistent unphased diplotypes so
# the most likely one is reported as primary — NOT a clinical frequency estimate.
# Rationale: *6 is a common haplotype (516+785 usually inherited in cis), so a
# {516,785}-het is most likely *6/*1 not *9/*4; but the 3-variant super-alleles
# (*7/*13/*34/*36) are rare, so a {516,785,1459}-het is most likely *6/*5 not *7/*1.
_ALLELE_FREQ_DEFAULT = 0.001
_ALLELE_FREQ = {
    "*1": 0.60, "*6": 0.25, "*5": 0.10, "*18": 0.07, "*4": 0.05, "*9": 0.05,
    "*2": 0.05, "*22": 0.03, "*17": 0.02, "*8": 0.01, "*7": 0.005, "*13": 0.005,
    "*36": 0.005, "*15": 0.005, "*29": 0.005, "*34": 0.003, "*30": 0.003,
}

# LD-rescue pair, keyed by rsID so it is build-agnostic: the genomic positions are
# resolved at runtime from the loaded, build-specific star table, so hg38/hg19/b37
# each get their own coordinates (the locus is a pure coordinate shift between
# builds). The rescued marker c.785A>G (rs2279343) sits in a PSV-sparse window where
# CYP2B7P reads (paralog base == ref A) dilute the real het below the ML caller's bar
# -> it is ML-rejected. But 785 and its anchor c.516G>T (rs3745274) are in strong cis
# LD (together they define *6), and 516 is in a dense PSV block that calls cleanly. So
# a kept-but-rejected 785 is re-activated only when its 516 anchor is a confident PASS
# call and 785 itself carries real ALT support. Contamination cannot fabricate 785
# (paralog is ref there), so this recovers *6 without inflating *9/*4. See genes.yaml
# keep_rejected_positions + gene.py _keep_var.
_LD_ANCHOR_RSID = {"rs2279343": "rs3745274"}  # rescued rsID -> cis anchor rsID
_LD_MIN_ALT = 4      # min ALT reads at the rescued position (guards against noise)
_LD_MIN_FRAC = 0.15  # min ALT fraction (dilution-aware; a het pushed down by ref leak)

# Max total span (bp) for an allele's markers to count as a "tight cluster" eligible
# for sequence-equivalence rescue (see _rescue_clustered_alleles). In the CYP2B6 table
# only *17 (4 adjacent exon-1 SNPs within 11 bp) qualifies; every other multi-marker
# allele spans >2 kb and is matched by exact keys exactly as before.
_CLUSTER_WINDOW_BP = 30


class CYP2B6StarAlleleCaller(StarAlleleCaller):
    """Star-allele caller with require-all-variants (cis) matching for CYP2B6.

    Reuses the CYP2D6 loaders/VCF-extraction/reference-encoding helpers but
    replaces the consumption loop: an allele is only called when EVERY one of
    its defining + associated variants is present, and (when phased) they lie on
    a single haplotype. Alleles are tried most-specific-first so combination
    alleles (e.g. *6 = c.516 + c.785) win over their single-variant subsets
    (*9 = c.516, *4 = c.785).
    """

    def _extract_variants_from_vcf(self, vcf_path: str) -> Dict[str, Dict]:
        """As the base extractor, but annotate each variant with its FILTER
        (low_conf = ML-rejected) and ALT depth, then LD-rescue the kept-rejected
        c.785 call (or drop it). Only positions in keep_rejected_positions can
        appear ML-rejected here — all other rejects were filtered upstream."""
        variants = super()._extract_variants_from_vcf(vcf_path)
        meta = self._read_filter_ad(vcf_path)
        for key, v in variants.items():
            filt, ad = meta.get(key, ("", None))
            v["low_conf"] = "MLrejected" in filt
            v["ad"] = ad
        self._ld_rescue(variants)
        return variants

    def _read_filter_ad(self, vcf_path: str) -> Dict[str, Tuple[str, Any]]:
        """Map each ALT record's variant key -> (FILTER string, (ref_ad, alt_ad))."""
        out: Dict[str, Tuple[str, Any]] = {}
        if not vcf_path or not os.path.exists(vcf_path):
            return out
        try:
            vcf = vcflib.VCF(vcf_path)
            for v in vcf:
                if not v.alt:
                    continue
                filt = v.filter
                filt_s = filt if isinstance(filt, str) else ",".join(filt or [])
                ad_raw = v.samples[0].get("AD") if v.samples else None
                ad_list = None
                if ad_raw is not None:
                    parts = ad_raw.split(",") if isinstance(ad_raw, str) else list(ad_raw)
                    try:
                        ad_list = [int(x) for x in parts]
                    except (ValueError, TypeError):
                        ad_list = None
                alts = v.alt if isinstance(v.alt, list) else [v.alt]
                for i, alt in enumerate(alts):
                    key = f"{self._canon_chrom(v.chrom)}:{v.pos + 1}:{v.ref}:{alt}"
                    pair = None
                    if ad_list and len(ad_list) >= i + 2:
                        pair = (ad_list[0], ad_list[i + 1])
                    out[key] = (filt_s, pair)
            vcf.close()
        except Exception:
            pass
        return out

    @property
    def _ld_anchor(self) -> Dict[int, int]:
        """{rescued_pos: anchor_pos} (1-based), built from the loaded build-specific
        star table via _LD_ANCHOR_RSID, so hg38/hg19/b37 each resolve their own
        coordinates. Cached per instance."""
        cached = getattr(self, "_ld_anchor_cache", None)
        if cached is not None:
            return cached
        pos_by_rsid: Dict[str, int] = {}
        for d in self.allele_defs:
            rsid, pos = d.get("rsid"), d.get("pos")
            if rsid and isinstance(pos, int):
                pos_by_rsid.setdefault(rsid, pos)
        anchor: Dict[int, int] = {}
        for rescued_rsid, anchor_rsid in _LD_ANCHOR_RSID.items():
            rp, ap = pos_by_rsid.get(rescued_rsid), pos_by_rsid.get(anchor_rsid)
            if rp is not None and ap is not None:
                anchor[rp] = ap
        self._ld_anchor_cache = anchor
        return anchor

    def _ld_rescue(self, variants: Dict[str, Dict]) -> None:
        """Promote a low-confidence (ML-rejected) LD-target variant to a real call
        when its cis LD anchor is confidently present and it carries real ALT
        support; drop any other low-confidence variant."""
        by_pos: Dict[int, List[Dict]] = {}
        for v in variants.values():
            by_pos.setdefault(v["pos"], []).append(v)
        for key in list(variants):
            v = variants[key]
            if not v.get("low_conf"):
                continue
            anchor_pos = self._ld_anchor.get(v["pos"])
            rescued = False
            if anchor_pos is not None:
                anchor_ok = any(
                    not a.get("low_conf", False) for a in by_pos.get(anchor_pos, [])
                )
                ad = v.get("ad")
                alt_ok = (
                    bool(ad)
                    and ad[1] >= _LD_MIN_ALT
                    and ad[1] / max(1, ad[0] + ad[1]) >= _LD_MIN_FRAC
                )
                rescued = anchor_ok and alt_ok
            if rescued:
                v["low_conf"] = False
                v["rescued"] = True
            else:
                del variants[key]

    def _required_keys(self, allele_def: Dict[str, Any]) -> List[str]:
        """All variant keys an allele requires: the defining SNV (unless the
        allele is structural) plus every associated variant."""
        keys = []
        if allele_def.get("chrom") != "STRUCTURAL":
            keys.append(self._make_variant_key(allele_def))
        keys.extend(allele_def.get("associated_variants", []))
        return keys

    @staticmethod
    def _num_required(allele_def: Dict[str, Any]) -> int:
        n = len(allele_def.get("associated_variants", []))
        if allele_def.get("chrom") != "STRUCTURAL":
            n += 1
        return n

    def _all_cis(self, variants: Dict[str, Dict], required: List[str]) -> bool:
        """Return False only when two required variants sit on opposite
        haplotypes within the same phase block (clear trans evidence). Unphased
        or single-copy ("both") variants are treated as compatible."""
        haps_by_ps: Dict[Any, set] = {}
        for k in required:
            v = variants.get(k)
            if not v or not v.get("phased") or v.get("ps") is None:
                continue
            hap = v.get("haplotype")
            if hap is None or hap == "both":
                continue
            haps_by_ps.setdefault(v["ps"], set()).add(hap)
        return all(len(haps) <= 1 for haps in haps_by_ps.values())

    def _ref_seq(self, chrom: str, start: int, end: int) -> Optional[str]:
        """Reference bases for 1-based inclusive [start, end] on the build genome,
        or None. Tries the contig with and without a 'chr' prefix (b37 vs hg38/hg19)."""
        if self.ref is None or end < start:
            return None
        for c in (chrom, f"chr{self._canon_chrom(chrom)}", self._canon_chrom(chrom)):
            try:
                s = self.ref.get(c, start - 1, end)
                if s:
                    return s.upper()
            except Exception:
                continue
        return None

    @staticmethod
    def _apply_edits(
        seq: str, wstart: int, edits: List[Tuple[int, str, str]]
    ) -> Optional[str]:
        """Splice (pos1, ref, alt) edits into a window sequence whose first base is
        1-based `wstart`. Applied right-to-left so earlier offsets stay valid. Returns
        None if any edit runs off the window or its REF does not match the sequence
        there (a representation we cannot interpret -> caller declines to rescue)."""
        s = seq
        for pos, ref, alt in sorted(edits, key=lambda e: e[0], reverse=True):
            off = pos - wstart
            if off < 0 or off + len(ref) > len(s):
                return None
            if s[off:off + len(ref)].upper() != ref.upper():
                return None
            s = s[:off] + alt.upper() + s[off + len(ref):]
        return s

    def _rescue_clustered_alleles(
        self, variants: Dict[str, Dict], variant_counts: Dict[str, int]
    ) -> None:
        """Rescue an allele whose markers are a tight cluster (total span
        <= _CLUSTER_WINDOW_BP) that the variant caller emitted in an equivalent but
        non-canonical form -- e.g. CYP2B6*17's three adjacent exon-1 SNPs
        (40991388/90/91) written as an insertion+deletion (40991387 G>GGC +
        40991389 CCG>C). Exact pos:ref:alt matching misses those, so reconstruct the
        sample's alt sequence over the cluster window; if it equals the allele's
        expected alt sequence, synthesize the canonical marker keys (phased cis with
        the backbone SNP) so the normal consumption logic can call the allele.

        Span-gated to clustered alleles ONLY (in the CYP2B6 table just *17; every
        other multi-marker allele spans >2 kb and is left byte-identical), so this
        can never change a non-clustered call. Fails safe: no reference, a REF
        mismatch, or a window that does not reconstruct exactly => no injection."""
        if self.ref is None:
            return
        for d in self.allele_defs:
            if d.get("chrom") == "STRUCTURAL":
                continue
            required = self._required_keys(d)
            if len(required) < 2:
                continue
            if all(variant_counts.get(k, 0) > 0 for k in required):
                continue  # already matchable by exact keys; never disturb it
            markers: List[Tuple[int, str, str]] = []
            for k in required:
                p = k.split(":")
                if len(p) != 4:
                    markers = []
                    break
                try:
                    markers.append((int(p[1]), p[2], p[3]))
                except ValueError:
                    markers = []
                    break
            if not markers:
                continue
            m_start = min(pos for pos, _, _ in markers)
            m_end = max(pos + len(ref) - 1 for pos, ref, _ in markers)
            if m_end - m_start + 1 > _CLUSTER_WINDOW_BP:
                continue  # not a tight cluster -> not this representation problem
            canon = self._canon_chrom(d["chrom"])
            win_vars = [
                v for v in variants.values()
                if self._canon_chrom(v["chrom"]) == canon
                and v["pos"] + len(v["ref"]) - 1 >= m_start - 1
                and v["pos"] <= m_end + 1
            ]
            if not win_vars:
                continue
            wstart = min([m_start] + [v["pos"] for v in win_vars])
            wend = max([m_end] + [v["pos"] + len(v["ref"]) - 1 for v in win_vars])
            ref_seq = self._ref_seq(d["chrom"], wstart, wend)
            if not ref_seq or len(ref_seq) != wend - wstart + 1:
                continue
            expected = self._apply_edits(ref_seq, wstart, markers)
            observed = self._apply_edits(
                ref_seq, wstart, [(v["pos"], v["ref"], v["alt"]) for v in win_vars]
            )
            if expected is None or observed is None or expected != observed:
                continue
            # Matched: take phase from the backbone (defining) SNP so the injected
            # markers are cis with it and the diplotype/phasing logic is satisfied.
            bb = variants.get(self._make_variant_key(d))
            hap = bb.get("haplotype") if bb else None
            ps = bb.get("ps") if bb else None
            phased = bool(bb.get("phased")) if bb else False
            cnt = bb.get("count", 1) if bb else 1
            for k in required:
                if variant_counts.get(k, 0) > 0:
                    continue
                p = k.split(":")
                variants[k] = {
                    "gt": "1/1" if cnt >= 2 else "0/1",
                    "count": cnt,
                    "ps": ps,
                    "phased": phased,
                    "haplotype": hap,
                    "chrom": d["chrom"],
                    "pos": int(p[1]),
                    "ref": p[2],
                    "alt": p[3],
                    "cluster_rescued": True,
                }
                variant_counts[k] = max(variant_counts.get(k, 0), cnt)

    def call_alleles(
        self, vcf_path: str, copy_number: int, conversions: List[Dict],
        intergenic_cn: int = 2,
    ) -> Dict[str, Any]:
        cn = copy_number if copy_number and copy_number > 0 else 2
        max_alleles = cn

        variants = self._extract_variants_from_vcf(vcf_path)
        ld_rescued = sorted(v["pos"] for v in variants.values() if v.get("rescued"))
        variant_counts = {k: v.get("count", 1) for k, v in variants.items()}
        # No-op on GRCh38 (no marker is reference-encoded); wired for a future
        # GRCh37 port where CYP2B7P-shared bases may be reference-encoded.
        self._apply_reference_encoding(variants, variant_counts, cn)

        # Rescue clustered-marker alleles the caller wrote in an equivalent but
        # non-canonical form (e.g. *17's adjacent SNPs emitted as an ins+del).
        # Only clustered alleles (*17) are eligible, so all other calls are unchanged.
        self._rescue_clustered_alleles(variants, variant_counts)

        called: List[str] = []

        # Tier 1: structural hybrids (*29/*30) from two orthogonal signals — the
        # conversion/fusion junction (PSV switch) and the single-copy intergenic CN
        # dosage (CN<2 => *29 partial deletion, CN>2 => *30 duplication). Take the max
        # count per allele so either signal can call the event; agreement = corroboration.
        conv_struct = self._get_structural_alleles(conversions)
        depth_struct = self._structural_from_intergenic(intergenic_cn)
        struct_counts: Dict[str, int] = dict(conv_struct)
        for name, cnt in depth_struct.items():
            struct_counts[name] = max(struct_counts.get(name, 0), cnt)
        struct_sources: Dict[str, List[str]] = {
            name: [
                s for s, present in
                (("depth", name in depth_struct), ("conversion", name in conv_struct))
                if present
            ]
            for name in struct_counts
        }
        for allele_name in sorted(struct_counts):
            for _ in range(struct_counts[allele_name]):
                if len(called) < max_alleles:
                    called.append(allele_name)

        # SNV alleles, most-specific-first (more required variants win ties).
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
                len(called) < max_alleles
                and all(variant_counts.get(k, 0) > 0 for k in required)
                and self._all_cis(variants, required)
            ):
                # A >1-variant allele called from variants that are all
                # heterozygous and unphased could instead be the trans
                # (separate-allele) configuration; record it for the alternative.
                if len(required) > 1 and not self._phased_evidence(variants, required):
                    flagged.append(d["allele"])
                called.append(d["allele"])
                for k in required:
                    variant_counts[k] -= 1

        while len(called) < max_alleles:
            called.append("*1")
        called = called[:max_alleles]

        # Phase unresolved: a combination allele called from unphased hets could
        # instead be the trans configuration (its variants split across
        # haplotypes). Enumerate both, rank by population-frequency likelihood,
        # and report the most likely as primary with the rest as alternatives.
        candidates = [list(called)]
        trans = self._trans_alternative(called, flagged)
        if trans and trans not in candidates:
            candidates.append(trans)
        candidates.sort(key=self._diplotype_freq_score, reverse=True)
        primary = candidates[0]
        alternatives = candidates[1:]

        result = {
            "alleles": primary,
            "method": "phased_consumption_require_all",
            "copy_number": copy_number,
        }
        if struct_counts:
            result["structural"] = dict(struct_counts)
            result["intergenic_cn"] = intergenic_cn
        notes = []
        for name in sorted(struct_counts):
            srcs = struct_sources.get(name, [])
            if len(srcs) >= 2:
                detail = "intergenic-depth dosage and conversion junction agree"
            elif srcs == ["depth"]:
                detail = (
                    "intergenic-depth dosage; no corroborating conversion junction "
                    "call (the PSV-sparse 5' junction can hide it)"
                )
            elif srcs == ["conversion"]:
                detail = "conversion junction; not corroborated by intergenic dosage"
            else:
                detail = "structural"
            notes.append(f"{name} structural allele - {detail}")
        if copy_number is not None and copy_number >= 3 and "*30" not in struct_counts:
            notes.append("Duplication present - allele assignment may be ambiguous")
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
        if ld_rescued:
            result["ld_rescued"] = ld_rescued
            notes.append(
                "c.785 (rs2279343) recovered by 516 LD: real het diluted by CYP2B7P "
                "reads and ML-rejected, re-activated because c.516 is a confident call"
            )
        if notes:
            result["note"] = "; ".join(notes)
        return result

    @staticmethod
    def _diplotype_freq_score(alleles: List[str]) -> float:
        """Relative likelihood of a diplotype = product of its allele
        frequencies. Used only to order equally-consistent unphased calls."""
        score = 1.0
        for a in alleles:
            score *= _ALLELE_FREQ.get(a, _ALLELE_FREQ_DEFAULT)
        return score

    def _trans_alternative(
        self, called: List[str], flagged: List[str]
    ) -> Optional[List[str]]:
        """Build the trans-configuration alternative diplotype for combination
        alleles that were called from unphased hets. Each flagged allele X paired
        with a `*1` slot is split into two named sub-alleles (subA/subB) whose
        required variants partition X's — e.g. *7 -> *6 + *5, *6 -> *9 + *4.
        Returns the alternative allele list, or None if it equals `called` or no
        decomposition exists."""
        alt = list(called)
        for xname in flagged:
            decomp = self._decompose_allele(xname)
            if not decomp:
                continue
            if xname not in alt or "*1" not in alt:
                continue
            sub_a, sub_b = decomp
            alt[alt.index(xname)] = sub_a
            alt[alt.index("*1")] = sub_b
        return alt if alt != called else None

    def _decompose_allele(self, name: str) -> Optional[tuple]:
        """If `name`'s required variants partition into exactly two other named
        alleles, return (subA, subB) with subA the larger/most-specific; else
        None. Used to enumerate the trans configuration of a combination allele."""
        xdef = next((d for d in self.allele_defs if d["allele"] == name), None)
        if not xdef:
            return None
        full = set(self._required_keys(xdef))
        if len(full) < 2:
            return None
        cands = []
        for d in self.allele_defs:
            if d["allele"] == name:
                continue
            rk = set(self._required_keys(d))
            if rk and rk < full:  # proper subset
                cands.append((d["allele"], d.get("tier", 9), rk))
        cands.sort(key=lambda c: (-len(c[2]), c[1]))  # largest sub-allele first
        for na, _, ra in cands:
            need = full - ra
            for nb, _, rb in cands:
                if nb != na and rb == need:
                    return (na, nb)
        return None

    @staticmethod
    def _phased_evidence(variants: Dict[str, Dict], required: List[str]) -> bool:
        """True if at least two required variants share a phase block on the same
        haplotype (positive cis evidence) or any is homozygous ('both')."""
        seen = {}
        for k in required:
            v = variants.get(k)
            if not v:
                continue
            if v.get("haplotype") == "both":
                return True
            if v.get("phased") and v.get("ps") is not None and v.get("haplotype") in (0, 1):
                key = (v["ps"], v["haplotype"])
                seen[key] = seen.get(key, 0) + 1
                if seen[key] >= 2:
                    return True
        return False

    def _get_structural_alleles(self, conversions: List[Dict]) -> Dict[str, int]:
        """Map detected CYP2B7P<->CYP2B6 hybrids to *29/*30.

        A 5' fusion (CYP2B7P sequence at the 5' end / a copy-number loss) is the
        *29 partial deletion; a copy-number gain is the *30 duplication hybrid.
        """
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
            if "*30" in blob:
                structural["*30"] = structural.get("*30", 0) + count
            elif "*29" in blob or (conv.get("is_fusion_candidate") and converted < 0):
                structural["*29"] = structural.get("*29", 0) + count
        return structural

    @staticmethod
    def _structural_from_intergenic(intergenic_cn: Optional[int]) -> Dict[str, int]:
        """Depth-based *29/*30 from the single-copy CYP2B6/CYP2B7P intergenic CN
        variable. That ~40 kb is uniquely mappable, so its copy number is a direct
        dosage readout of the reciprocal intron4-junction NAHR events that the
        whole-gene CN averages away: CN<2 => *29 (partial deletion, alleles lost),
        CN>2 => *30 (duplication). Returns {} at the CN=2 baseline, so normal samples
        are unaffected."""
        out: Dict[str, int] = {}
        if intergenic_cn is None:
            return out
        if intergenic_cn < 2:
            out["*29"] = 2 - intergenic_cn
        elif intergenic_cn > 2:
            out["*30"] = intergenic_cn - 2
        return out


# CPIC allele-function -> label used by the phenotype rules.
_FUNC_NORMAL = "Normal"
_FUNC_INCREASED = "Increased"
_FUNC_DECREASED = "Decreased"
_FUNC_NONE = "No Function"
_FUNC_UNCERTAIN = "Uncertain Function"


class CYP2B6(Gene):
    """CYP2B6 pharmacogene with star-allele calling and CPIC phenotype."""

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()
        lines: List[str] = []

        cyp2b6_cn = result_data["Copy numbers"].get("CYP2B6", 2)
        cyp2b7p_cn = result_data["Copy numbers"].get("CYP2B7P", 2)
        intergenic_cn = result_data["Copy numbers"].get("CYP2B6_intergenic", 2)
        lines.append(f"CYP2B6: CN={cyp2b6_cn}")
        lines.append(f"CYP2B7P: CN={cyp2b7p_cn}")
        if intergenic_cn != 2:
            ev = "*29 partial deletion" if intergenic_cn < 2 else "*30 duplication"
            lines.append(
                f"CYP2B6/CYP2B7P intergenic CN={intergenic_cn} -> {ev} "
                f"({abs(2 - intergenic_cn)} allele(s); single-copy dosage probe over the "
                "intron4-junction NAHR interval)"
            )

        conversions = result_data.get("gene_conversions", []) or []
        if conversions:
            lines.append(f"Structural/conversion events detected: {len(conversions)}")
            for i, conv in enumerate(conversions, 1):
                region = conv.get("region", "")
                converted = conv.get("converted_alleles", 0)
                hybrid = self._interpret_hybrid(conv)
                conv["star_allele"] = hybrid
                lines.append(f"  Event {i}: {hybrid}")
                lines.append(f"    Region: {region}")
                lines.append(f"    Allele change: {converted:+d}")
                sites = conv.get("conversion_sites", [])
                if sites:
                    lines.append(f"    Signature sites: {', '.join(sites)}")
        else:
            lines.append("No structural/conversion events detected")

        star_caller = CYP2B6StarAlleleCaller(self.config, ref=self.ref)
        star_result = star_caller.call_alleles(
            result_data.get("Variants", ""), cyp2b6_cn, conversions,
            intergenic_cn=intergenic_cn,
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
                    f"Phenotype: {phenotype} (primary); "
                    f"alternative configuration(s): {', '.join(alt_phenos)}"
                )
            else:
                suffix = " (all configurations agree)" if alt_diplotypes else ""
                lines.append(f"Phenotype: {phenotype}{suffix}")
            efavirenz = self._efavirenz_note(phenotype)
            if efavirenz:
                lines.append(f"  Efavirenz: {efavirenz}")

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data

    def _interpret_hybrid(self, conv: Dict[str, Any]) -> str:
        """Label a detected CYP2B7P<->CYP2B6 hybrid event."""
        interp = conv.get("interpretation")
        if isinstance(interp, list):
            for m in interp:
                name = str(m.get("event_name", ""))
                if "*29" in name or "*30" in name or "hybrid" in name.lower():
                    return name
        converted = conv.get("converted_alleles", 0)
        is_fusion = conv.get("is_fusion_candidate", False)
        if is_fusion and converted < 0:
            return "*29 (CYP2B7P-CYP2B6 hybrid / partial deletion)"
        if is_fusion and converted > 0:
            return "*30 (CYP2B7P/CYP2B6 hybrid duplication)"
        if converted > 0:
            return "Conversion (CYP2B7P->CYP2B6)"
        return "Conversion (CYP2B6->CYP2B7P)"

    def _allele_function_map(self, star_caller: CYP2B6StarAlleleCaller) -> Dict[str, str]:
        m = {"*1": _FUNC_NORMAL}
        for d in star_caller.allele_defs:
            m[d["allele"]] = d.get("function", _FUNC_UNCERTAIN)
        return m

    def _predict_phenotype(
        self, alleles: List[str], star_caller: CYP2B6StarAlleleCaller
    ) -> str:
        """CPIC CYP2B6 diplotype -> metabolizer phenotype (efavirenz guideline)."""
        if not alleles:
            return "Indeterminate"
        fmap = self._allele_function_map(star_caller)
        funcs = [fmap.get(a.split("/")[0].strip(), _FUNC_UNCERTAIN) for a in alleles]

        inc = funcs.count(_FUNC_INCREASED)
        nor = funcs.count(_FUNC_NORMAL)
        dec = funcs.count(_FUNC_DECREASED)
        nof = funcs.count(_FUNC_NONE)
        unc = funcs.count(_FUNC_UNCERTAIN)
        total = len(funcs)

        if unc > 0:
            return "Indeterminate (uncertain-function allele present)"

        if total != 2:
            # Non-diploid (hybrid deletion/duplication): summarize by function.
            if nof + dec == total and total > 0:
                return "Poor metabolizer (no functional alleles; non-diploid)"
            return f"Indeterminate (non-diploid: {total} alleles)"

        if inc == 2:
            return "Ultrarapid metabolizer"
        if inc == 1 and nor == 1:
            return "Rapid metabolizer"
        if nor == 2:
            return "Normal metabolizer"
        if nof == 2 or dec == 2 or (dec == 1 and nof == 1):
            return "Poor metabolizer"
        # remaining diploid mixes: normal/increased with one decreased/no-fn
        if (nor + inc) == 1 and (dec + nof) == 1:
            return "Intermediate metabolizer"
        return "Indeterminate"

    @staticmethod
    def _efavirenz_note(phenotype: str) -> str:
        if phenotype.startswith("Poor"):
            return "consider 400 or 200 mg/day (elevated exposure; higher adverse-effect risk)"
        if phenotype.startswith("Intermediate"):
            return "consider 400 mg/day"
        if phenotype.startswith(("Normal", "Rapid", "Ultrarapid")):
            return "standard 600 mg/day"
        return ""
