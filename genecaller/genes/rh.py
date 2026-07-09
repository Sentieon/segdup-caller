"""RH gene-specific class for the Rh blood group locus (RHD / RHCE).

RHD and RHCE are highly homologous (~97-98%) paralogs on chr1p36.11 in *inverted*
(tail-to-tail) orientation. The dominant clinical event is the whole-gene RHD
deletion (RhD-negative), reported here as RHD copy number / predicted RhD zygosity.

Antigen layer (this class), driven by the curated ISBT allele table
data/<build>/rh_alleles.tsv (see segdup/rh/isbt/gen_rh_alleles.py):

  * RhD refinement — **weak D** and **DEL** alleles are typed by require-ALL-variants
    matching against the ISBT definitions (each allele = a defining SNV set). The
    most-specific matching allele is reported (a sub-allele suppresses its parent).
    This is genotype/zygosity-level and needs NO phasing.
  * RhCE **E/e** antigen is typed at its determinant site (c.676 Ala226Pro). Antigens are
    codominant, so presence/zygosity is read directly from the site's genotype (phase-free;
    validated cohort E+ ~27%, matching population).

  NOT typed here (see doc/RH.md):
  * **C/c** — the C determinant is the RHCE exon-2 block (c.307 Pro103Ser + linked), which
    arose from an RHD->RHCE conversion and falls in the caller's exon-2 blind spot
    (validated 0-1/226 cohort calls vs pop AF ~0.42). Not readable from the result VCF.
  * **Partial-D** (RHD-CE-D hybrids) — only *flagged* as candidates from gene-conversion
    events (need cis-phase).
  * Exact RHCE haplotype naming (Ce/cE/CE/ce).
  * RHD*Psi and indel/MNV-defined alleles (the SNV-only require-all table excludes them).

The inverted orientation is handled entirely by the base GeneMapping (reverse map),
so this subclass carries NO orientation-specific logic (cf. genes/ikbkg.py).
"""

from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.util import get_data_file
from genecaller.genes.cyp2d6 import StarAlleleCaller


class RHAntigenCaller(StarAlleleCaller):
    """Marker-based RH antigen typing from the curated ISBT allele table.

    Reuses the base VCF extractor + variant-key helpers only; the CYP diplotype
    machinery is not used. RH has no reference-encoded backbone to reconstruct, so
    backbone loading is skipped and a purpose-built table loader is used.
    """

    def __init__(self, config: Dict[str, Any], ref: Any = None):
        self.ref = ref
        self.alleles = self._load_rh_alleles(config)
        self.backbone_vars: List[Dict[str, Any]] = []

    @staticmethod
    def _load_rh_alleles(config: Dict[str, Any]) -> List[Dict[str, Any]]:
        path = get_data_file(config["star_allele_file"])
        assert path is not None, f"RH allele file not found: {config['star_allele_file']}"
        alleles = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if parts[0] == "allele":
                    continue
                parts = (parts + [""] * 6)[:6]
                allele, gene, category, antigen, variants, protein = parts
                varkeys = set()
                for v in variants.split(","):
                    toks = v.strip().split(":")
                    if len(toks) == 3:                       # pos:ref:alt; skip malformed
                        varkeys.add(f"1:{toks[0]}:{toks[1]}:{toks[2]}")
                if not varkeys:
                    continue
                varkeys = frozenset(varkeys)
                alleles.append(
                    dict(allele=allele, gene=gene, category=category, antigen=antigen,
                         protein=protein, varkeys=varkeys)
                )
        return alleles

    def type_rhd(self, variants: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Require-ALL-variants match of RHD weak/DEL alleles, keeping only the most
        specific (an allele whose SNV set is a proper subset of another matched allele's
        set is dropped). Alleles that share an *identical* defining-SNV set are collapsed
        into one hit (they are SNV-indistinguishable, e.g. a weak-D and a DEL that differ
        only by a feature not in the table) — reported as `category='ambiguous'` rather
        than a spurious compound. zygosity = homozygous iff every defining SNV is
        homozygous. Each hit carries `categories` (the set of ISBT categories it spans)."""
        matched = [a for a in self.alleles
                   if a["gene"] == "RHD" and all(k in variants for k in a["varkeys"])]
        groups: Dict[Any, List[Dict[str, Any]]] = {}
        for a in matched:
            groups.setdefault(a["varkeys"], []).append(a)
        kept = []
        for vk, group in groups.items():
            if any(vk < other for other in groups):   # a strictly more-specific set matched
                continue
            hom = all(variants[k].get("count", 1) >= 2 for k in vk)
            cats = sorted({a["category"] for a in group})
            kept.append({
                "allele": " / ".join(a["allele"] for a in group),
                "category": cats[0] if len(cats) == 1 else "ambiguous",
                "categories": cats,
                "antigen": " / ".join(dict.fromkeys(a["antigen"] for a in group)),
                "protein": group[0]["protein"],
                "zygosity": "homozygous" if hom else "heterozygous",
            })
        return kept

    def type_rhce(self, variants: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        """Per-site codominant C/c and E/e typing. For each determinant site the alt
        base is the named (upper-case) antigen and the reference base is the other; a
        called alt (het/hom) means the named antigen is present, its absence means the
        reference antigen. Returns {system: {present, absent, genotype, called}}."""
        profile = {}
        for a in self.alleles:
            if a["category"] != "antigen_site":
                continue
            try:
                up, lo = a["antigen"].split(" ")[0].split("/")  # "E/e" -> ("E","e")
                key = next(iter(a["varkeys"]))                  # single site
            except (ValueError, StopIteration):                 # malformed row -> skip, don't crash
                continue
            v = variants.get(key)
            if v is None:
                pres, absent, geno = lo, up, f"{lo}{lo}"     # ref-hom -> reference antigen
            elif v.get("count", 1) >= 2:
                pres, absent, geno = up, lo, f"{up}{up}"     # alt-hom
            else:
                pres, absent, geno = f"{up}{lo}", "", f"{up}{lo}"  # het -> both codominant
            profile[f"{up}/{lo}"] = {"present": pres, "absent": absent,
                                     "genotype": geno, "called": v is not None}
        return profile

    def call(self, vcf_path: str) -> Dict[str, Any]:
        variants = self._extract_variants_from_vcf(vcf_path)
        return {"rhd_hits": self.type_rhd(variants), "rhce": self.type_rhce(variants),
                "n_variants": len(variants)}


class RH(Gene):
    """RH (RHD/RHCE): RhD zygosity from RHD CN, ISBT marker-based weak-D/DEL refinement,
    RhCE C/c + E/e antigen profile, and gene-conversion (partial-D candidate) summary."""

    # Exact base-type match only: types 1/2/3 are the established "manage as RhD-positive"
    # weak-D types. Sub-alleles (1.1, 1.2, …) and higher types deliberately fall through to
    # the conservative caution message — a prefix match would wrongly catch "type 10".."19".
    _MANAGE_POSITIVE = ("weak D type 1", "weak D type 2", "weak D type 3")

    @staticmethod
    def rhd_zygosity(rhd_cn: int) -> str:
        """Map RHD copy number to a predicted RhD zygosity statement."""
        if rhd_cn <= 0:
            return "RhD-negative — no intact RHD (homozygous RHD deletion/absence)"
        if rhd_cn == 1:
            return ("RhD-positive, hemizygous — single RHD copy "
                    "(one RhD-negative allele in trans)")
        if rhd_cn == 2:
            return "RhD-positive, homozygous — two RHD copies"
        return f"atypical RHD copy number ({rhd_cn}) — possible duplication"

    def _rhd_refine_lines(self, hits: List[Dict[str, Any]], rhd_cn: int) -> List[str]:
        if not hits:
            return ["  RhD antigen markers: none of the curated weak-D/DEL alleles "
                    "detected — conventional RHD."]
        lines = ["  RhD antigen alleles detected:"]
        for h in hits:
            lines.append(f"    {h['allele']} ({h['antigen']}) — {h['zygosity']}")
        cats_of = lambda h: h.get("categories") or [h.get("category")]
        weak = [h for h in hits if "weak" in cats_of(h)]
        dele = [h for h in hits if "DEL" in cats_of(h)]
        if any(h.get("category") == "ambiguous" for h in hits):
            lines.append("    => an allele above is SNV-indistinguishable between the "
                         "listed types (weak-D vs DEL); serology required to resolve.")
        expressed_alone = (rhd_cn <= 1
                           or any(h["zygosity"] == "homozygous" for h in hits)
                           or len(hits) >= 2)
        if len(hits) >= 2:
            lines.append("    => two reduced-D alleles: compound reduced-D phenotype "
                         "likely (exact cis-haplotypes need phasing).")
        elif not expressed_alone:
            lines.append("    => single allele heterozygous at CN2: D is likely expressed "
                         "from the apparently-normal RHD in trans (carrier).")
        if weak and expressed_alone:
            managed = all(h["antigen"] in self._MANAGE_POSITIVE for h in weak)
            if managed:
                lines.append("    => Weak D (types 1/2/3): reduced D antigen, "
                             "conventionally managed as RhD-positive (RhIg not indicated) "
                             "— confirm per transfusion policy.")
            else:
                lines.append("    => Weak D (type not in 1/2/3, e.g. type 4/DAR): reduced D; "
                             "some types can form anti-D — do NOT assume RhD-positive "
                             "management; confirm serologically.")
        if dele and expressed_alone:
            lines.append("    => DEL: RHD present but D undetectable by routine serology; "
                         "donor vs recipient management differs — confirm per policy.")
        lines.append("    NOTE: genotype-based prediction (require-all-variants); confirm "
                     "serologically. Partial-D hybrids are flagged separately below.")
        return lines

    @staticmethod
    def _rhce_lines(profile: Dict[str, Dict[str, Any]], rhce_cn: int) -> List[str]:
        p = profile.get("E/e") if profile else None
        if not p or rhce_cn <= 0:
            return []
        ee = (f"E{'+' if 'E' in p['present'] else '-'}"
              f"e{'+' if 'e' in p['present'] else '-'}")
        # An absent E-determinant (called=False) is inferred hom-ref: a hom-ref call and an
        # uncovered site are indistinguishable in the result VCF, so say so rather than
        # asserting a positively-observed genotype.
        if p["called"]:
            head = (f"  Predicted RhCE E/e antigen: {ee}  (genotype {p['genotype']}; "
                    "c.676 Ala226Pro, phase-free)")
        else:
            head = (f"  Predicted RhCE E/e antigen: {ee}  (no c.676 E-determinant variant "
                    "called — E-negative only if the site is covered)")
        out = [head]
        if rhce_cn != 2:
            out.append(f"    NOTE: E/e +/- zygosity assumes a diploid RHCE, but RHCE CN="
                       f"{rhce_cn} here — the antigen call may be approximate.")
        out.append("    NOTE: C/c is NOT typed — the RHCE exon-2 C determinant lies in the "
                   "RHD->RHCE gene-conversion region, a caller blind spot. Exact RHCE "
                   "haplotype and RHCE hybrids are not resolved here.")
        return out

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()
        cns = result_data.get("Copy numbers", {})
        rhd_cn = cns.get("RHD", self.baseline_cn)
        rhce_cn = cns.get("RHCE", self.baseline_cn)

        lines = ["RH (Rh blood group) Analysis:"]
        lines.append(f"  RHD:  CN={rhd_cn}")
        lines.append(f"  RHCE: CN={rhce_cn}")
        lines.append(f"  Predicted RhD zygosity: {self.rhd_zygosity(rhd_cn)}")

        cfg = getattr(self, "config", None)
        if cfg and cfg.get("star_allele_file"):
            typed = None
            try:
                caller = RHAntigenCaller(cfg, ref=getattr(self, "ref", None))
                typed = caller.call(result_data.get("Variants", ""))
            except Exception:
                typed = None
            if typed is not None and typed.get("n_variants", 0) > 0:
                if rhd_cn > 0:
                    lines.extend(self._rhd_refine_lines(typed["rhd_hits"], rhd_cn))
                lines.extend(self._rhce_lines(typed["rhce"], rhce_cn))
            else:
                # Distinguish "typing did not run / no data" from a real negative call —
                # a missing/empty result VCF or parse error must not read as conventional.
                lines.append("  RH antigen typing unavailable — no variants read from the "
                             "result VCF (missing/empty VCF or parse error); CN calls above "
                             "are unaffected.")

        # Gene conversion events (RHD-CE-D hybrids / partial-D candidates).
        gene_conversions = result_data.get("gene_conversions", [])
        if gene_conversions:
            lines.append("")
            lines.append(f"Gene conversion events: {len(gene_conversions)}")
            for i, conv in enumerate(gene_conversions, 1):
                n = conv.get("converted_alleles", 0)
                lines.append(
                    f"  Event {i}: {conv.get('region', '?')}, converted_alleles={n}"
                )
                if n:
                    lines.append(
                        "    RHCE-like sequence within RHD — RHD-CE-D hybrid "
                        "(partial-D) candidate; confirm serologically."
                    )
                sites = conv.get("conversion_sites", [])
                if sites:
                    lines.append(f"    Sites: {', '.join(sites)}")
                interp = conv.get("interpretation")
                if isinstance(interp, list):
                    for match in interp:
                        score = (
                            f" (score {match['match_score']:.1%})"
                            if "match_score" in match
                            else ""
                        )
                        lines.append(f"    Known event: {match['event_name']}{score}")

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data
