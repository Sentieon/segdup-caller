"""STRC gene-specific class for DFNB16 autosomal-recessive hearing loss (STRC / STRCP1).

STRC (chr15q15.3, - strand) and its pseudogene STRCP1 (~80 kb distal, - strand,
~99.6% identical coding) form a direct tandem segmental duplication. DFNB16 is the second
most common cause of autosomal-recessive hearing loss after GJB2 (GJB2/DFNB1). Pathogenic
STRC alleles arise mainly three ways:
  1. a recurrent ~100 kb NAHR deletion that removes the whole STRC copy (dominant mechanism);
     77-88% of these co-delete the neighbouring unique gene CATSPER2 (which sits between the
     STRC and STRCP1 repeat copies), adding male infertility => Deafness-Infertility Syndrome
     (OMIM #611102);
  2. STRCP1 -> STRC gene conversion transferring pseudogene-inactivating sequence (uncommon);
  3. independent point mutations.

STRC and STRCP1 are ~100% identical over exons 1-15; only exons 16-29 (the 3' half) carry
PSVs, so paralog assignment - and therefore which copy is deleted - rests on those 3' PSVs.
The per-paralog copy number (STRC vs STRCP1) is what distinguishes a pathogenic STRC deletion
from a benign STRCP1 pseudogene deletion; depth over the block alone cannot. CATSPER2 is unique
sequence, so its copy number is a clean depth probe for the co-deletion extent.

v1 scope: STRC/STRCP1 copy number (recurrent deletion via the STRC paralog CN), the CATSPER2
co-deletion extent (deafness-infertility), paralog-differentiated small variants, and - where
configured - STRCP1 -> STRC conversion PSVs typed by zygosity with a corroboration gate
(mirroring genes/otoa.py and genes/sbds.py). Two pathogenic STRC alleles => affected (DFNB16);
one => carrier. Benign 3-copy gains of STRC or STRCP1 (~1.85% of the population) are not called
pathogenic.
"""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.conversion_detector import GeneConversionDetector
import vcflib


class STRC(Gene):
    """STRC gene: recurrent-deletion + CATSPER2-extent + STRCP1->STRC conversion for DFNB16."""

    def __init__(self, cfg: dict, ref_file: str, **kwargs) -> None:
        super().__init__(cfg, ref_file, **kwargs)
        # Load conversion variants from the shared known-conversion DB (recipient == STRC).
        self.conversion_variants = self._load_conversion_variants_from_db()

    def _load_conversion_variants_from_db(self) -> Dict[str, Dict[str, Any]]:
        """Load STRC-specific single-position conversion variants from the DB.

        Returns {event_name: {chrom, pos, ref, alt}} for STRC-recipient rows that
        carry concrete ref/alt alleles (skips complex/multi-variant events). Empty
        until a STRCP1->STRC row is added to known_gene_conversions.tsv (Phase D).
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

        variants: Dict[str, Dict[str, Any]] = {}
        for kc in known_conversions:
            if kc.recipient_gene != "STRC":
                continue
            if kc.ref == "NA" or kc.alt == "NA":
                continue  # skip complex events without a single ref/alt

            try:
                locus_parts = kc.recipient_locus.split(":")
                if len(locus_parts) != 2:
                    continue
                chrom = locus_parts[0]
                pos_str = locus_parts[1]
                pos = int(pos_str.split("-")[0]) if "-" in pos_str else int(pos_str)
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
            f"Loaded {len(variants)} conversion variants from database for STRC"
        )
        return variants

    def _check_conversion_variants(self, vcf_path: str) -> List[Dict[str, Any]]:
        """Scan the output VCF for the known STRC conversion variants.

        Returns a list of {variant, type (mono-/bi-allelic), alleles, pos, ref, alt}.
        """
        found_conversions: List[Dict[str, Any]] = []

        if not os.path.exists(vcf_path) or not self.conversion_variants:
            return found_conversions

        try:
            vcf = vcflib.VCF(vcf_path)
            for v in vcf:
                if not v.alt:
                    continue

                chrom = v.chrom if v.chrom.startswith("chr") else f"chr{v.chrom}"
                pos = v.pos + 1  # vcflib v.pos is 0-based

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
                                "alleles": alt_count,
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
        lines = ["STRC (DFNB16 autosomal-recessive hearing loss) Analysis:"]

        cns = result_data.get("Copy numbers", {})
        strc_cn = cns.get("STRC", self.baseline_cn)
        strcp1_cn = cns.get("STRCP1", self.baseline_cn)
        catsper2_cn = cns.get("CATSPER2")  # present only when the CATSPER2 probe is configured
        lines.append(f"  STRC:     CN={strc_cn}")
        lines.append(f"  STRCP1:   CN={strcp1_cn}")
        if catsper2_cn is not None:
            lines.append(f"  CATSPER2: CN={catsper2_cn}")

        # ---- Pathogenic mechanism 1: recurrent whole-gene STRC deletion ----
        # STRC copy loss = pathogenic DFNB16 allele(s), read from the STRC paralog CN (the 3'
        # PSVs split the STRC/STRCP1 block). This is the copy that matters: a STRCP1-side
        # deletion of the same block is benign (handled below).
        del_alleles = max(0, self.baseline_cn - strc_cn)
        if del_alleles > 0:
            lines.append("")
            lines.append(
                f"STRC copy loss: {del_alleles} allele(s) deleted — recurrent DFNB16 "
                "~100 kb NAHR deletion of the functional STRC copy."
            )
            # CATSPER2-extent: CATSPER2 sits between the STRC and STRCP1 repeats and is
            # co-deleted in 77-88% of these events. A reduced CATSPER2 CN => the deletion
            # extends over it (deafness-infertility); an intact CATSPER2 => STRC-limited.
            if catsper2_cn is not None:
                if catsper2_cn < self.baseline_cn:
                    lines.append(
                        f"  Extent: co-deletes CATSPER2 (CN={catsper2_cn}) — Deafness-"
                        "Infertility Syndrome (OMIM #611102): the contiguous deletion also "
                        "removes CATSPER2, adding male infertility to the hearing loss."
                    )
                else:
                    lines.append(
                        f"  Extent: CATSPER2 intact (CN={catsper2_cn}) — STRC-limited "
                        "deletion (hearing loss without the CATSPER2 infertility component)."
                    )

        # ---- Benign STRCP1-side deletion of the block ----
        # A copy lost on the STRCP1 (pseudogene) side is not DFNB16: functional STRC is intact.
        # Report it (with any CATSPER2 co-deletion, which still carries infertility relevance)
        # so a block copy loss is never silently read as pathogenic.
        strcp1_del = max(0, self.baseline_cn - strcp1_cn)
        if strcp1_del > 0 and del_alleles == 0:
            lines.append("")
            msg = (
                f"STRCP1 copy loss: {strcp1_del} allele(s) — benign pseudogene deletion; "
                f"NOT DFNB16 (functional STRC intact, CN={strc_cn})."
            )
            if catsper2_cn is not None and catsper2_cn < self.baseline_cn:
                msg += (
                    f" NB CATSPER2 is co-deleted (CN={catsper2_cn}): CATSPER2 loss causes "
                    "male infertility (OMIM #611102) independently of the hearing-loss status."
                )
            lines.append(msg)

        # ---- Benign copy gains ----
        # STRC/STRCP1 3-copy gains are a common benign polymorphism (~1.85%); flag as benign
        # so a CN>2 is not mistaken for anything pathogenic.
        if strc_cn > self.baseline_cn or strcp1_cn > self.baseline_cn:
            lines.append("")
            lines.append(
                f"Copy gain (STRC CN={strc_cn}, STRCP1 CN={strcp1_cn}) — benign polymorphism "
                "(~1.85% of the population carry a 3-copy gain); not a pathogenic allele."
            )

        # ---- Pathogenic mechanism 2: STRCP1 -> STRC conversion PSVs (corroboration-gated) ----
        # STRCP1->STRC conversion copies pseudogene-inactivating sequence into the functional
        # gene. Direct PSV genotypes are prone to STRCP1 mismapping (paralog reads MAPQ ~6-15),
        # so a converted allele counts as pathogenic only when a positive-direction conversion
        # tract SPANS the PSV; an uncorroborated direct call is surfaced as LOW CONFIDENCE.
        # Inert until a STRCP1->STRC row is added to the conversion DB (Phase D).
        conversion_variants = self._check_conversion_variants(
            result_data.get("Variants", "")
        )
        pos_tracts = [
            c
            for c in result_data.get("gene_conversions", [])
            if c.get("converted_alleles", 0) > 0
        ]

        def _spanning_tract(pos: int) -> bool:
            for t in pos_tracts:
                region = t.get("region", "")
                if ":" not in region:
                    continue
                try:
                    a, b = (int(x) for x in region.split(":")[1].split("-"))
                except (ValueError, IndexError):
                    continue
                if a <= pos <= b:
                    return True
            return False

        corroborated, uncorroborated = [], []
        for c in conversion_variants:
            (corroborated if _spanning_tract(c["pos"]) else uncorroborated).append(c)
        conv_alleles = sum(c["alleles"] for c in corroborated)

        if corroborated:
            lines.append("")
            lines.append("STRCP1->STRC conversion variant(s) — tract-corroborated:")
            for c in corroborated:
                lines.append(
                    f"  {c['type']} {c['variant']} "
                    f"({self.chr}:{c['pos']} {c['ref']}>{c['alt']}; pseudogene-derived, "
                    "conversion tract spans the PSV)"
                )
        if uncorroborated:
            lines.append("")
            lines.append(
                "LOW-CONFIDENCE conversion PSV call(s) — NOT counted as a pathogenic allele: "
                "the PSV was directly typed but no conversion tract spans the site, the "
                "signature of STRCP1 mismapping rather than a real conversion. Confirm "
                "orthogonally:"
            )
            for c in uncorroborated:
                lines.append(
                    f"  {c['type']} {c['variant']} "
                    f"({self.chr}:{c['pos']} {c['ref']}>{c['alt']}; direct PSV call only)"
                )

        # Segmenter-detected positive-direction conversion tracts (STRCP1 -> STRC), shown as
        # context. STRC has no crisp canonical conversion allele (unlike OTOA p.Glu787*), so a
        # tract is ratio-based and non-diagnostic on its own (allelic dropout of a divergent STRC
        # haplotype, or STRCP1 mismapping); it drives a pathogenic call only where it spans a
        # directly-typed PSV (handled above). Surfaced so a real conversion signal is not hidden.
        if pos_tracts:
            lines.append("")
            lines.append(
                f"Positive-direction conversion tract(s): {len(pos_tracts)} "
                "— ratio-based context, non-diagnostic alone (allelic dropout / STRCP1 "
                "mismapping); pathogenic only where a tract spans a typed PSV (see above)"
            )
            for i, conv in enumerate(pos_tracts, 1):
                lines.append(
                    f"  Tract {i}: {conv.get('region', '?')}, "
                    f"converted_alleles={conv.get('converted_alleles', 0)}"
                )

        # ---- Recessive-disease summary: combine STRC deletion + corroborated conversion ----
        # DFNB16 is autosomal recessive; a deleted STRC allele and a tract-corroborated
        # converted allele are distinct pathogenic alleles.
        pathogenic_alleles = del_alleles + conv_alleles
        lines.append("")
        if pathogenic_alleles >= 2:
            lines.append(
                f"SUMMARY: {pathogenic_alleles} pathogenic STRC allele(s) "
                f"(deletion={del_alleles}, conversion={conv_alleles}) — consistent with "
                "DFNB16 hearing loss (bi-allelic). Confirm; DFNB16 is autosomal recessive."
            )
        elif pathogenic_alleles == 1:
            kind = "deletion" if del_alleles else "conversion"
            lines.append(
                f"SUMMARY: 1 pathogenic STRC allele ({kind}) — carrier (heterozygous)."
            )
        else:
            msg = "SUMMARY: no pathogenic STRC deletion or conversion allele detected."
            if uncorroborated:
                msg += (
                    " (An uncorroborated conversion PSV was detected — low-confidence, "
                    "likely STRCP1 mismapping; see above.)"
                )
            lines.append(msg)
        lines.append(
            "  NOTE: only the recurrent deletion, its CATSPER2 extent, and configured "
            "STRCP1->STRC conversion PSVs are interpreted here; other STRC pathogenic point "
            "mutations are reported via the base variant calls. STRC and STRCP1 are ~100% "
            "identical over exons 1-15, so variants there are low-confidence (no discriminating "
            "PSVs) — confirm orthogonally."
        )

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data
