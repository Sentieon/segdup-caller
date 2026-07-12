"""OTOA gene-specific class for DFNB22 nonsyndromic hearing loss (OTOA / OTOAP1).

OTOA (chr16p12.2, + strand) and its pseudogene OTOAP1 (~0.8 Mb proximal, + strand,
~99% identical over OTOA's 3' exons) form a direct segmental duplication. DFNB22 is
autosomal recessive; pathogenic OTOA alleles arise three ways:
  1. a recurrent ~500 kb NAHR deletion that removes the whole gene (dominant mechanism);
  2. OTOAP1 -> OTOA gene conversion transferring the pseudogene stop p.Glu787*
     (c.2359G>T, exon 22) into the functional gene (Laurent 2021, PMID 33492714);
  3. independent point mutations.

The same p.Glu787* is a notorious short-read false positive from OTOAP1 mismapping
(Kim 2025, PMID 39943967), so it is typed from the paralog-aware genotype and gated by
corroboration, mirroring genes/sbds.py. Two pathogenic alleles (any combination of
deletion + conversion) => affected (DFNB22); one => carrier.

v1 scope: OTOA/OTOAP1 copy number (recurrent deletion via the OTOA paralog CN and the
optional longdel probe), paralog-differentiated small variants, the recurrent p.Glu787*
conversion PSV typed by zygosity with a corroboration gate, and segmenter-detected
conversion tracts. NB a ~15 kb mid-block PSV desert (chr16:21,739-21,754k, 5 coding
exons) is a dark region where paralog assignment is unreliable — variants there are
low-confidence and should be confirmed orthogonally.
"""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.conversion_detector import GeneConversionDetector
import vcflib


class OTOA(Gene):
    """OTOA gene: OTOAP1->OTOA conversion + recurrent-deletion interpretation for DFNB22."""

    def __init__(self, cfg: dict, ref_file: str, **kwargs) -> None:
        super().__init__(cfg, ref_file, **kwargs)
        # Load conversion variants from the shared known-conversion DB (recipient == OTOA).
        self.conversion_variants = self._load_conversion_variants_from_db()

    def _load_conversion_variants_from_db(self) -> Dict[str, Dict[str, Any]]:
        """Load OTOA-specific single-position conversion variants from the DB.

        Returns {event_name: {chrom, pos, ref, alt}} for OTOA-recipient rows that
        carry concrete ref/alt alleles (skips complex/multi-variant events).
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
            if kc.recipient_gene != "OTOA":
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
            f"Loaded {len(variants)} conversion variants from database for OTOA"
        )
        return variants

    def _check_conversion_variants(self, vcf_path: str) -> List[Dict[str, Any]]:
        """Scan the output VCF for the known OTOA conversion variants.

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
        lines = ["OTOA (DFNB22 nonsyndromic hearing loss) Analysis:"]

        cns = result_data.get("Copy numbers", {})
        otoa_cn = cns.get("OTOA", self.baseline_cn)
        otoap1_cn = cns.get("OTOAP1", self.baseline_cn)
        lines.append(f"  OTOA:   CN={otoa_cn}")
        lines.append(f"  OTOAP1: CN={otoap1_cn}")

        # ---- Pathogenic mechanism 1: recurrent whole-gene deletion ----
        # OTOA copy loss = pathogenic DFNB22 allele(s). Primary signal is the OTOA
        # paralog CN (3' block allele imbalance); the longdel probe over the unique 5'
        # corroborates when configured (reported as OTOA_recurrentDEL in Copy numbers).
        del_alleles = max(0, self.baseline_cn - otoa_cn)
        longdel_cn = cns.get("OTOA_recurrentDEL")
        if longdel_cn:  # present and non-zero only when the deletion fires
            del_alleles = max(del_alleles, abs(longdel_cn))
        if del_alleles > 0:
            lines.append("")
            probe = (
                " (5'-unique longdel probe corroborates)"
                if longdel_cn
                else ""
            )
            lines.append(
                f"OTOA copy loss: {del_alleles} allele(s) deleted{probe} — "
                "recurrent DFNB22 ~500 kb NAHR deletion of the whole gene."
            )

        # ---- Pathogenic mechanism 2: recurrent conversion PSV p.Glu787* ----
        # The direct PSV genotype at p.Glu787* is prone to OTOAP1 mismapping — a minority
        # of pseudogene reads leak into the OTOA-assigned pileup, yielding spurious het
        # calls with variable VAF (~0.2-0.5) that a VAF filter cannot separate (cohort
        # Phase G-b: ~4% false carriers). A *real* OTOAP1->OTOA conversion transferring
        # p.Glu787* also converts the neighbouring PSVs, so it is CORROBORATED by a
        # positive-direction conversion tract SPANNING the PSV position. A conversion
        # allele therefore counts as pathogenic only when so corroborated; an
        # uncorroborated direct call is surfaced as LOW CONFIDENCE (likely mismapping) —
        # never silently dropped, never a definitive carrier.
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
            lines.append("Recurrent OTOAP1->OTOA conversion variant(s) — tract-corroborated:")
            for c in corroborated:
                lines.append(
                    f"  {c['type']} {c['variant']} "
                    f"({self.chr}:{c['pos']} {c['ref']}>{c['alt']}; pseudogene-derived stop, "
                    "conversion tract spans the PSV)"
                )
        if uncorroborated:
            lines.append("")
            lines.append(
                "LOW-CONFIDENCE p.Glu787*-type PSV call(s) — NOT counted as a pathogenic "
                "allele: the PSV was directly typed but no conversion tract spans the site, "
                "the signature of OTOAP1 mismapping rather than a real conversion. Confirm "
                "orthogonally:"
            )
            for c in uncorroborated:
                lines.append(
                    f"  {c['type']} {c['variant']} "
                    f"({self.chr}:{c['pos']} {c['ref']}>{c['alt']}; direct PSV call only)"
                )

        # Segmenter-detected positive-direction conversion tracts (OTOAP1 -> OTOA), shown
        # as context. Ratio-based and non-diagnostic on their own (allelic dropout of a
        # divergent OTOA haplotype, or OTOAP1 mismapping); a tract drives a pathogenic call
        # only where it spans a directly-typed PSV (handled above).
        if pos_tracts:
            lines.append("")
            lines.append(
                f"Positive-direction conversion tract(s): {len(pos_tracts)} "
                "— ratio-based context, non-diagnostic alone (allelic dropout / OTOAP1 "
                "mismapping); pathogenic only where a tract spans a typed PSV (see above)"
            )
            for i, conv in enumerate(pos_tracts, 1):
                lines.append(
                    f"  Tract {i}: {conv.get('region', '?')}, "
                    f"converted_alleles={conv.get('converted_alleles', 0)}"
                )

        # ---- Recessive-disease summary: combine deletion + corroborated conversion ----
        # DFNB22 is autosomal recessive; a deleted allele and a tract-corroborated
        # converted allele are distinct alleles, so pathogenic = del_alleles + conv_alleles
        # (conv_alleles counts only corroborated p.Glu787* calls; see mechanism 2).
        pathogenic_alleles = del_alleles + conv_alleles
        lines.append("")
        if pathogenic_alleles >= 2:
            lines.append(
                f"SUMMARY: {pathogenic_alleles} pathogenic OTOA allele(s) "
                f"(deletion={del_alleles}, conversion={conv_alleles}) — consistent with "
                "DFNB22 hearing loss (bi-allelic). Confirm; DFNB22 is autosomal recessive."
            )
        elif pathogenic_alleles == 1:
            kind = "deletion" if del_alleles else "conversion"
            lines.append(
                f"SUMMARY: 1 pathogenic OTOA allele ({kind}) — carrier (heterozygous)."
            )
        else:
            msg = "SUMMARY: no recurrent OTOA deletion or conversion allele detected."
            if uncorroborated:
                msg += (
                    " (An uncorroborated p.Glu787*-type PSV was detected — low-confidence, "
                    "likely OTOAP1 mismapping; see above.)"
                )
            lines.append(msg)
        lines.append(
            "  NOTE: only the recurrent deletion and the p.Glu787* conversion PSV are typed "
            "here; other OTOA pathogenic point mutations are reported via the base variant "
            "calls. Variants in the ~15 kb PSV desert mid-block of the OTOA 3' duplication "
            "(5 coding exons) are low-confidence (dark paralog region) — confirm orthogonally."
        )

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data
