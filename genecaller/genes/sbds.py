"""SBDS gene-specific class for Shwachman-Diamond syndrome (SBDS / SBDSP1).

SBDS (chr7q11.21, - strand) and its pseudogene SBDSP1 (~5.8 Mb distal, + strand,
~97% identical) recombine by gene conversion. >90% of SDS patients carry an exon-2
conversion variant, and two recurrent pseudogene-derived PSVs — c.258+2T>C (splice
donor of intron 2) and c.183_184delinsCT (p.Lys62* / frameshift) — account for ~80%
of disease alleles. Standard short-read pipelines miscall these via SBDSP1
mismapping; this class reports the conversion status plus copy number, mirroring
genes/gba.py.

Direction: the pathogenic conversion is SBDSP1 -> SBDS (pseudogene sequence copied
into the functional gene). SDS is autosomal recessive, so two pathogenic alleles
(homozygous, or compound heterozygous across the two PSVs) => affected; one => carrier.

v1 scope: SBDS/SBDSP1 copy number, paralog-differentiated small variants (base
machinery), the two recurrent conversion PSVs typed by zygosity, and segmenter-detected
conversion tracts. Full pathogenic-variant interpretation beyond the two recurrent PSVs
is a deferred layer; this class is its extension point.
"""

import os
from typing import Dict, Any, List
from genecaller.gene import Gene
from genecaller.conversion_detector import GeneConversionDetector
import vcflib


class SBDS(Gene):
    """SBDS gene: SBDSP1->SBDS gene-conversion interpretation for Shwachman-Diamond syndrome."""

    def __init__(self, cfg: dict, ref_file: str, **kwargs) -> None:
        super().__init__(cfg, ref_file, **kwargs)
        # Load conversion variants from the shared known-conversion DB (recipient == SBDS).
        self.conversion_variants = self._load_conversion_variants_from_db()

    def _load_conversion_variants_from_db(self) -> Dict[str, Dict[str, Any]]:
        """Load SBDS-specific single-position conversion variants from the DB.

        Returns {event_name: {chrom, pos, ref, alt}} for SBDS-recipient rows that
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
            if kc.recipient_gene != "SBDS":
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
            f"Loaded {len(variants)} conversion variants from database for SBDS"
        )
        return variants

    def _check_conversion_variants(self, vcf_path: str) -> List[Dict[str, Any]]:
        """Scan the output VCF for the known SBDS conversion variants.

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
        lines = ["SBDS (Shwachman-Diamond syndrome) Analysis:"]

        cns = result_data.get("Copy numbers", {})
        sbds_cn = cns.get("SBDS", self.baseline_cn)
        sbdsp1_cn = cns.get("SBDSP1", self.baseline_cn)
        lines.append(f"  SBDS:   CN={sbds_cn}")
        lines.append(f"  SBDSP1: CN={sbdsp1_cn}")

        # Recurrent pathogenic conversion PSVs (SBDSP1 -> SBDS direction).
        conversion_variants = self._check_conversion_variants(
            result_data.get("Variants", "")
        )
        pathogenic_alleles = sum(c["alleles"] for c in conversion_variants)

        if conversion_variants:
            lines.append("")
            lines.append("Recurrent SBDSP1->SBDS conversion variants detected:")
            for c in conversion_variants:
                lines.append(
                    f"  {c['type']} {c['variant']} "
                    f"(chr7:{c['pos']} {c['ref']}>{c['alt']}; pseudogene-derived)"
                )

        # Segmenter-detected conversion tracts (SBDSP1 -> SBDS = positive direction).
        # These are RATIO-based and NOT diagnostic on their own: a positive-direction
        # tract can arise from allelic dropout of a divergent SBDS haplotype (one gene
        # allele's reads fail to lift, inflating the pseudogene ratio to ~2/3), which
        # mimics a heterozygous conversion (observed on healthy HG005). Pathogenic calls
        # therefore come from the directly typed PSV genotypes above; a tract's
        # known-event match is surfaced as pathogenic ONLY when CORROBORATED by such a
        # direct call, otherwise it is flagged non-diagnostic.
        typed_events = {c["variant"] for c in conversion_variants}
        gene_conversions = [
            c for c in result_data.get("gene_conversions", []) if c["converted_alleles"] > 0
        ]
        if gene_conversions:
            lines.append("")
            lines.append(
                f"Positive-direction conversion tract(s): {len(gene_conversions)} "
                "— ratio-based, non-diagnostic alone (can reflect allelic dropout of a "
                "divergent haplotype; pathogenic status is set by the direct PSV calls above)"
            )
            for i, conv in enumerate(gene_conversions, 1):
                lines.append(
                    f"  Tract {i}: {conv.get('region', '?')}, "
                    f"converted_alleles={conv.get('converted_alleles', 0)}"
                )
                interp = conv.get("interpretation")
                if isinstance(interp, list):
                    for match in interp:
                        ev = match.get("event_name", "?")
                        score = (
                            f" (score {match['match_score']:.1%})"
                            if "match_score" in match
                            else ""
                        )
                        if ev in typed_events:
                            lines.append(
                                f"    matches {ev}{score} — CORROBORATED by direct PSV call"
                            )
                        else:
                            lines.append(
                                f"    matches {ev}{score} — NOT corroborated by a direct PSV "
                                "call (non-diagnostic; likely allelic-dropout artifact)"
                            )

        # Recessive-disease summary: 2 pathogenic alleles (homozygous or compound
        # het across the two PSVs) => affected; 1 => carrier.
        lines.append("")
        if pathogenic_alleles >= 2:
            lines.append(
                f"SUMMARY: {pathogenic_alleles} pathogenic SBDS conversion alleles "
                "detected — consistent with Shwachman-Diamond syndrome (bi-allelic). "
                "Confirm; SDS is autosomal recessive."
            )
        elif pathogenic_alleles == 1:
            lines.append(
                "SUMMARY: 1 pathogenic SBDS conversion allele detected — carrier "
                "(heterozygous)."
            )
        else:
            lines.append(
                "SUMMARY: no recurrent SBDS conversion variant detected."
            )
        lines.append(
            "  NOTE: only the two recurrent exon-2 PSVs (c.258+2T>C, c.183_184delinsCT) "
            "are typed here; other SBDS pathogenic variants are reported via the base "
            "variant calls."
        )

        result_data["cn_interpretation"] = "\n".join(lines)
        return result_data
