"""LPA (KIV-2 VNTR) caller.

Unlike the paralog genes, LPA has no pseudogene: the clinical signal is the copy
number of the KIV-2 coding tandem repeat (~5.5 kb unit, ~2-80 total copies), whose
count inversely sets plasma Lp(a) / ASCVD risk. GRCh38 collapses the array to ~6
reference copies at chr6:160,611,568-160,646,868.

This class subclasses ``Gene`` only for the registry's ``issubclass`` contract; it does
NOT run the paralog CN solver (capped at CN 0-6) and defines its own ``__init__`` with
no diff_vcf / map / GeneMapping. It reuses the framework's existing primitives:

  (1) ``Bam.liftover`` collapses ALL array reads onto a single unit (no ``.dat``:
      ``gene_mapping=None`` -> the failed-region step is skipped). Forcing every unit
      onto one copy rescues the multimappers (~90% become MAPQ>=60), so depth and AD
      work on the collapsed pileup.
  (2) Total count = slope * (unit_depth / flank_depth), the unit depth read off the
      collapsed BAM with ``Bam.get_depth`` (MAPQ0), the flank the diploid CN2 baseline.
      Validated r=0.994 vs HPRC assembly truth (docs/LPA_BUILD_PLAN.md Probe #2).
  (3) Allele split = a reference-allele-pair deconvolution over a 50-marker unit
      panel. The pooled alt-fraction vector (read via ``Bam.get_pileup`` on the
      collapsed pileup) is the copy-count-weighted mixture of the two alleles'
      per-copy alt-fractions; we find the pair of catalogued KIV-2 allele-haplotypes
      whose best mixture matches it, and the mix fraction gives the split. A
      two-signal confidence gate abstains to total-only when the split is unreliable:
      the two-allele fit must beat the best single allele (gain) AND be identifiable
      -- few near-optimal pairs may disagree on the count (spread). This catches the
      novel-allele failure mode (a sample whose alleles are absent from the library
      matches a wrong pair at low residual, but with a wide spread) instead of emitting
      a confident wrong split. ~54% of samples confidently split (86% of those within 2
      copies of assembly truth; 1.1% leak) vs ~42% for single fixed markers; the rest
      report total-only. Validated leave-one-out on 181 HPRC samples, the caller's own
      collapse. (The confident fraction is honest, not maximal: a denser library finds more
      plausible competing pairs and correctly abstains on genuinely ambiguous samples.)

No small variants are called (``skip_variant_call``); the acceptance gate is count
concordance vs assembly truth, not F1.
"""

import copy
from typing import Any, Dict, List, Optional, Tuple

from genecaller.gene import Gene
from genecaller.bam_process import Bam
from genecaller.util import get_data_file
from genecaller.logging import get_logger

logger = get_logger("Gene.LPA")


class LPA(Gene):
    # Opt out of the orchestrator's paralog liftover; we collapse the array ourselves.
    needs_liftover = False

    def __init__(self, cfg: dict, ref_file: str, sex: Optional[str] = None) -> None:
        # NOTE: intentionally does not call Gene.__init__ -- LPA has no paralog
        # diff_vcf/map/GeneMapping and uses none of the CN solver. Only the
        # orchestrator contract is implemented.
        self.cfg = cfg
        self.ref_file = ref_file
        self.ref = ref_file  # Bam() takes the reference path
        self.sample_sex = sex
        self.gene_names = cfg.get("gene_names", ["LPA"])
        self.config = dict(cfg.get("config", {}))
        self.config.setdefault("skip_variant_call", True)
        self.baseline_cn = cfg.get("baseline_cn", 2)
        self.variant_call_threads = 4
        self.converted_segments = []

        c = self.config
        self.kiv2_array = c["kiv2_array"]        # array of ref copies to collapse
        self.kiv2_unit = c["kiv2_unit"]          # single-copy target unit
        self.flank_regions = c["flank_regions"]  # unique L+R flanks = CN2 baseline
        self.ref_copies = c.get("ref_copies", 6)
        self.min_marker_depth = c.get("min_marker_depth", 20)
        self.min_markers = c.get("min_markers", 8)   # min panel markers covered to deconvolve
        self.gain_gate = c.get("gain_gate", 1.25)    # 2-allele-vs-best-single fit-gain confidence gate
        self.spread_gate = c.get("spread_gate", 4)   # max near-optimal minor-count spread (identifiability)
        self.spread_margin = c.get("spread_margin", 1.15)  # residual factor defining "near-optimal"
        cf = c.get("catalog_file")
        self.catalog_file = get_data_file(cf) if cf else None
        self._lib: Optional[Dict[str, Any]] = None   # lazily-loaded {POS, ALT, LIB}
        self._split_diag: Optional[Dict[str, Any]] = None
        self._result: Optional[Dict[str, Any]] = None

    # --- orchestrator no-ops (LPA has no paralog CN merge / gene conversion) ---
    def merge_regions_by_cn(self) -> None:  # noqa: D401
        pass

    def detect_gene_conversions(self) -> None:  # noqa: D401
        pass

    def call_small_var(self) -> None:  # skipped via skip_variant_call; stub for safety
        pass

    def resolve_phased(self, params: dict) -> None:
        pass

    # --- the estimator ---
    def call_cn(self, read_data: dict, params: dict) -> None:
        sr = read_data["short_read"]
        src_bam = sr["bam"]
        lo_params = copy.deepcopy(sr["params"])
        lo_params["prefix"] = lo_params["prefix"] + ".lpa_collapse"
        collapsed_fname = lo_params["prefix"] + ".bam"

        # (1) collapse the whole KIV-2 array onto one unit (reuse liftover, no .dat)
        src_bam.liftover(
            extract=self.kiv2_array,
            target=self.kiv2_unit,
            out_bam=collapsed_fname,
            gene_mapping=None,
            failed_margin=0,
        )
        collapsed = Bam(collapsed_fname, src_bam.ref, lo_params)

        # (2) total count from the collapse depth ratio. RAW coverage (MAPQ0), NOT the
        # CNVscope-normalized get_depth: the collapsed BAM has no genome-wide coverage to
        # normalize against, so CNVscope emits no windows for it. Raw unit/flank coverage
        # is the validated primitive (Probe #1/#2) and both are on the same scale.
        unit_dp = self._raw_depth(src_bam, collapsed_fname, self.kiv2_unit, "lpa.unit")
        flank_dp = self._raw_depth(src_bam, src_bam.bam, self.flank_regions, "lpa.flank")
        ratio = (unit_dp / flank_dp) if flank_dp else 0.0
        # The flank is the diploid (2-copy) baseline; the collapsed unit carries all N
        # haploid KIV-2 copies -> N = 2 * ratio. dp_norm corrects the collapse's block
        # depth inflation -- the same segdup-depth knob the paralog genes use (STRC 1.08,
        # RCCX/HBA ~1.1) -- calibrated vs assembly truth.
        dp_norm = float(self.config.get("dp_norm", 1.0)) or 1.0
        total = 2.0 * ratio / dp_norm
        total_units = max(0, int(round(total)))

        # (3) allele split from marker-panel AD on the collapsed pileup
        alleles, phaseable, marker_detail = self._allele_split(collapsed, total_units)

        logger.info(
            f"KIV-2 total={total_units} (ratio={ratio:.3f}, unit_dp={unit_dp:.1f}, "
            f"flank_dp={flank_dp:.1f}); phaseable={phaseable}"
            + (f" split={alleles[0]}+{alleles[1]}" if alleles else "")
        )
        self._result = {
            "total_units": total_units,
            "total_units_raw": round(total, 2),
            "ratio": round(ratio, 4),
            "unit_depth": round(unit_dp, 2),
            "flank_depth": round(flank_dp, 2),
            "phaseable": phaseable,
            "alleles": alleles,
            "markers": marker_detail,
            "split_diag": self._split_diag,
        }

    def _raw_depth(self, ref_bam: Bam, bam_path: str, regions: str, tag: str) -> float:
        """Length-weighted mean RAW coverage (MAPQ0) over `regions` via CoverageMetrics.

        Uses the same `sentieon driver ... CoverageMetrics --min_map_qual 0` the caller
        already relies on. Raw (un-normalized) so the collapsed-unit and flank depths are
        directly comparable, which the CNVscope path is not for a sparse collapsed BAM.
        """
        import subprocess
        import os
        from genecaller import throttle

        tmpdir = ref_bam.tmpdir
        bed = os.path.join(tmpdir, f"{os.path.basename(bam_path)}.{tag}.bed")
        with open(bed, "w") as fh:
            for reg in regions.split(","):
                chrom, span = reg.split(":")
                start, end = span.split("-")
                fh.write(f"{chrom}\t{int(start) - 1}\t{int(end)}\n")
        out = os.path.join(tmpdir, f"{os.path.basename(bam_path)}.{tag}.cov")
        cmd = (
            f"{ref_bam.sentieon} driver -i {bam_path} -r {ref_bam.ref} "
            f"--interval {bed} --algo CoverageMetrics --min_map_qual 0 {out}"
        )
        with throttle.slot():
            subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        num = den = 0.0
        with open(out + ".sample_interval_summary") as fh:
            next(fh)  # header
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                span = p[0].split(":")[1]
                start, end = span.split("-")
                length = int(end) - int(start) + 1
                num += float(p[2]) * length  # col 3 = average_coverage
                den += length
        return num / den if den else 0.0

    def _load_library(self) -> None:
        """Lazily load the shipped unit-haplotype library (POS, ALT, LIB)."""
        import numpy as np

        if self._lib is None:
            d = np.load(self.catalog_file)
            self._lib = {"POS": d["POS"], "ALT": d["ALT"], "LIB": d["LIB"]}

    def _allele_split(
        self, collapsed: Bam, total_units: int
    ) -> Tuple[Optional[List[int]], bool, List[Dict[str, Any]]]:
        """Reference-allele-pair deconvolution -> minor allele copy count (if confident).

        The pooled alt-fraction over the 50-marker unit panel is the copy-count-weighted
        mixture of the two alleles' per-copy alt-fractions. We find the pair of catalogued
        reference alleles whose optimal mixture best matches the observation; the mixing
        fraction maps to the split. Each marker's AD is read off the collapsed pileup with
        the Bam.get_pileup primitive; unit-local panel positions map to genomic coordinates
        via the (build-specific) kiv2_unit start, so the library itself is build-independent.

        Gate: the two-allele fit must beat the best SINGLE reference allele by a factor of
        >= gain_gate. When the sample carries an allele absent from the library, no pair
        fits well, the gain collapses toward 1, and we abstain (report total-only) rather
        than emit a confident wrong split.
        """
        import numpy as np

        if self.catalog_file is None or total_units <= 0:
            return None, False, []
        self._load_library()
        pos, alt, lib = self._lib["POS"], self._lib["ALT"], self._lib["LIB"]
        chrom = self.kiv2_unit.split(":")[0]
        unit_start = int(self.kiv2_unit.split(":")[1].split("-")[0])

        f = np.full(len(pos), np.nan)
        detail: List[Dict[str, Any]] = []
        for i, p in enumerate(pos):
            # unit-local (1-based) -> genomic 0-based; get_pileup expects a 0-based pos
            # (it queries region pos+1 and matches pysam reference_pos == pos).
            g = unit_start + int(p) - 2
            rec = collapsed.get_pileup(chrom, g)
            if rec is None:
                continue
            ad = rec.AD  # {base: count}
            tot = sum(ad.values())
            if tot < self.min_marker_depth:
                continue
            altc = int(ad.get(str(alt[i]), 0))
            frac = altc / tot
            f[i] = frac
            detail.append(
                {"pos": int(p), "alt": str(alt[i]),
                 "alt_reads": altc, "depth": int(tot), "alt_frac": round(frac, 3)}
            )
        ncov = int(np.sum(~np.isnan(f)))
        if ncov < self.min_markers:
            self._split_diag = {"markers_covered": ncov, "gated": "too_few_markers"}
            return None, False, detail

        minor, gain, resid, spread = self._deconvolve(
            f, total_units, lib, self.spread_margin
        )
        # Two-signal gate: the split must both (a) improve on a single allele (gain) and
        # (b) be identifiable -- few near-optimal pairs disagree on it (spread). Novel
        # alleles absent from the library slip past gain alone (they fit a wrong pair at low
        # residual) but betray themselves with a wide spread.
        phaseable = gain >= self.gain_gate and spread <= self.spread_gate
        self._split_diag = {
            "markers_covered": ncov, "minor_pred": int(minor),
            "gain": round(gain, 3), "gain_gate": self.gain_gate,
            "spread": int(spread), "spread_gate": self.spread_gate,
            "resid": round(resid, 4),
            "gated": None if phaseable else ("low_gain" if gain < self.gain_gate else "wide_spread"),
        }
        if not phaseable:
            return None, False, detail
        return [total_units - minor, minor], True, detail

    @staticmethod
    def _deconvolve(
        f, total_units: int, lib, spread_margin: float = 1.3
    ) -> Tuple[int, float, float, int]:
        """Best two-allele mixture over the library. Returns (minor, gain, resid, spread).

        For each ordered pair of library alleles (gi, gj) solve the mixing fraction x that
        minimizes ||x*gi + (1-x)*gj - f||^2 over the covered markers; keep the best-fitting
        pair. Two orthogonal confidence signals are returned alongside the split:

          gain  = (best single-allele residual) / (best pair residual): how much two distinct
                  alleles improve on one. Low gain -> the sample looks homogeneous.
          spread = range of the implied minor-copy count across ALL near-optimal pairs (those
                  within spread_margin x the best residual): the split's identifiability. A
                  novel-allele sample matches several wrong pairs about equally well and gives
                  a wide spread even at low residual -- which a residual/gain test alone misses.
        """
        import numpy as np

        mask = ~np.isnan(f)
        ff = f[mask]
        G = np.asarray(lib)[:, mask]
        k = ff.shape[0]
        pairs: List[Tuple[float, float]] = []  # (residual, x)
        M = G.shape[0]
        for i in range(M):
            gi = G[i]
            for j in range(i, M):
                d = gi - G[j]
                denom = float(d @ d)
                if denom < 1e-9:
                    x = 0.5
                else:
                    x = min(1.0, max(0.0, float((ff - G[j]) @ d / denom)))
                pred = x * gi + (1.0 - x) * G[j]
                res = float(np.sum((pred - ff) ** 2)) / k
                pairs.append((res, x))
        pairs.sort(key=lambda t: t[0])
        best_res, best_x = pairs[0]
        single = float(np.min(np.mean((G - ff) ** 2, axis=1)))
        n1 = int(round(best_x * total_units))
        minor = max(0, min(n1, total_units - n1))
        # gain: floor the denominator so a near-perfect pair fit reads as high confidence (not
        # 0); a homogeneous sample (single fits as well as the pair) still lands at gain ~1.
        gain = single / max(best_res, 1e-6)

        def _minor(x: float) -> int:
            a = int(round(x * total_units))
            return min(a, total_units - a)

        thr = best_res * spread_margin
        near = [_minor(x) for res, x in pairs if res <= thr]
        spread = (max(near) - min(near)) if near else 0
        return minor, gain, best_res, spread

    def prepare_output(self) -> Dict[str, Any]:
        r = self._result or {}
        total = r.get("total_units", 0)
        cns: Dict[str, Any] = {"LPA_KIV2_total": total}
        alleles = r.get("alleles")
        if alleles:
            cns["LPA_KIV2_allele1"] = alleles[0]
            cns["LPA_KIV2_allele2"] = alleles[1]

        lines = [
            "LPA (Lp(a) / KIV-2 VNTR) analysis:",
            f"  KIV-2 total copies (both alleles): {total}",
        ]
        if r.get("phaseable") and alleles:
            lines.append(f"  allele split: {alleles[0]} + {alleles[1]} copies")
            lines.append(
                "  (smaller allele drives plasma Lp(a): fewer KIV-2 copies "
                "-> higher Lp(a) -> higher ASCVD risk)"
            )
        else:
            lines.append(
                "  allele split: not resolvable (allele absent from reference library or "
                "subtype-homozygous; ~46% of samples report total-only)"
            )
        lines.append(
            "  Qualitative Lp(a)/ASCVD association only (RUO); no quantitative Lp(a) value."
        )

        return {
            "Copy numbers": cns,
            "Variants": "",
            "kiv2_total_units": total,
            "kiv2_alleles": alleles,
            "phaseable": bool(r.get("phaseable", False)),
            "kiv2_detail": {
                k: r.get(k)
                for k in ("ratio", "unit_depth", "flank_depth", "total_units_raw")
            },
            "split_diag": r.get("split_diag"),
            "markers": r.get("markers", []),
            "interpretation": "\n".join(lines),
        }
