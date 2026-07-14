# LPA (KIV-2 VNTR) — Lp(a) / ASCVD copy number

## Overview

`LPA` encodes apolipoprotein(a), the protein component of lipoprotein(a) [Lp(a)]. Its
signature feature is KIV-2, a ~5.5 kb coding tandem repeat (kringle IV type-2) present in
~2–80 copies per allele. Copy number is inversely correlated with plasma Lp(a): fewer KIV-2
copies → higher Lp(a) → higher risk of atherosclerotic cardiovascular disease (ASCVD) and
calcific aortic valve disease. The smaller (low-copy) allele dominates the circulating Lp(a) level.

Unlike the other segdup-caller genes, LPA has no pseudogene/paralog — the clinical signal is the
count of a within-gene VNTR, not a gene-vs-paralog copy-number split. GRCh38 (and GRCh37) collapse
the array to 6 reference copies (`chr6:160,611,568-160,646,868` in hg38), so short reads over the
array multi-map and are discarded by ordinary genotyping. This module recovers the count by
collapsing the whole array onto a single unit.

This is a research-use-only (RUO), qualitative Lp(a)/ASCVD association; the module reports KIV-2
copy number, not a quantitative Lp(a) concentration.

## What this module does

- KIV-2 total copy number (both alleles summed) — the primary output, always reported.
- Allele split (the two per-allele KIV-2 counts) — reported when it can be resolved confidently;
  otherwise the module reports total-only.
- No small variants are called (`skip_variant_call`); the acceptance gate is copy-number
  concordance vs assembly truth, not variant F1.

## Genomic structure (GRCh38)

| Feature | Coordinates (hg38) | Role |
|---|---|---|
| KIV-2 array (collapsed) | chr6:160,611,568-160,646,868 | the ~6 reference copies; all reads here are collapsed |
| KIV-2 unit | chr6:160,612,109-160,617,655 (5,547 bp) | single-copy collapse target + marker frame |
| Left / right unique flanks | chr6:160,607,000-160,611,000 / chr6:160,647,500-160,651,500 | diploid CN=2 depth baseline |

## Output examples

### Confident allele split (HPRC sample NA18945)

```yaml
LPA:
  Copy numbers:
    LPA_KIV2_total: 42
    LPA_KIV2_allele1: 22
    LPA_KIV2_allele2: 20
  interpretation: |-
    LPA (Lp(a) / KIV-2 VNTR) analysis:
      KIV-2 total copies (both alleles): 42
      allele split: 22 + 20 copies
      (smaller allele drives plasma Lp(a): fewer KIV-2 copies -> higher Lp(a) -> higher ASCVD risk)
      Qualitative Lp(a)/ASCVD association only (RUO); no quantitative Lp(a) value.
```

### Total-only, split not resolvable (HPRC sample HG02523)

The two alleles are too similar, or not represented in the reference library, so the module
reports the total but abstains on the split:

```yaml
LPA:
  Copy numbers:
    LPA_KIV2_total: 39
  interpretation: |-
    LPA (Lp(a) / KIV-2 VNTR) analysis:
      KIV-2 total copies (both alleles): 39
      allele split: not resolvable (allele absent from reference library or subtype-homozygous)
      Qualitative Lp(a)/ASCVD association only (RUO); no quantitative Lp(a) value.
```

## Validation

All numbers are held-out against HPRC assembly truth (per-haplotype KIV-2 copy counts from the
phased assemblies). The caller runs on the real short-read CRAM; nothing about the sample's own
assembly is used at call time. Cohort: 183 HPRC samples with a resolved KIV-2 total; 181 of them
are also represented in the allele-split reference library.

### Total copy number (183-sample HPRC cohort)

The KIV-2 array averages ~38 copies, so absolute copy errors are best read against that scale:

| metric | value |
|---|---|
| mean copy count (truth) | 38.5 (range 16–57) |
| median \|error\| | 1 copy |
| mean \|error\| | 1.23 copies (3.2% of the mean count) |
| max \|error\| | 7 copies |
| within 2 copies of truth | 165/183 = 90.2% |
| within 4 copies of truth | 180/183 = 98.4% |

(Pearson r = 0.977 vs truth; the copy-error distribution above is the more actionable view — a
"small" 3% mean error is still ~1 copy.) The total count is the primary Lp(a)/ASCVD determinant
and is always reported.

### Allele split (leave-one-out, 181-sample library)

Each sample is scored against a library that excludes it (leave-one-out) and uses the caller's own
total count, so the numbers estimate the accuracy of the emitted `allele1`/`allele2` on a novel
patient. The split is reported only when the confidence gate passes; otherwise the caller reports
total-only. Errors below are for the minor allele (the smaller count, which drives Lp(a)):

| metric | value |
|---|---|
| split reported (confident) | 99/181 = 54% |
| mean minor-allele copy count (truth) | 16.4 (range 7–25) |
| median \|error\| | 1 copy |
| mean \|error\| | 1.45 copies |
| max \|error\| | 10 copies |
| within 2 copies of truth | 81/99 = 81.8% |
| within 4 copies of truth | 97/99 = 98.0% |
| confident-but-wrong (leak) | 2/181 = 1.1% |

The abstaining ~46% are subtype-homozygous or carry an allele absent from the reference library;
the gate reports total-only for them rather than emit a confident wrong split. 

### Cross-build reproducibility (hg38 / hg19 / b37)

The KIV-2 unit is colinear across builds: the hg38 unit sequence aligns to both the hg19 and b37
references as `5547M` (100% identity, no indels), and the marker frame is preserved base-for-base.
The unit-local reference library is shared across builds unchanged; only the `kiv2_unit` genomic
coordinate differs per build. hg19 collapses ~7 reference KIV-2 copies (hg38 has 6), which does not
affect the copy-number estimate. Cross-build caller runs on GIAB agree: HG002 gives an identical
total (30) and split (19 + 11) on hg38 and b37, and HG005 agrees within one copy (40 vs 41).

### Platform concordance (Ultima)

The caller was checked on Ultima Genomics reads (single-end, distinct error profile) by concordance
with Illumina on GIAB HG002–HG005 (Illumina 30×, Ultima 40×, hg38). The depth-ratio estimator is
self-normalizing, so the Illumina-calibrated `dp_norm` transfers unchanged despite Ultima's higher
absolute depth.

| sample | Illumina total | Ultima total | Illumina split | Ultima split |
|---|---|---|---|---|
| HG002 | 30 | 31 | 19 + 11 | 20 + 11 |
| HG003 | 40 | 39 | total-only | total-only |
| HG004 | 27 | 29 | 16 + 11 | 18 + 11 |
| HG005 | 40 | 38 | total-only | 37 + 1 |

The KIV-2 total agrees within 2 copies on all four samples with no systematic offset
(Δ = +1, −1, +2, −2) — no worse than the method's own Illumina-vs-truth spread. The allele split
agrees on the minor allele wherever both platforms split confidently (HG002 and HG004 both give
minor = 11) and both abstain together on HG003. HG005 is a boundary case: Illumina abstains
(total-only) while Ultima reports a confident 37 + 1, illustrating that the split — near-threshold
for some samples — is more platform-sensitive than the total, which is why it is gated and secondary.

### Scope and known limitations

- RUO / qualitative. Reports KIV-2 copy number, not a quantitative Lp(a) concentration.
- Allele split is not always resolvable. ~46% of samples are subtype-homozygous or carry an allele
  absent from the reference library; these correctly report total-only. This is a property of pooled
  short-read data, not a bug — the total count (the dominant risk determinant) is always reported.
- Residual novel-allele leaks (~1.1%). A small fraction of samples whose true alleles are absent
  from the library can match a wrong library pair confidently; a fuller catalog reduces this.
  Long-read input removes this failure mode (a HiFi read spans a full unit, so copies are genotyped
  and clustered directly, without a reference library).

## References

- Kraft HB, Utermann G, et al. Apolipoprotein(a) KIV-2 repeat and Lp(a) levels.
- Coassin S, Kronenberg F. The Lp(a) / KIV-2 copy-number determinant of cardiovascular risk.
- Illumina DRAGEN LPA caller (short-read KIV-2 copy number + linked-marker allele split).
