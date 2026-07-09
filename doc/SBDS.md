# SBDS (SBDS / SBDSP1) — Shwachman-Diamond syndrome locus

## Overview

Shwachman-Diamond syndrome (SDS) is an autosomal-recessive ribosomopathy — exocrine
pancreatic insufficiency, bone-marrow failure (neutropenia) with a predisposition to
myelodysplastic syndrome / acute myeloid leukemia, and skeletal abnormalities. About 90% of
cases are caused by biallelic variants in SBDS on chromosome 7q11.21.

SBDS has a nearly identical pseudogene, SBDSP1, ~5.8 Mb distal on the same arm. The two
share ~97% nucleotide identity, and the great majority of pathogenic SBDS alleles arise not
by point mutation but by gene conversion — a block of pseudogene sequence copied into the
functional gene (SBDSP1 → SBDS). Because SBDSP1 carries sequence that is inactivating in the
SBDS context, more than 90% of SDS patients carry one of a small number of recurrent exon-2
conversion variants.

Because SBDSP1 is so similar, standard short-read pipelines routinely mismap the two
paralogs and miscall the conversion variants. SBDS is a canonical case for paralog-aware
genotyping, which is what this caller provides for SDS diagnosis and carrier screening.

## What This Module Does

The SBDS module performs:

- Copy-number calling: SBDS and SBDSP1 copy number.
- Conversion-variant typing: types the two recurrent exon-2 pathogenic variants (c.258+2T>C
  and c.183_184delinsCT) by zygosity and summarizes a recessive-disease call — two pathogenic
  alleles → consistent with Shwachman-Diamond syndrome (bi-allelic), one allele → carrier,
  zero → negative.
- Small variant calling: paralog-differentiated small variants in SBDS and SBDSP1.
- Conversion-tract reporting: reports SBDSP1 → SBDS conversion tracts as supporting context,
  flagged non-diagnostic unless a direct variant call corroborates them — so a tract on its
  own never produces a carrier or affected call.
- Validated on 226 HPRC samples (copy number and conversion specificity) and against curated
  GIAB truth (small variants), reproducible across GRCh38/GRCh37 and Illumina/Ultima.

The SBDS module does not:

- Type SDS pathogenic alleles beyond the two recurrent exon-2 variants into a named disease
  call. Rarer conversion variants and private point mutations still appear in the standard
  small-variant calls, but are not summarized into the carrier/affected call.

## Genomic structure (GRCh38)

| Gene | GRCh38 region | Strand | Role |
|---|---|---|---|
| SBDS   | chr7:66,987,679–66,995,693 | − | functional gene (recipient) |
| SBDSP1 | chr7:72,829,656–72,836,701 | + | pseudogene (donor) |

The two paralogs are in opposite (inverted) orientation and ~5.8 Mb apart. Unlike the
tandem SMN1/SMN2 pair, this distal, inverted arrangement does not form simple fusion genes;
pathogenic exchange happens by internal gene conversion.

## Pathogenic conversion variants

Two recurrent pseudogene-derived variants in exon 2 account for the large majority of SDS
disease alleles:

| Variant (HGVS) | Effect | GRCh38 position(s) |
|---|---|---|
| c.258+2T>C | intron-2 splice-donor disruption | chr7:66,994,210 A>G |
| c.183_184delinsCT (p.Lys62\*) | nonsense / frameshift (K62X) | chr7:66,994,286 T>A + 66,994,287 A>G |

Both are the SBDSP1 (pseudogene) allele copied into SBDS. Because SDS is autosomal recessive,
disease requires two pathogenic alleles: homozygous for one variant, or compound
heterozygous across the two (the most common genotype, since a single conversion tract
spanning both sites is frequent). A single pathogenic allele is a carrier.

Copy-number changes at SBDS are rare — it is a conserved, dosage-sensitive gene — so the
clinically actionable signal is the conversion genotype, not copy number.

## Conversion tracts and the non-diagnostic flag

Alongside the directly typed variants, the caller reports ratio-based conversion tracts. On
some healthy samples a divergent SBDS haplotype can make a tract appear where there is no
real conversion — about 27% of individuals show at least one such tract. To keep the disease
call reliable, pathogenic status is set only by the directly typed variant genotypes: a
tract that matches a known pathogenic event is labeled pathogenic only when a direct call
corroborates it, and is otherwise flagged non-diagnostic (see the HG005 example below). A
tract on its own never produces a carrier or affected call.

## Output Examples

All examples are actual caller output. Copy number and a SUMMARY line are always reported;
the recurrent variants and any conversion tracts appear when present.

### Normal, no conversion (GIAB HG002; also HG001, HG003–HG005, and every HPRC sample tested)

```yaml
SBDS:
  Copy numbers:
    SBDS: 2
    SBDSP1: 2
  cn_interpretation: |-
    SBDS (Shwachman-Diamond syndrome) Analysis:
      SBDS:   CN=2
      SBDSP1: CN=2

    SUMMARY: no recurrent SBDS conversion variant detected.
      NOTE: only the two recurrent exon-2 PSVs (c.258+2T>C, c.183_184delinsCT) are typed here; other SBDS pathogenic variants are reported via the base variant calls.
```

### Carrier, one pathogenic allele (simulated positive control)

```yaml
SBDS:
  Copy numbers:
    SBDS: 2
    SBDSP1: 2
  cn_interpretation: |-
    SBDS (Shwachman-Diamond syndrome) Analysis:
      SBDS:   CN=2
      SBDSP1: CN=2

    Recurrent SBDSP1->SBDS conversion variants detected:
      mono-allelic SBDS_c.258+2T>C (chr7:66994210 A>G; pseudogene-derived)

    SUMMARY: 1 pathogenic SBDS conversion allele detected — carrier (heterozygous).
      NOTE: only the two recurrent exon-2 PSVs (c.258+2T>C, c.183_184delinsCT) are typed here; other SBDS pathogenic variants are reported via the base variant calls.
```

### Affected, compound heterozygous (simulated positive control)

Both recurrent variants are present, one allele each, and the conversion tract that spans
them is corroborated by the direct calls:

```yaml
SBDS:
  Copy numbers:
    SBDS: 2
    SBDSP1: 2
  cn_interpretation: |-
    SBDS (Shwachman-Diamond syndrome) Analysis:
      SBDS:   CN=2
      SBDSP1: CN=2

    Recurrent SBDSP1->SBDS conversion variants detected:
      mono-allelic SBDS_c.258+2T>C (chr7:66994210 A>G; pseudogene-derived)
      mono-allelic SBDS_c.183_184delinsCT (chr7:66994286 T>A; pseudogene-derived)

    Positive-direction conversion tract(s): 1 — ratio-based, non-diagnostic alone (can reflect allelic dropout of a divergent haplotype; pathogenic status is set by the direct PSV calls above)
      Tract 1: chr7:66993897-66994328, converted_alleles=1
        matches SBDS_c.258+2T>C (score 100.0%) — CORROBORATED by direct PSV call
        matches SBDS_c.183_184delinsCT (score 100.0%) — CORROBORATED by direct PSV call

    SUMMARY: 2 pathogenic SBDS conversion alleles detected — consistent with Shwachman-Diamond syndrome (bi-allelic). Confirm; SDS is autosomal recessive.
      NOTE: only the two recurrent exon-2 PSVs (c.258+2T>C, c.183_184delinsCT) are typed here; other SBDS pathogenic variants are reported via the base variant calls.
```

### Non-diagnostic ratio tract on a healthy sample (GIAB HG005)

HG005 has a divergent SBDS haplotype that produces ratio tracts matching the pathogenic
events, but no direct variant call supports them, so the tracts are flagged non-diagnostic
and the SUMMARY stays negative:

```yaml
SBDS:
  Copy numbers:
    SBDS: 2
    SBDSP1: 2
  cn_interpretation: |-
    SBDS (Shwachman-Diamond syndrome) Analysis:
      SBDS:   CN=2
      SBDSP1: CN=2

    Positive-direction conversion tract(s): 3 — ratio-based, non-diagnostic alone (can reflect allelic dropout of a divergent haplotype; pathogenic status is set by the direct PSV calls above)
      Tract 1: chr7:66989236-66989756, converted_alleles=1
      Tract 2: chr7:66991539-66991741, converted_alleles=1
      Tract 3: chr7:66993959-66994286, converted_alleles=1
        matches SBDS_c.258+2T>C (score 100.0%) — NOT corroborated by a direct PSV call (non-diagnostic; likely allelic-dropout artifact)
        matches SBDS_c.183_184delinsCT (score 100.0%) — NOT corroborated by a direct PSV call (non-diagnostic; likely allelic-dropout artifact)

    SUMMARY: no recurrent SBDS conversion variant detected.
      NOTE: only the two recurrent exon-2 PSVs (c.258+2T>C, c.183_184delinsCT) are typed here; other SBDS pathogenic variants are reported via the base variant calls.
```

## Validation

### Copy number and conversion specificity (HPRC assembly-truth cohort)

Validated on the 226-sample HPRC assembly cohort (GRCh38), with copy-number and conversion
truth taken from the phased assemblies. 224 samples have usable truth; copy number is scored
on the 216 with an unambiguous phase at the locus (the other 8 are excluded from the CN
count), and all 224 are used for the conversion check.

| metric | result |
|---|---|
| SBDS copy number (exact) | 216 / 216 (100%) |
| SBDSP1 copy number (exact) | 215 / 216 (99.5%) |
| Pathogenic-conversion specificity | 1.000 (0 / 224 false positive) |

No sample in this unaffected cohort carries a pathogenic conversion, so the cohort tests
specificity — and it is perfect, including the 60 / 224 (27%) of samples that show the
non-diagnostic ratio-tract artifact: not one becomes a false carrier or affected call. The
single copy-number discordance (HG03195) is a soft one-copy SBDSP1 (pseudogene) deletion
reported as 2 rather than 1 — clinically inconsequential, since SBDSP1 dosage carries no SDS
meaning (an identical deletion in HG03834 was caught).

### Small-variant accuracy (GIAB, curated truth)

Benchmarked against curated GIAB truth within each gene's high-confidence regions — HG002
against v5.0q (which covers the segmental-duplication core that the older v4.2.1 benchmark
excludes), and HG001/HG003/HG004/HG005 against NIST v4.2.1.

| Gene | Samples | SNV F1 | indel F1 | overall F1 |
|---|---|---|---|---|
| SBDS | HG001–HG005 | 1.000 | 1.000 | 1.000 |
| SBDSP1 | HG001–HG005 | 1.000 | 1.000 | 1.000 |

Perfect (0 false positive / 0 false negative) for both paralogs across all five samples, SNV
and indel — 34 SBDS + 22 SBDSP1 truth variants in total, every one called correctly.

### Sensitivity (simulated positive control)

No SDS-positive samples were available, so the pathogenic genotypes were simulated by
injecting the recurrent variant alleles into GIAB reads at the SBDS locus. The caller typed
each correctly:

| simulated genotype | caller call |
|---|---|
| carrier (1 × c.258+2T>C) | 1 pathogenic allele — carrier |
| affected, homozygous (2 × c.258+2T>C) | 2 alleles — Shwachman-Diamond syndrome (bi-allelic) |
| affected, compound-het (c.258+2T>C / c.183_184delinsCT) | 2 alleles — SDS (bi-allelic), phased in trans |

### Cross-build reproducibility (hg38 / hg19 / b37)

The SBDS gene is identical between GRCh38 and GRCh37, and SBDSP1 differs at only six
positions, none at a pathogenic site — so the GRCh37 (hg19 and b37) support is a clean
coordinate lift of the GRCh38 model. Runs on HG001–HG005 reproduce the GRCh38 result:
identical copy number, identical conversion status, and no genotype disagreements at shared
sites.

### Cross-platform (Ultima vs Illumina)

On HG002, HG003 and HG004, Ultima and Illumina short reads give identical copy number and
conversion status, and both reach F1 = 1.000 on small variants (SBDS and SBDSP1, SNV and
indel). The indel-recall gap Ultima shows at some homopolymer-rich loci does not affect SBDS.

### Validation scope and known gaps

Validated (hg38 Illumina unless noted):

- SBDS / SBDSP1 copy number — 100% / 99.5% exact on the HPRC cohort.
- The two recurrent pathogenic variants (carrier / affected typing) — specificity 1.000
  across 224 unaffected samples; sensitivity confirmed on simulated positives.
- Small variants (SNV/indel) — F1 = 1.000 vs curated GIAB truth (HG001–HG005).
- Cross-build (hg38 / hg19 / b37) and cross-platform (Ultima) reproducibility.

Not yet validated (deliberately out of v1 scope):

- SDS pathogenic alleles beyond the two recurrent exon-2 variants — surfaced via the
  paralog-aware small-variant calls, but not typed into a named disease call.
- Conversion-tract sensitivity — the tract detector is tuned for specificity and gated
  non-diagnostic; it is not benchmarked against conversion-tract truth. Pathogenic calls rely
  on the direct variant genotypes.
- Real SDS patient material — sensitivity is established by simulation, not by an affected
  clinical sample.

## References

- Boocock GRB et al. *Mutations in SBDS are associated with Shwachman-Diamond syndrome.*
  Nat Genet 2003. (SBDS gene, gene conversion, the two recurrent exon-2 variants.)
- Nelson AS, Myers KC. *Shwachman-Diamond syndrome.* (clinical review.)
