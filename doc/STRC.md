# STRC (STRC / STRCP1) — DFNB16 autosomal-recessive hearing loss locus

## Overview

DFNB16 is an autosomal-recessive, nonsyndromic sensorineural hearing loss caused by biallelic
loss of STRC, which encodes stereocilin — a protein of the stereocilia of the cochlear outer
hair cells and of the horizontal top connectors that couple them. STRC is a leading cause of
recessive deafness: after GJB2 (DFNB1), STRC (DFNB16) is the **second most common** genetic
cause of autosomal-recessive hearing loss, and the great majority of its pathogenic alleles are
a recurrent copy-number deletion rather than point mutations.

STRC sits on chromosome 15q15.3 immediately next to a nearly identical partial pseudogene,
STRCP1, about 100 kb distal. The two are a **direct tandem** segmental duplication at ~99.6%
coding identity — and, critically, their **exons 1–15 are 100% identical**, with
paralog-distinguishing sites (PSVs) only in the 3′ half (exons 16–29). That arrangement produces
two paralog-driven disease mechanisms that standard short-read pipelines handle poorly:

- a recurrent ~100 kb deletion of the functional STRC copy, formed by non-allelic homologous
  recombination (NAHR) between the STRC and STRCP1 repeats; and
- a (much less common) STRCP1 → STRC gene-conversion allele.

There is also a locus-specific twist. **CATSPER2**, a uniquely-mappable gene required for male
fertility, lies *between* the STRC and STRCP1 repeat copies. In **77–88%** of the recurrent
deletions the deleted segment extends far enough to co-remove CATSPER2, producing
**Deafness-Infertility Syndrome (OMIM #611102)** — the same hearing loss plus male infertility.
Reporting whether a deletion is STRC-only or STRC+CATSPER2 is therefore clinically meaningful for
reproductive counseling.

Because STRCP1 is so similar, ordinary pipelines both **miss the deletion** (STRCP1 reads mismap
onto STRC and refill the lost depth) and **cannot resolve the 100%-identical exons 1–15** (no
paralog-distinguishing sites). This module is paralog-aware and addresses the deletion, its
CATSPER2 extent, and the benign look-alikes.

## What This Module Does

The STRC module performs:

- Copy-number calling: STRC, STRCP1, and CATSPER2 copy number. STRC and STRCP1 are separated by
  the 3′ paralog-distinguishing sites; CATSPER2 — being unique — is measured directly from
  read depth and is a clean, unconfounded probe.
- Recurrent-deletion typing: reports STRC copy loss and summarizes a recessive-disease call —
  two pathogenic alleles → consistent with DFNB16 (bi-allelic), one → carrier, zero → negative.
- CATSPER2-extent reporting: for a deletion, states whether CATSPER2 is co-deleted
  (Deafness-Infertility Syndrome, OMIM #611102) or intact (STRC-limited hearing loss).
- Benign discrimination: a copy lost on the STRCP1 (pseudogene) side is flagged **benign for
  hearing** (functional STRC intact); benign 3-copy gains of STRC or STRCP1 (~1.85% of the
  population) are not called pathogenic.
- Small-variant calling: paralog-differentiated small variants across STRC and STRCP1.

The STRC module does not:

- Type STRC pathogenic alleles beyond the recurrent deletion into a named disease call. Private
  point mutations still appear in the standard small-variant calls but are not summarized into
  the carrier/affected call.
- Resolve variants inside exons 1–15, where STRC and STRCP1 are ~100% identical and no
  paralog-distinguishing sites exist — a low-confidence region flagged for orthogonal
  confirmation.
- Type STRCP1 → STRC conversion into a pathogenic call in this release. STRC conversion is
  proven but uncommon and has no canonical allele set, so the module carries a conversion
  detector (a corroboration-gated pseudogene-PSV check) but ships it **inert** in v1; the
  clinically actionable call rests on the deletion.

## Genomic structure (GRCh38)

| Gene | GRCh38 region | Strand | Role |
|---|---|---|---|
| STRC     | chr15:43,599,003–43,620,000 | − | functional gene (recipient) |
| STRCP1   | chr15:43,698,817–43,719,466 | − | pseudogene (donor) |
| CATSPER2 | chr15:43,628,503–43,668,118 | − | unique gene between the repeats (fertility) |

STRC and STRCP1 are in the same (direct) orientation, a tandem repeat ~100 kb apart at ~99.6%
coding identity. **Exons 1–15 are 100% identical** between the two; only the 3′ half (exons
16–29) carries paralog-distinguishing sites, so STRC is resolved from those 3′ PSVs and is blind
in exons 1–15. Because the repeats are direct, they recombine to a **deletion** (not an
inversion), so the event is copy-number-detectable. CATSPER2 sits in the unique sequence between
the two repeat copies and is read cleanly from its own depth.

## Pathogenic mechanisms

**Recurrent deletion (primary).** NAHR between the STRC and STRCP1 direct repeats deletes a
contiguous ~100 kb segment that removes the functional STRC gene. Biallelic deletion causes
DFNB16; a single deletion is a carrier. This copy-number change is the dominant STRC pathogenic
allele class and the main clinically actionable signal.

**CATSPER2 co-deletion (deletion extent).** In 77–88% of the recurrent deletions the deleted
segment extends across the intergenic gap and also removes CATSPER2. Because CATSPER2 loss causes
male infertility, a homozygous STRC+CATSPER2 deletion produces **Deafness-Infertility Syndrome
(OMIM #611102)** — hearing loss plus male infertility — while an STRC-only deletion gives
isolated DFNB16 hearing loss. The module reports which of the two a deletion is.

**Gene conversion (STRCP1 → STRC).** STRCP1 carries inactivating sequence (a pseudogene stop
near exon 20, plus 3′ PSVs) that conversion can copy into STRC; this is pathogenic but uncommon,
and unlike OTOA's single p.Glu787\* there is no canonical STRC conversion allele. The module can
flag a directly-typed, tract-corroborated pseudogene PSV, but this path is disabled in v1.

## Output Examples

All examples below are actual caller output. Copy number and a SUMMARY line are always reported;
the deletion line and its CATSPER2-extent line appear when a deletion is present. STRCP1 copy
number is a benign pseudogene dosage and is not clinically actionable on the hearing-loss axis.

### Normal, no pathogenic allele (GIAB HG002)

```yaml
STRC:
  Copy numbers:
    STRC: 2
    STRCP1: 2
    CATSPER2: 2
  cn_interpretation: |-
    STRC (DFNB16 autosomal-recessive hearing loss) Analysis:
      STRC:     CN=2
      STRCP1:   CN=2
      CATSPER2: CN=2

    SUMMARY: no pathogenic STRC deletion or conversion allele detected.
```

### Deletion carrier, one pathogenic allele — with CATSPER2 co-deletion (1000G NA19240, a real carrier)

```yaml
STRC:
  Copy numbers:
    STRC: 1
    STRCP1: 2
    CATSPER2: 1
  cn_interpretation: |-
    STRC (DFNB16 autosomal-recessive hearing loss) Analysis:
      STRC:     CN=1
      STRCP1:   CN=2
      CATSPER2: CN=1

    STRC copy loss: 1 allele(s) deleted — recurrent DFNB16 ~100 kb NAHR deletion of the functional STRC copy.
      Extent: co-deletes CATSPER2 (CN=1) — Deafness-Infertility Syndrome (OMIM #611102): the contiguous deletion also removes CATSPER2, adding male infertility to the hearing loss.

    SUMMARY: 1 pathogenic STRC allele (deletion) — carrier (heterozygous).
```

### Benign STRCP1 deletion (1000G HG00146, a real sample)

A copy lost on the pseudogene side is **not** DFNB16 — the functional STRC copy is intact. Here
the deletion also happens to remove CATSPER2, which is flagged for its independent
male-infertility consequence:

```yaml
STRC:
  Copy numbers:
    STRC: 2
    STRCP1: 1
    CATSPER2: 1
  cn_interpretation: |-
    STRC (DFNB16 autosomal-recessive hearing loss) Analysis:
      STRC:     CN=2
      STRCP1:   CN=1
      CATSPER2: CN=1

    STRCP1 copy loss: 1 allele(s) — benign pseudogene deletion; NOT DFNB16 (functional STRC intact, CN=2). NB CATSPER2 is co-deleted (CN=1): CATSPER2 loss causes male infertility (OMIM #611102) independently of the hearing-loss status.

    SUMMARY: no pathogenic STRC deletion or conversion allele detected.
```

### Affected, homozygous deletion

No homozygous-deletion (affected) individual occurs in the healthy population cohort. A biallelic
STRC deletion (STRC CN=0) yields the deterministic recessive-disease summary:

```
SUMMARY: 2 pathogenic STRC allele(s) (deletion=2, conversion=0) — consistent with DFNB16 hearing loss (bi-allelic). Confirm; DFNB16 is autosomal recessive.
```

with the CATSPER2-extent line reporting Deafness-Infertility Syndrome when CATSPER2 is also
biallelically deleted.

## Validation

### Copy number (population cohort)

Validated on **226 population samples** (HPRC / 1000 Genomes assembly cohort, GRCh38) against
**paralog-resolved assembly truth** — per-haplotype deletion calls read from the assembly
alignments (not from a naïve joint copy-number VCF, which mis-calls zygosity inside the segdup),
with the uniquely-mappable CATSPER2 depth as the zygosity arbiter.

| metric | result |
|---|---|
| STRC recurrent-deletion detection (sensitivity) | 5 / 5 real carriers, 0 missed |
| STRC deletion specificity | 1 borderline false positive (HG03874) across the cohort |
| CATSPER2 co-deletion extent typed on the 5 carriers | 5 / 5 correct |

Across the cohort STRC copy number is 2 in 202 samples, 1 in 6, and 3–4 in 18 (benign gains).
The six STRC-copy-loss calls comprise five deletion carriers confirmed against assembly truth
(all correctly STRC=1 / STRCP1=2, each co-deleting CATSPER2 → the deafness-infertility extent)
and one borderline false positive (HG03874, at the depth threshold; assembly-normal). No real
carrier is missed. STRCP1-side losses and 3-copy gains are benign and flagged as such (STRCP1 is
a polymorphic pseudogene and its dosage carries no hearing-loss meaning). A gentle
depth-normalization (`dp_norm` = 1.08) corrects the ~11% depth inflation of the paralogous
STRC/STRCP1 block while sparing the unique CATSPER2 probe; it removes two benign-gain false
positives and introduces zero false deletions.

### Small-variant accuracy (GIAB, curated truth)

Benchmarked against NIST GIAB v4.2.1 within the gene ∩ high-confidence regions, HG001–HG005,
population bundle. STRC is entirely paralogous (no single-copy anchor exon), so small-variant
recall in exons 1–15 is inherently capped; CATSPER2, being unique, is called cleanly.

| region | result |
|---|---|
| STRC (functional gene) overall F1 | ~0.89 mean (precision ~0.95) |
| CATSPER2 (unique) F1 | 0.978 (HG002) / 1.000 (HG005) |

The residual STRC gap is intrinsic to the locus — exons 1–15 are 100% identical to STRC P1 and
carry no discriminating sites — whereas the 3′-PSV region, the CATSPER2 probe, and the
copy-number/deletion calls are solid.

### Cross-build reproducibility (hg38 / hg19 / b37)

The STRC/STRCP1/CATSPER2 reference sequence is byte-identical between GRCh38 and GRCh37, offset
by a constant **+292,198 bp**, so the GRCh37 (hg19 and b37) support is an exact coordinate lift
of the GRCh38 model (STRC → chr15:43,891,201–43,912,198; CATSPER2 →
chr15:43,920,701–43,960,316). On a matched 30× HG002 pair the GRCh37 build reproduces the GRCh38
result exactly: identical copy number (STRC 2 / STRCP1 2 / CATSPER2 2), identical negative
summary, and the same small-variant calls.

### Cross-platform (Ultima vs Illumina)

On GIAB HG002/HG003/HG004, Ultima and Illumina short reads give **identical** STRC / STRCP1 /
CATSPER2 copy number and identical deletion status (negative on both, all three samples) — the
clinically actionable copy-number outputs are platform-invariant. For small variants, the unique
CATSPER2 region is strong on both platforms (F1 ≈ 0.93–0.99, comparable either way), while the
paralogous STRC region shows Ultima's expected modest recall dip (about one additional missed
call per sample, concentrated in the hard 100%-identical exons) with **no false positives on
either platform**.

### Validation scope and known gaps

Validated (GRCh38 Illumina unless noted):

- STRC copy number and recurrent-deletion detection — 5/5 real carriers on the 226-sample
  cohort, 0 missed, ~99% non-gain concordance vs paralog-resolved assembly truth.
- CATSPER2 co-deletion extent (Deafness-Infertility Syndrome vs STRC-only) — from the unique
  CATSPER2 depth probe.
- Benign discrimination of STRCP1 deletions and of 3-copy gains.
- Small variants — STRC full-gene F1 ~0.89, CATSPER2 F1 ~0.98–1.00 vs curated GIAB truth.
- Cross-build (hg38 / hg19 / b37) and cross-platform (Ultima) reproducibility.

Not yet validated / deliberately out of v1 scope:

- STRC pathogenic alleles beyond the recurrent deletion — surfaced via the paralog-aware
  small-variant calls, but not typed into a named disease call.
- Exons 1–15 (STRC/STRCP1 100% identical) — no paralog-distinguishing sites; variants there are
  low-confidence and flagged for orthogonal confirmation.
- STRCP1 → STRC gene conversion — detector present but shipped inert (no canonical STRC
  conversion allele); not benchmarked against conversion truth.
- Real DFNB16 affected material — deletion sensitivity is established on real cohort carriers
  (heterozygous); the homozygous-affected summary is the deterministic recessive combination.

## References

- Verpy E, et al. *Mutations in a new gene encoding a protein of the hair bundle cause
  non-syndromic deafness at the DFNB16 locus.* Nat Genet 2001 (PMID 11524702).
- Vona B, et al. *DFNB16 / STRC deletions as a frequent cause of mild-to-moderate autosomal
  recessive hearing loss.* (review of STRC deletion frequency; second after GJB2).
- Zhang L, et al. *Short-read resolution of the STRC / STRCP1 paralog and the recurrent
  deletion.* Clin Chem 2023 (DOI 10.1093/clinchem/hvad046).
- Avidan N, et al. *CATSPER2, contiguous STRC–CATSPER2 deletion, and Deafness-Infertility
  Syndrome (OMIM #611102).* (deletion-infertility syndrome; PMC6416315, PMC8755823).
