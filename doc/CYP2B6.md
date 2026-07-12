# CYP2B6 Gene Caller

> See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.

## Overview

CYP2B6 (cytochrome P450 2B6) is a clinically important pharmacogene on chromosome
19q13.2. It metabolizes efavirenz, nevirapine, bupropion, methadone,
cyclophosphamide, ketamine, artemisinin, and other drugs. The CYP2B6 locus sits
in the *CYP2ABFGST* gene cluster next to its highly similar non-functional
paralog CYP2B7P (~96% exonic identity), which causes the short-read misalignment
and gene-conversion events that confound standard variant callers. The clinically
important signal is the star diplotype and its metabolizer phenotype; whole-gene
copy-number change is rare, so the normal result is CN=2 for both paralogs.

## What This Module Does

The CYP2B6 module performs:

- Star-allele calling — translates the genotype into PharmVar/CPIC nomenclature
  (`*1`, `*4`, `*6`, `*9`, `*18`, …).
- Combination-allele resolution — the clinically dominant `*6` is defined by two
  variants in *cis* (`c.516G>T` + `c.785A>G`); an allele is assigned only when all of
  its variants are present.
- CPIC efavirenz phenotype — poor / intermediate / normal / rapid / ultrarapid
  metabolizer, with dosing implications.
- Copy-number calling — CYP2B6 and CYP2B7P copy number (normally 2 each), and the
  copy-number dosage that flags the structural alleles below.
- Structural allele detection — the rare `*29` (partial deletion) and `*30` (duplication),
  reciprocal hybrids formed by recombination in the intron 4 / exon 5 region. Because each
  changes the copy number of only part of the gene, they are read from copy-number dosage;
  a whole-gene CYP2B6 copy number of 2 does not exclude them.
- Small-variant calling — paralog-aware variant calls across the locus.

Supported on GRCh38 (hg38), GRCh37 (hg19), and b37; the build is auto-detected from the
input. Coordinates in this document are GRCh38, and calls are equivalent across builds
(validated identical on HG001–HG005).

The CYP2B6 module does not:

- Assign rare PharmVar alleles outside the core CPIC-actionable set (e.g. `*17` is not
  resolved from short reads).
- Guarantee `*29`/`*30` on real carrier data — copy-number detection is validated on
  constructed positive controls in both directions (deletion and duplication) and is
  specific across the reference cohort (no false calls), but no natural carrier was
  available to confirm sensitivity, so a `*29`/`*30` call should be confirmed orthogonally.

## Gene Cluster Structure (GRCh38)

| Gene | GRCh38 region | Strand | Description |
|------|---------------|--------|-------------|
| CYP2B6 | chr19:40,991,282–41,018,398 | + | Functional drug-metabolizing enzyme (9 exons, 491 aa) |
| CYP2B7P | chr19:40,924,265–40,950,660 | + | Non-functional pseudogene, ~41 kb upstream, head-to-tail, ~96% exonic identity |

Reference transcript NM_000767.5 (MANE Select); RefSeqGene NG_007929.1. The two
paralogs are ~96% identical across the gene body; short-read mismapping concentrates
in exon 5 (only ~2 nt divergent from CYP2B7P) and exon 9.

### Recombination / gene conversion

A high-homology stretch spanning intron 4 into exon 5 is the recombination-prone zone.
The reciprocal products are two rare hybrid alleles:

| Star Allele | Structure | Function | Frequency |
|-------------|-----------|----------|-----------|
| CYP2B6*29 | CYP2B7P–CYP2B6 partial hybrid (exons 1–4 from CYP2B7P, 5–9 from CYP2B6); a single collapsed copy | Decreased | ~0.5% (African-American) |
| CYP2B6*30 | Reciprocal CYP2B6/CYP2B7P hybrid duplication (extra copy carries a `p.R378X` stop) | Likely nonfunctional | ~0.25% (Asian) |

Unlike CYP2D6, CYP2B6 has no clean whole-gene deletion (no `*5`-style allele) and no
expanded hybrid family (`*13/*36/*68` analogs). Whole-gene copy number is normally CN=2.

Both hybrids delete or duplicate a shared segment between the two paralogs, so they show
up as a copy-number change over that segment even though each gene body still looks like
CN=2. The caller reads this dosage directly: a `*29` drops the affected copy number to 1
(or 0 if homozygous), a `*30` raises it to 3. A flagged allele is reported with the
affected copy number, the partner star allele, and the resulting metabolizer phenotype
(a `*29` is decreased-function, contributing toward a poor/intermediate result).

## Star Allele Calling

An allele is assigned only when all of its defining and associated variants are present,
and — when the reads phase them — lie on the same haplotype. More specific combination
alleles take precedence over their single-variant subsets: a haplotype carrying both
c.516 and c.785 is called `*6`, not `*9` or `*4`.

### Core alleles (CPIC function)

| Allele | Defining variant(s) | Protein | GRCh38 chr19 | Function |
|--------|---------------------|---------|--------------|----------|
| `*1` | reference | — | — | Normal |
| `*2` | rs8192709 (c.64C>T) | R22C | 40,991,369 | Normal |
| `*4` | rs2279343 (c.785A>G) | K262R | 41,009,358 | Increased |
| `*5` | rs3211371 (c.1459C>T) | R487C | 41,016,810 | Normal |
| `*6` | rs3745274 + rs2279343 | Q172H + K262R | 41,006,936 + 41,009,358 | Decreased |
| `*7` | `*6` + c.1459C>T | Q172H+K262R+R487C | +41,016,810 | Decreased |
| `*8` | rs12721655 (c.415A>G) | K139E | 41,004,377 | No function |
| `*9` | rs3745274 (c.516G>T) | Q172H | 41,006,936 | Decreased |
| `*13` | `*8` + `*6` | K139E+Q172H+K262R | 41,004,377 + … | No function |
| `*15` | rs35979566 (c.1172T>A) | I391N | 41,012,693 | Uncertain |
| `*17` | c.76/83/85/86 | T26S+D28G+R29T | 40,991,381–391 | Normal |
| `*18` | rs28399499 (c.983T>C) | I328T | 41,012,316 | No function |
| `*22` | rs34223104 (c.-82T>C) | promoter ↑expr | 40,991,224 | Uncertain |
| `*29` | CYP2B7P/CYP2B6 hybrid | — | structural | Decreased |
| `*34` | `*22` + `*7` | — | 40,991,224 + … | Decreased |
| `*36` | `*22` + `*6` | Q172H+K262R | 40,991,224 + … | Decreased |

> `*6` (the most common decreased-function allele) = `*9` (c.516) + `*4` (c.785).
> A haplotype with only c.516 is `*9`; only c.785 is `*4`. The caller distinguishes
> these by phase.

### Metabolizer phenotype (CPIC efavirenz guideline)

CPIC uses a diplotype → phenotype lookup (not an activity score):

| Phenotype | Diplotype rule | Efavirenz |
|-----------|----------------|-----------|
| Ultrarapid (UM) | two increased-function alleles | standard 600 mg/day |
| Rapid (RM) | normal + increased | standard 600 mg/day |
| Normal (NM) | two normal-function alleles | standard 600 mg/day |
| Intermediate (IM) | normal/increased + decreased/no-function | consider 400 mg/day |
| Poor (PM) | two decreased / two no-function / decreased+no-function | consider 400 or 200 mg/day |

Poor metabolizers have elevated efavirenz plasma exposure and higher risk of CNS
adverse effects and treatment discontinuation.

## Output Examples

Each `<sample>.yaml` carries a `CYP2B6` block with copy numbers and a human-readable
`cn_interpretation`; the structured `star_alleles` object additionally exposes
`alleles` (primary diplotype), `alt_alleles` (other configurations consistent with
unphased data), `phase_unresolved`, `ld_rescued` (positions recovered by the c.785
LD-rescue), and `note`. The examples below are reproducible on the named public sample.

### Normal metabolizer, `*1/*1` (Sample HG001)

```yaml
CYP2B6:
  Copy numbers:
    CYP2B6: 2
    CYP2B7P: 2
  cn_interpretation: |-
    CYP2B6: CN=2
    CYP2B7P: CN=2
    No structural/conversion events detected
    Star Alleles: *1 / *1
    Phenotype: Normal metabolizer
      Efavirenz: standard 600 mg/day
```

### Normal metabolizer, `*2/*5` (Sample HG002)

```yaml
CYP2B6:
  Copy numbers:
    CYP2B6: 2
    CYP2B7P: 2
  cn_interpretation: |-
    CYP2B6: CN=2
    CYP2B7P: CN=2
    No structural/conversion events detected
    Star Alleles: *2 / *5
    Phenotype: Normal metabolizer
      Efavirenz: standard 600 mg/day
```

### Intermediate metabolizer, phase unresolved `*6/*5` (Sample HG003)

Three heterozygous markers (c.516 + c.785 in exon 4/5 and c.1459 in exon 9) are too far
apart to phase from short reads, so the caller reports the most likely diplotype and the
alternative it cannot exclude. Both map to the same phenotype.

```yaml
CYP2B6:
  Copy numbers:
    CYP2B6: 2
    CYP2B7P: 2
  cn_interpretation: |-
    CYP2B6: CN=2
    CYP2B7P: CN=2
    Star Alleles: *6 / *5 (most likely; also possible: *7 / *1)
      Note: Phase unresolved from short reads - primary is the most likely diplotype;
      also consistent with: *7 / *1 (compound-heterozygote configuration cannot be excluded)
    Phenotype: Intermediate metabolizer (all configurations agree)
      Efavirenz: consider 400 mg/day
```

### Intermediate metabolizer, `*6/*1` with LD-rescued c.785 (Sample HG00126)

The `*6`-defining c.785 is diluted by CYP2B7P reads and would be missed on its own; it is
recovered because its *cis* partner c.516 is a confident call. The `ld_rescued` field flags
this, and the *trans* alternative (`*9/*4`) is reported since the two are linked rather than
read-phased.

```yaml
CYP2B6:
  Copy numbers:
    CYP2B6: 2
    CYP2B7P: 2
  cn_interpretation: |-
    CYP2B6: CN=2
    CYP2B7P: CN=2
    No structural/conversion events detected
    Star Alleles: *6 / *1 (most likely; also possible: *9 / *4)
      Note: Phase unresolved from short reads - primary is the most likely diplotype;
      also consistent with: *9 / *4 (compound-heterozygote configuration cannot be
      excluded); c.785 (rs2279343) recovered by 516 LD: real het diluted by CYP2B7P
      reads and ML-rejected, re-activated because c.516 is a confident call
    Phenotype: Intermediate metabolizer (all configurations agree)
      Efavirenz: consider 400 mg/day
```

## Validation

Ground truth for CYP2B6 comes from:

- HPRC phased diploid assemblies (226 samples) — the primary per-sample truth. Star
  diplotypes are read directly from each sample's phased assembly genotypes at the
  star-marker positions (phasing makes *cis*/*trans* unambiguous). Copy number is derived
  from how each haplotype's assembly aligns to the CYP2B6 vs CYP2B7P gene bodies, not from
  an assembly CNV track (which is unreliable at segmental duplications).
- Long-read / assembly-based star callers for GIAB HG001–HG005.
- GeT-RM consensus diplotypes (Gaedigk et al.) as an orthogonal cross-check where
  available. The assembly truth is independently corroborated by population genetics:
  across the cohort the derived allele frequencies (`*6` 33%, `*5` 6.5%, `*2` 5.1%,
  `*4` 3.8%, `*18` 2.5%) match published CYP2B6 frequencies.

GIAB reference samples (GRCh38, short-read):

| Sample | CYP2B6 CN | Star | Phenotype | vs truth |
|--------|-----------|------|-----------|----------|
| HG001 (NA12878) | 2 | *1/*1 | Normal | reference at all 8 core markers; CYP2B7P reads at c.785 correctly rejected, not miscalled *4 |
| HG002 (son) | 2 | *2/*5 | Normal | matches Q100 assembly (hap1 *2/c.64, hap2 *5/c.1459) |
| HG003 (father) | 2 | *6/*5 (alt *7/*1) | Intermediate | Mendelian truth *6/*5 |
| HG004 (mother) | 2 | *2/*5 | Normal | consistent |
| HG005 (NA24631) | 2 | *1/*1 | Normal | reference at all core markers; consistent |

HG001 and HG003 together show both directions of the paralog-aware value: HG001 rejects
CYP2B7P-derived reads at the c.785 hotspot (avoiding a false *4), while HG003 detects a
genuine *6 (c.516+c.785) at the same locus.

The markers used for calling were curated across the HPRC assembly cohort to remove
spurious gene-conversion artifacts — for example, on HG002 this reduced 7 spurious
conversion events to 0 while preserving the correct *2/*5 call. Markers outside the
high-homology exon 4/5 block (such as *2/c.64 and *5/c.1459) lie in uniquely-mappable
sequence and are called correctly regardless.

### HPRC cohort concordance (226 samples, GRCh38, short-read)

Short-read calls were compared against the phased-assembly truth for the full 226-sample
HPRC cohort:

| Metric | Concordance |
|--------|-------------|
| Star diplotype | ~93% (recovers exactly the *6-vs-*9 cases the gene targets) |
| Copy number (CYP2B6 & CYP2B7P) | ~99% — CYP2B6's own copy number correct in every sample; the single cohort-wide CN miss is one CYP2B7P pseudogene duplication (HG02622, true 2/3) |

The dominant short-read failure mode is at c.785 (rs2279343): it is not a
paralog-differentiating site (CYP2B7P carries the reference base there) and sits in a
marker-sparse window, so CYP2B7P reads can dilute a genuine heterozygous *6 call below the
variant caller's threshold, dropping *6 to *9. Because the contamination can only add
reference reads (never fabricate the c.785 alt), the caller LD-rescues a diluted c.785
when its *cis* partner c.516 is a confident call, which is what *6 requires. This recovers
the lost *6 diplotypes with no false gains (+5 diplotypes, 0 regressions on a clean
before/after comparison, lifting star concordance from ~88% to ~93%).

> Coverage caveat. Every sample in this cohort has CYP2B6 at CN=2 (whole-gene CYP2B6 CNV
> is very rare), so the ~99% CN figure is essentially a specificity result — the caller
> does not over-call copy-number changes. `*29`/`*30` sensitivity could not be measured on
> real data because no carriers were present (consistent with the rarity noted in the
> hybrid-allele table above). Their copy-number detection is instead validated on constructed
> positive controls — a simulated `*29` heterozygote and homozygote and a `*30` duplication
> are each recovered with the correct dosage and star call — and a copy-number scan across
> the reference cohort raised no false `*29`/`*30` (clean specificity). The one real CN
> event in the cohort (HG02622's CYP2B7P pseudogene duplication) is a CYP2B7P gain, not a
> CYP2B6 structural allele, and is still missed — see below.

Reproducible per-case examples. Each sample below is a public HPRC/1000-Genomes ID;
running the caller on its GRCh38 short-read alignment reproduces the listed call, which
matches the phased-assembly truth:

| Case | Example samples | CN | Star (= truth = call) | Phenotype |
|------|-----------------|----|-----------------------|-----------|
| Reference | HG00097, HG00099, HG00133 | 2/2 | `*1/*1` | Normal |
| `*6` het — LD-rescued at c.785 | HG00126, HG00350, HG02293 | 2/2 | `*6/*1` | Intermediate |
| `*6` het (directly phased) | HG00140, HG00235 | 2/2 | `*6/*1` | Intermediate |
| `*6` homozygous | HG00642, HG00738, HG01167 | 2/2 | `*6/*6` | Poor |
| `*6/*5` — c.785 LD-rescued | HG01175, HG01261 | 2/2 | `*6/*5` | Intermediate |
| `*5` het | HG00272, HG00280, HG00321 | 2/2 | `*5/*1` | Normal |
| `*2` het | HG00673, HG01099, HG01243 | 2/2 | `*2/*1` | Normal |
| `*4` het | HG00128, HG00408, HG00609 | 2/2 | `*4/*1` | Rapid |
| `*18` het | HG01150, HG02922, HG03050 | 2/2 | `*18/*1` | Intermediate |
| `*2` homozygous | HG01081 | 2/2 | `*2/*2` | Normal |

When the c.785 call is LD-rescued the caller flags it (`ld_rescued` in the star output)
and, since c.516 and c.785 are inferred *cis* by linkage rather than read-phased, reports
the *trans* alternative too (e.g. `*6/*1` primary with `*9/*4` also possible) — both map
to the same Intermediate phenotype.

Known-hard cases in the cohort (short-read limits, not specific to this caller):

| Sample | Truth | Call | Why |
|--------|-------|------|-----|
| HG00658, HG01960 | `*4/*1`, `*2/*4` | misses `*4` | `*4` is c.785 alone — no c.516 anchor, so the LD-rescue cannot fire (a lone diluted c.785 is indistinguishable from contamination) |
| HG01993 | `*6/*1` | `*6/*6` | one haplotype carries a deletion spanning exon 4/5; the other haplotype's `*6` reads read as homozygous — unresolvable from short reads |
| HG02622 | `*6/*2`, CN 2/3 | `*7/*2`, CN 2/2 | a CYP2B7P duplication; the extra pseudogene copy is not detected and perturbs phasing |
| HG02145 | `*17/*6` | `*1/*6` | `*17` (a 4-SNV cluster at c.76–86) not resolved from short reads |

Long-read input resolves the phase-unresolved and diluted-c.785 cases outright.

## Known Limitations

- The star-allele database covers the core CPIC-actionable alleles; rare PharmVar alleles
  not in the table receive no designation (e.g. `*17`, see HG02145 above).
- `*29`/`*30` hybrids are rare and are detected from copy-number dosage, corroborated by
  the gene-conversion signal when present. Detection is validated on constructed positive
  controls (both deletion and duplication) and is specific across the reference cohort, but
  no natural carrier was available, so real-data sensitivity is unconfirmed and a call
  should be confirmed orthogonally. A `*30` duplication is reported, but its metabolizer
  impact is left uncertain because its functional consequence is not established.
- c.785 (rs2279343) at low allele fraction. The `*6`-defining c.785 is diluted by CYP2B7P
  reads and is LD-rescued only when its c.516 partner is a confident call. A c.785 that
  appears alone (i.e. `*4`) has no such anchor, so it can be missed when its alt-read
  fraction is low. In the two observed cases (HG00658, HG01960) the locus is actually
  paralog-clean — the flanking PSV shows no CYP2B7P leak — and the deficit is
  reference/mapping bias against the alt allele (VAF ~0.22-0.30), leaving no dilution
  or cis-phasing signal available to rescue it. Long-read input resolves it directly.
- CYP2B7P copy-number gains (extra pseudogene copies) are not detected (HG02622); CYP2B6's
  own copy number is called reliably.
- Compound-heterozygote vs. combination-allele ambiguity: when a combination allele's
  defining variants are heterozygous but unphased, the caller reports the most likely
  diplotype as primary and lists the alternative(s) in `alt_alleles` (e.g. c.516+c.785+c.1459
  unphased → *6/*5 primary, *7/*1 alternative). The "most likely" ranking uses approximate
  population allele frequencies as a tiebreak only, not a clinical frequency estimate.
  Supply long reads to resolve the phase outright.

## Resources

- PharmVar CYP2B6: https://www.pharmvar.org/gene/CYP2B6
- CPIC efavirenz guideline: Desta et al., *Clin Pharmacol Ther* 2019;106:726–733 (PMID 31006110)
- PharmVar GeneFocus: CYP2B6: Desta et al., *Clin Pharmacol Ther* 2021;110:82–97 (PMC8693800)
