# CYP2A6 Gene Caller

> See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.

## Overview

CYP2A6 (cytochrome P450 2A6) is a pharmacogene on chromosome 19q13.2. It is the
principal enzyme of nicotine C-oxidation (nicotine → cotinine, ~70–80% of nicotine
clearance) and also metabolizes coumarin, letrozole, tegafur, valproate, and others.
The CYP2A6 locus sits in the *CYP2ABFGST* gene cluster next to its highly similar,
non-functional paralog CYP2A7 (~95% identity), a head-to-tail tandem pair on the minus
strand — the same architecture as CYP2D6/CYP2D7, and the source of the short-read
misalignment and gene-conversion events that confound standard callers.

Three signals matter for CYP2A6: (1) the star diplotype and its nicotine-metabolism
rate; (2) the common whole-gene deletion `*4` (East Asian allele frequency ~15–20%,
so real deletion sensitivity is measurable in cohorts); and (3) the `*1×2` duplication.
Unlike CYP2D6 and CYP2B6, **CYP2A6 has no CPIC drug-dosing guideline** — the module
reports a research nicotine-metabolism grouping (activity-score based), not a
prescribing recommendation.

## What This Module Does

The CYP2A6 module performs:

- Star-allele calling — translates the genotype into PharmVar nomenclature (`*1`, `*2`,
  `*4`, `*7`, `*9`, `*10`, `*17`, `*35`, …).
- Copy-number calling — CYP2A6 and CYP2A7 copy number (normally 2 each), on an
  asymmetric (gamma) prior because both a real duplication (`*1×2`) and a real deletion
  (`*4`) occur.
- Structural allele detection — the whole-gene deletion `*4` (drops CYP2A6 below CN 2),
  the `*1×2` duplication (raises it above 2), and the `*12` / `*34` CYP2A7::CYP2A6 5′
  hybrids (read by fusion detection).
- Gene-conversion detection — the common benign `*1B`/`*46`-type 3′-flank
  CYP2A7↔CYP2A6 conversion, reported as a normal-function event so it is not mistaken
  for a pathogenic change.
- Combination-allele resolution — require-all-*cis* matching (`*10` = `*7` + `*8`;
  `*13` = `*9` + a second promoter variant).
- Small-variant calling — paralog-aware variant calls across the locus, using the
  population short-read model.
- Nicotine-metabolism phenotype — an activity-score grouping (normal / intermediate /
  slow / poor metabolizer), reported as research context, with no dosing note.

Supported on GRCh38 (hg38), GRCh37 (hg19), and b37; the build is auto-detected from the
input. Coordinates in this document are GRCh38, and calls are equivalent across builds.

The CYP2A6 module does not:

- Provide CPIC dosing guidance — no CYP2A6 guideline exists; the phenotype grouping is a
  nicotine-metabolism research annotation only.
- Assign rare PharmVar alleles outside the implemented core set.
- Resolve every 3′ terminal-exon paralog-cluster SNV from short reads (`*35`/`*10`/`*9` are
  the residual miss — see Validation).

## Gene Cluster Structure (GRCh38)

| Gene | GRCh38 region | Strand | Description |
|------|---------------|--------|-------------|
| CYP2A6 | chr19:40,841,915–40,850,493 | − | Functional nicotine C-oxidase (9 exons, 494 aa) |
| CYP2A7 | chr19:40,873,525–40,882,752 | − | Non-functional paralog, ~23 kb upstream (5′ side), ~95% identity |

Canonical transcript ENST00000301141 (RefSeq NM_000762). The ~95% paralog identity
concentrates short-read mismapping in the exon-9 / 3′ region, where several of the
reduced-function star markers sit.

### Structural and conversion events

| Star | Structure | Function | How detected |
|------|-----------|----------|--------------|
| `*4` | Whole-gene deletion — the crossover removes the exon body but retains the 3′-flank | None | Copy-number dosage (see the longdel note) |
| `*1×2` | Whole-gene duplication | Increased | Copy-number dosage (gamma prior) |
| `*12` | CYP2A7::CYP2A6 5′ hybrid (intron-2 junction) | Decreased | Fusion detection |
| `*34` | CYP2A7::CYP2A6 5′ hybrid (intron-4 junction) | Decreased | Fusion detection |
| `*1B` / `*46` | Benign 3′-flank CYP2A7↔CYP2A6 gene conversion | Normal | Conversion detection |

> **The `*4` deletion and the longdel model.** The `*4` crossover deletes the CYP2A6
> exon body (chr19:40,843,800–40,850,493) but **retains the 3′-flank**
> (40,841,651–40,843,795). A single whole-gene copy-number average would blend the
> deleted body (CN 0) with the retained flank (CN 2) and read a spurious CN 1. The caller
> instead carries a `*4` exon-body deletion degree-of-freedom, so a homozygous `*4` floors
> the functional exon copy number to 0 while the retained flank is reported separately
> (`whole-locus CN`). This is why a `*4/*4` sample prints `CYP2A6: CN=0` with a
> `whole-locus CN=1` annotation.

## Star Allele Calling

An allele is assigned only when all of its defining variants are present and — when the
reads phase them — lie on the same haplotype. More specific combination alleles take
precedence over their single-variant subsets (`*10` = `*7` + `*8`).

### Core alleles

| Allele | Defining variant(s) | Protein | GRCh38 chr19 | Function |
|--------|---------------------|---------|--------------|----------|
| `*1` | reference | — | — | Normal |
| `*2` | rs1801272 | L160H | 40,848,628 A>T | No function |
| `*4` | whole-gene deletion | — | structural | No function |
| `*5` | rs5031017 | G479V | 40,843,845 C>A | No function |
| `*7` | rs5031016 | I471T | 40,843,869 A>G | Decreased |
| `*8` | rs28399468 | R485L | 40,843,827 C>A | Normal |
| `*9` | rs28399433 | −48T>G (TATA box) | 40,850,474 A>C | Decreased |
| `*10` | `*7` + `*8` | I471T + R485L | 40,843,869 + 40,843,827 | Decreased |
| `*12` | CYP2A7::CYP2A6 5′ hybrid | — | structural | Decreased |
| `*13` | `*9` + rs28399434 | −48T>G + promoter | 40,850,474 + 40,850,414 | Decreased |
| `*17` | rs28399454 | V365M | 40,845,362 C>T | Decreased |
| `*20` | rs568811809 | frameshift (2-bp del) | 40,848,284 CTT>C | No function |
| `*34` | CYP2A7::CYP2A6 5′ hybrid | — | structural | Decreased |
| `*35` | rs143731390 | N438Y | 40,843,969 T>A | Decreased |
| `*1×2` | whole-gene duplication | — | structural | Increased |

> `*7` and `*8` share the exon-9 3′ region; `*10` is their *cis* combination. `*9` and
> `*13` share the −48T>G promoter (TATA-box) variant. The caller distinguishes the
> combinations by requiring all *cis* variants of the more-specific allele.

### Nicotine-metabolism phenotype (research grouping — no CPIC guideline)

The grouping sums an activity score over the two alleles (normal = 1, decreased = 0.5,
none = 0):

| Group | Activity score | Example diplotypes |
|-------|----------------|--------------------|
| Normal metabolizer | ~2.0 | `*1/*1` |
| Intermediate metabolizer | ~1.0–1.5 | `*9/*1`, `*1/*4` |
| Slow metabolizer | ~0.5 | `*9/*4` |
| Poor metabolizer | 0 | `*4/*4` |

CYP2A6 activity is associated with nicotine-metabolism rate, smoking behavior and
cessation, and coumarin/letrozole/tegafur handling, but **no CPIC dosing guideline
exists**, so the module attaches no prescribing note — only the research grouping.

## Output Examples

Each `<sample>.yaml` carries a `CYP2A6` block with copy numbers and a human-readable
`cn_interpretation`; the structured `star_alleles` object exposes `alleles` (primary
diplotype), `alt_alleles`, and `copy_number`. The examples below are reproducible on the
named public sample.

### Normal metabolizer, `*1/*1` (Sample HG001)

```yaml
CYP2A6:
  Copy numbers:
    CYP2A6: 2
    CYP2A7: 2
  cn_interpretation: |-
    CYP2A6: CN=2
    CYP2A7: CN=2
    No structural/conversion events detected
    Star Alleles: *1 / *1
    Nicotine-metabolism phenotype: Normal metabolizer
      Note: research nicotine-metabolism grouping (activity score); CYP2A6 has no CPIC dosing guideline
```

### Intermediate metabolizer, `*9/*1` with a benign 3′-flank conversion (Sample HG005)

The `*9` (−48T>G) promoter allele is a reduced-expression variant. HG005 also carries the
common benign `*1B`/`*46`-type 3′-flank conversion, reported as a normal-function event.

```yaml
CYP2A6:
  Copy numbers:
    CYP2A6: 2
    CYP2A7: 2
  cn_interpretation: |-
    CYP2A6: CN=2
    CYP2A7: CN=2
    Structural/conversion events detected: 1
      Event 1: 3'-flank conversion (CYP2A7<->CYP2A6; *1B/*46-type, benign/normal function)
        Region: chr19:40846605-40847165
        Allele change: +1
    Star Alleles: *9 / *1
    Nicotine-metabolism phenotype: Intermediate metabolizer
      Note: research nicotine-metabolism grouping (activity score); CYP2A6 has no CPIC dosing guideline
```

### Slow metabolizer, `*9/*4` het whole-gene deletion (Sample NA18945)

One haplotype carries the `*4` deletion; the exon-body copy number drops to 1 while the
3′-flank is retained (`whole-locus CN=2`).

```yaml
CYP2A6:
  Copy numbers:
    CYP2A6: 1
    CYP2A6_4del: -1
    CYP2A7: 2
  cn_interpretation: |-
    CYP2A6: CN=1
      (*4 exon-body deletion -1; 3'-flank retained, whole-locus CN=2)
    CYP2A7: CN=2
    Star Alleles: *9 / *4
      Note: *4 whole-gene deletion x1 (CYP2A6 CN=1)
    Nicotine-metabolism phenotype: Slow metabolizer
      Note: research nicotine-metabolism grouping (activity score); CYP2A6 has no CPIC dosing guideline
```

### Poor metabolizer, `*4/*4` homozygous deletion (Sample HG02523)

```yaml
CYP2A6:
  Copy numbers:
    CYP2A6: 0
    CYP2A6_4del: -1
    CYP2A7: 2
  cn_interpretation: |-
    CYP2A6: CN=0
      (*4 exon-body deletion -1; 3'-flank retained, whole-locus CN=1)
    CYP2A7: CN=2
    Star Alleles: *4 / *4
    Nicotine-metabolism phenotype: Poor metabolizer
      Note: research nicotine-metabolism grouping (activity score); CYP2A6 has no CPIC dosing guideline
```

## Validation

Ground truth for CYP2A6 comes from:

- HPRC phased diploid assemblies (223 samples with usable truth) — the primary per-sample
  truth. Each sample's two assembled haplotypes are aligned against the CYP2A6 and CYP2A7
  gene sequences (minimap2, secondary alignments retained); the hits are clustered into
  distinct genomic loci, and each locus is assigned to whichever paralog it matches best. A
  locus is counted as one clean gene copy when its alignment identity is ≥ 0.92 over at
  least half the gene body, and copy number is the number of clean copies summed over the
  two haplotypes. Because copies are counted directly rather than inferred from a single
  spanning alignment, this **represents `*1×2` duplications as CN 3 and `*4` whole-gene
  deletions as CN < 2** — a `*4` crossover drops the exon-body identity to ~0.83, so the
  fused copy is not counted as clean. An assembly CNV track is not used (it is unreliable at
  segmental duplications).
- GIAB HG001–HG005 (short-read), with the Ashkenazi trio (HG002/3/4) providing a
  Mendelian cross-check.
- PharmVar allele definitions. The assembly truth is corroborated by population genetics:
  the cohort-derived allele frequencies (`*9` 13.1%, `*4` 6.5%, `*17` 3.8%, `*35` 2.0%)
  match published CYP2A6 spectra for a multi-ethnic cohort.

### GIAB reference samples (GRCh38, short-read)

| Sample | CYP2A6 CN | Star | Phenotype | Note |
|--------|-----------|------|-----------|------|
| HG001 (NA12878) | 2 | `*1/*1` | Normal | reference at all markers |
| HG002 (son) | 2 | `*1/*1` | Normal | trio-consistent; a PSV gate rejects an Illumina exon-body depth dip that would otherwise mis-call `*4` (see below) |
| HG003 (father) | 2 | `*1/*1` | Normal | |
| HG004 (mother) | 2 | `*1/*1` | Normal | trio Mendelian-consistent |
| HG005 (NA24631) | 2 | `*9/*1` | Intermediate | `*9` (−48T>G) heterozygote; assembly-confirmed CN 2 |

### HPRC cohort concordance (223 samples, GRCh38, short-read, population bundle)

Copy number and star diplotype are compared per sample against the assembly truth described
above. Because that truth counts gene copies directly, duplications and deletions are both
represented, and no sample needs to be excluded as an assembly-truth artifact. One sample
(HG03492) has no usable assembly alignment and is not scored.

| Metric | Concordance | Comment |
|--------|-------------|---------|
| Copy number (CYP2A6 & CYP2A7) | 223/223 (100%) | every CN state matches, including three CYP2A6 `*1×2` duplications (CN 3: HG01891/HG02055/HG02451), one CYP2A7 duplication (HG02080), and the 23 heterozygous + 3 homozygous `*4` deletions. |
| `*4` whole-gene-deletion sensitivity | 26/26 (100%) | every sample with truth CYP2A6 CN < 2 is called as a deletion. |
| Star diplotype (primary-exact) | 212/223 (95.1%) | the 11 discordances are single-SNV differences in the 3′ terminal-exon marker cluster (`*35` N438Y, `*7`/`*10` I471T±R485L, promoter `*9`), where CYP2A6 and CYP2A7 are ~95% identical and a heterozygous alt can be diluted or mimicked by paralog reads. Misses (7) and over-calls (6) are balanced — no directional bias. |

**The PSV-consistency gate.** A `*4` deletion call is vetoed when the CYP2A6/CYP2A7
paralog-sequence-variant (PSV) allele fraction over the exon body is diploid-level
(≥ 0.42 over ≥ 5 sites) — a more reliable copy-ratio than mappability-fragile unique depth.
This rejects an Illumina-only false `*4` on HG002 (a coverage dip that Ultima and the
diploid assembly both contradict) while changing none of the 28 cohort samples called
CYP2A6 CN < 2: it fixes HG002 without touching a single true deletion.

### Cross-platform (Ultima vs Illumina)

HG002/HG003/HG004/HG005 were each called on both Illumina and Ultima short reads on GRCh38.
Post-gate, all eight runs are platform-concordant and truth-correct (HG002 `*1/*1`, HG003
`*1/*1`, HG004 `*1/*1`, HG005 `*9/*1`). The one platform-specific effect — an Illumina HG002
exon-body depth dip — is handled by the PSV gate above.

### Small-variant accuracy (GIAB reference samples)

Base-level accuracy was measured with `aardvark` against GIAB truth, gated to the CYP2A6
gene ∩ each sample's high-confidence region, using the population short-read model.
**Recall is 1.000 on every arm — no truth variant is missed;** the only residual is a few
Illumina false positives (CYP2A7 paralog leakage), which Ultima does not show.

| Platform / build | Truth | SNV recall | SNV F1 | Note |
|---|---|---|---|---|
| Illumina, GRCh38 (HG001–005) | v4.2.1 | 1.000 | 0.90–1.00 (mean ≈0.955) | 0 FN; per-sample FPs 2/2/0/1/0 |
| Ultima, GRCh38 (HG001–005) | v4.2.1 | 1.000 | 1.000 (×5) | 0 FP, 0 FN |
| Illumina, GRCh38 (HG002) | v5.0q\* | 0.917 | 0.880 SNV; 0.78 with indels | stricter/complete truth |

\* v5.0q (HG002 only) is a stricter, more complete truth including exon-region indels; the
ALL-variant F1 dips to 0.78 on 4 indels — the same short-read segmental-duplication indel
limit seen on CYP2B6/STRC. On the standard v4.2.1 benchmark, SNV calling is essentially
perfect (Illumina mean F1 ≈0.955, Ultima a clean 1.000 across all five samples), with the
`*4`/`*9`/`*17` star markers all in the called set.

## Known Limitations

- **No CPIC guideline.** CYP2A6 has no CPIC drug-dosing guideline; the metabolizer
  grouping is a research nicotine-metabolism annotation, not a prescribing recommendation.
- **3′ terminal-exon paralog-cluster SNVs.** The reduced-function markers `*35` (N438Y),
  `*7`/`*10` (I471T ± R485L), and promoter `*9` sit in the ~95%-identity 3′ region, where a
  genuine heterozygous alt can be diluted — or mimicked — by CYP2A7 reads; these account for
  all 11 single-SNV star discordances in the cohort (7 misses, 6 over-calls — balanced). A
  population model bundle (recommended) mitigates the dilution, as on the other
  pharmacogenes. In the same region, Illumina short reads also produce a small number of
  CYP2A7-leakage false positives at the base level (0–2 per GIAB sample; recall stays
  1.000), which Ultima does not show.
- **`*12`/`*34` hybrids.** These structural alleles are rare and no natural carrier was
  available in the cohort to confirm sensitivity; a call should be confirmed orthogonally.
  (`*1×2` duplications, by contrast, are confirmed: three carriers in the cohort are called
  at CN 3, matching the assembly truth.)

## Resources

- PharmVar CYP2A6: https://www.pharmvar.org/gene/CYP2A6
- PharmGKB CYP2A6: https://www.pharmgkb.org/gene/PA121
- CYP2A6 is the principal nicotine-metabolizing enzyme; the nicotine-metabolism ratio
  (3′-hydroxycotinine / cotinine) is a genotype-associated biomarker of CYP2A6 activity
  used in smoking-cessation research. CYP2A6 has no CPIC prescribing guideline.
