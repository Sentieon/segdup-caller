# IKBKG Gene Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## What This Module Does

The IKBKG module performs:

- **Copy number calling**: Determines copy numbers for IKBKG flanking regions, IKBKGP1 flanking regions, and the shared 11.7kb exon 4-10 segment (consolidated total). Sex-specific baselines (male CN=1, female CN=2) are applied automatically. Validated at >99% accuracy overall (100% for males) against HPRC v2.0 assembly truth.
- **Gene conversion detection in the 5' region (Block 1)**: Uses 42 high-quality PSV markers in the 5' flank of IKBKG to detect pathogenic (IKBKGP1-to-IKBKG) and benign (IKBKG-to-IKBKGP1) gene conversion events via HMM-based segmentation of allele ratios.

The IKBKG module does **not** perform:

- **CNV gene attribution in the 11.7kb IKBKGdel region**: The 11.7kb segment (exons 4-10) is fully homogenized between IKBKG and IKBKGP1 with no gene-differentiating sites. CNV events (deletions or duplications) in this region cannot be attributed to gene vs pseudogene from short reads. These events are reported as consolidated totals with an `[UNRESOLVED]` flag.
- **Small variant calling**: The liftover region outside Block 1 is heavily homogenized between IKBKG and IKBKGP1. Gene-differentiating sites are absent in these regions, making variant assignment to gene vs pseudogene unreliable from phased liftover BAMs. Variant calling is therefore skipped for this gene.

## Overview

IKBKG (Inhibitor of Kappa light polypeptide gene enhancer in B-cells, Kinase Gamma), also known as NEMO (NF-kB Essential Modulator), is an X-linked gene on Xq28 encoding a regulatory subunit of the IKK complex essential for NF-kB signaling. Its pseudogene IKBKGP1 is an **inverted** duplicate located ~76kb downstream, sharing an 11.7kb segmental duplication that contains exons 4-10.

**Unique features of IKBKG analysis**:
- **Pervasive gene conversion**: Exons 4-10 and surrounding regions are heavily homogenized between IKBKG and IKBKGP1, making standard variant calling unreliable
- **Sex-specific ploidy**: Males are hemizygous (baseline CN=1), females are diploid (baseline CN=2)
- **Bidirectional CNV detection**: The 11.7kb segment can be deleted OR duplicated
- **Mosaic detection**: Living males with pathogenic IKBKG events must be mosaic
- **Gene conversion detection**: Identifies pathogenic IKBKGP1-to-IKBKG conversions using 42 population-trained markers in the 5' region

## IKBKG Gene Cluster Structure

### Genomic Organization

| Gene | Genomic Coordinates | Description |
|------|---------------------|-------------|
| **IKBKG** | chrX:154,555,612-154,575,000 | Functional NF-kB essential modulator (10 exons) |
| **IKBKGP1** | chrX:154,632,000-154,648,856 | Non-functional inverted pseudogene ~76kb downstream |

### Inverted Pseudogene Architecture

Unlike most gene/pseudogene pairs in segmental duplications (e.g., GBA1/GBAP1, PMS2/PMS2CL) where the two copies are arranged in tandem, IKBKGP1 is in **reverse complement orientation** relative to IKBKG. This inverted arrangement has two important consequences:

1. **NAHR between inverted repeats produces inversions**, not simple deletions or fusions. The common exon 4-10 deletion is mediated by this mechanism.
2. **Gene conversion still occurs regardless of orientation**. Non-reciprocal sequence transfer between the paralogs does not require same-strand alignment.

### Segmental Duplication and Copy Number Regions

The 11.7kb segmental duplication encompasses exons 4-10 of IKBKG and the corresponding inverted region in IKBKGP1. Each gene is modeled as three sub-regions:

```
chrX (Xq28)
                                         ~76kb gap
IKBKG (+ strand)                     |--------------|          IKBKGP1 (inverted)
 ---------+-----------+----------     --------+-----------+----------
| 5' flank | 11.7kb segment | 3' flank |    | 5' flank | 11.7kb segment | 3' flank |
| 2.4kb    | (exons 4-10)   | 2.7kb    |    | 2.7kb    | (exons 4-10)   | 2.4kb    |
 ---------+-----------+----------     --------+-----------+----------
|42 PSVs   | NO PSVs (0 sites) |            |            | NO PSVs (0 sites) |          |
|diff ~100%| fully homogenized |            |            | fully homogenized |          |
 ---------+-----------+----------     --------+-----------+----------
 154555612  154558020   154569680 154572430  154632000  154634748   154646416  154648856

           -------- gene -> -------->              <-------- pseudo (inverted) --------

 distinguishable    indistinguishable                       indistinguishable
  (unique seq)      (identical in both)                     (identical in both)
```

The tool reports three copy number values by collapsing the internal sub-regions:
- **IKBKG**: Copy number of the IKBKG flanking regions (distinguishable from pseudogene)
- **IKBKGP1**: Copy number of the IKBKGP1 flanking regions (distinguishable from gene)
- **IKBKGdel_region**: Total copy number of the 11.7kb IKBKGdel segment across both IKBKG and IKBKGP1 (indistinguishable by short reads)

> **Why a single consolidated segment CN?** The 11.7kb segment contains no gene-differentiating variants due to complete sequence homogenization. Short reads cannot resolve whether a segment-level CNV affects IKBKG or IKBKGP1. The tool reports only the total segment CN to avoid implying a gene-vs-pseudogene distinction that the data cannot support. Flank regions ARE distinguishable and reported separately.

### Sex-Specific Ploidy

IKBKG is located in the non-PAR region of the X chromosome. Copy number baselines differ by sex:

| Sex | IKBKG Baseline | IKBKGP1 Baseline | X Copies |
|-----|---------------|------------------|----------|
| **Male** | CN=1 | CN=1 | 1 |
| **Female** | CN=2 | CN=2 | 2 |

Sex information is required for IKBKG analysis (provided via `--sex` or inferred from the BAM file).

## Sequence Homogenization Between IKBKG and IKBKGP1

### Why This Locus Is Challenging

The IKBKG/IKBKGP1 locus presents a unique analytical challenge: **prevalent gene conversion has homogenized large portions of the two paralog sequences**, making it difficult to distinguish the functional gene from the pseudogene using standard approaches.

Analysis of ~460 Human Pangenome Reference Consortium (HPRC) haplotypes reveals the following pattern:

- **5' region (Block 1)**: 42 SNP positions reliably differentiate IKBKG from IKBKGP1 across ~100% of samples. This region is well-conserved and provides strong gene-vs-pseudogene distinction. Gene conversion detection is performed using these markers.
- **Exons 4-10 region (11.7kb segment)**: The entire 11.7kb segmental duplication contains **zero reliable differentiating positions** due to complete sequence homogenization. This is the fundamental reason short reads cannot determine whether a segment-level CNV affects IKBKG or IKBKGP1.
- **3' region**: The 3' flank is largely homogenized across the population, with only a handful of positions retaining weak differentiation in a minority of individuals. This region does not provide sufficient signal for reliable gene conversion detection from short reads.

This means that across the population, the majority of the IKBKG/IKBKGP1 locus appears nearly identical at the sequence level. Standard approaches that rely on fixed reference differences between gene and pseudogene fail at this locus because those "differences" have been erased by recurrent gene conversion in most people.

### PSV Discovery Process

Paralog-Specific Variant (PSV) sites were identified by analyzing ~460 HPRC long-read assembly haplotypes. Each haplotype was aligned to the reference and the gene and pseudogene tracks were compared at every position across the full 16.8kb IKBKG region to find sites where the two paralogs carry different alleles.

Unlike genes such as PMS2 where most differentiating sites are fixed across the population, the pervasive gene conversion at IKBKG means a **paired per-sample approach** is required — comparing gene vs. pseudogene alleles within each individual rather than relying on population-level consensus. The discovery process also used a two-pass training strategy: an initial pass to discover candidate blocks, followed by strict relabeling against the best anchor block to eliminate false positives.

The result is **42 SNP-only PSV markers** in a single high-quality linkage block in the 5' region:

| Block | Region | SNP Sites | Differentiation | Description |
|-------|--------|-----------|-----------------|-------------|
| **Block 1 (5')** | chrX:154,555,619-154,555,882 | 42 | ~100% of samples | Reliable anchor for CN calling and gene conversion detection |

## Pathogenic Events

### 11.7kb Exon 4-10 Deletion

The most common pathogenic event is deletion of the 11.7kb segment containing exons 4-10. This deletion removes most of the IKBKG coding region, abolishing NF-kB signaling.

**Segment attribution is unresolved**: Because the 11.7kb segment is fully homogenized between IKBKG and IKBKGP1, short reads cannot determine which gene carries the deletion. The tool reports segment-level CNVs with an `[UNRESOLVED]` flag and recommends long-read confirmation for definitive gene assignment. The mechanism is inherently NAHR between the inverted segmental duplications.

### Whole-Gene Deletion

When the IKBKG flanking regions show reduced CN **and** the total segment CN is reduced, a whole-gene deletion is reported. Because flank CN IS distinguishable between gene and pseudogene, whole-gene deletions can be confidently attributed to IKBKG. This represents a larger chromosomal event beyond the segmental duplication and is more severe than segment-only deletion.

### Gene Conversion (IKBKGP1 to IKBKG)

Pseudogene sequence from IKBKGP1 can overwrite functional IKBKG sequence through gene conversion, introducing pathogenic variants. The tool detects conversion events in the 5' region using Block 1's 42 PSV markers:

- **5' conversion** (IKBKGP1 to IKBKG): Pseudogene sequence replacing the 5' region of IKBKG, involving up to 42 signature variant positions

**Direction determines pathogenicity**:

| Direction | Allele Change | Clinical Significance |
|-----------|---------------|----------------------|
| **IKBKGP1 to IKBKG** | Positive (+1, +2) | **Pathogenic**: pseudogene variants disrupt functional gene |
| **IKBKG to IKBKGP1** | Negative (-1, -2) | **Benign**: functional sequence replacing pseudogene |

### 11.7kb Segment Duplication

Duplication of the 11.7kb segment is detected but is usually benign regardless of which gene carries the extra copy. Like segment deletions, duplications are reported with an `[UNRESOLVED]` flag since gene attribution cannot be determined from short reads.

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Normal Female

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 2
    IKBKGP1: 2
    IKBKGdel_region (IKBKG+IKBKGP1 total): 4
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (female):
      IKBKG flanks: CN=2
      IKBKGP1 flanks: CN=2
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=4 (expected 4)

    No structural events detected.
```

### Example 2: Normal Male

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 1
    IKBKGP1: 1
    IKBKGdel_region (IKBKG+IKBKGP1 total): 2
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (male):
      IKBKG flanks: CN=1
      IKBKGP1 flanks: CN=1
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=2 (expected 2)

    No structural events detected.
```

### Example 3: Female Heterozygous 11.7kb Deletion

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 2
    IKBKGP1: 2
    IKBKGdel_region (IKBKG+IKBKGP1 total): 3
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (female):
      IKBKG flanks: CN=2
      IKBKGP1 flanks: CN=2
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=3 (expected 4)

    [UNRESOLVED] 11.7kb segment deletion detected (total CN=3, expected=4)
      The 11.7kb segment is identical between IKBKG and IKBKGP1.
      Short reads cannot resolve which gene carries this event.
      Mechanism: NAHR between segmental duplications.
      If IKBKG: pathogenic (exon 4-10 loss). If IKBKGP1: benign.
      Recommend long-read confirmation for definitive assignment.
```

### Example 4: Female Whole-Gene IKBKG Deletion

When flanks are also reduced, gene attribution IS possible:

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 1
    IKBKGP1: 2
    IKBKGdel_region (IKBKG+IKBKGP1 total): 3
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (female):
      IKBKG flanks: CN=1
      IKBKGP1 flanks: CN=2
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=3 (expected 4)

    IKBKG whole-gene deletion detected (flank CN=1, segment total CN=3)
      Flank CN reduction confirms IKBKG involvement.
      Clinical: Heterozygous whole-gene IKBKG deletion —
      Incontinentia Pigmenti carrier.
```

### Example 5: Male with Mosaic Conversion

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 1
    IKBKGP1: 1
    IKBKGdel_region (IKBKG+IKBKGP1 total): 2
  mosaic:
    event: IKBKG_conversion
    type: conversion
    fraction: 0.35
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (male):
      IKBKG flanks: CN=1
      IKBKGP1 flanks: CN=1
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=2 (expected 2)

    Mosaic Event Detected: IKBKG_conversion (conversion), fraction=35.0%
      Mosaicism explains viability in affected male —
      35% of cells carry the pathogenic event.
```

### Example 6: Gene Conversion (Pathogenic, Female)

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 2
    IKBKGP1: 2
    IKBKGdel_region (IKBKG+IKBKGP1 total): 4
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (female):
      IKBKG flanks: CN=2
      IKBKGP1 flanks: CN=2
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=4 (expected 4)

    Gene Conversion Events: 1
      Event 1: IKBKGP1->IKBKG (pathogenic direction), alleles=1
        Region: chrX:154555619-154555882
        Clinical: Pseudogene sequence in functional IKBKG —
        Incontinentia Pigmenti carrier.
  gene_conversions:
  - converted_alleles: 1
    region: chrX:154555619-154555882
```

## Interpretation Guide

### Copy Number Patterns

| IKBKG | IKBKGP1 | IKBKGdel_region total | Interpretation | Attribution | Clinical |
|-------------|---------------|--------------|----------------|-------------|----------|
| 2 | 2 | 4 | Normal female | — | Normal |
| 1 | 1 | 2 | Normal male | — | Normal |
| 2 | 2 | 3 | Het segment deletion | **[UNRESOLVED]** | If IKBKG: IP carrier. If IKBKGP1: benign. |
| 1 | 1 | 1 | Hemi segment deletion | **[UNRESOLVED]** | If IKBKG: non-viable unless mosaic. If IKBKGP1: benign. |
| 1 | 2 | 3 | Het whole-gene IKBKG deletion | **Resolved** (flanks) | IP carrier |
| 2 | 0 | 2 | Hom whole-gene IKBKGP1 deletion | **Resolved** (flanks) | Benign |
| 2 | 2 | 5 | Segment duplication | **[UNRESOLVED]** | Usually benign regardless of gene |

> **Attribution column**: "Resolved" means flank CN provides gene-level evidence. "[UNRESOLVED]" means only total segment CN is known — gene vs pseudogene assignment requires long-read confirmation.

### IKBKGP1 Events

IKBKGP1 is a non-functional pseudogene. Because the 11.7kb segment is reported as a consolidated total CN, segment-level events are not attributed to IKBKG or IKBKGP1 individually. IKBKGP1 flanking CN changes (if detected) are reported separately and have no direct clinical significance.

## Validation

The IKBKG caller was validated by comparing short-read results against long-read HPRC v2.0 assembly truth for 223 overlapping samples (112 male, 111 female). Assembly truth was derived from phased haplotype assemblies of 1000 Genomes samples.

### Copy Number Accuracy

| Metric | Overall | Male (n=112) | Female (n=111) |
|--------|---------|--------------|----------------|
| IKBKG flank CN | 221/223 (99.1%) | 112/112 (100%) | 109/111 (98.2%) |
| IKBKGP1 flank CN | 222/223 (99.6%) | 112/112 (100%) | 110/111 (99.1%) |
| IKBKGdel_region CN | 223/223 (100%) | 112/112 (100%) | 111/111 (100%) |

Male copy number is 100% accurate across all three metrics. The 3 female CN mismatches are explained by specific structural events:

| Sample | Gene | Assembly | Caller | Explanation |
|--------|------|----------|--------|-------------|
| HG01891 | IKBKG | 2 | 3 | Tandem duplication in shared segment creates skewed allele balance |
| NA19909 | IKBKG | 2 | 3 | Large-scale CNV gain (segment CN=8) |
| NA18565 | IKBKGP1 | 3 | 2 | Complete Block 1 conversion (42/42 PSVs) makes extra IKBKGP1 copy indistinguishable from IKBKG in short reads |

### Gene Conversion — Block 1 (5' region)

Block 1 gene conversion is rare. In the HPRC v2.0 dataset (223 samples), only 3 samples show Block 1 conversion in the assembly, all female, all in the benign direction (IKBKG to IKBKGP1):

| Sample | Detail |
|--------|--------|
| HG00133 | hap2: Conv_in_IKBKGP1 (38/42 PSVs) |
| HG01496 | mat: Conv_in_IKBKGP1 (38/42 PSVs) |
| NA18565 | hap2: Conv_in_IKBKGP1 (42/42 PSVs, complete conversion) |

No pathogenic (IKBKGP1 to IKBKG) Block 1 conversions were observed in the HPRC v2.0 dataset.

### Validation Summary

- Copy number accuracy is excellent (>99% overall, 100% for males)
- The rare CN mismatches are attributable to unusual structural events (tandem duplications, large CNV gains, complete Block 1 gene conversion)
- No pathogenic Block 1 gene conversions were observed in the 223-sample HPRC v2.0 validation set, so pathogenic conversion detection sensitivity cannot be assessed from this dataset

## Known Limitations

- **11.7kb segment gene attribution**: The most significant limitation. The 11.7kb segment is fully homogenized between IKBKG and IKBKGP1 with no differentiating variants. Short reads cannot determine whether a segment-level CNV (deletion or duplication) affects the functional gene or the pseudogene. These events are reported with an `[UNRESOLVED]` flag. Long-read sequencing or other orthogonal methods are required for definitive gene assignment. Flanking regions ARE distinguishable, so whole-gene deletions (where flanks are also affected) CAN be attributed to a specific gene.
- **No small variant calling**: Due to pervasive homogenization outside Block 1, variant assignment from phased liftover BAMs is unreliable. Small variant calling is skipped for this gene.
- **Inverted pseudogene**: IKBKGP1 is inverted relative to IKBKG, which precludes classical fusion detection. Structural events are detected via copy number changes rather than fusion breakpoints.
- **Balanced reciprocal exchanges**: Balanced sequence exchanges between IKBKG and IKBKGP1 that leave no net allele imbalance are undetectable by short-read sequencing.
- **Novel conversions**: Events outside the Block 1 gene conversion detection region will not be reported due to lack of reliable gene-differentiating sites.

## IKBKG-Specific Resources

### Clinical Databases
- **ClinVar IKBKG**: https://www.ncbi.nlm.nih.gov/clinvar/?term=IKBKG
- **OMIM Incontinentia Pigmenti**: https://www.omim.org/entry/308300
- **OMIM IKBKG**: https://www.omim.org/entry/300248

### Literature

1. Smahi A, Courtois G, Vabres P, et al. (2000) "Genomic rearrangement in NEMO impairs NF-kappaB activation and is a cause of incontinentia pigmenti." *Nature* 405:466-472. https://doi.org/10.1038/35013114
2. Fusco F, Paciolla M, Conte MI, et al. (2008) "Genomic architecture at the Incontinentia Pigmenti locus favours de novo pathological alleles through different mechanisms." *Human Molecular Genetics* 17(10):1521-1530. https://doi.org/10.1093/hmg/ddn040
3. Conte MI, Pescatore A, Paciolla M, et al. (2014) "Insight into IKBKG/NEMO locus: report of new mutations and complex genomic rearrangements leading to incontinentia pigmenti disease." *Human Mutation* 35(2):165-177. https://doi.org/10.1002/humu.22483
4. Liao WW, Asber M, Didion JP, et al. (2023) "A draft human pangenome reference" *Nature* 617:312-324. https://doi.org/10.1038/s41586-023-05896-x
