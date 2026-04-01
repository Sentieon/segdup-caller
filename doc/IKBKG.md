# IKBKG Gene Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## Overview

IKBKG (Inhibitor of Kappa light polypeptide gene enhancer in B-cells, Kinase Gamma), also known as NEMO (NF-kB Essential Modulator), is an X-linked gene on Xq28 encoding a regulatory subunit of the IKK complex essential for NF-kB signaling. Its pseudogene IKBKGP1 is an **inverted** duplicate located ~76kb downstream, sharing an 11.7kb segmental duplication that contains exons 4-10.

**Unique features of IKBKG analysis**:
- **Pervasive gene conversion**: Exons 4-10 and surrounding regions are heavily homogenized between IKBKG and IKBKGP1, making standard variant calling unreliable
- **Sex-specific ploidy**: Males are hemizygous (baseline CN=1), females are diploid (baseline CN=2)
- **Bidirectional CNV detection**: The 11.7kb segment can be deleted OR duplicated
- **Mosaic detection**: Living males with pathogenic IKBKG events must be mosaic

## What This Module Does

The IKBKG module performs:

- **Copy number calling**: Determines copy numbers for IKBKG flanking regions, IKBKGP1 flanking regions, and the shared 11.7kb exon 4-10 segment (consolidated total). Sex-specific baselines (male CN=1, female CN=2) are applied automatically. Validated at 99.6% accuracy overall (100% for males) against HPRC v2.0 assembly truth.

The IKBKG module does **not** perform:

- **CNV gene attribution in the 11.7kb IKBKGdel region**: The 11.7kb segment (exons 4-10) is fully homogenized between IKBKG and IKBKGP1 with no gene-differentiating sites. CNV events (deletions or duplications) in this region cannot be attributed to gene vs pseudogene from short reads. These events are reported as consolidated totals with an `[UNRESOLVED]` flag.
- **Small variant calling**: The liftover region outside Block 1 is heavily homogenized between IKBKG and IKBKGP1. Gene-differentiating sites are absent in these regions, making variant assignment to gene vs pseudogene unreliable from phased liftover BAMs. Variant calling is therefore skipped for this gene.
- **Gene conversion detection**: The Block 1 PSVs used for CN calling are the same sites that differentiate IKBKG from IKBKGP1 — a converted region would simply be classified as the other allele, making gene conversion detection circular. The 3' region lacks sufficient differentiating sites for reliable detection.

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

> **Why a single consolidated segment CN?** The 11.7kb segment contains no confident gene-differentiating variants due to complete sequence homogenization. Short reads cannot resolve whether a segment-level CNV affects IKBKG or IKBKGP1. The tool reports only the total segment CN to avoid implying a gene-vs-pseudogene distinction that the data cannot support. Flank regions ARE distinguishable and reported separately.

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

- **5' region (Block 1)**: 42 SNP positions reliably differentiate IKBKG from IKBKGP1 across ~100% of samples. This region is well-conserved and provides strong gene-vs-pseudogene distinction for copy number calling.
- **Exons 4-10 region (11.7kb segment)**: The entire 11.7kb segmental duplication contains **zero reliable differentiating positions** due to complete sequence homogenization. This is the fundamental reason short reads cannot determine whether a segment-level CNV affects IKBKG or IKBKGP1.
- **3' region**: The 3' flank is largely homogenized across the population, with only a handful of positions retaining weak differentiation in a minority of individuals. This region does not provide sufficient signal for reliable gene conversion detection from short reads.

This means that across the population, the majority of the IKBKG/IKBKGP1 locus appears nearly identical at the sequence level. Standard approaches that rely on fixed reference differences between gene and pseudogene fail at this locus because those "differences" have been erased by recurrent gene conversion in most people.

## Pathogenic Events

### 11.7kb Exon 4-10 Deletion

The most common pathogenic event is deletion of the 11.7kb segment containing exons 4-10. This deletion removes most of the IKBKG coding region, abolishing NF-kB signaling.

**Segment attribution is unresolved**: Because the 11.7kb segment is fully homogenized between IKBKG and IKBKGP1, short reads cannot determine which gene carries the deletion. The tool reports segment-level CNVs with an `[UNRESOLVED]` flag and recommends long-read confirmation for definitive gene assignment. The mechanism is inherently NAHR between the inverted segmental duplications.

### Whole-Gene Deletion

When the IKBKG flanking regions show reduced CN **and** the total segment CN is reduced, a whole-gene deletion is reported. Because flank CN IS distinguishable between gene and pseudogene, whole-gene deletions can be confidently attributed to IKBKG. This represents a larger chromosomal event beyond the segmental duplication and is more severe than segment-only deletion.

### 11.7kb Segment Duplication

Duplication of the 11.7kb segment is detected but is usually benign regardless of which gene carries the extra copy. Like segment deletions, duplications are reported with an `[UNRESOLVED]` flag since gene attribution cannot be determined from short reads.

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Normal Female (Example: HG00235)

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

### Example 2: Normal Male (Example: HG00140)

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

### Example 3: Female 11.7kb Deletion (Example: HG00733)

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

### Example 4: Male Whole-Gene IKBKGP1 Deletion (Example: HG03225)

When flanks are also reduced, gene attribution IS possible:

```yaml
IKBKG:
  Copy numbers:
    IKBKG: 1
    IKBKGP1: 0
    IKBKGdel_region (IKBKG+IKBKGP1 total): 1
  cn_interpretation: |-
    IKBKG (NEMO) Analysis (male):
      IKBKG flanks: CN=1
      IKBKGP1 flanks: CN=0
      IKBKGdel_region — 11.7kb shared segment (exons 4-10): total CN=1 (expected 2)

    [UNRESOLVED] 11.7kb segment deletion detected (total CN=1, expected=2)
      The 11.7kb segment is identical between IKBKG and IKBKGP1.
      Short reads cannot resolve which gene carries this event.
      Mechanism: NAHR between segmental duplications.
      If IKBKG: pathogenic (exon 4-10 loss). If IKBKGP1: benign.
      Recommend long-read confirmation for definitive assignment.
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

## Validation

The IKBKG caller was validated by comparing short-read results against long-read HPRC v2.0 assembly truth for 223 overlapping samples (112 male, 111 female). Assembly truth was derived from phased haplotype assemblies of 1000 Genomes samples.

### Copy Number Accuracy

| Metric | Overall | Male (n=112) | Female (n=111) |
|--------|---------|--------------|----------------|
| IKBKG CN | 222/223 (99.6%) | 112/112 (100%) | 110/111 (99.1%) |
| IKBKGP1 CN | 222/223 (99.6%) | 112/112 (100%) | 110/111 (99.1%) |
| IKBKGdel region | 222/223 (99.6%) | 112/112 (100%) | 110/111 (99.1%) |

Male copy number is 100% accurate across all three metrics. The only errors are in two female samples with extreme structural events or high copy-number gain events.

### Validation Summary

- Copy number accuracy is 99.6% overall, 100% for males across all metrics
- The only mismatches are in two female samples with extreme structural events (full Block 1 conversion, high-gain duplication)

## Known Limitations

- **11.7kb segment gene attribution**: The most significant limitation. The 11.7kb segment is fully homogenized between IKBKG and IKBKGP1 with no differentiating variants. Short reads cannot determine whether a segment-level CNV (deletion or duplication) affects the functional gene or the pseudogene. These events are reported with an `[UNRESOLVED]` flag. Long-read sequencing or other orthogonal methods are required for definitive gene assignment. Flanking regions ARE distinguishable, so whole-gene deletions (where flanks are also affected) CAN be attributed to a specific gene.
- **No gene conversion detection**: The Block 1 PSVs used for CN calling are the same sites that differentiate IKBKG from IKBKGP1 — a converted region would simply be classified as the other allele, making gene conversion detection circular. The 3' region lacks sufficient differentiating sites for reliable detection.
- **No small variant calling**: Due to pervasive homogenization outside Block 1, variant assignment from phased liftover BAMs is unreliable. Small variant calling is skipped for this gene.
- **Full Block 1 conversion**: When IKBKGP1 has complete Block 1 gene conversion (42/42 PSVs), it becomes indistinguishable from IKBKG in short-read pileup analysis, causing CN misattribution (NA18565 case).
- **High-gain events**: Samples with many tandem duplications in the shared exon 4-10 region can create skewed allele balance that causes false IKBKG gain (NA19909 case).
- **Inverted pseudogene**: IKBKGP1 is inverted relative to IKBKG, which precludes classical fusion detection. Structural events are detected via copy number changes rather than fusion breakpoints.

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
