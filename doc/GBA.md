# GBA Gene Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## Overview

GBA (glucocerebrosidase) is a critical lysosomal enzyme gene associated with Gaucher disease and Parkinson's disease risk. The GBA gene cluster on chromosome 1 (chr1:155,210,382-155,241,868) contains the functional GBA1 gene and its highly homologous pseudogene GBAP1 with ~96% sequence identity.

**Unique features of GBA analysis**:
- **Gene conversion detection**: Identifies RecNciI, RecTL, and other conversion events
- **Fusion allele detection**: Detects GBA1::GBAP1 fusion creating chimeric genes
- **Conversion hotspot**: Exons 9-11 region is a major recombination site

## What This Module Does

The GBA module performs:

- **Copy number calling**: Determines copy numbers for GBA1 and GBAP1 regions. There is no assembly truth for GBA1/GBAP1 copy number yet.
- **Gene conversion detection**: Identifies RecNciI, RecTL, and individual conversion variants (L483P, A495P, D409H (legacy name of D448H), c.1263del) via HMM-based segmentation. Compared with Tayebi et al. (2025) on 11 1000 Genomes high-coverage samples: A495P, RecTL, N409S and c.115+1G>A are found as in the study; L483P and c.1263del are each missed in one sample; one RecNciI is reported as an unnamed GBA1::GBAP1 fusion (see [1000 Genomes Validation](#1000-genomes-validation)).
- **Fusion allele detection**: Detects GBA1::GBAP1 fusion events where recombination extends to the 3' gene terminus.
- **Small variant calling**: Calls variants in the GBA1/GBAP1 locus using liftover and phasing.

## GBA Gene Cluster Structure

### Genomic Organization

| Gene | Genomic Coordinates | Description |
|------|---------------------|-------------|
| **GBA1** | chr1:155,230,976-155,241,868 | Functional glucocerebrosidase gene (11 exons) |
| **GBAP1** | chr1:155,210,382-155,219,657 | Non-functional pseudogene with ~96% sequence identity |

### Gene Conversion Hotspot

**Recombination region**: chr1:155,234,450-155,236,000 (exons 9-11 of GBA1)

This region is a hotspot for homologous recombination that transfers GBAP1 sequence into GBA1, creating pathogenic conversion alleles.

## Common Gene Conversion Events

> **Nomenclature note**: This document uses HGVS protein names, which count the 39-residue signal peptide. The one exception is D409H, a legacy name (HGVS D448H), which this page marks "(legacy)" because the caller reports it under that name. Legacy names still common in the literature: L444P → L483P, N370S → N409S, A456P → A495P, V460V → V499V, D409H → D448H.

### RecNciI (RecTL-NciI)

**Structure**: ~50 bp conversion event spanning three variants

**Signature variants**:
- L483P (c.1448T>C): Pathogenic missense
- A495P (c.1483G>C): Pathogenic missense
- V499V (c.1497G>C): Synonymous

**Clinical significance**:
- Gaucher disease type 1 when biallelic
- Strong Parkinson's disease risk factor (5-10× increased risk)
- One of the most common GBA1 conversion alleles

### RecTL (RecA)

**Structure**: Larger conversion event spanning exons 9-11

**Signature variants**:
- c.1263del (55-bp deletion, c.1263_1317del): Frameshift
- D409H (legacy; HGVS D448H, c.1342G>C): Pathogenic missense
- L483P (c.1448T>C): Pathogenic missense
- A495P (c.1483G>C): Pathogenic missense
- V499V (c.1497G>C): Synonymous

**Clinical significance**:
- Gaucher disease type 1 when biallelic
- Parkinson's disease risk factor
- Common in certain populations

### Individual Conversion Variants

| Variant | Type | Clinical Significance |
|---------|------|----------------------|
| **c.1263del** | 55-bp deletion (c.1263_1317del), frameshift | Gaucher disease; can occur alone or as part of RecTL |
| **D409H** (legacy) | Missense | Gaucher disease; isolated conversion or part of RecTL |
| **L483P** | Missense | Gaucher disease; isolated conversion or part of RecNciI/RecTL |
| **A495P** | Missense | Rare isolated conversion; part of RecNciI/RecTL |
| **V499V** | Synonymous | Usually part of larger RecNciI event |

## Gene Fusion Events

### GBA1::GBAP1_Fusion

**Structure**: Recombination at 3' end of GBA1 (exons 9-11) creating chimeric gene

**Detection method**:
- Position-based (last converted segment toward 3' terminus at chr1:155,234,452)
- Extends to gene terminus (distinguishing feature from internal conversions)

**Breakpoint location**: chr1:155,230,000-155,238,000 (recombination hotspot region)

**Clinical significance**:
- Creates chimeric GBA1::GBAP1 protein with altered C-terminal structure
- May cause Gaucher disease or increase Parkinson's disease risk
- Frequency: Rare (<1%) but clinically significant

**Key distinction from conversions**: Fusion events replace the 3' end of GBA1 with GBAP1 sequence extending to the gene terminus, while gene conversions transfer discrete internal segments.

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Gene Conversion (RecTL)

```yaml
GBA1:
  Copy numbers:
    GBA1: 2
    GBAP1: 2
  Variants: HG00115.GBA1.result.vcf.gz
  cn_interpretation: |-
    GBA1: CN=2
    GBAP1: CN=2
    Gene conversion/fusion events detected: 1
      Event 1 (CONVERSION): chr1:155234902-155235917
        Converted alleles: 1
        Conversion sites: rs708606, V499V, A495P, L483P, rs426516, D409H, c.1263del
        Interpretation: RecTL (match score: 100.0%, variants: D409H, V499V, A495P, c.1263del, L483P)
  gene_conversions:
  - conversion_sites:
    - rs708606
    - V499V
    - A495P
    - L483P
    - rs426516
    - D409H
    - c.1263del
    converted_alleles: 1
    interpretation:
    - event_name: RecTL
      match_score: 1.0
      matched_variants:
      - D409H
      - V499V
      - A495P
      - c.1263del
      - L483P
    is_fusion_candidate: false
    region: chr1:155234902-155235917
```

### Example 2: Gene Fusion (GBA1::GBAP1)

```yaml
GBA1:
  Copy numbers:
    GBA1: 2
    GBAP1: 1
  Variants: HG00422.GBA1.result.vcf.gz
  cn_interpretation: |-
    GBA1: CN=2
    GBAP1: CN=1
    Gene conversion/fusion events detected: 1
      Event 1 (FUSION): chr1:155233045-155235726
        Converted alleles: 1
        Conversion sites: None
        Interpretation: GBA1::GBAP1_Fusion (position-based match)
        Clinical note: Gene fusion event (GBA1::GBAP1) detected. This may result in altered protein structure and function.
  gene_conversions:
  - converted_alleles: 1
    fusion_type: 3_PRIME
    interpretation:
    - event_name: GBA1::GBAP1_Fusion
    is_fusion_candidate: true
    region: chr1:155233045-155235726
```

## Interpretation Guide

### Copy Number Patterns

| GBA1 | GBAP1 | Interpretation | Clinical Significance |
|------|-------|----------------|----------------------|
| 2 | 2 | Normal/Reference | No copy number abnormality |
| 2 | 1 | Heterozygous GBAP1 deletion | Common (~20%), benign |
| 2 | 0 | Homozygous GBAP1 deletion | Uncommon (~5%), benign |
| 2 | 3 | GBAP1 duplication | Rare, generally benign |
| 1 | 2 | Heterozygous GBA1 deletion | Rare, requires confirmation |
| 0 | 2 | Homozygous GBA1 deletion | Extremely rare, likely artifact |

### Gene Conversion/Fusion Patterns

| Event Type | Event Name | Alleles | Clinical Significance |
|------------|-----------|---------|----------------------|
| CONVERSION | RecNciI (L483P+A495P+V499V) | 1 | Gaucher carrier; PD risk factor |
| CONVERSION | RecNciI (L483P+A495P+V499V) | 2 | Gaucher type 1; high PD risk |
| CONVERSION | RecTL (c.1263del+D409H (legacy)+L483P+A495P+V499V) | 1 | Gaucher carrier; PD risk factor |
| CONVERSION | RecTL (c.1263del+D409H (legacy)+L483P+A495P+V499V) | 2 | Gaucher type 1; high PD risk |
| CONVERSION | c.1263del only | 1 | Gaucher carrier (frameshift) |
| CONVERSION | c.1263del only | 2 | Gaucher disease |
| CONVERSION | D409H (legacy) only | 1 | Gaucher carrier; PD risk |
| CONVERSION | D409H (legacy) only | 2 | Gaucher type 1; high PD risk |
| CONVERSION | L483P only | 1 | Gaucher carrier; PD risk |
| CONVERSION | L483P only | 2 | Gaucher type 1; high PD risk |
| FUSION | GBA1::GBAP1_Fusion | 1 | Gaucher carrier; potential PD risk |
| FUSION | GBA1::GBAP1_Fusion | 2 | Gaucher disease; high PD risk |

## Gene-Specific Implementation Details

### Copy Number Prior

GBA uses **Categorical priors** with separate distributions for each region, reflecting their distinct population frequencies:

**GBA1**:
- Very tight prior around CN=2 (99.75% expected diploid)
- Categorical distribution: `{0: 0.0005, 1: 0.001, 2: 0.9975, 3: 0.0005, 4: 0.0005}`
- Reflects extreme rarity of GBA1 deletions

**GBAP1**:
- Broader prior allowing deletions
- Categorical distribution: `{0: 0.05, 1: 0.20, 2: 0.74, 3: 0.005, 4: 0.005}` (74% CN=2, 20% CN=1, 5% CN=0)
- Reflects common GBAP1 deletions in population

### Conversion Detection Parameters

HMM statistical segmentation:
- Uses CBS (Circular Binary Segmentation), HMM, or PELT algorithms
- Minimum signal threshold: ≥2 conversion signals (configurable via `min_num_signals`)
- Strand-aware detection distinguishes 3' vs 5' events

### Fusion Detection

- **3' terminus position**: chr1:155,234,452 for GBA1
- **Recombination region**: chr1:155,230,000-155,238,000
- **Method**: Last converted segment toward 3' terminus identifies fusion candidates

## Validation Results

### 1000 Genomes Validation

Validated against **Tayebi et al. (2025)** findings on 1000 Genomes high-coverage samples. For this validation, we only used Illumina short reads downloaded from s3://human-pangenomics/ or from 1000-genome data repository. 

The study column gives each sample's finding from Tayebi et al.; the caller columns give what the caller reports on the 1000 Genomes high-coverage short reads.

| Sample | Study finding (Tayebi et al. 2025) | Caller total CN (GBA1 / GBAP1) | Caller-reported event | Agreement |
|--------|------------------------------------|--------------------------------|---------------------------------|-----------|
| **CN Loss (GBA1 Intact)** | | | | |
| NA21133 | CN loss, GBA1 intact; WT,WT | 3 (2 / 1) | GBA1::GBAP1 fusion (position-based match) | CN loss agrees; the caller also reports a fusion |
| NA21099 | CN loss, GBA1 intact; WT,WT | 3 (2 / 1) | none | agrees |
| NA19922 | CN loss, GBA1 intact; WT,WT | 3 (2 / 1) | mono-allelic rs708606 (possible GBAP1→GBA1 conversion) | CN loss agrees; the caller also flags rs708606 |
| **Monoallelic Conversion** | | | | |
| HG02615 | A495P,WT | 4 (2 / 2) | mono-allelic A495P (possible GBAP1→GBA1 conversion) | agrees |
| NA20815 | L483P,WT | 4 (2 / 2) | mono-allelic rs708606; L483P not called | L483P missed |
| HG02439 | c.1263del,WT | 4 (2 / 2) | none; c.1263del not called | c.1263del missed |
| **Specific Recombinants** | | | | |
| HG00422 | RecNciI,WT | 3 (2 / 1) | GBA1::GBAP1 fusion (position-based match); mono-allelic rs708606 | reported as a fusion, not named RecNciI |
| HG00119 | RecTL,WT | 4 (2 / 2) | RecTL (match score 80%, c.1263del not matched); bi-allelic rs708606 | agrees |
| HG00115 | RecTL,WT | 4 (2 / 2) | RecTL (match score 80%, c.1263del not matched); bi-allelic rs708606 | agrees |
| **Other (Heterozygous SNV)** | | | | |
| NA20503 | N409S,WT | 4 (2 / 2) | heterozygous chr1:155,235,843 T>C (N409S; no VCF name) | agrees; SNV called without a name |
| HG02514 | c.115+1G>A,WT | 4 (2 / 2) | heterozygous chr1:155,240,629 C>T (c.115+1G>A; VCF ID 155240628_C_T) | agrees; SNV called without a name |

**Key findings**:
- **Copy number:** the three samples the study lists as CN loss read total 3 (GBA1 2, GBAP1 1). HG00422 also reads total 3, and the other seven read 4.
- **Named as in the study:** A495P (HG02615) and RecTL (HG00119, HG00115).
- **Called without a name:** N409S (NA20503) and c.115+1G>A (HG02514) are heterozygous SNVs in the VCF.
- **Not reproduced:** L483P (NA20815) and c.1263del (HG02439) are not called. HG00422's RecNciI is reported as a position-based GBA1::GBAP1 fusion, without the RecNciI name.
- **Also reported:** a fusion in NA21133 and a possible conversion at rs708606 in NA19922. The study lists neither.

### HPRC Validation

There is no assembly truth for GBA1/GBAP1 copy number yet.

## Known Limitations

- **Novel conversions**: Events not in the database will be reported but may lack specific interpretation
- **Fusion breakpoint precision**: Position-based detection identifies candidates; exact breakpoint may vary within recombination region
- **Compound heterozygotes**: Phasing may be ambiguous for samples with multiple conversion events
- **Low coverage**: Conversion detection sensitivity decreases below 20× coverage (≥30× recommended for fusion detection)

## GBA-Specific Resources

### Clinical Databases
- **ClinVar GBA**: https://www.ncbi.nlm.nih.gov/clinvar/?term=GBA
- **NCBI Gaucher Disease GeneReview**: https://www.ncbi.nlm.nih.gov/books/NBK1269/

### Reference Data
- **Human Pangenome Reference Consortium (HPRC)**: https://humanpangenome.org/ - Population-scale haplotype data used for validation

### Literature

1. Tayebi N, Lichtenberg J, Hertz E, et al. (2025) "Is Gauchian genotyping of GBA1 variants reliable?" *Communications Biology* 8:718. https://doi.org/10.1038/s42003-025-08059-y
2. Liao WW, Asber M, Didion JP, et al. (2023) "A draft human pangenome reference" *Nature* 617:312-324. https://doi.org/10.1038/s41586-023-05896-x
