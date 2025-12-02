# HBA Gene Cluster Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## Overview

The HBA gene cluster on chromosome 16p13.3 (chr16:169,454-177,522) contains two tandemly duplicated alpha-globin genes (HBA1 and HBA2) with ~97% sequence identity. Deletions in these genes cause alpha-thalassemia, one of the most common inherited blood disorders worldwide.

**Unique features of HBA analysis**:
- **Tandem duplication structure**: Two highly similar genes arranged in tandem
- **Common deletion detection**: Identifies -α3.7 and -α4.2 deletions (account for ~90% of α-thalassemia)
- **Large deletion detection**: Detects --SEA, --MED, --FIL, --THAI deletions via non-duplication region
- **Clinical phenotype prediction**: Maps genotype to α-thalassemia clinical severity

## HBA Gene Cluster Structure

### Genomic Organization

The HBA caller analyzes three copy number regions:

| Region | Genomic Coordinates | Description |
|--------|---------------------|-------------|
| **-α3.7 deletion** | chr16:172,871-176,674 | 3.7 kb deletion region (rightward deletion) |
| **-α4.2 deletion** | chr16:169,818-174,075 | 4.2 kb deletion region (leftward deletion) |
| **Non-duplication** | chr16:169,455-177,521 | Unique flanking regions (detects large deletions) |

### Tandem Gene Arrangement

| Gene | Locations | Description |
|------|-----------|-------------|
| **HBA2** (α2) | chr16:169,455-170,791 and chr16:175,796-177,521 | Upstream gene, higher expression |
| **HBA1** (α1) | chr16:173,712-175,046 and chr16:171,990-173,710 | Downstream gene, lower expression |

Each gene exists in highly homologous duplication regions that facilitate unequal crossing-over during meiosis, leading to recurrent deletions.

### Common Deletion Types

**Small deletions** (~90% of cases):
- **-α3.7**: 3.7 kb rightward deletion creating functional α2α1 fusion gene. 10-30% heterozygous frequency (ethnicity-dependent)
- **-α4.2**  4.2 kb leftward deletion removing HBA2 gene. 1-5% heterozygous frequency

**Large deletions** (--SEA, --MED, --FIL, --THAI):
- **Mechanism**: Large chromosomal deletions (>20 kb)
- **Result**: Complete removal of entire alpha-globin gene cluster
- **Effect**: Removes both alpha-globin genes (4 genes → 2 genes when heterozygous)

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Heterozygous -α3.7 Deletion (Silent Carrier)

```yaml
HBA:
  Copy numbers:
    a3.7: 1
    a4.2: 2
    non-duplication: 2
  Variants: out.gba/NA19240/NA19240.HBA.result.vcf.gz
  cn_interpretation: |-
    non-duplication: Normal copy number (CN=2)
    -α3.7: Heterozygous deletion (CN=1, chr16:172871-176674)
    -α4.2: Normal copy number (CN=2)
    Clinical interpretation: Heterozygous -α3.7 deletion (αα/-α3.7) - Silent carrier
```

### Example 2: Homozygous -α3.7 Deletion (Alpha-Thalassemia Trait)

```yaml
HBA:
  Copy numbers:
    a3.7: 0
    a4.2: 2
    non-duplication: 2
  Variants: out.hba/HG03453/HG03453.HBA.result.vcf.gz
  cn_interpretation: |-
    -α3.7: Homozygous deletion (CN=0, chr16:172871-176674)
    -α4.2: Normal copy number (CN=2)
    non-duplication: Normal copy number (CN=2)
    Clinical interpretation: Homozygous -α3.7 deletion (-α3.7/-α3.7) - Alpha-thalassemia trait
```

### Example 3: Heterozygous Large Deletion (Silent Carrier for --SEA)

```yaml
HBA:
  Copy numbers:
    a3.7: 0
    a4.2: 0
    non-duplication: 1
  Variants: sample003.HBA.result.vcf.gz
  cn_interpretation: |-
    -α3.7: Deletion (CN=0, chr16:172871-176674)
    -α4.2: Deletion (CN=0, chr16:169818-174075)
    non-duplication: Heterozygous deletion (CN=1)
    Clinical interpretation: Heterozygous large deletion (αα/--), likely --SEA or similar - Silent carrier
```

### Example 4: -α3.7 duplication

```yaml
HBA:
  Copy numbers:
    a3.7: 3
    a4.2: 2
    non-duplication: 2
  Variants: out.hba/HG00735/HG00735.HBA.result.vcf.gz
  cn_interpretation: |-
    -α3.7: Triplication (CN=3, chr16:172871-176674)
    -α4.2: Normal copy number (CN=2)
    non-duplication: Normal copy number (CN=2)
    Clinical interpretation: α-globin triplication involving -α3.7 region - Usually normal
```

## Interpretation Guide

### Copy Number to Phenotype Mapping

| a3.7 | a4.2 | non-dup | Interpretation | α-Globin Genes | Clinical Phenotype |
|------|------|---------|----------------|----------------|-------------------|
| 2 | 2 | 2 | Normal (αα/αα) | 4 | Normal |
| 1 | 2 | 2 | Heterozygous -α3.7 (αα/-α3.7) | 3 | Silent carrier |
| 2 | 1 | 2 | Heterozygous -α4.2 (αα/-α4.2) | 3 | Silent carrier |
| 0 | 2 | 2 | Homozygous -α3.7 (-α3.7/-α3.7) | 2 | Alpha-thalassemia trait |
| 2 | 0 | 2 | Homozygous -α4.2 (-α4.2/-α4.2) | 2 | Alpha-thalassemia trait |
| 1 | 1 | 2 | Compound heterozygous (-α3.7/-α4.2) | 2 | Alpha-thalassemia trait |
| 0 | 0 | 1 | Heterozygous large deletion (αα/--) | 2 | Silent carrier |
| 0 | 0 | 0 | Homozygous large deletion (--/--) | 0 | Hemoglobin Bart's (lethal) |
| 1 | 2 | 1 | Compound -α3.7 + large deletion (-α3.7/--) | 1 | Hemoglobin H disease |
| 2 | 1 | 1 | Compound -α4.2 + large deletion (-α4.2/--) | 1 | Hemoglobin H disease |
| 3 | 2 | 2 | α-globin triplication | 5 | Usually normal |
| 2 | 3 | 2 | α-globin triplication | 5 | Usually normal |

**Note**: Large deletions (--SEA, --MED, --FIL, --THAI) are detected via **non-duplication region CN**:
- **non-dup CN=1**: Heterozygous large deletion (αα/--)
- **non-dup CN=0**: Homozygous large deletion (--/--)

This provides comprehensive detection of both common small deletions and rare large deletions.

## Gene-Specific Implementation Details

### Copy Number Prior

HBA uses **Gaussian priors** reflecting population-specific deletion frequencies:

**Implementation**: Priors can be adjusted for ethnicity where deletion frequencies vary significantly (e.g., Southeast Asian vs European populations).

### Depth Normalization Factor

HBA analysis uses a **depth normalization factor of 1.12** to account for the tandem duplication structure and improve CN estimation accuracy.

### Copy Number State Detection

The caller determines CN state (0, 1, 2, 3, 4) for each of the three regions using maximum likelihood optimization across all regions simultaneously.

## Validation Results

### HPRC v2.0 Benchmark

Validated on **215 diverse samples** from HPRC v2.0 with truth CN calls derived from assemblies. For this validation, we only used Illumina short reads downloaded from s3://human-pangenomics/

| Metric | Value |
|--------|-------|
| **Total samples** | 215 |
| **Total CN calls** | 645 (215 samples × 3 regions) |
| **Concordant calls** | 645 |
| **Discrepancies** | 0 |
| **Accuracy** | **100%** (0% discordance rate) |

### Performance by Region

| Region | Total Calls | Discrepancies | Accuracy |
|--------|-------------|---------------|----------|
| **-α3.7** | 215 | 0 | 100% |
| **-α4.2** | 215 | 0 | 100% |
| **non-duplication** | 215 | 0 | 100% |

**CN distribution in HPRC v2.0**:

| CN State | -α3.7 | -α4.2 | non-dup |
|----------|-------|-------|---------|
| CN=0 | 10 (4.7%) | 1 (0.5%) | 0 (0%) |
| CN=1 | 45 (20.9%) | 7 (3.3%) | 0 (0%) |
| CN=2 | 160 (74.4%) | 207 (96.3%) | 215 (100%) |
| CN=3 | 0 (0%) | 0 (0%) | 0 (0%) |
| CN=4 | 0 (0%) | 0 (0%) | 0 (0%) |

**Key findings**:
- **Perfect concordance** across all 645 copy number calls
- **Common -α3.7 deletions** successfully detected (25.6% combined CN=0+CN=1)
- **Rare -α4.2 deletions** successfully detected (3.8% combined CN=0+CN=1)
- **Non-duplication region** perfectly stable (100% CN=2), validating large deletion detection approach

## Known Limitations

- **Small variant detection**: Focus is on common deletions; rare point mutations in HBA1/HBA2 may require manual review
- **Complex rearrangements**: Unusual rearrangements (e.g., HBA1/HBA2 gene conversion, anti-Lepore variants) may need additional analysis
- **Duplication phasing**: For triplications (CN=3), specific gene assignment may be ambiguous
- **Non-deletion α-thalassemia**: Rare non-deletion mutations (e.g., Constant Spring, Hb Quong Sze) require small variant analysis

## HBA-Specific Resources

### Clinical Guidelines
- **ACMG Guidelines**: Recommendations for α-thalassemia screening and genetic counseling
- **Regional screening programs**: Many regions with high α-thalassemia prevalence have specific guidelines

### Clinical Databases
- **HbVar Database**: http://globin.cse.psu.edu/hbvar/ - Hemoglobin variant database
- **IthaGenes**: http://www.ithagenes.org/ - Hemoglobinopathy mutation database
