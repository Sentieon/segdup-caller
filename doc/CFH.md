# CFH Gene Cluster Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## Overview

The CFH gene cluster on chromosome 1 (chr1:196,652,043-196,948,803) spans approximately 297 kilobases containing the complement factor H (CFH) gene and four related genes (CFHR1-4) with high sequence similarity (>99% in some regions). This gene family regulates the alternative complement pathway.

**Unique features of CFH analysis**:
- **Four copy number regions**: CFH, CFHR31, CFHR142, CFHR2 with distinct CNV patterns
- **Common deletions**: CFHR3+CFHR1 deletions very common (~30-40% heterozygous, ~20% homozygous)
- **Region-specific priors**: Different Bayesian priors for each region reflecting population frequencies
- **High homology resolution**: Liftover and phasing resolve >99% identical regions

## CFH Gene Cluster Structure

### Genomic Organization

The CFH caller analyzes four distinct copy number regions:

| Region | Genomic Coordinates | Description | Typical CN |
|--------|---------------------|-------------|------------|
| **CFH** | chr1:196,742,575-196,747,504 | Core CFH gene region | Highly conserved (>99.5% CN=2) |
| **CFHR31** | chr1:196,774,840-196,832,189 | CFHR3 and part of CFHR1 | High variability (~52% CN≠2) |
| **CFHR142** | chr1:196,814,300-196,936,300 | CFHR1, CFHR4, part of CFHR2 | Low-moderate variability (~4% CN≠2) |
| **CFHR2** | chr1:196,936,300-196,948,803 | CFHR2 terminal region | Highly conserved (>99.5% CN=2) |

### Homology and Complexity

The cluster contains multiple pairs of high-homology regions with >99% sequence identity, making standard variant calling challenging. The Sentieon Segdup Caller resolves these regions through:
- **Liftover**: Coordinate transformation between homologous segments
- **Global CN optimization**: Simultaneous optimization across all four regions
- **Phasing**: Distinguishing variants between similar genes

### Common Copy Number Variations

1. **CFHR3+CFHR1 deletion** (affects CFHR31 region)
   - **Heterozygous**: ~30-40% in general populations
   - **Homozygous**: ~20% in general populations
   - Partially or completely deletes CFHR3 and CFHR1 genes

2. **CFHR1+CFHR4 deletion** (affects CFHR142 region)
   - Rare in most populations
   - Deletes portions of CFHR1 and CFHR4

3. **CFHR duplications**
   - **CFHR31 duplications**: Rare (<1%)
   - **CFHR142 duplications**: Rare (<1%)

4. **CFH and CFHR2 terminal region**
   - Highly conserved
   - CNVs extremely rare (>99.5% have CN=2)

### CNV Frequency and Clinical Impact

| CNV Event | Frequency | Clinical Significance |
|-----------|-----------|----------------------|
| **CFHR3+CFHR1 deletion (het)** | ~30-40% | Common; protective for AMD in some contexts |
| **CFHR3+CFHR1 deletion (hom)** | ~20% | Common; benign or protective |
| **CFHR1+CFHR4 deletion** | Rare | May affect complement regulation |
| **CFHR duplications** | <1% | Clinical significance unclear |
| **CFH deletions/CNVs** | Very rare | Clinically significant; aHUS risk |

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Common CFHR3+CFHR1 Deletions

```yaml
CFH:
  Copy numbers:
    CFH: 2
    CFHR142: 1
    CFHR2: 2
    CFHR31: 1
  Variants: out.cfh/HG03239/HG03239.CFH.result.vcf.gz
  cn_interpretation: |-
    CFH: Normal copy number (CN=2)
    CFHR31: Heterozygous deletion (CN=1, chr1:196774840-196832189,)
    CFHR142: Heterozygous deletion (CN=1, chr1:196814300-196936300,)
    CFHR2: Normal copy number (CN=2)
    CNV events detected: 2
```

### Example 2: Homozygous CFHR3+CFHR1 Deletion

```yaml
CFH:
  Copy numbers:
    CFH: 2
    CFHR142: 2
    CFHR2: 2
    CFHR31: 0
  Variants: out.cfh/HG01109/HG01109.CFH.result.vcf.gz
  cn_interpretation: |-
    CFH: Normal copy number (CN=2)
    CFHR31: Homozygous deletion (CN=0, chr1:196774840-196832189,)
    CFHR142: Normal copy number (CN=2)
    CFHR2: Normal copy number (CN=2)
    CNV events detected: 1
```

### Example 3: Normal (No Deletions)

```yaml
CFH:
  Copy numbers:
    CFH: 2
    CFHR142: 2
    CFHR2: 2
    CFHR31: 2
  Variants: out.cfh/HG02273/HG02273.CFH.result.vcf.gz
  cn_interpretation: |-
    CFH: Normal copy number (CN=2)
    CFHR31: Normal copy number (CN=2)
    CFHR142: Normal copy number (CN=2)
    CFHR2: Normal copy number (CN=2)
```

## Interpretation Guide

### Copy Number Patterns

| CFH | CFHR31 | CFHR142 | CFHR2 | Interpretation | Clinical Notes |
|-----|--------|---------|-------|----------------|----------------|
| 2 | 2 | 2 | 2 | Normal/Reference | No CNVs detected |
| 2 | 1 | 2 | 2 | Heterozygous CFHR3+CFHR1 deletion | Common (~30-40%); may be protective for AMD |
| 2 | 0 | 2 | 2 | Homozygous CFHR3+CFHR1 deletion | Common (~20%); generally benign |
| 2 | 2 | 1 | 2 | Heterozygous CFHR1+CFHR4 deletion | Rare; clinical significance uncertain |
| 2 | 2 | 0 | 2 | Homozygous CFHR1+CFHR4 deletion | Very rare; may affect complement function |
| 2 | 1 | 1 | 2 | Compound heterozygous deletions | Uncommon (~1.4%); clinical significance varies |
| 2 | 3 | 2 | 2 | CFHR3 duplication | Rare (<1%); clinical significance unclear |
| 2 | 2 | 3 | 2 | CFHR1/CFHR4 duplication | Rare (<1%); clinical significance unclear |
| 1 | 2 | 2 | 2 | CFH heterozygous deletion | **Very rare; clinically significant; aHUS risk** |

## Gene-Specific Implementation Details

### Region-Specific Copy Number Priors

CFH uses **Categorical priors** with distinct distributions for each region reflecting population frequencies:

**CFH region**:
- Very tight prior around CN=2 (>99.5% diploid)
- Reflects extreme rarity of CFH deletions

**CFHR31 region** (CFHR3+CFHR1):
- Broad prior allowing common deletions
- Approximate distribution based on HPRC data:
  - CN=0: ~20% (homozygous deletion)
  - CN=1: ~32% (heterozygous deletion)
  - CN=2: ~48% (normal)

**CFHR142 region** (CFHR1+CFHR4+CFHR2):
- Moderate prior allowing rare deletions/duplications
- Approximate distribution:
  - CN=0: ~0.5%
  - CN=1: ~2.8%
  - CN=2: ~95.8%
  - CN=3: ~0.9%

**CFHR2 terminal region**:
- Very tight prior around CN=2 (>99.5% diploid)
- Highly conserved like CFH

### Global Copy Number Optimization

The caller performs **simultaneous optimization** across all four regions, which:
- Accounts for overlapping deletion regions (CFHR31 and CFHR142 share CFHR1)
- Improves accuracy by leveraging correlations between regions
- Resolves ambiguous cases through joint likelihood

## Validation Results

### HPRC v2.0 Benchmark

Validated on **213 diverse samples** from HPRC v2.0 with truth CN calls derived from assemblies. For this validation, we only used Illumina short reads downloaded from s3://human-pangenomics/. 

| Metric | Value |
|--------|-------|
| **Total samples** | 213 |
| **Total CN calls** | 852 (213 samples × 4 regions) |
| **Concordant calls** | 848 |
| **Discrepancies** | 4 |
| **Nominal accuracy** | **99.5%** |

### Discrepancy Resolution

All 4 discrepancies were manually reviewed by comparing to raw PacBio HiFi long-read alignments:

| Sample | Region | HPRC Call | Sentieon Call | HiFi Read Evidence | Resolution |
|--------|--------|-----------|---------------|-------------------|------------|
| HG00272 | CFHR142 | 2 | 3 | Clear duplication (1 extra copy) | **Sentieon correct** |
| HG01934 | CFHR142 | 2 | 3 | Clear duplication (1 extra copy) | **Sentieon correct** |
| HG04204 | CFHR31 | 1 | 2 | No deletion evidence | **Sentieon correct** |
| HG02273 | CFHR31 | 1 | 2 | No deletion evidence | **Sentieon correct** |

**Key finding**: In all 4 cases, Sentieon calls were **consistent with PacBio HiFi long-read evidence**, suggesting HPRC v2.0 annotations may contain errors. This indicates **effective 100% accuracy** when accounting for likely annotation errors.

### Performance by Region

| Region | Total Calls | Discrepancies | Nominal Accuracy | Validated Accuracy |
|--------|-------------|---------------|------------------|-------------------|
| CFH | 213 | 0 | 100% | 100% |
| CFHR31 | 213 | 2 | 99.1% | **100%** (both HPRC errors) |
| CFHR142 | 213 | 2 | 99.1% | **100%** (both HPRC errors) |
| CFHR2 | 213 | 0 | 100% | 100% |

### Population CNV Frequencies (HPRC v2.0)

Based on 213 samples, consistent with literature:

**CFHR31 region**:
- CN=0 (homozygous deletion): 19.8%
- CN=1 (heterozygous deletion): 32.3%
- CN=2 (normal): 47.9%

**CFHR142 region**:
- CN=0 (homozygous deletion): 0.5%
- CN=1 (heterozygous deletion): 2.8%
- CN=2 (normal): 95.8%
- CN=3 (duplication): 0.9%

**Compound deletions** (both CFHR31 and CFHR142): 1.4%

These frequencies validate the region-specific Bayesian priors used in the model.

## Known Limitations

- **Rare complex rearrangements**: Unusual structural variants beyond simple CNVs may require manual review
- **Small variants in homologous regions**: Variants in >99% identical segments may be challenging to phase
- **Novel fusion genes**: CFH::CFHR or CFHR::CFHR fusion genes from non-allelic homologous recombination may not be fully characterized

## CFH-Specific Resources

### Clinical Guidelines
- **aHUS**: Genetic testing and complement studies recommended for patients with suspected aHUS
- **AMD**: CFH genotyping (Y402H and other variants) used for risk stratification in some contexts

### Clinical Databases
- **ClinVar CFH**: https://www.ncbi.nlm.nih.gov/clinvar/?term=CFH
- **aHUS Database**: Information on complement-mediated aHUS genetic causes
