# CFH Gene Cluster Caller

## Overview

The CFH Gene Cluster Caller is capable of genotyping the complement factor H (CFH) gene cluster from whole-genome sequencing (WGS) data. This specialized caller can accurately determine copy numbers and call small variants in the CFH gene family, which consists of CFH, CFHR1, CFHR2, CFHR3, and CFHR4.

The CFH gene cluster is located on chromosome 1 (chr1:196,652,043-196,948,803) and spans approximately 297 kilobases containing multiple paralogous genes with high sequence similarity that pose significant challenges for standard variant callers.

## Clinical Significance

Structural variations in the CFH gene cluster are associated with important clinical conditions:

- **Atypical hemolytic uremic syndrome (aHUS)**: A rare disease characterized by abnormal blood clotting in small vessels
- **Age-related macular degeneration (AMD)**: A leading cause of vision loss in older adults

Copy number variations (CNVs) in this region, particularly deletions involving CFHR3, CFHR1, and CFHR4, are common in the population and have been implicated in both protective and risk effects for these diseases.

## CFH Gene Cluster Structure

### Genomic Organization

The CFH gene cluster contains four distinct copy number regions analyzed by the Sentieon Segdup Caller:

| Region | Genomic Coordinates | Description |
|--------|---------------------|-------------|
| **CFH** | chr1:196,742,575-196,747,504 | Core CFH gene region (highly conserved) |
| **CFHR3 - CFHR1** | chr1:196,774,840-196,832,189 | CFHR3 and part of CFHR1 (high variability) |
| **CFHR1 - CFHR4 - CFHR2** | chr1:196,814,300-196,936,300 | CFHR1, CFHR4, and part of CFHR2 |
| **CFHR2** | chr1:196,936,300-196,948,803 | CFHR2 terminal region (highly conserved) |

### Homology Regions

The cluster contains multiple pairs of high-homology regions (>99% sequence identity). Sentieon Segdup caller is able to resolve these regions through liftover, global copy-number optimization and phasing.

### Common Copy Number Variations

The following CNV events are commonly observed in the population:

1. **CFHR3+CFHR1 deletion** (heterozygous ~30-40%, homozygous ~20%): Complete or partial deletion of CFHR3 and CFHR1 genes
2. **CFHR1+CFHR4 deletion**: Deletion spanning CFHR1 and CFHR4 genes
3. **CFHR3 duplication**: Duplication events in the CFHR31 region
4. **CFHR1 duplication**: Duplication events affecting CFHR142 region

The CFH and CFHR2 terminal regions are highly conserved with CNVs being extremely rare (>99.5% of individuals have CN=2).

## What Sentieon Segdup Caller Does

The Sentieon Segdup Caller performs a comprehensive three-phase analysis of the CFH gene cluster:

### Phase 1: Liftover and Read Mapping

1. Maps sequencing reads between high-homology regions
2. Uses sophisticated coordinate transformation based on known sequence alignments
3. Generates lifted-over BAM files for improved copy number estimation
4. Handles complex rearrangements and partial overlaps

### Phase 2: Copy Number Determination

1. Calculates normalized read depth across all CFH regions
2. Uses maximum likelihood optimization
3. Incorporates region-specific Bayesian priors reflecting population frequencies:
   - **CFH**: Very tight prior around CN=2 
   - **CFHR31**: Broad prior allowing deletions 
   - **CFHR142**: Moderate deletion frequency
   - **CFHR2**: Very tight prior around CN=2
4. Determines the most likely copy number state for each region (0, 1, 2, 3, or 4 copies)
5. Outputs copy number calls per region with high confidence

### Phase 3: Small Variant Calling

1. Performs ploidy-aware variant calling using Sentieon DNAscope
2. Adjusts expected ploidy based on determined copy numbers:
   - CN=0: No variant calling (region deleted)
   - CN=1: Haploid calling
   - CN=2: Diploid calling (standard)
   - CN=3+: Polyploid calling with appropriate parameters
3. Phases variants
4. Generates final VCF file with accurate genotypes for the CFH gene

## Output Files

The Sentieon Segdup Caller generates the following outputs for CFH analysis:

### 1. Summary YAML File

**Filename**: `{outdir}/{sample}_CFH_summary.yaml`

Contains copy number calls and analysis metadata:

```yaml
sample: HG00253
gene: CFH
Copy numbers:
    CFH: 2
    CFHR31: 1
    CFHR142: 1
    CFHR2: 2
cn_interpretation: |
  CFH: Normal copy number (CN=2)
  CFHR31: Heterozygous deletion (CN=1, chr1:196774840-196832189, 57.3 kb)
  CFHR142: Heterozygous deletion (CN=1, chr1:196814300-196936300, 122.0 kb)
  CFHR2: Normal copy number (CN=2)

  CNV events detected: 2
  Total affected region: 179.3 kb
Variants: output/HG03239_CFH.vcf.gz
```

### 2. VCF File

**Filename**: `{outdir}/{sample}_CFH.vcf.gz`

Contains phased small variants (SNVs and indels) called in the CFH gene cluster. Variants are reported with:
- Accurate genotypes adjusted for local copy number
- Phase information (when available)
- Standard VCF quality metrics (GQ, DP, AD, etc.)
- Appropriate ploidy representation

### Interpretation Guide

| CFH | CFHR31 | CFHR142 | CFHR2 | Interpretation |
|-----|--------|---------|-------|----------------|
| 2   | 2      | 2       | 2     | Normal/Reference |
| 2   | 1      | 2       | 2     | Heterozygous CFHR3+CFHR1 deletion |
| 2   | 0      | 2       | 2     | Homozygous CFHR3+CFHR1 deletion |
| 2   | 2      | 1       | 2     | Heterozygous CFHR1+CFHR4 deletion |
| 2   | 2      | 0       | 2     | Homozygous CFHR1+CFHR4 deletion |
| 2   | 1      | 1       | 2     | Compound heterozygous deletions |
| 2   | 3      | 2       | 2     | CFHR3 duplication |

## Accuracy Benchmark

The Sentieon Segdup Caller was benchmarked against the **HPRC v2.0 dataset** (Human Pangenome Reference Consortium) comprising **213 diverse samples** with truth copy number calls derived from PacBio HiFi long-read sequencing assemblies.

### Benchmark Results

- **Total samples analyzed**: 213
- **Total copy number calls**: 852 (217 samples � 4 regions)
- **Concordant calls**: 848
- **Discrepancies**: 4 (0.47% discordance rate)

### Discrepancy Analysis

All 5 discrepancies were manually reviewed by comparing to raw PacBio HiFi long-read alignments:

| Sample | Region | HPRC v2.0 Call | Sentieon Call | Resolution |
|--------|--------|----------------|---------------|------------|
| HG00272 | CFHR142 | 2 | 3 | **Sentieon correct** - HiFi reads show clear duplication of 1 extra copy |
| HG01934 | CFHR142 | 2 | 3 | **Sentieon correct** - HiFi reads show clear duplication of 1 extra copy|
| HG04204 | CFHR31 | 1 | 2 | **Sentieon correct** - HiFi reads do not show evidence of deletion |
| HG02273 | CFHR31 | 1 | 2 | **Sentieon correct** - HiFi reads do not show evidence of deletion |

**Key Finding**: In all 4 discrepancy cases, the Sentieon Segdup Caller calls were **consistent with PacBio HiFi long-read evidence**, suggesting that the HPRC v2.0 annotations may contain errors in these specific samples. This indicates that the Sentieon caller achieves **>99.5% accuracy** with potential **100% accuracy** when accounting for likely annotation errors.

### Performance by Region

| Region | Total Calls | Discrepancies | Accuracy |
|--------|-------------|---------------|----------|
| CFH | 213 | 0 | 100% |
| CFHR31 | 213 | 2 | 99.1% (100% after validation) |
| CFHR142 | 213 | 2 | 99.1% (100% after validation) |
| CFHR2 | 213 | 0 | 100% |

### Population Variant Frequencies (from HPRC dataset)

Based on the 213-sample benchmark:

- **CFHR31 homozygous deletion (CN=0)**: 19.8%
- **CFHR31 heterozygous deletion (CN=1)**: 32.3% 
- **CFHR142 homozygous deletion (CN=0)**: 0.5% 
- **CFHR142 heterozygous deletion (CN=1)**: 2.8% 
- **CFHR142 duplication (CN=3)**: 0.9% 
- **Compound deletions of CFHR31 and CFHR142**: 1.4% 

These frequencies align well with previous population studies and validate the region-specific prior distributions used in the Bayesian model.

## Usage Example

```bash
# Run CFH gene cluster analysis
segdup-caller \
  --short sample.bam \
  --reference GRCh38.fa \
  --sr_model SentieonGeneCaller0.1.bundle \
  --genes CFH \
  --outdir results/

```

## See Also

- [Sentieon DNAscope Variant Caller](https://www.sentieon.com/products/)
- [HPRC Data Portal](https://humanpangenome.org/)
