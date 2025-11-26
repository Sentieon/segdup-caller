# GBA Gene Caller

## Overview

The GBA Gene Caller is capable of genotyping the glucocerebrosidase gene (GBA) from whole-genome sequencing (WGS) data. This specialized caller can accurately determine copy numbers, detect gene conversion events, and call small variants in the GBA gene and its highly homologous pseudogene GBAP1.

The GBA gene cluster is located on chromosome 1 (chr1:155,210,382-155,241,868) and spans approximately 31 kilobases containing the functional GBA1 gene and the non-functional pseudogene GBAP1 with ~96% sequence identity that poses significant challenges for standard variant callers.

## Clinical Significance

Variants in the GBA gene are associated with important clinical conditions:

- **Gaucher disease**: A lysosomal storage disorder caused by biallelic pathogenic variants in GBA1, characterized by accumulation of glucocerebroside in macrophages
  - **Type 1** (non-neuronopathic): Most common form with hepatosplenomegaly, anemia, thrombocytopenia, and bone disease
  - **Type 2** (acute neuronopathic): Severe form with early-onset neurological involvement
  - **Type 3** (chronic neuronopathic): Intermediate form with later-onset neurological symptoms

- **Parkinson's disease risk**: Heterozygous GBA1 variants (particularly gene conversion alleles) are the most common genetic risk factor for Parkinson's disease, increasing lifetime risk 5-10 fold

Gene conversion events between GBA1 and GBAP1 create complex alleles that are difficult to detect with standard variant calling approaches. These conversion alleles account for approximately 15-30% of pathogenic GBA1 variants in Ashkenazi Jewish populations and 5-10% in other populations.

## GBA Gene Cluster Structure

### Genomic Organization

The GBA gene cluster contains two highly similar genes analyzed by the Sentieon Segdup Caller:

| Gene | Genomic Coordinates | Description |
|------|---------------------|-------------|
| **GBA1** | chr1:155,230,976-155,241,868 | Functional glucocerebrosidase gene (11 exons) |
| **GBAP1** | chr1:155,210,382-155,219,657 | Non-functional pseudogene with ~96% sequence identity |

### Gene Conversion Hotspot

The gene conversion region (chr1:155,234,450-155,236,000) corresponds to exons 9-11 of GBA1, a hotspot for recombination events that transfer variants from GBAP1 into GBA1.

### Common Gene Conversion Events

The following gene conversion events are commonly observed:

> **Note on nomenclature**: This document uses HGVS standard nomenclature. Legacy names are still common in literature: L444P → L483P, N370S → N409S, D409H → D448H, A495P → A534P. The 39-residue difference reflects the signal peptide.

1. **RecNciI** (also called RecTL-NciI or complex recombinant allele): A ~50 bp conversion event spanning three variants:
   - L483P (c.1448T>C): Pathogenic missense variant
   - A495P (c.1483G>C): Pathogenic missense variant
   - V499V (c.1497G>C): Synonymous variant
   - **Clinical significance**: Causes Gaucher disease type 1; strong Parkinson's disease risk factor

2. **RecTL** (also called RecA or terminal locus recombinant): A larger conversion event spanning exons 9-11:
   - c.1263del (84GG deletion): Frameshift variant
   - D409H (c.1226A>C): Pathogenic missense variant
   - L483P (c.1448T>C): Pathogenic missense variant
   - A495P (c.1483G>C): Pathogenic missense variant
   - V499V (c.1497G>C): Synonymous variant
   - **Clinical significance**: Causes Gaucher disease type 1; Parkinson's disease risk factor

3. **Individual conversion variants**:
   - **c.1263del alone**: 84GG deletion (~55 bp), frameshift variant
   - **D409H alone**: Can occur as isolated conversion event
   - **L483P alone**: Can occur as isolated conversion event
   - **A495P alone**: Rare isolated conversion event
   - **V499V alone**: Usually part of larger RecNciI event

### Gene Fusion Events

In addition to gene conversions, the caller can detect **GBA1::GBAP1 fusion events** where recombination occurs at the 3' end of GBA1 (exons 9-11), creating a chimeric gene structure:

**GBA1::GBAP1_Fusion**: A recombination event extending to the 3' terminus of GBA1:
   - **Breakpoint location**: chr1:155,230,000-155,238,000 (recombination hotspot region)
   - **Detection method**: Position-based (last converted segment toward 3' terminus)
   - **Clinical significance**: Creates chimeric GBA1::GBAP1 protein with altered C-terminal structure; may cause Gaucher disease or increase Parkinson's disease risk
   - **Frequency**: Rare (<1% of GBA variants), but clinically significant when present

**Key distinction from conversions**: While gene conversions transfer discrete segments from GBAP1 into GBA1, fusion events result in a chimeric gene where the 3' end of GBA1 is replaced by the corresponding GBAP1 sequence extending to the gene terminus.


### GBAP1 Copy Number Variations

GBAP1 deletions are relatively common in the population:

- **Heterozygous GBAP1 deletion** (~20% in some populations): Generally benign but may complicate genetic counseling
- **Homozygous GBAP1 deletion** (~5%): Benign; may affect gene conversion frequency
- **GBAP1 duplications** (rare): Generally benign

The GBA1 gene is highly conserved with deletions being extremely rare (<0.1% of individuals have CN≠2).

## What Sentieon Segdup Caller Does

The Sentieon Segdup Caller performs a comprehensive three-phase analysis of the GBA gene cluster:

### Phase 1: Liftover and Read Mapping

1. Maps sequencing reads between GBA1 and GBAP1 regions
2. Uses sophisticated coordinate transformation based on known sequence alignments
3. Generates lifted-over BAM files for improved variant detection
4. Handles the 96% sequence identity between GBA1 and GBAP1

### Phase 2: Copy Number Determination

1. Calculates normalized read depth across GBA1 and GBAP1 regions
2. Uses maximum likelihood optimization with gene-specific Bayesian priors:
   - **GBA1**: Very tight prior around CN=2 (99.8% expected diploid)
   - **GBAP1**: Broader prior allowing deletions (74% diploid, 20% heterozygous deletion, 5% homozygous deletion)
3. Determines the most likely copy number state for each gene (0, 1, 2, 3, or 4 copies)
4. Outputs copy number calls per gene with high confidence

### Phase 3: Gene Conversion and Fusion Detection

1. Identifies variants that match GBAP1 sequence in the GBA1 locus
2. Detects regions with clustered conversion signals using statistical segmentation (HMM, CBS, or PELT)
3. Performs strand-aware fusion detection:
   - Identifies the last converted segment toward the 3' terminus (chr1:155,234,452 for GBA1)
   - Optionally validates breakpoint location in recombination region (chr1:155,230,000-155,238,000)
   - Distinguishes fusion events from internal gene conversions
4. Matches detected events against known conversion/fusion allele database
5. Reports:
   - Event type (CONVERSION or FUSION)
   - Event coordinates
   - Number of converted alleles (heterozygous or homozygous conversion)
   - Specific conversion sites (variant IDs)
   - Matched known events (RecNciI, RecTL, GBA1::GBAP1_Fusion, etc.)
   - Match quality scores

### Phase 4: Small Variant Calling

1. Performs ploidy-aware variant calling using Sentieon DNAscope
2. Adjusts expected ploidy based on determined copy numbers:
   - CN=0: No variant calling (gene deleted)
   - CN=1: Haploid calling
   - CN=2: Diploid calling (standard)
   - CN=3+: Polyploid calling with appropriate parameters
3. Phases variants to distinguish GBA1 from GBAP1 alleles
4. Generates final VCF file with accurate genotypes for both genes

## Output Files

The Sentieon Segdup Caller generates the following outputs for GBA analysis:

### 1. Summary YAML File

**Filename**: `{outdir}/{sample}_GBA_summary.yaml`

Contains copy number calls, gene conversion events, and analysis metadata:

**Example 1: Gene conversion event (RecNciI)**
```yaml
sample: HG00253
gene: GBA
Copy numbers:
  GBA1: 2
  GBAP1: 1
gene_conversions:
  - region: "chr1:155235203-155235252"
    converted_alleles: 1
    is_fusion_candidate: false
    conversion_sites:
      - L483P
      - A495P
      - V499V
    interpretation:
      - event_name: RecNciI
        match_score: 1.0
        matched_variants:
          - L483P
          - A495P
          - V499V
cn_interpretation: |
  GBAP1: Heterozygous deletion detected (CN=1)

  Gene conversion/fusion events detected: 1
    Event 1 (CONVERSION): chr1:155235203-155235252
      Converted alleles: 1
      Conversion sites: L483P, A495P, V499V
      Interpretation: RecNciI (match score: 100.0%, variants: L483P, A495P, V499V)
Variants: output/HG00253_GBA.vcf.gz
```

**Example 2: Gene fusion event (GBA1::GBAP1_Fusion)**
```yaml
sample: HG00254
gene: GBA
Copy numbers:
  GBA1: 2
  GBAP1: 2
gene_conversions:
  - region: "chr1:155234000-155238000"
    converted_alleles: 1
    is_fusion_candidate: true
    conversion_sites: []
    interpretation:
      - event_name: GBA1::GBAP1_Fusion
        match_score: 1.0
        matched_variants: []
cn_interpretation: |
  GBAP1: Normal copy number (CN=2)

  Gene conversion/fusion events detected: 1
    Event 1 (FUSION): chr1:155234000-155238000
      Converted alleles: 1
      Interpretation: GBA1::GBAP1_Fusion (position-based detection)
      Clinical note: Gene fusion event (GBA1::GBAP1) detected. This may result in altered protein structure and function.
Variants: output/HG00254_GBA.vcf.gz
```

### 2. VCF File

**Filename**: `{outdir}/{sample}_GBA.vcf.gz`

Contains phased small variants (SNVs and indels) called in the GBA gene cluster. Variants are reported with:
- Genotypes adjusted for local copy number
- Phase information distinguishing GBA1 from GBAP1 alleles
- Standard VCF quality metrics (GQ, DP, AD, etc.)
- Appropriate ploidy representation

### Interpretation Guide

#### Copy Number Patterns

| GBA1 | GBAP1 | Interpretation |
|------|-------|----------------|
| 2    | 2     | Normal/Reference |
| 2    | 1     | Heterozygous GBAP1 deletion (common, benign) |
| 2    | 0     | Homozygous GBAP1 deletion (benign) |
| 2    | 3     | GBAP1 duplication (rare, benign) |
| 1    | 2     | Heterozygous GBA1 deletion (rare, requires confirmation) |
| 0    | 2     | Homozygous GBA1 deletion (extremely rare, likely artifact) |

#### Gene Conversion and Fusion Patterns

| Event Type | Event Name | Converted Alleles | Clinical Significance |
|------------|-----------|-------------------|----------------------|
| CONVERSION | RecNciI (L483P+A495P+V499V) | 1 | Gaucher disease carrier; Parkinson's risk factor |
| CONVERSION | RecNciI (L483P+A495P+V499V) | 2 | Gaucher disease type 1; high Parkinson's risk |
| CONVERSION | RecTL (c.1263del+D409H+L483P+A495P+V499V) | 1 | Gaucher disease carrier; Parkinson's risk factor |
| CONVERSION | RecTL (c.1263del+D409H+L483P+A495P+V499V) | 2 | Gaucher disease type 1; high Parkinson's risk |
| CONVERSION | c.1263del only | 1 | Gaucher disease carrier (frameshift) |
| CONVERSION | c.1263del only | 2 | Gaucher disease |
| CONVERSION | D409H only | 1 | Gaucher disease carrier; Parkinson's risk factor |
| CONVERSION | D409H only | 2 | Gaucher disease type 1; high Parkinson's risk |
| CONVERSION | L483P only | 1 | Gaucher disease carrier; Parkinson's risk factor |
| CONVERSION | L483P only | 2 | Gaucher disease type 1; high Parkinson's risk |
| CONVERSION | A495P only | 1 | Gaucher disease carrier |
| CONVERSION | A495P only | 2 | Gaucher disease |
| CONVERSION | V499V only | 1 | Synonymous; usually part of larger event |
| FUSION | GBA1::GBAP1_Fusion | 1 | Gaucher disease carrier; potential Parkinson's risk; altered protein structure |
| FUSION | GBA1::GBAP1_Fusion | 2 | Gaucher disease; high Parkinson's risk; chimeric protein |
| - | No events detected | - | No known conversion/fusion alleles |

## Accuracy and Validation

The Sentieon Segdup Caller for GBA has been validated on diverse sample cohorts in HRPC and 1000-genome samples:

- **GBA1 and GBA copy number**: 100% accuracy on validation cohorts of HPRC v1.0 release samples

### Gene conversion detection validation

Validated the findings by Tayebi N, et al. Below is a selected summary of samples with varying types of structural events. The call is fully consistent with the paper's finding.

| Sample Name | Coverage | Total CN | Recombinant Variant Genotype | Complete Genotype |
|-------------|----------|----------|------------------------------|-------------------|
| **CN Loss (GBA1 Intact)** | | | | |
| NA21133.final | 1000 Genomes (High Coverage) | 3 | - | WT,WT |
| NA21099.final | 1000 Genomes (High Coverage) | 3 | - | WT,WT |
| NA19922.final | 1000 Genomes (High Coverage) | 3 | - | WT,WT |
| **Monoallelic Conversion** | | | | |
| HG02615.final | 1000 Genomes (High Coverage) | 4 | A495P/ | A495P,WT |
| NA20815.final | 1000 Genomes (High Coverage) | 4 | L483P/ | L483P,WT |
| HG02439.final | 1000 Genomes (High Coverage) | 4 | c.1263del/ | c.1263del,WT |
| **Specific Recombinants** | | | | |
| HG00422.final | 1000 Genomes (High Coverage) | 3 | RecNciI/ | RecNciI,WT |
| HG00119.final | 1000 Genomes (High Coverage) | 4 | c.1263del+RecTL/ | c.1263del+RecTL,WT |
| HG00115.final | 1000 Genomes (High Coverage) | 4 | c.1263del+RecTL/ | c.1263del+RecTL,WT |
| **Other (Heterozygous SNV)** | | | | |
| NA20503.final | 1000 Genomes (High Coverage) | 4 | - | c.1226A>G(p.Asn409Ser),WT (g.155235843T>C) |
| HG02514.final | 1000 Genomes (High Coverage) | 4 | - | c.115+1G>A,WT (g.155240629C>T) |

**Key findings:**
- **CN Loss samples** (NA21133, NA21099, NA19922): Total CN=3 indicates GBAP1 heterozygous deletion with intact GBA1 (CN=2)
- **Monoallelic conversions** successfully detected for L483P, A495P, and c.1263del variants
- **Specific recombinants** correctly identified:
  - HG00422: RecNciI conversion event detected
  - HG00119 and HG00115: Both c.1263del and RecTL conversion events detected
- **Non-conversion variants** properly distinguished from conversions (NA20503, HG02514)
- All samples validated against Tayebi et al. (2025) findings with 100% concordance


## Usage Example

```bash
# Run GBA gene analysis
segdup-caller \
  --short sample.bam \
  --reference GRCh38.fa \
  --sr_model SentieonGeneCaller0.1.bundle \
  --genes GBA \
  --outdir results/
```

## Technical Notes

- **Reference genome**: GRCh38/hg38 (coordinates not compatible with hg19)
- **Recommended coverage**: ≥30× for optimal copy number and conversion/fusion detection 
- **Gene conversion detection threshold**: ≥2 conversion signals (configurable via `min_num_signals`)

## Interpretation Caveats

- **GBAP1 sequence**: Conversion/fusion detection depends on accurate GBAP1 sequence in the reference genome
- **Novel conversions/fusions**: Events other than mentioned above will be reported but not interpreted
- **Fusion breakpoint precision**: Position-based detection identifies fusion candidates; precise breakpoint may vary within recombination region
- **Compound heterozygotes**: Phasing may be ambiguous for samples with multiple conversion/fusion events
- **Low coverage**: Conversion/fusion detection sensitivity decreases below 20× coverage; ≥30× recommended for fusion detection

## References

1. Tayebi N, Lichtenberg J, Hertz E, et al. (2025) "Is Gauchian genotyping of GBA1 variants reliable?" *Communications Biology* 8:718. https://doi.org/10.1038/s42003-025-08059-y

## See Also

- [Sentieon DNAscope Variant Caller](https://www.sentieon.com/products/)
- [NCBI Gaucher Disease GeneReview](https://www.ncbi.nlm.nih.gov/books/NBK1269/)
- [GBA Gene Conversion Database](https://www.ncbi.nlm.nih.gov/clinvar/?term=GBA)
