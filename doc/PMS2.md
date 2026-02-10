# PMS2 Gene Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## Overview

PMS2 (PMS1 homolog 2, mismatch repair system component) is a DNA mismatch repair gene associated with Lynch syndrome. The PMS2 gene cluster on chromosome 7p22.1 contains the functional PMS2 gene and its highly homologous pseudogene PMS2CL with ~98% sequence identity in exons 11-15.

**Unique features of PMS2 analysis**:
- **Haplotype-aware PSV identification**: Uses population-scale haplotype data rather than reference genome alone
- **Block-based conversion detection**: Identifies clinically significant Exon 13-15 conversion events
- **Direction-aware interpretation**: Distinguishes pathogenic (PMS2CL→PMS2) from benign (PMS2→PMS2CL) conversions
- **Fusion detection**: Identifies recombinant alleles extending to 3' terminus

## PMS2 Gene Cluster Structure

### Genomic Organization

| Gene | Genomic Coordinates | Description |
|------|---------------------|-------------|
| **PMS2** | chr7:5,965,000-5,993,410 | Functional mismatch repair gene (15 exons, reverse strand) |
| **PMS2CL** | chr7:6,733,563-6,759,796 | Non-functional pseudogene (~98% identity in exons 11-15) |

### Reverse Strand Orientation

Since *PMS2* is on the **reverse strand**, "Upstream" (5') is at a *higher* genomic coordinate than "Downstream" (3'). The following diagram clarifies the relationship between genomic coordinates, gene structure, and block architecture:

```
Genomic Coordinates (GRCh38):
  5,990k  <-- (Upstream / 5') ----------------------- (Downstream / 3') -->  5,970k

  Gene:      [ Exons 1-10 ] ------- [ Exons 11-12 ] ------- [ Exon 13 ] ------- [ Exon 14-15 ]
  Role:        Upstream          Stable Anchor in PMS2        Hotspot               Tail
                               Mosaic conversion in PMS2CL
                                    (Reference)              (Trigger)            (Extension)
```

### Homology Region

**Active recombination region**: chr7:5,977,708-5,979,621 (exons 13-15 of PMS2)

This region is a hotspot for gene conversion events that can transfer PMS2CL sequence into PMS2, creating pathogenic alleles associated with Lynch syndrome. Because PMS2 is located on the **reverse strand**, the block coordinates map from the 5' (stable) to 3' (active) end of the gene.

### Long Deletion Events

The caller also detects structural deletions affecting PMS2 function:

| Deletion | Region | Description |
|----------|--------|-------------|
| **PMS2_exon13-15.DEL** | chr7:5,970,925-5,979,200 | Deletion spanning exons 13-15 |
| **PMS2CL_exon13-14.DEL** | chr7:6,748,000-6,750,000 | Pseudogene deletion (benign) |

## Common Gene Conversion Events

### Classic Exon 13-15 Conversion

**Structure**: Conversion event spanning the active hotspot region 

**Signature**: PMS2CL sequence replaces PMS2 in exons 13-15

**Clinical significance**:
- Loss of mismatch repair function
- Lynch syndrome carrier when monoallelic
- Lynch syndrome when biallelic
- One of the most common pathogenic PMS2 structural variants

### Recombinant Alleles (3' Extension)

**Structure**: Conversion extending to/beyond the 3' terminus boundary (Block 1)

**Detection method**:
- Position-based: event extends to chr7:5,977,708 (Block 1 boundary)
- Creates fusion-like alleles with PMS2CL sequence at 3' end
- Segments must be contiguous (gap < 500 bp). This tolerance bridges small drops in coverage or short non-informative regions to reconstruct the full event.

**Clinical significance**:
- Similar to classic conversion: loss of MMR function
- May indicate more extensive recombination event

### Reciprocal Exchanges

**Structure**: Simultaneous bidirectional exchange where PMS2 sequence transfers to PMS2CL on one chromosome while PMS2CL sequence transfers to PMS2 on the homologous chromosome.

**Population frequency**: Based on HPRC v2.0 assembly analysis, balanced reciprocal exchanges occur in approximately 3-5% of the population. These events are mechanistically distinct from unidirectional gene conversions and likely arise from meiotic crossover events with gene conversion at the breakpoints.

**Detection challenge**: Balanced reciprocal exchanges are inherently undetectable using short-read sequencing because:
- **No copy number change**: Total PMS2 + PMS2CL copy number remains 4
- **Balanced PSV frequencies**: Each PSV site shows 50% gene allele and 50% pseudogene allele, identical to the normal diploid state
- **No distinguishing signal**: The aggregate read data from both chromosomes appears indistinguishable from a sample with no conversion events

**Clinical implications**: When a balanced reciprocal exchange occurs, one allele carries a pathogenic PMS2CL→PMS2 conversion (loss of MMR function) while the other carries a benign PMS2→PMS2CL conversion. The individual is effectively a Lynch syndrome carrier, but this cannot be detected from short-read data alone.

**Detection strategies**:
- **Haplotype-aware detection**: By leveraging distinctive patterns of the PMS2->PMS2CL and PMS2CL->PMS2 conversions, we are able to discover most of the unequal reciprocal swap with short-read sequencing only. Validation with HPRC v2.0 data shows high sensitivity and specificity with this approach. 

- **Long-read sequencing**: Phase-resolved long reads can directly observe haplotype-specific PSV patterns, revealing the reciprocal nature of the exchange
- **Trio-based phasing**: Parental samples enable computational phasing to assign alleles to specific haplotypes
- **Orthogonal confirmation**: Clinical testing may include MLPA or targeted long-range PCR to resolve ambiguous cases

**Current approach**: The caller reports these cases as "Normal (No Conversion)" since no unambiguous signal is present. When clinical suspicion is high but short-read results are negative, the interpretation notes recommend consideration of long-read or trio-based follow-up testing.

## Haplotype-Aware PSV Identification

### Algorithm Overview

A critical limitation of relying solely on standard reference genomes (such as GRCh38) is that **the reference sequence itself is a composite that may contain historical or fixed gene conversion events.** This means the reference genome might represent a "hybrid" allele at certain positions rather than the pure, ancestral gene sequence. Consequently, using the reference as the sole ground truth for defining Paralog-Specific Variants (PSVs) is unreliable; a site might appear to differ between the gene and pseudogene in the reference, but be identical in the actual population, or vice-versa.

Our solution is a **Haplotype-Aware PSV Identification Algorithm**. Instead of relying on a static reference genome, this algorithm leverages population-scale haplotype data to identify reliable PSVs—genomic sites that statistically differentiate the functional *PMS2* gene from the *PMS2CL* pseudogene across diverse human populations.

The core principle is **Population Consistency**: A valid marker must not only differ between the gene and pseudogene in the reference genome but must reliably differentiate them across hundreds of diverse human haplotypes. This approach effectively filters out unreliable sites caused by reference artifacts or common population polymorphisms (SNPs), isolating stable "Blocks" of sequence that define the gene's functional identity.

### Training Process

The algorithm was trained using the **Human Pangenome Reference Consortium (HPRC) v2.0** dataset. This dataset provides high-quality, fully phased diploid assemblies for over 230 diverse individuals (460+ haplotypes), offering a comprehensive view of human genetic variation.

The training process involved three key steps:

1. **Haplotype Alignment:** Each haplotype from the HPRC assembly was aligned to the GRCh38 reference for both the *PMS2* and *PMS2CL* loci.
2. **Consistency Scoring:** Every genomic position was evaluated for its ability to distinguish the gene from the pseudogene. A "Consistency Score" was calculated for each site, measuring how frequently the site matched the canonical gene allele versus the pseudogene allele across the entire population.
3. **Block Definition:** Sites with high consistency (>95-99%) were aggregated into contiguous "Blocks." Isolated sites were discarded as noise, while clusters of consistent PSVs were defined as high-confidence structural blocks.

### Training Validation

The resulting block structure was rigorously validated using a "Leave-One-Out" strategy and cross-referenced against known biological samples.

* **Mosaic Analysis:** A visualization technique (Mosaic Plot) was used to map the allelic state of every training sample against the defined blocks. This confirmed that the vast majority of the population shows a clean, binary separation between the gene and pseudogene profiles.
* **Candidate Scanning:** The algorithm was tested against samples with known gene conversion events. The block logic successfully identified distinct "signatures" for biological events—such as the classic Exon 13-15 conversion—distinguishing them from random alignment noise or benign polymorphisms.

### Confirmation with Literature

The learned structure identifies haplotypes which confirm the "3' Active / 5' Stable" model described in Lynch Syndrome literature. The algorithm correctly identified that the sequence homology—and thus the risk of gene conversion—is strongest in the **Exon 13-15 region**, while the upstream region of PMS2 remains protected and stable in the general population. These boundaries allow the caller to specifically target pathogenic events while ignoring upstream noise.

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Normal (HG00097)

```yaml
PMS2:
  Copy numbers:
    PMS2: 2
    PMS2CL: 2
  Variants: []
  cn_interpretation: |-
    PMS2: CN=2
    PMS2CL: CN=2
    No gene conversion events detected
  gene_conversions: []
```

### Example 2: Classic Exon 13-15 (or only Exon 13) Conversion (Pathogenic) (HG01346)

```yaml
PMS2:
  Copy numbers:
    PMS2: 2
    PMS2CL: 2
  Variants: []
  cn_interpretation: |-
    PMS2: CN=2
    PMS2CL: CN=2
    Gene Conversion Events: 1
      Event 1: Classic Conversion (Exon 13) (PMS2CL->PMS2)
        Region: chr7:5978823-5979012
        Allele change: +1
        Clinical: Potential loss of MMR function - relevant for Lynch syndrome
        Evidence: B2=0.90
  gene_conversions:
  - clinical_significance: Potential loss of MMR function - relevant for Lynch syndrome
    converted_alleles: 1
    interpretation: Classic Conversion (Exon 13-15) (PMS2CL->PMS2)
    is_fusion_candidate: false
    is_inferred: false
    region: chr7:5978823-5979012
```

### Example 3: Recombinant Allele (3' Extension) (HG00738)

```yaml
PMS2:
  Copy numbers:
    PMS2: 2
    PMS2CL: 2
  Variants: []
  cn_interpretation: |-
    PMS2: CN=2
    PMS2CL: CN=2
    Gene Conversion Events: 1
      Event 1: Recombinant Conversion (3' Extension) (PMS2CL->PMS2)
        Region: chr7:5977708-5980942
        Allele change: +1
        Clinical: Potential loss of MMR function - relevant for Lynch syndrome
        Evidence: B1=1.00, B2=1.00
        Recombinant: Extends to 3' terminus (Exon 14/15)
  gene_conversions:
  - clinical_significance: Potential loss of MMR function - relevant for Lynch syndrome
    converted_alleles: 1
    interpretation: Recombinant Conversion (3' Extension) (PMS2CL->PMS2)
    is_fusion_candidate: true
    is_inferred: false
    region: chr7:5977708-5980942
```

### Example 4: Benign Gain Event (PMS2→PMS2CL) (HG00099)

```yaml
PMS2:
  Copy numbers:
    PMS2: 2
    PMS2CL: 2
  Variants: []
  cn_interpretation: |-
    PMS2: CN=2
    PMS2CL: CN=2
    Gene Conversion Events: 1
      Event 1: Recombinant Extended Conversion (Exon 11-15 + 3' Ext) (PMS2->PMS2CL)
        Region: chr7:6736904-6747083
        PMS2 Reference: chr7:5977708-5987790
        Allele change: -1
        Clinical: Functional gain in pseudogene - likely benign
        Evidence: B1=1.00, B2=1.00, B3=0.69, B4=1.00
        Recombinant: Extends to 3' terminus (Exon 14/15)
  gene_conversions:
  - clinical_significance: Functional gain in pseudogene - likely benign
    converted_alleles: -1
    interpretation: Recombinant Extended Conversion (Exon 11-15 + 3' Ext) (PMS2->PMS2CL)
    is_fusion_candidate: true
    is_inferred: false
    pms2_reference_region: chr7:5977708-5987790
    region: chr7:6736904-6747083
```

### Example 4: Reciprocal swap (HG01496)

```yaml
PMS2:
  Copy numbers:
    PMS2: 2
    PMS2CL: 2
  Variants: []
  cn_interpretation: |-
    PMS2: CN=2
    PMS2CL: CN=2
    Gene Conversion Events: 2
      Event 1: Recombinant Extended Conversion [Allele Mismatch] (PMS2->PMS2CL) [Fragmented]
        Region: chr7:6736904-6745959
        PMS2 Reference: chr7:5978823-5987790
        Allele change: -2
        Clinical: Functional gain in pseudogene - likely benign
        Evidence: B1=1.00, B2=1.00, B3=0.69, B4=1.00
        Inferred: Reciprocal Swap (Allele Mismatch) in chr7:6745959-6747083
        PMS2 Reference (swap): chr7:5978823-5977708
        Recombinant: Extends to 3' terminus (Exon 14/15)
        Note: Event reconstructed from 2 fragments
      Event 2 (Inferred): Inferred Reciprocal Swap (Exon 13) (PMS2CL->PMS2)
        Region: chr7:5978823-5977708
        Allele change: +1
        Clinical: Pathogenic (Reciprocal Crossover) - PMS2 promoter drives pseudogene
        Note: Derived from artifactual gap pattern
  gene_conversions:
  - clinical_significance: Functional gain in pseudogene - likely benign        
    converted_alleles: -2
    inferred_mechanism: Reciprocal Swap (Allele Mismatch)
    inferred_swap_region: chr7:6745959-6747083
    interpretation: Recombinant Extended Conversion [Allele Mismatch] (PMS2->PMS2CL)
      [Fragmented]
    is_fusion_candidate: true
    is_inferred: false
    pms2_reference_region: chr7:5978823-5987790
    pms2_reference_swap_region: chr7:5978823-5977708
    region: chr7:6736904-6745959
  - clinical_significance: Pathogenic (Reciprocal Crossover) - PMS2 promoter drives
      pseudogene
    conversion_sites: []
    converted_alleles: 1
    inferred_mechanism: Reciprocal Swap (Allele Mismatch)
    inferred_swap_region: chr7:5978823-5977708
    interpretation: Inferred Reciprocal Swap (Exon 13) (PMS2CL->PMS2)
    is_fusion_candidate: false
    is_inferred: true
    region: chr7:5978823-5977708
```

## Interpretation Guide

### Conversion Direction

| Allele Change | Direction | Interpretation | Clinical Significance |
|---------------|-----------|----------------|----------------------|
| **Positive (+1, +2)** | PMS2CL→PMS2 | Pseudogene sequence replacing gene | **Pathogenic**: Loss of MMR function |
| **Negative (-1, -2)** | PMS2→PMS2CL | Gene sequence replacing pseudogene | **Benign**: Functional gain in pseudogene |

### Event Types

| Event Type | Block Coverage | Clinical Notes |
|------------|----------------|----------------|
| **Classic Conversion (Exon 13-15)** | Block 2 (Hotspot) ≥80% | Standard pathogenic conversion event |
| **Recombinant Conversion (3' Extension)** | Extends to Block 1 boundary | Fusion-like event at 3' terminus |
| **No significant event** | Below thresholds | Raw signals below block criteria (likely noise) |

### Copy Number Patterns

| PMS2 | PMS2CL | Interpretation | Clinical Notes |
|------|--------|----------------|----------------|
| 2 | 2 | Normal/Reference | Standard diploid |
| 2 | 1 | Heterozygous PMS2CL deletion | Generally benign |
| 2 | 0 | Homozygous PMS2CL deletion | Generally benign |
| 1 | 2 | Heterozygous PMS2 deletion | **Clinically significant** |
| 0 | 2 | Homozygous PMS2 deletion | **Clinically significant** |

### Liftover Regions

| Region | Coordinates | Description |
|--------|-------------|-------------|
| **PMS2** | chr7:5,965,000-5,989,044, chr7:5,991,653-5,992,796, chr7:5,992,797-5,993,410 | Exons 11-15 + tail, Intron 9, Exon 9 upstream |
| **PMS2CL** | chr7:6,735,644-6,759,796, chr7:6,734,551-6,735,547, chr7:6,733,563-6,734,235 | Corresponding pseudogene regions |

### Fusion Detection

- **3' terminus boundary**: chr7:5,977,708 (Block 1 boundary)
- **Method**: Events extending to/beyond this boundary are flagged as recombinant/fusion candidates
- **Requirement**: Segments must be contiguous (gap < 500 bp) to qualify as recombinant

## Validation

### Validation Procedure

The PMS2 gene conversion caller was validated using the Human Pangenome Reference Consortium (HPRC) v2.0 dataset:

1. **Dataset**: All 226 HPRC v2.0 samples with Illumina short-read BAM files
2. **Caller execution**: Ran segdup-caller on the PMS2 gene for each sample
3. **Ground truth generation**: Aligned diploid assemblies to the PMS2/PMS2CL loci and verified gene/pseudogene assignment using anchor PSV sites
4. **Quality filtering**: Excluded 50+ samples with alignment artifacts (broken or artifactual conversion patterns) from the truth set. 175 samples with high-quality assemblies remain for validation. 
5. **Comparison**: Evaluated segdup-caller gene conversion calls against the assembly-derived truth set

### Result Summary

| Gene   | Normal | Conversion events | TP  | FN | FP | Precision | Recall |
|--------|--------|-------------------|-----|----|----|-----------|--------|
| PMS2   | 145    | 30                | 170 | 4  | 1  | 99.4%     | 97.7%  |
| PMS2CL | 89     | 86                | 172 | 2  | 1  | 99.4%     | 98.9%  |

### Observations

**Conversion frequency**: The observed conversion rate in this cohort (17% for PMS2, 49% for PMS2CL) is consistent with published population-level estimates for PMS2/PMS2CL gene conversion events. The higher rate of PMS2CL conversions (PMS2→PMS2CL direction) compared to PMS2 conversions (PMS2CL→PMS2 direction) suggests that gene-to-pseudogene conversion events are more frequent, which aligns with the expectation that selection pressure preserves functional PMS2 alleles.

**False negatives**: Segdup caller for PMS2 gene correctly identified 80% of the reciprocal gene conversions. The few missed cases of false negatives are primarily due to **balanced reciprocal exchanges**—where equal amounts of sequence are swapped between PMS2 and PMS2CL on homologous chromosomes. These reciprocal events result in no net copy number change and identical PSV allele frequencies (50% gene, 50% pseudogene), making them indistinguishable from the normal diploid state using short-read sequencing alone. Detection of such balanced exchanges would require long-read sequencing or trio-based phasing to resolve haplotype-specific assignments.

## Known Limitations

- **Balanced reciprocal exchanges**: Balanced swaps are limited by short-reads sequencing platform to distinguish 
- **Novel conversions**: Events outside the defined block regions may not be fully characterized
- **Upstream events**: Conversions in the stable anchor region (Block 3) are rare and may need manual review
- **Complex rearrangements**: Multiple overlapping conversion events may be challenging to interpret
- **Low coverage**: Conversion detection sensitivity decreases below 20× coverage
- **Fragmented events**: Disconnected segments are flagged but may require confirmation

## PMS2-Specific Resources

### Clinical Guidelines
- **NCCN Guidelines**: Lynch syndrome genetic testing recommendations
- **ACMG/AMP Guidelines**: Variant interpretation for mismatch repair genes

### Clinical Databases
- **ClinVar PMS2**: https://www.ncbi.nlm.nih.gov/clinvar/?term=PMS2
- **InSiGHT Database**: https://www.insight-group.org/ - Lynch syndrome variant database

### Reference Data
- **Human Pangenome Reference Consortium (HPRC)**: https://humanpangenome.org/ - Population-scale haplotype data used for training and validation

### Literature

1. Liao WW, Asber M, Didion JP, et al. (2023) "A draft human pangenome reference" *Nature* 617:312-324. https://doi.org/10.1038/s41586-023-05896-x
2. Vaughn CP, et al. (2010) "Clinical analysis of PMS2: mutation detection and avoidance of pseudogenes" *Human Mutation* 31(5):588-593
