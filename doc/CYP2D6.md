# CYP2D6 Gene Caller

> **See [general.md](general.md) for overall methodology, output file formats, validation approach, and usage instructions.**

## Overview

CYP2D6 (cytochrome P450 2D6) is one of the most clinically important pharmacogenes, metabolizing approximately 25% of all commonly prescribed medications. The CYP2D6/CYP2D7 gene cluster on chromosome 22q13.2 (chr22:42,123,192-42,145,873) contains the functional CYP2D6 gene and its non-functional paralog CYP2D7 with ~93% sequence identity.

**Unique features of CYP2D6 analysis**:
- **Star allele calling**: Translates genotype into pharmacogenomic nomenclature (*1, *2, *3, *4, *5, etc.)
- **Phenotype prediction**: Predicts metabolizer status (ultrarapid, normal, intermediate, poor)
- **Hybrid allele detection**: Identifies structural variants *13, *36, *68 (non-functional hybrid genes)
- **Duplication-aware calling**: Handles CYP2D6xN alleles common in some populations (1-5%)

## CYP2D6 Gene Cluster Structure

### Genomic Organization

| Gene | Analysis region | Description |
|------|-------------|-------------|
| **CYP2D6** | chr22:42,123,192-42,132,193 | Functional drug-metabolizing enzyme (9 exons) |
| **CYP2D7** | chr22:42,135,344-42,145,873 | Non-functional pseudogene (~93% identical to CYP2D6) |

### Recombination Hotspots

Gene conversions and fusions occur at three major hotspots:

1. **Exon 9** (chr22:42,126,500-42,126,800)
   - Associated with **CYP2D6*36** (exon 9 conversion → frameshift)
   - Creates non-functional allele

2. **Intron 7** (chr22:42,126,800-42,128,000)
   - Major recombination hotspot
   - Site of many internal gene conversions

3. **Intron 1** (chr22:42,129,700-42,130,900)
   - Major recombination hotspot
   - Breakpoint region for **CYP2D6*13** and **CYP2D6*68** hybrid alleles

### Structural Variants (Hybrid Alleles)

| Star Allele | Structure | Frequency | Function | Clinical Impact |
|-------------|-----------|-----------|----------|-----------------|
| **CYP2D6*5** | Complete gene deletion | 2-7% | No function | Poor metabolizer when homozygous/compound het |
| **CYP2D6*13** | CYP2D6-CYP2D7 hybrid (D6 5' + D7 3') | 1-2% | No function | Poor metabolizer |
| **CYP2D6*68** | CYP2D7-CYP2D6 hybrid (D7 5' + D6 3') | <1% | No function | Poor metabolizer |
| **CYP2D6*36** | Exon 9 conversion from CYP2D7 | <1% | No function (frameshift) | Poor metabolizer |
| **CYP2D6xN** | Gene duplication/multiplication | 1-5% | Variable (depends on which allele duplicated) | Can cause ultrarapid metabolism |

### Copy Number Variations

**CYP2D6 deletions** (CYP2D6*5):
- Heterozygous: 2-7% (varies by ethnicity)
- Homozygous: Rare, results in poor metabolizer phenotype

**CYP2D6 duplications/multiplications**:
- 1-5% frequency (higher in some populations)
- Designated as CYP2D6xN (e.g., *1x2, *2x3)
- Functional allele duplications → ultrarapid metabolizer
- Non-functional allele duplications → no phenotype change

**CYP2D7 copy number**: Usually CN=2, variations are typically benign

## Star Allele Calling

### Methodology

The caller uses a **phased consumption algorithm** with tiered matching:

**Tier 1 - Structural variants**: *5 (deletion), *13/*36/*68 (hybrids)
**Tier 2 - No function**: *3, *4, *6, *7, *8, *11, *12, *15, etc. (frameshifts, stop codons, splicing defects)
**Tier 3 - Reduced function**: *9, *10, *14, *17, *29, *41, *54, *55, *59
**Tier 4 - Uncertain function**: *25, *26, *28, *71, *72, *73, *82, *84, *85
**Tier 5 - Normal function**: *1 (wild-type), *2, *35

### Backbone Variants

Some alleles require specific "backbone" variants in cis:

| Star Allele | Defining Variant | Required Backbone | Function |
|-------------|------------------|-------------------|----------|
| **CYP2D6*41** | rs28371725 (intron 6 splicing) | rs16947 (R296C) | Reduced |
| **CYP2D6*29** | rs59421388 (V338M) | *41 backbone | Reduced |
| **CYP2D6*35** | rs769258 (splicing) | *2 backbone | Normal |

The caller verifies these variants are phased together before assigning the allele.

### Copy Number Integration

- **CN < 2**: Automatically assigns *5 (deletion) for missing copies
- **CN = 2**: Standard diploid calling (two alleles)
- **CN ≥ 3**: Calls alleles for all copies, adds note: "Duplication present - allele assignment may be ambiguous"

### Duplication Handling

For CYP2D6xN alleles (CN ≥ 3):
- Variant dosage helps assign alleles (0/1 = one copy has variant; 1/1 = multiple copies)
- Note added about ambiguity in allele-to-copy assignment
- Long-read sequencing or targeted methods may be needed for definitive assignment

## Example Outputs

> **Note**: Clinical significance and interpretation information in the output examples below is for **Educational and Informational purposes only**.

### Example 1: Simple Heterozygous Deletion (*1/*5)

```yaml
CYP2D6:
  Copy numbers:
    CYP2D6: 1
    CYP2D7: 2
  Variants: NA18945.CYP2D6.result.vcf.gz
  cn_interpretation: |-
    CYP2D6: CN=1
    CYP2D7: CN=2
    Structural variants detected: 1
      Event 1: Conversion (CYP2D6→CYP2D7)
        Region: chr22:42129265-42129730
        Allele change: -2
        Clinical: May affect enzyme function
    Star Alleles: *1 / *5
    Phenotype: Intermediate metabolizer (mix of functional and impaired alleles)
  star_alleles:
    alleles: ['*1', '*5']
    copy_number: 1
    method: phased_consumption
```

### Example 2: Hybrid Allele with Duplication (*68/*3/*4)

```yaml
CYP2D6:
  Copy numbers:
    CYP2D6: 2
    CYP2D7: 2
  Variants: NA12878.CYP2D6.result.vcf.gz
  cn_interpretation: |-
    CYP2D6: CN=2
    CYP2D7: CN=2
    Structural variants detected: 2
      Event 1: *68-like (5' D7-D6 fusion)
        Region: chr22:42129756-42131695
        Allele change: -1
        Clinical: Non-functional allele, poor metabolizer
      Event 2: Conversion (CYP2D6→CYP2D7)
        Region: chr22:42132059-42132064
        Allele change: -2
        Clinical: May affect enzyme function
    Star Alleles: *68 / *3 / *4
      Note: Duplication present - allele assignment may be ambiguous
    Phenotype: Poor metabolizer (no functional alleles)
  gene_conversions:
  - converted_alleles: -1
    fusion_type: 5_PRIME
    is_fusion_candidate: true
    region: chr22:42129756-42131695
    star_allele: '*68-like (5'' D7-D6 fusion)'
  star_alleles:
    alleles: ['*68', '*3', '*4']
    copy_number: 3
    method: phased_consumption
    note: Duplication present - allele assignment may be ambiguous
```

### Example 3: Normal Diploid (*2/*2)

```yaml
CYP2D6:
  Copy numbers:
    CYP2D6: 2
    CYP2D7: 2
  Variants: NA18952.CYP2D6.result.vcf.gz
  cn_interpretation: |-
    CYP2D6: CN=2
    CYP2D7: CN=2
    Structural variants detected: 1
      Event 1: Conversion (CYP2D6→CYP2D7)
        Region: chr22:42129192-42129756
        Allele change: -1
        Signature sites: rs267608305, rs775232205, rs373668214, rs528764348
        Clinical: May affect enzyme function
    Star Alleles: *2 / *2
    Phenotype: Normal metabolizer
  star_alleles:
    alleles: ['*2', '*2']
    copy_number: 2
    method: phased_consumption
```

### Example 4: Poor Metabolizer (*56/*5)

```yaml
CYP2D6:
  Copy numbers:
    CYP2D6: 1
    CYP2D7: 2
  Variants: HG03225.CYP2D6.result.vcf.gz
  cn_interpretation: |-
    CYP2D6: CN=1
    CYP2D7: CN=2
    No structural variants detected
    Star Alleles: *56 / *5
    Phenotype: Poor metabolizer (no functional alleles)
  star_alleles:
    alleles: ['*56', '*5']
    copy_number: 1
    method: phased_consumption
```

## Gene-Specific Implementation Details

### Copy Number Prior

CYP2D6 uses **Gamma distribution** priors to allow for common duplications:
- **Shape** (α): 3.0
- **Scale** (β): 0.67
- **Mean**: α×β = 2.0 (centered at diploid)
- **Asymmetry**: Allows for higher CN states with appropriate probability

### Conversion Detection Parameters

Specialized HMM transition rates:
- **Background regions**: 1e-8 (very sticky in stable exonic regions)
- **Recombination hotspots**: 0.001 (loose in known hotspot regions)

### Fusion Detection

- **5' terminus**: Position 42,130,865
- **3' terminus**: Position 42,126,499
- **Method**: Identifies conversions extending to gene termini

## Validation Status

**Note**: CYP2D6 validation is ongoing with limited samples to date.

### Current Validation
- Copy number calling: High concordance with orthogonal methods (qPCR, MLPA) in preliminary testing
- Structural variant detection: Successfully detects *13, *36, *68 hybrid alleles
- Star allele calling: Concordance with reference methods being evaluated

### Reference Materials
- GeT-RM (Genetic Testing Reference Materials) characterized samples
- CPIC and PharmVar reference datasets

## Known Limitations

- **Star allele database**: Covers major alleles but may not include all rare variants
- **Novel variants**: Variants not in PharmVar database will not receive star allele designation
- **Duplication phasing**: Specific allele-to-copy assignment ambiguous for CN≥3 without long-read data
- **Rare structural variants**: Some complex SVs may not be fully characterized

## CYP2D6-Specific Resources

### Star Allele Database
- **PharmVar CYP2D6**: https://www.pharmvar.org/gene/CYP2D6
  - Comprehensive star allele definitions
  - Updated regularly with new variants

### Literature

1. Gaedigk A, Turner A, Everts RE, et al. (2019) "Characterization of Reference Materials for Genetic Testing of CYP2D6 Alleles: A GeT-RM Collaborative Project" *Journal of Molecular Diagnostics* 21(6):1034-1052. https://doi.org/10.1016/j.jmoldx.2019.06.007

2. Sentieon DNAscope: https://www.sentieon.com/products/
