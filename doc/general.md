# Sentieon Segdup Caller: General Overview

## Introduction

The **Sentieon Segdup Caller** (also known as **segdup-caller**) is a specialized variant caller designed for genes located in segmental duplication regions that contain highly similar homologs or pseudogenes. These regions pose significant challenges for standard variant callers due to read mapping ambiguity and copy number variations.

The tool processes both short-read and long-read whole-genome sequencing (WGS) data to provide:
- **Accurate copy number determination** for each gene/region
- **Gene conversion and fusion detection** (for applicable genes)
- **Ploidy-aware variant calling** adjusted for local copy number
- **Star allele nomenclature** (for pharmacogenes like CYP2D6)

## Supported Genes

Currently supported gene families include:

| Gene Family | Clinical Application | Key Features |
|-------------|---------------------|--------------|
| **CYP2D6/CYP2D7** | Pharmacogenomics (drug metabolism) | Copy number, conversions, fusions, star alleles |
| **GBA/GBAP1** | Gaucher disease, Parkinson's risk | Gene conversions (RecNciI, RecTL), fusions |
| **HBA1/HBA2** | Alpha-thalassemia | Common deletions (-α3.7, -α4.2) |
| **SMN1/SMN2** | Spinal muscular atrophy | Copy number, known deletions |
| **CFH/CFHR1-4** | Atypical HUS, AMD | Complement factor H family |
| **PMS2/PMS2CL** | Lynch syndrome | DNA mismatch repair |
| **STRC/STRCP1** | Hearing loss | Stereocilin genes |
| **NCF1/NCF1B** | Chronic granulomatous disease | NADPH oxidase |
| **CYP11B1/CYP11B2** | Aldosterone disorders | Aldosterone synthase |

## General Methodology

The Sentieon Segdup Caller employs a multi-phase analysis pipeline optimized for each gene family:

### Phase 1: Liftover and Read Mapping

**Purpose**: Improve read assignment between highly similar paralogous regions.

1. **Coordinate transformation**: Maps reads between gene paralogs (e.g., GBA1 ↔ GBAP1, CYP2D6 ↔ CYP2D7) using known sequence alignments
2. **Liftover BAM generation**: Creates lifted-over BAM files that reassign reads based on sequence similarity
3. **Improved depth estimation**: Enables more accurate read depth calculation for copy number calling
4. **Variant detection enhancement**: Improves ability to distinguish true variants from paralog differences

**Technical details**:
- Handles complex rearrangements and partial overlaps
- Bidirectional mapping ensures comprehensive coverage

### Phase 2: Copy Number Determination

**Purpose**: Determine the copy number (CN) state for each gene or genomic region.

1. **Read depth normalization**: Calculates normalized read depth across all target regions
2. **Maximum likelihood optimization**: Uses sophisticated hill-climbing algorithm to find optimal CN state
3. **Bayesian prior integration**: Incorporates gene-specific priors reflecting population frequencies:
   - **Gaussian priors**: For genes with tight CN distribution around diploid (e.g., GBA1)
   - **Categorical priors**: For genes with known discrete CN states (e.g., SMN1)
   - **Gamma priors**: For genes with asymmetric CN distributions (e.g., CYP2D6 duplications)
4. **Region-specific priors**: Some genes use different priors for different regions (e.g., CFH cluster)
5. **CN state output**: Reports most likely copy number (0, 1, 2, 3, 4+) for each region

### Phase 3: Small Variant Calling

**Purpose**: Call SNVs and indels with appropriate ploidy for each region.

1. **Ploidy-aware calling**: Uses Sentieon DNAscope with ploidy adjusted by copy number:
   - **CN=0**: No variant calling (region deleted)
   - **CN=1**: Haploid calling mode
   - **CN=2**: Standard diploid calling
   - **CN=3+**: Polyploid calling with appropriate parameters
2. **Phasing**: Uses whatshap to phase variants and distinguish paralogs
3. **VCF generation**: Produces high-quality VCF files with genotypes reflecting local ploidy
4. **Quality metrics**: Standard VCF fields (GQ, DP, AD, etc.) for downstream filtering

### Phase 4: Gene Conversion and Fusion Detection

**Purpose**: Identify structural rearrangements between paralogs (for applicable genes).

**Applicable to**: GBA, CYP2D6, and other genes with known conversion events.

1. **Variant signature extraction**: Identifies variants matching paralog sequence (e.g., GBAP1 sequence in GBA1)
2. **Statistical segmentation**: Uses Hidden Markov Model (HMM) with specialized transition rates:
   - **Background regions**: Very low transition rate (e.g., 1e-8) for stable exonic regions
   - **Recombination hotspots**: Higher transition rate (e.g., 0.001) in known hotspot regions
3. **Fusion detection**: Strand-aware detection of terminal fusions:
   - **5' fusions**: Conversions extending to 5' gene terminus
   - **3' fusions**: Conversions extending to 3' gene terminus
4. **Database matching**: Matches detected events against known structural variants:
   - GBA: RecNciI, RecTL, GBA1::GBAP1_Fusion
   - CYP2D6: *13, *36, *68 hybrids
5. **Quality scoring**: Reports match quality and confidence metrics

### Phase 5: Star Allele Calling (Pharmacogenes)

**Purpose**: Translate genotype into pharmacogenomic nomenclature (CYP2D6 currently).

**Applicable to**: CYP2D6 (extensible to other pharmacogenes).

1. **Variant extraction**: Extracts phased variants from VCF with genotype information
2. **Structural allele integration**: Incorporates copy number and conversion data (*5, *13, *68)
3. **Tiered allele matching**: Matches variants to star allele definitions in priority order:
   - **Tier 1**: Structural alleles from CN/conversion analysis
   - **Tier 2**: No-function alleles (frameshifts, stop codons, splicing defects)
   - **Tier 3**: Reduced-function alleles
   - **Tier 4**: Uncertain-function alleles
   - **Tier 5**: Normal-function reference alleles
4. **Backbone handling**: For haplotype-dependent alleles, verifies backbone variants are in cis
5. **Phasing-aware consumption**: Uses phase information to group variants correctly
6. **Copy number integration**: Accounts for duplications and deletions
7. **Default assignment**: Assigns reference allele (*1) for unmatched haplotypes

### Phase 6: Phenotype Prediction (Pharmacogenes)

**Purpose**: Predict functional phenotype based on star alleles.

**Applicable to**: CYP2D6 metabolizer phenotypes (extensible to other pharmacogenes).

- **Ultrarapid metabolizer (UM)**: Functional allele duplications
- **Normal metabolizer (NM)**: Two functional alleles
- **Intermediate metabolizer (IM)**: Mix of functional and reduced/non-functional alleles
- **Poor metabolizer (PM)**: No functional alleles

## Output Files

The Segdup Caller generates consistent output files across all genes:

### 1. Summary YAML File

**Filename**: `{outdir}/{sample}.yaml`

**Structure**: Single YAML file containing results for **all genes** analyzed in the run.

**Contents**:
- Sample identifier
- Copy number calls for each gene/region
- Gene conversion/fusion events (if applicable)
- Star alleles and phenotype predictions (for pharmacogenes)
- Clinical interpretation summary
- Path to gene-specific VCF files

**Example structure** (multiple genes):
```yaml
sample: HG00253

GBA1:
  Copy numbers:
    GBA1: 2
    GBAP1: 1
  Variants: output/HG00253.GBA1.result.vcf.gz
  cn_interpretation: |-
    GBA1: CN=2
    GBAP1: CN=1 (heterozygous deletion)
    Gene conversion/fusion events detected: 1
      Event 1 (CONVERSION): chr1:155235203-155235252
        Converted alleles: 1
        Interpretation: RecNciI
  gene_conversions:
  - region: chr1:155235203-155235252
    converted_alleles: 1
    interpretation:
    - event_name: RecNciI
      match_score: 1.0
CYP2D6:
  Copy numbers:
    CYP2D6: 2
    CYP2D7: 2
  Variants: output/HG00253.CYP2D6.result.vcf.gz
  cn_interpretation: |-
    CYP2D6: CN=2
    CYP2D7: CN=2
    Star Alleles: *1 / *4
    Phenotype: Intermediate metabolizer
  star_alleles:
    alleles:
    - '*1'
    - '*4'
    copy_number: 2
    method: phased_consumption
HBA:
  Copy numbers:
    a3.7: 1
    a4.2: 2
    non-duplication: 2
  Variants: output/HG00253.HBA.result.vcf.gz
  cn_interpretation: |-
    -α3.7: Heterozygous deletion (CN=1)
    -α4.2: Normal copy number (CN=2)
    Clinical interpretation: Heterozygous -α3.7 deletion (αα/-α3.7) - Silent carrier
```

### 2. VCF Files (Per Gene)

**Filename**: `{outdir}/{sample}.{GENE}.result.vcf.gz`

**One VCF per gene**: Each analyzed gene produces a separate indexed VCF file.

**Contents**:
- Phased small variants (SNVs and indels)
- Genotypes adjusted for local copy number
- Phase information distinguishing paralogs (e.g., HBA1 vs HBA2, CYP2D6 vs CYP2D7)
- Standard VCF quality metrics (GQ, DP, AD, AF, etc.)
- Appropriate ploidy representation (haploid, diploid, polyploid)

**Example**: For a sample analyzed for GBA, CYP2D6, and HBA:
- `{sample}.GBA1.result.vcf.gz` - GBA1/GBAP1 variants
- `{sample}.CYP2D6.result.vcf.gz` - CYP2D6/CYP2D7 variants
- `{sample}.HBA.result.vcf.gz` - HBA1/HBA2 variants

**VCF features**:
- **Phasing**: Variants are phased using whatshap
- **Ploidy handling**: Genotypes reflect local copy number (e.g., 0/1/2 for CN=3)
- **Quality filtering**: Can be filtered using standard VCF quality fields

## Validation Methodology

The Sentieon Segdup Caller has been validated using well-characterized reference datasets:

### Validation Cohorts

1. **HPRC (Human Pangenome Reference Consortium)**
   - **Version**: v1.0 and v2.0 releases
   - **Samples**: 215+ diverse samples
   - **Truth data**: Assembled haplotypes provide ground truth for copy number
   - **Coverage**: High-coverage WGS (>30×)
   - **Use**: Primary benchmark for copy number accuracy

2. **1000 Genomes Project**
   - **Version**: High coverage release (30×)
   - **Samples**: Diverse population samples
   - **Truth data**: Orthogonal validation methods, published variant calls
   - **Use**: Validation of gene conversions, star alleles, and small variants

3. **Coriell Cell Repository Reference Materials**
   - **Samples**: GeT-RM characterized reference materials (e.g., for CYP2D6)
   - **Truth data**: Characterized by multiple orthogonal methods
   - **Use**: Star allele calling validation

### Validation Metrics

- **Copy number accuracy**: Concordance with truth calls from assemblies
- **Variant calling concordance**: Agreement with truth VCF calls
- **Gene conversion detection**: Sensitivity and specificity for known conversion events
- **Star allele concordance**: Agreement with reference genotyping methods

## Usage

### Basic Command

```bash
segdup-caller \
  --short <sample.bam> \
  --reference <reference.fa> \
  --sr_model <model_bundle> \
  --genes <gene1,gene2,...> \
  --outdir <output_directory>
```

### Parameters

- `--short`: Input BAM file (short-read WGS, ≥30× recommended)
- `--reference`: Reference genome FASTA (GRCh38/hg38)
- `--sr_model`: Sentieon model bundle (SentieonGeneCaller0.1.bundle)
- `--genes`: Comma-separated list of genes to analyze (or "all")
- `--outdir`: Output directory for results

### Optional Parameters

- `--log_level`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `--keep_temp`: Retain intermediate files for debugging

### Example: Multi-Gene Analysis

```bash
segdup-caller \
  --short NA12878.bam \
  --reference GRCh38.fa \
  --sr_model SentieonGeneCaller0.1.bundle \
  --genes GBA,CYP2D6,HBA,SMN1,PMS2 \
  --outdir results/NA12878/
```

**Output**:
- `results/NA12878/NA12878.yaml` - Combined summary for all genes
- `results/NA12878/NA12878.GBA1.result.vcf.gz` - GBA variants
- `results/NA12878/NA12878.CYP2D6.result.vcf.gz` - CYP2D6 variants
- `results/NA12878/NA12878.HBA.result.vcf.gz` - HBA variants
- `results/NA12878/NA12878.SMN1.result.vcf.gz` - SMN1 variants
- `results/NA12878/NA12878.PMS2.result.vcf.gz` - PMS2 variants

## Technical Requirements

- **Reference genome**: GRCh38/hg38 (coordinates not compatible with hg19/GRCh37)
- **Recommended coverage**: ≥30× for optimal performance
  - Minimum: 20× (reduced accuracy for CNV and conversions)
  - Optimal: 30-50×
  - Higher coverage (>50×) provides marginal improvements
- **Read length**: ≥100 bp (150 bp preferred)
- **Sequencing platform**: Illumina short-read sequencing

### Software Dependencies

- **Sentieon Genomics Software** (v202503+): Commercial variant calling platform
- **samtools** (v1.16+): BAM file manipulation
- **bcftools** (v1.10+): VCF file processing
- **whatshap** (v2.3+): Haplotype phasing
- **Python** (v3.10+) with packages: pysam, pyyaml, colorlog, scipy, pandas, vcflib

## Caveats and Limitations

### General Limitations

- **Novel structural variants**: Events not in the database will be reported but may lack interpretation
- **Low coverage**: Performance degrades below 20× coverage
- **Reference genome dependency**: Requires accurate reference sequence for paralogs
- **Phase ambiguity**: Long-range phasing may be incomplete without long-read data
- **Complex rearrangements**: Some rare complex structural variants may not be fully characterized

### Gene-Specific Limitations

- **Duplication phasing**: For CN≥3, specific allele-to-copy assignment may be ambiguous
- **Fusion breakpoints**: Exact breakpoint positions may vary within recombination regions
- **Novel variants**: Variants not in PharmVar (CYP2D6) or clinical databases may lack annotation
- **Compound heterozygotes**: Phasing of distant variants may be incomplete

## Gene-Specific Documentation

For detailed information about each gene including:
- Clinical significance
- Gene cluster structure
- Common variants and their interpretation
- Gene-specific validation results

See the gene-specific documentation files:
- [GBA.md](GBA.md) - Gaucher disease, Parkinson's risk
- [CYP2D6.md](CYP2D6.md) - Pharmacogenomics, drug metabolism
- [HBA.md](HBA.md) - Alpha-thalassemia
- [CFH.md](CFH.md) - Complement factor H family
- Additional genes: See `/doc` directory

## References

1. Sentieon DNAscope: https://www.sentieon.com/products/

## See Also

- [Sentieon Genomics](https://www.sentieon.com/)
- [PharmGKB Pharmacogenomics Database](https://www.pharmgkb.org/)
- [CPIC Clinical Pharmacogenetics Guidelines](https://cpicpgx.org/)
- [HPRC Human Pangenome Reference Consortium](https://humanpangenome.org/)
