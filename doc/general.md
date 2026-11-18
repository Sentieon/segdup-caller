# Sentieon Segdup Caller: General Overview

## Introduction

The **Sentieon Segdup Caller** (also known as **segdup-caller**) is a specialized variant caller designed for genes located in segmental duplication regions that contain highly similar homologs or pseudogenes. These regions pose significant challenges for standard variant callers due to read mapping ambiguity and copy number variations.

The tool processes short-read whole-genome sequencing (WGS) data, optionally together with long reads (PacBio HiFi or Oxford Nanopore), to provide:
- **Accurate copy number determination** for each gene/region
- **Gene conversion and fusion detection** (for applicable genes)
- **Ploidy-aware variant calling** adjusted for local copy number
- **Star alleles and metabolizer phenotypes** for the pharmacogenes CYP2D6, CYP2B6 and CYP2A6
- **RhD zygosity and RH antigen typing**: weak D, DEL and the E/e antigen
- **LPA KIV-2 copy number**: a VNTR count for a gene with no paralog

### Research Use Only Disclaimer

**THIS TOOL IS FOR RESEARCH USE ONLY AND IS NOT A VALIDATED CLINICAL TEST.**

The Sentieon Segdup Gene Caller has not been validated for clinical diagnostic purposes. Clinical notes and significance information provided in the output YAML files are for **Educational and Informational purposes only** and should not be used as the sole basis for clinical decision-making.

## Supported Genes

The caller supports the following regions. The region name is the name that `--genes` (`-g`) takes.

| Region | Genes Encoded | Clinical Application | Key Outputs |
|--------|---------------|----------------------|-------------|
| **SMN1** | SMN1, SMN2 | Spinal muscular atrophy | CN, variants |
| **PMS2** | PMS2, PMS2CL | Lynch syndrome | CN, variants, gene conversion |
| **CYP2D6** | CYP2D6, CYP2D7 | Pharmacogenomics (drug metabolism) | CN, variants, gene conversion, star alleles |
| **CYP2B6** | CYP2B6, CYP2B7P | Pharmacogenomics (efavirenz, nevirapine, bupropion, methadone) | CN, variants, gene conversion, star alleles |
| **CYP2A6** | CYP2A6, CYP2A7 | Pharmacogenomics (nicotine metabolism, research grouping) | CN, variants, gene conversion, star alleles (*4 deletion, nicotine metabolism) |
| **GBA** | GBA1, GBAP1 | Gaucher disease, Parkinson's risk | CN, variants, gene conversion |
| **STRC** | STRC, STRCP1 | Hearing loss (DFNB16) | CN, variants, deletion/carrier (DFNB16), CATSPER2 extent |
| **NCF1** | NCF1, NCF1B | Chronic granulomatous disease | CN, variants |
| **CFH** | CFH, CFHR1, CFHR2, CFHR3, CFHR4 | Atypical HUS, AMD | CN, variants |
| **CYP11B1** | CYP11B1, CYP11B2 | Aldosterone disorders | CN, variants, gene conversion |
| **HBA** | HBA1, HBA2 | Alpha-thalassemia | CN, variants |
| **RCCX** | RCCX1 (CYP21A2), RCCX2 (CYP21A1P) | Congenital adrenal hyperplasia | CN, variants, gene conversion |
| **IKBKG** | IKBKG (NEMO), IKBKGP1 | Incontinentia pigmenti | CN, gene conversion (5' region only) |
| **RH** | RHD, RHCE | Rh blood group (RhD matching, hemolytic disease of the fetus and newborn) | CN (RhD zygosity), variants, gene conversion |
| **SBDS** | SBDS, SBDSP1 | Shwachman-Diamond syndrome | CN, variants, gene conversion (SDS carrier/affected) |
| **OTOA** | OTOA, OTOAP1 | Hearing loss (DFNB22) | CN, variants, gene conversion (DFNB22 deletion/carrier) |
| **LPA** | LPA (KIV-2 VNTR) | Lp(a) / ASCVD risk (qualitative) | KIV-2 copy number: total + allele split |

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
- Duplicate-flagged reads are excluded from both the direct depth probes and the lifted (collapsed) reads. LPA is the exception: it counts duplicates in both its unit and its flank depths
- Each gene has a direct depth constant (`dp_norm`) and may set a separate collapsed-frame constant (`dp_norm_liftover`) in `genes.yaml`

### Phase 2: Copy Number Determination

**Purpose**: Determine the copy number (CN) state for each gene or genomic region.

1. **Read depth normalization**: Calculates normalized read depth across all target regions
2. **Maximum a posteriori (MAP) optimization**: A hill-climbing algorithm finds the CN state maximizing `log P(data | CN) + log P(CN)`, combining the data likelihood (read depth and PSV allele balance) with a gene-specific prior:
   - **Gaussian priors**: For genes with tight CN distribution around diploid (e.g., PMS2)
   - **Categorical priors**: For genes with discrete, well-characterized CN states, including region-specific distributions (e.g., CFH cluster, GBA1/GBAP1)
   - **Gamma priors**: For genes with asymmetric CN distributions (e.g., CYP2D6 duplications)
3. **CN state output**: Reports most likely copy number (0, 1, 2, 3, 4+) for each region

LPA has no paralog. Its KIV-2 total comes from the depth ratio of the collapsed repeat unit to the unique flanks, and its allele split from a 50-marker panel matched against a reference library of KIV-2 alleles (see [LPA.md](LPA.md)).

### Phase 3: Small Variant Calling

**Purpose**: Call SNVs and indels with appropriate ploidy for each region.

1. **Ploidy-aware calling**: Uses Sentieon DNAscope with ploidy adjusted by copy number:
   - **CN=0**: No variant calling (region deleted)
   - **CN=1**: Haploid calling mode
   - **CN=2**: Standard diploid calling
   - **CN=3+**: Polyploid calling with appropriate parameters
2. **Phasing**: The caller phases with its own phasers and picks the method automatically. Haploid regions need no phasing; short reads are phased by read clustering (k-modes), and long reads by an anchor-based phaser. Phasing adds phase only and keeps the called genotypes.
3. **VCF generation**: Produces high-quality VCF files with genotypes reflecting local ploidy
4. **Quality metrics**: Standard VCF fields (GQ, DP, AD, etc.) for downstream filtering

IKBKG and LPA have no small-variant calling. With a long-read BAM (`--long`), the caller also calls variants from the long reads and combines them with the short-read calls in the same per-gene VCF; copy number comes from the short reads.

### Phase 4: Gene Conversion and Fusion Detection

**Purpose**: Identify structural rearrangements between paralogs (for applicable genes).

**Applicable to**: PMS2, CYP2D6, CYP2B6, CYP2A6, GBA, CYP11B1, RCCX, IKBKG (5' region only), RH, SBDS and OTOA.

1. **Variant signature extraction**: Identifies variants matching paralog sequence (e.g., GBAP1 sequence in GBA1)
2. **Statistical segmentation**: Uses Hidden Markov Model (HMM) with specialized transition rates:
   - **Background regions**: Very low transition rate (e.g., 1e-8) for stable exonic regions
   - **Recombination hotspots**: Higher transition rate (e.g., 0.001) in known hotspot regions
3. **Fusion detection**: Strand-aware detection of terminal fusions:
   - **5' fusions**: Conversions extending to 5' gene terminus
   - **3' fusions**: Conversions extending to 3' gene terminus
4. **Database matching**: Matches detected events against a table of known conversions and fusions, for example:
   - GBA: RecNciI, RecTL, GBA1::GBAP1_Fusion
   - CYP2D6: *36 (exon 9 conversion)
   - CYP2B6: *29 hybrid; CYP2A6: *12 and *34 hybrids
   - SBDS: c.183_184delinsCT and c.258+2T>C; OTOA: c.2359G>T
5. **Quality scoring**: Reports match quality and confidence metrics

### Phase 5: Star Allele Calling (Pharmacogenes)

**Purpose**: Translate genotype into pharmacogenomic nomenclature.

**Applicable to**: CYP2D6, CYP2B6 and CYP2A6.

CYP2D6 uses tiered matching:

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

CYP2B6 and CYP2A6 use a stricter matcher, because many of their alleles are defined by two or more variants together. A star allele matches only when all of its defining variants are present, and in cis when the variants are phased; the most specific matching allele is reported.

RH antigen typing uses the same require-all-variants matching against an allele table curated from the ISBT RHD allele definitions, to type weak D and DEL alleles. The E/e antigen is read at its determinant site (c.676).

### Phase 6: Phenotype Prediction (Pharmacogenes)

**Purpose**: Predict functional phenotype based on star alleles, and RH antigen status.

**Applicable to**: CYP2D6, CYP2B6, CYP2A6 and RH.

- **CYP2D6** metabolizer phenotypes:
  - **Ultrarapid metabolizer (UM)**: Functional allele duplications
  - **Normal metabolizer (NM)**: Two functional alleles
  - **Intermediate metabolizer (IM)**: Mix of functional and reduced/non-functional alleles
  - **Poor metabolizer (PM)**: No functional alleles
- **CYP2B6** metabolizer phenotypes from the CPIC efavirenz diplotype-to-phenotype lookup: ultrarapid, rapid, normal, intermediate and poor metabolizer.
- **CYP2A6** nicotine-metabolism grouping from an activity score: ultrarapid, rapid, normal, intermediate, slow and poor metabolizer. It is a research grouping; CYP2A6 has no CPIC guideline.
- **RH**: predicted RhD zygosity from the RHD copy number, notes for weak D and DEL alleles, and the E/e antigen.

## Output Files

The Segdup Caller generates consistent output files across all genes:

### 1. Summary YAML File

**Filename**: `{outdir}/{sample}.yaml`

**Structure**: Single YAML file containing results for **all genes** analyzed in the run.

**Contents**: one section per region, named after the region's first gene (for example `GBA1` for GBA, `RCCX1` for RCCX and `RHD` for RH), with:
- Copy number calls for each gene/region
- Gene conversion/fusion events (if applicable)
- Star alleles and phenotype predictions (for pharmacogenes)
- Clinical interpretation summary
- Path to the region's VCF file

The file also has a `config` section (run settings, software versions, resource usage and run date) and, with `--merge_vcf`, a `merged_vcf` entry.

> **Note**: Clinical significance and interpretation information in the YAML output is for **Educational and Informational purposes only**.

**Example structure** (an excerpt of HG002's output with `--merge_vcf`, with paths shown for `-o output`; other sections and some fields are left out):
```yaml
CYP2D6:
  Copy numbers:
    CYP2D6: 2
    CYP2D7: 2
  Variants: output/HG002.CYP2D6.result.vcf.gz
  cn_interpretation: |-
    CYP2D6: CN=2
    CYP2D7: CN=2
    No structural variants detected
    Star Alleles: *4 / *2
    Phenotype: Intermediate metabolizer (mix of functional and impaired alleles)
  gene_conversions: []
  star_alleles:
    alleles:
    - '*4'
    - '*2'
    copy_number: 2
    method: phased_consumption
HBA:
  Copy numbers:
    a3.7: 2
    a4.2: 2
    non-duplication: 2
  Variants: output/HG002.HBA.result.vcf.gz
  cn_interpretation: |-
    non-duplication: Normal copy number (CN=2)
    -α3.7: Normal copy number (CN=2)
    -α4.2: Normal copy number (CN=2)
    Clinical interpretation: Normal/Reference (αα/αα) - Normal
LPA:
  Copy numbers:
    LPA_KIV2_allele1: 21
    LPA_KIV2_allele2: 12
    LPA_KIV2_total: 33
  Variants: ''
  interpretation: |-
    LPA (Lp(a) / KIV-2 VNTR) analysis:
      KIV-2 total copies (both alleles): 33
      allele split: 21 + 12 copies
      (smaller allele drives plasma Lp(a): fewer KIV-2 copies -> higher Lp(a) -> higher ASCVD risk)
      Qualitative Lp(a)/ASCVD association only (RUO); no quantitative Lp(a) value.
  kiv2_total_units: 33
config:
  job_options:        # the run's command-line settings
  resource_usage:     # memory, CPU and wall time
  run_date:
  software_versions:  # segdup-caller, Sentieon, samtools, bcftools
merged_vcf: output/HG002.merged.vcf.gz
```

### 2. VCF Files (Per Gene)

**Filename**: `{outdir}/{sample}.{GENE}.result.vcf.gz`, where `{GENE}` is the region's first gene, as in the YAML section names.

**One VCF per gene**: Each analyzed region produces a separate indexed VCF file, except IKBKG and LPA, which call no small variants. With `--merge_vcf`, the caller also writes `{outdir}/{sample}.merged.vcf.gz`, which concatenates the per-gene VCFs and tags each record with `INFO/GENE` for its source gene.

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
- **Phasing**: Variants are phased by the caller's own phasers (k-modes for short reads, anchor-based for long reads); haploid regions are not phased
- **Ploidy handling**: Genotypes reflect local copy number (e.g., 0/1/2 for CN=3)
- **Quality filtering**: Can be filtered using standard VCF quality fields

## Validation Methodology

The Sentieon Segdup Caller has been validated using well-characterized reference datasets:

### Validation Cohorts

1. **HPRC and 1000 Genomes samples with phased assemblies**
   - **Samples**: up to 278 samples; the number with usable truth varies by gene (see each gene page)
   - **Truth data**: copy numbers counted from each sample's phased assemblies
   - **Coverage**: Illumina 30× short reads
   - **Use**: Primary benchmark for copy number accuracy

2. **GIAB HG001–HG005**
   - **Truth data**: NIST v4.2.1 small-variant benchmarks
   - **Use**: small-variant accuracy; the HG002–HG004 trio also checks the calls for Mendelian consistency

3. **1000 Genomes Project**
   - **Version**: High coverage release (30×)
   - **Use**: GBA gene conversions, compared with the findings of Tayebi et al. (2025) (see [GBA.md](GBA.md))

### Validation Metrics

- **Copy number accuracy**: Concordance with copy numbers counted from the assemblies
- **Variant calling concordance**: Agreement with the GIAB truth VCFs
- **Gene conversion detection**: Agreement with published findings for GBA, and specificity on the HPRC cohort where the gene page reports it (for example SBDS and OTOA)
- **Star allele concordance**: Agreement with star alleles read from the assemblies, where the gene page reports it (CYP2B6, CYP2A6)

## Usage

### Basic Command

```bash
segdup-caller -s <sample.bam or .cram> -r <reference.fa> --sr_model <short-read model bundle> -o <output_directory>
```

Without `-g`, the caller runs every supported region. To add long reads, pass `-l <long-read BAM> --lr_model <long-read model bundle>`. The [README](../README.md) lists every option and the Sentieon model bundles.

### Example: Multi-Gene Analysis

```bash
segdup-caller -s NA12878.bam -r GRCh38.fa --sr_model <short-read model bundle> \
  -g GBA,CYP2D6,HBA,SMN1,PMS2 -o results/NA12878
```

**Output** (named after the sample's `SM` tag, or `--sample_name`):
- `results/NA12878/NA12878.yaml` - Combined summary for all genes
- `results/NA12878/NA12878.GBA1.result.vcf.gz` - GBA variants
- `results/NA12878/NA12878.CYP2D6.result.vcf.gz` - CYP2D6 variants
- `results/NA12878/NA12878.HBA.result.vcf.gz` - HBA variants
- `results/NA12878/NA12878.SMN1.result.vcf.gz` - SMN1 variants
- `results/NA12878/NA12878.PMS2.result.vcf.gz` - PMS2 variants

## Technical Requirements

- **Reference genome**: GRCh38 (hg38) or GRCh37 (hg19 or b37). The build is detected automatically from the input. NCF1 is available on GRCh38 only.
- **Coverage**: at least 20×; 30× is ideal
- **Short-read sequencing platform**: short-read WGS such as Illumina, Element Biosciences, MGI and others, with the matching Sentieon short-read model
- **Long reads (optional)**: aligned PacBio HiFi or Oxford Nanopore BAM files, with the matching Sentieon long-read model

### Software Dependencies

- **Sentieon Genomics Software** (202503 or later): Commercial variant calling platform
- **samtools** (1.16 or later): BAM file manipulation
- **bcftools** (1.16 or later): VCF file processing
- **Python** (3.9 or later) with packages: pysam, pyyaml, colorlog, packaging, scipy, pandas, vcflib

## Caveats and Limitations

### General Limitations

- **Novel structural variants**: Events not in the database will be reported but may lack interpretation
- **Low coverage**: Performance degrades below 20× coverage
- **Reference genome dependency**: Requires accurate reference sequence for paralogs
- **Phase ambiguity**: Long-range phasing may be incomplete without long-read data
- **Complex rearrangements**: Some rare complex structural variants may not be fully characterized
- **Library preparation**: The depth constants are fitted on HPRC/1000 Genomes Illumina data. On other library preparations, an extra CYP2D6 copy can be missed (see [CYP2D6.md](CYP2D6.md)), and one person's LPA totals differ by up to 8 copies between preparations (see [LPA.md](LPA.md))

### Gene-Specific Limitations

- **Duplication phasing**: For CN≥3, specific allele-to-copy assignment may be ambiguous
- **Fusion breakpoints**: Exact breakpoint positions may vary within recombination regions
- **Novel variants**: Variants not in PharmVar (CYP2D6, CYP2B6, CYP2A6) or clinical databases may lack annotation
- **Compound heterozygotes**: Phasing of distant variants may be incomplete
- **SMN1/SMN2 split**: The split between SMN1 and SMN2 rests on the few sites that tell the paralogs apart.
  - When one copy carries a short SMN2-like tract, the data can fit two splits, and a small depth change can move the call between them.
  - In NA19240 (assembly SMN1 4 / SMN2 0), the caller calls SMN1 3 / SMN2 1 with an SMN2 6.3 kb deletion. The total of 4 copies is correct.
  - On the HPRC assembly cohort, the SMN1 + SMN2 total is right in 201 of 204 samples and the split in 153 of 156.
- **NCF1**: NCF1C, the third NCF1-family locus, is not modelled. Its reads can therefore push NCF1B, or NCF1, to 3 copies.
  - Against assembly truth, the caller's NCF1 and NCF1B match the assemblies' counts in 168 of 214 HPRC samples.
  - Most assemblies (201 of 214) carry two copies each of NCF1, NCF1B and NCF1C; the caller calls 158 of these NCF1 2 / NCF1B 2.
  - Of the 46 misses, 39 are NCF1B 3 calls on that common layout.

## Gene-Specific Documentation

For detailed information about each gene including:
- Clinical significance
- Gene cluster structure
- Common variants and their interpretation
- Gene-specific validation results

See the gene-specific documentation files:
- [CFH.md](CFH.md) - Complement factor H family (atypical HUS, AMD)
- [CYP2A6.md](CYP2A6.md) - Pharmacogenomics: star alleles, *4 deletion, nicotine metabolism
- [CYP2B6.md](CYP2B6.md) - Pharmacogenomics: star alleles and metabolizer phenotype
- [CYP2D6.md](CYP2D6.md) - Pharmacogenomics, drug metabolism
- [GBA.md](GBA.md) - Gaucher disease, Parkinson's risk
- [HBA.md](HBA.md) - Alpha-thalassemia
- [IKBKG.md](IKBKG.md) - Incontinentia pigmenti
- [LPA.md](LPA.md) - KIV-2 copy number, Lp(a) / ASCVD
- [OTOA.md](OTOA.md) - DFNB22 hearing loss
- [PMS2.md](PMS2.md) - Lynch syndrome
- [RCCX.md](RCCX.md) - RCCX modules, congenital adrenal hyperplasia
- [RH.md](RH.md) - Rh blood group: RhD zygosity and antigen typing
- [SBDS.md](SBDS.md) - Shwachman-Diamond syndrome
- [STRC.md](STRC.md) - DFNB16 hearing loss

## References

1. Sentieon DNAscope: https://www.sentieon.com/products/

## See Also

- [Sentieon Genomics](https://www.sentieon.com/)
- [PharmVar Star Allele Definitions](https://www.pharmvar.org/)
- [PharmGKB Pharmacogenomics Database](https://www.pharmgkb.org/)
- [CPIC Clinical Pharmacogenetics Guidelines](https://cpicpgx.org/)
- [HPRC Human Pangenome Reference Consortium](https://humanpangenome.org/)
