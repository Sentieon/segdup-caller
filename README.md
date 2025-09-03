# Sentieon Segdup Gene Caller

**Segdup Gene caller** is a variant caller for genes located in segmental duplication regions that contain highly similar homologs or pseudogenes. It can leverage available long-read sequencing data to improve variant calling accuracy.

Segdup caller currently supports whole-genome sequencing (WGS) paired-end short reads of 150 bp length, such as those from Illumina, Element Biosciences, MGI, and others. Sentieon provides platform-specific short-read variant calling models. For accurate copy number calling and phasing, Segdup caller recommends a minimum of 20× coverage, with 30× being ideal.

Segdup caller supports additional long-read data from **PacBio HiFi** or **Oxford Nanopore**, aligned BAM files, using the corresponding Sentieon models. 

### Supported Genes

Segdup caller currently supports the following genes and segmental duplication regions:

| Region Name | Genes Encoded                        |
|-------------|--------------------------------------|
| SMN1        | SMN1, SMN2                           |
| PMS2        | PMS2, PMS2CL                         |
| CYP2D6      | CYP2D6, CYP2D7                       |
| GBA         | GBA1, GBAP1                          |
| STRC        | STRC, STRCP1                         |
| NCF1        | NCF1, NCF1B                          |
| CFH         | CFH, CFHR1, CFHR2, CFHR3, CFHR4      |
| CYP11B1     | CYP11B1, CYP11B2                     |
| HBA         | HBA1, HBA2                           |

> **Note:** Segdup caller is expanding the list of supported genes. If there are particular genes you hope to support, or you encounter any issues, please [file an issue](https://github.com/Sentieon/segdup-caller/issues).

---

# Installation

To install the package from Git:

```bash
git clone https://github.com/Sentieon/segdup-caller.git
pip install ./segdup-caller
```

---

## Dependencies

- [Sentieon Genomics software (version 202503 or later)](https://www.sentieon.com/free-trial/)
- [samtools (version 1.16 or later)](http://www.htslib.org/)
- [bcftools (version 1.10 or later)](http://www.htslib.org/)

---

## Sentieon Models

You can find a list of Sentieon models in our [Sentieon models repository](https://github.com/Sentieon/sentieon-models).

---

# Usage

```bash
usage: segdup-caller [-h] --short SHORT [--long LONG] --sr_model SR_MODEL
                     [--lr_model LR_MODEL] --reference REFERENCE --genes GENES
                     [--sample_name SAMPLE_NAME] [--sr_prefix SR_PREFIX]
                     [--lr_prefix LR_PREFIX] [--config CONFIG] --outdir OUTDIR
                     [--version] [--keep_temp]
```

Targeted variant caller for genes with highly similar paralogs.

### Options:

```
  -h, --help                Show this help message and exit

  --short SHORT, -s SHORT   Input short-read BAM or CRAM (required)

  --long LONG, -l LONG      Input long-read BAM or CRAM (optional)

  --sr_model SR_MODEL       Short read model bundle (required)

  --lr_model LR_MODEL       Long read model bundle (required if --long is provided)

  --reference REFERENCE, -r REFERENCE
                            Reference file

  --genes GENES, -g GENES   List of genes to be called (comma separated).
                            Supported genes: CFH, CFHR3, CYP11B1, CYP2D6, GBA,
                            NCF1, PMS2, SMN1, STRC, HBA

  --sample_name SAMPLE_NAME Sample name (default: SM tag in the input short-read
                            BAM file will be used)

  --sr_prefix SR_PREFIX     Short read result prefix

  --lr_prefix LR_PREFIX     Long read result prefix

  --config CONFIG           Custom gene configuration file (advanced users only)

  --outdir OUTDIR, -o OUTDIR
                            Output directory

  --version                 Show program's version number and exit

  --keep_temp               Keep temporary files for debugging

  --threads THREADS, -t THREADS
                            Number of parallel threads for gene processing
                            (default: number of available CPU cores)
```

### Examples

**Check version:**
```bash
segdup-caller --version
```

**Calling SMN1/SMN2 with a short-read input BAM:**
```bash
segdup-caller -s $short_read_bam -r $hg38_REF --sr_model $short_read_model_bundle -g SMN1 -o $outdir 
```

---

# Output

Sentieon Segdup caller produces the following output:

- VCF files containing variant calls for each specified gene region
- A YAML file containing additional information, such as copy number states

