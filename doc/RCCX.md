# RCCX Module Interpretation

## Overview

The RCCX (Repeat C4-CYP21-TNX) module is a complex tandem repeat on chromosome 6p21.3 within the MHC class III region. Each module contains three gene units:

- **C4** (Complement Component 4): Can be C4A or C4B isotype; can have HERV-K insertion (Long) or deletion (Short)
- **CYP21** (Steroid 21-Hydroxylase): Functional CYP21A2 or pseudogene CYP21A1P
- **TNX** (Tenascin X): Functional TNXB or pseudogene TNXA

The reference genome contains two RCCX modules per haplotype (bimodular), arranged as:

```
5' ── C4A ── CYP21A1P ── TNXA ── RP2 ── C4B ── CYP21A2 ── TNXB ── 3'
      └─ Module 1 (pseudogene) ──┘    └── Module 2 (functional) ──┘
```

Non-Allelic Homologous Recombination (NAHR) between modules causes unequal crossovers that create monomodular (1 module) or expanded (3+ modules) haplotypes, along with chimeric genes relevant for Congenital Adrenal Hyperplasia (CAH) and Ehlers-Danlos Syndrome (CAH-X).

## What This Module Does

The RCCX module performs:

- **RCCX module counting**: Determines total RCCX module count and individual copy numbers for C4A, C4B, CYP21A2, CYP21A1P, TNXB, and TNXA.
- **C4 structural typing**: Resolves C4 Long/Short (HERV-K) status and detects C4A↔C4B isotype switching events.
- **CYP21 chimera detection**: Classifies CYP21A1P/A2 conversion events into 3 tiers (microconversion, boundary conversion, 30kb deletion chimera) for CAH diagnosis.
- **TNX/CAH-X detection**: Identifies TNXA→TNXB chimeric events associated with Ehlers-Danlos Syndrome.
- **Small variant calling**: Calls variants in the RCCX locus using liftover and phasing.
- Validated via population genetics consistency against established epidemiological expectations (HWE, carrier frequencies) on 212 HPRC samples. Assembly-based ground truth is unreliable for RCCX due to systematic collapse in current assemblies.

The RCCX module does **not** perform:

- **Gene conversion detection for C4**: C4A/C4B isotype differences are handled as structural typing rather than gene conversion events.

## Module Count and Haplotype

Module count is derived from total CYP21 copies (CYP21A2 + CYP21A1P), since each RCCX module contains exactly one CYP21 gene.

| Total Modules | Haplotype Description | Frequency |
|---------------|----------------------|-----------|
| 4 | Bimodular/Bimodular (Standard) | ~70% European |
| 3 | Bimodular/Monomodular | ~25% European |
| 2 | Monomodular/Monomodular | Common |
| 5+ | Expanded | Rare |
| 1 | Monomodular/Null | Rare |
| 0 | Null/Null | Very rare |

**Key principle**: Module count is a structural property independent of gene fusions. A 3-module (Bimodular/Monomodular) sample is a common benign variant, NOT evidence of a 30kb deletion chimera.

## C4 Interpretation

### Isotype and Size

C4 isotype (C4A vs C4B) and size (Long vs Short) are two independent dimensions:

- **Isotype**: Determined by exon 26 sequence differences (PCPVLD vs LSPVIH). 
- **Size**: Determined by HERV-K retroviral insertion status. Long form = HERV-K present; Short form = HERV-K deleted.

### HERV-K (Long/Short) Detection

The HERV-K sequences are **identical** between C4A and C4B copies. Short reads cannot determine which specific C4 copy carries the HERV-K deletion — only the **total** number of Short (HERV-K deleted) and Long (HERV-K present) copies is resolved.

The output reports a single `C4_HERV-K` value: negative cn_diff indicates that many copies are Short; the remainder are Long.

**Example**: If `C4A=2`, `C4B=1`, `C4_HERV-K=-1`:
- Total C4 = 3
- Total Short = 1, Long = 2
- We know 2 C4A + 1 C4B and 2 Long + 1 Short, but NOT which specific copy is Short with long read phasing

### Isotype Switching Detection

C4A and C4B differ by ~5 nucleotides in exon 26 (the isotypic region). Gene conversion at this small target can switch a C4A copy to C4B or vice versa — a well-documented and common event within the RCCX locus.

**Detection method**: In standard RCCX modules, C4A pairs with CYP21A1P (module 1) and C4B pairs with CYP21A2 (module 2). CYP21 identity is structurally stable (differences distributed across the entire gene, making wholesale conversion unlikely), so CYP21 copy numbers serve as the structural ground truth for module counts. A mismatch between C4 isotype and its expected CYP21 pairing indicates isotype switching.

| Mismatch | Interpretation |
|----------|---------------|
| 0 | No isotype switching — C4 isotypes match module structure |
| +N | N copies of C4B → C4A conversion (net C4A gain) |
| -N | N copies of C4A → C4B conversion (net C4B gain) |

**Interaction with CYP21 structural chimera**: When a CYP21 Class 3 chimera (30kb deletion) is also detected, the isotype mismatch may reflect the NAHR crossover point (between C4 exon 26 and CYP21) rather than a true isotype switching gene conversion. In this case, the output is annotated as "may be structural".

**Clinical relevance**:
- **C4A deficiency** (C4A→C4B switching): Associated with SLE (systemic lupus erythematosus) risk
- **C4B deficiency** (C4B→C4A switching): Associated with increased infection susceptibility

## NAHR Chimera Events

Two distinct NAHR events produce chimeric genes in the RCCX locus, each removing different intervening genes depending on which homologous pair undergoes unequal crossover.

### CYP21 Chimera (CYP21A1P ↔ CYP21A2 NAHR)

NAHR between CYP21A1P (module 1) and CYP21A2 (module 2) deletes the intervening C4B and TNXA, producing a CYP21A1P/A2 chimeric gene:

```
Normal bimodular haplotype:
  C4A ── CYP21A1P ── TNXA ── C4B ── CYP21A2 ── TNXB

NAHR crossover between CYP21A1P and CYP21A2:
  C4A ── CYP21A1P ── TNXA ── C4B ── CYP21A2 ── TNXB
              ╲____________________╱
                  deleted (~30kb)

Resulting monomodular haplotype:
  C4A ── CYP21A1P/A2 chimera ── TNXB
         ▲
         5' from A1P, 3' from A2
         (promoter edge = fusion junction)
```

**Genes lost**: C4B, TNXA

### TNX Chimera / CAH-X (TNXA ↔ TNXB NAHR)

NAHR between TNXA (module 1) and TNXB (module 2) deletes the intervening C4B and CYP21A2, producing a TNXA/B chimeric gene:

```
Normal bimodular haplotype:
  C4A ── CYP21A1P ── TNXA ── C4B ── CYP21A2 ── TNXB

NAHR crossover between TNXA and TNXB:
  C4A ── CYP21A1P ── TNXA ── C4B ── CYP21A2 ── TNXB
                        ╲____________________╱
                            deleted (~30kb)

Resulting monomodular haplotype:
  C4A ── CYP21A1P ── TNXA/B chimera
                     ▲
                     5' from TNXA, 3' from TNXB
                     (5' TNXB boundary = fusion junction)
```

**Genes lost**: C4B, CYP21A2
**Clinical**: CAH-X syndrome — loss of functional TNXB causes EDS; loss of CYP21A2 causes concurrent CAH

## CYP21 Conversion Interpretation

### 3-Tier Classification

CYP21 conversion events are classified using a 3-tier system based on two independent criteria: **boundary anchoring** and **CN deletion evidence**.

| Class | Name | Boundary | CN Evidence | Clinical Significance |
|-------|------|----------|-------------|----------------------|
| **1** | Microconversion (Internal) | No — entirely within gene body | N/A | Benign: small internal sequence patch |
| **2** | Gene Conversion (Boundary) | Yes — touches 5' promoter edge | No — total modules == baseline | Non-deletional: gene conversion without chromosome breakage |
| **3** | Unequal Crossover (30kb Deletion) | Yes — touches 5' promoter edge | Yes — total modules < baseline | Pathogenic: NAHR-mediated chimeric gene |

### Boundary Detection

Only the **5' (promoter) edge** of CYP21 genes is checked for boundary anchoring. The 3' edge faces TNXB/TNXA and conversions there are handled separately by the TNX interpretation pathway.

The 5' edge is determined dynamically as the edge **farthest from TNXB** (since CYP21 sits immediately upstream of TNX in each module).

### Direction Filter

Only **pathogenic direction** conversions (`converted_alleles > 0`, pseudogene → functional gene) are processed for classification. Benign direction events are skipped.

## TNX / CAH-X Interpretation

### 3-Tier Classification

TNX conversion events use a parallel 3-tier system:

| Class | Name | Boundary | CN Evidence | Clinical Significance |
|-------|------|----------|-------------|----------------------|
| **1** | TNXB Microconversion (Internal) | No — internal patch | N/A | Benign: small sequence patch within TNXB |
| **2** | TNXB Boundary Conversion | Yes — touches 5' TNXB | No — total modules == baseline | Non-deletional: chromosome structurally intact |
| **3** | CAH-X Chimeric Event | Yes — touches 5' TNXB | Yes — total modules < baseline | Pathogenic: EDS risk; evaluate concurrent CAH |

### Boundary Detection

The 5' TNXB boundary (facing CYP21A2) is the fusion junction for TNXA→TNXB NAHR events. The boundary edge is determined as the TNXB coordinate **closest to CYP21A2**.

### Clinical Output

| Status | Trigger | Clinical Note |
|--------|---------|---------------|
| CAH-X Chimeric Event | Class 3 | High risk for EDS; evaluate CYP21A2 for concurrent CAH |
| TNXB Boundary Conversion | Class 2 | TNXA sequence overwrote TNXB 5' boundary; chromosome intact |
| TNXB Haploinsufficiency | TNXB CN < 2, no conversion | Reduced copy number without gene conversion signal |
| Normal | None | No events detected |

### Edge Rule

Only conversions anchored to gene boundaries are fusion candidates. Internal patches — even large ones — are classified as microconversions (Class 1) and not escalated to structural events. This eliminates false chimera calls from benign sequence exchange.

## Output Examples

### Normal Sample (Bimodular/Bimodular) (Sample HG00126, HG00253, HG00423, HG00621, etc)

```yaml
RCCX1:
  Copy numbers:
    C4A: 2
    C4B: 2
    C4_HERV-K: -1
    CYP21A1P: 2
    CYP21A2: 2
    TNXA: 2
    TNXB: 2
  cn_interpretation: |-
    RCCX Module Analysis:
      Total Modules: 4
      Haplotype: Bimodular/Bimodular (Standard)

    C4 (Complement Component 4):
      Total C4 Copies: 4
      Isotypes: 2 C4A, 2 C4B
      Sizes: 3 Long (with HERV-K), 1 Short (HERV-K deleted)

    CYP21 (Steroid 21-Hydroxylase):
      Total CYP21 Copies: 4
      Functional (CYP21A2): 2
      Pseudogene (CYP21A1P): 2
      Chimera Classification: None
      Conversion Events: None

    TNX (Tenascin X):
      Functional TNXB: 2
      CAH-X Classification: None
      Conversion Events: None
```

### Monomodular Sample with CYP21 Microconversion (Sample HG00280)

```yaml
RCCX1:
  Copy numbers:
    C4A: 1
    C4B: 2
    C4_HERV-K: -1
    CYP21A1P: 1
    CYP21A2: 2
    TNXA: 1
    TNXB: 2
  cn_interpretation: |-
    RCCX Module Analysis:
      Total Modules: 3
      Haplotype: Bimodular/Monomodular

    C4 (Complement Component 4):
      Total C4 Copies: 3
      Isotypes: 1 C4A, 2 C4B
      Sizes: 2 Long (with HERV-K), 1 Short (HERV-K deleted)

    CYP21 (Steroid 21-Hydroxylase):
      Total CYP21 Copies: 3
      Functional (CYP21A2): 2
      Pseudogene (CYP21A1P): 1
      Chimera Classification: None
      Microconversions:
        chr6:32039537-32039801 (alleles: +1)

    TNX (Tenascin X):
      Functional TNXB: 2
      CAH-X Classification: None
      Conversion Events: None
```

### C4 Isotype Switching (HG00133)

```yaml
RCCX1:
  Copy numbers:
    C4A: 3
    C4B: 1
    C4_HERV-K: -1
    CYP21A1P: 2
    CYP21A2: 2
    TNXA: 2
    TNXB: 2
  cn_interpretation: |-
    RCCX Module Analysis:
      Total Modules: 4
      Haplotype: Bimodular/Bimodular (Standard)

    C4 (Complement Component 4):
      Total C4 Copies: 4
      Isotypes: 3 C4A, 1 C4B
      Sizes: 3 Long (with HERV-K), 1 Short (HERV-K deleted)
      Isotype Switching: 1 C4B->C4A isotype switching

    CYP21 (Steroid 21-Hydroxylase):
      Total CYP21 Copies: 4
      Functional (CYP21A2): 2
      Pseudogene (CYP21A1P): 2
      Chimera Classification: None
      Conversion Events: None

    TNX (Tenascin X):
      Functional TNXB: 2
      CAH-X Classification: None
      Conversion Events: None
```

### CYP21 Chimera (30kb Deletion)

```yaml
RCCX:
  Copy numbers:
    C4A: 2
    C4B: 0
    CYP21A2: 1
    CYP21A1P: 1
    TNXA: 0
    TNXB: 2
  cn_interpretation: |-
    RCCX Module Analysis:
      Total Modules: 2
      Haplotype: Monomodular/Monomodular

    C4 (Complement Component 4):
      Total C4 Copies: 2
      Isotypes: 2 C4A, 0 C4B
      Sizes: 2 Long (with HERV-K), 0 Short (HERV-K deleted)

    CYP21 (Steroid 21-Hydroxylase):
      Total CYP21 Copies: 2
      Functional (CYP21A2): 1
      Pseudogene (CYP21A1P): 1
      Chimera Classification: Class 3: Unequal Crossover (30kb Deletion)
      Boundary Conversions:
        chr6:32038000-32041153 (alleles: +1)

    TNX (Tenascin X):
      Functional TNXB: 2
      CAH-X Classification: None
      Conversion Events: None
```

### CAH-X Chimeric Event

```yaml
RCCX:
  Copy numbers:
    C4A: 2
    C4B: 1
    CYP21A2: 1
    CYP21A1P: 2
    TNXA: 2
    TNXB: 2
    C4_HERV-K: -1
  cn_interpretation: |-
    RCCX Module Analysis:
      Total Modules: 3
      Haplotype: Bimodular/Monomodular

    C4 (Complement Component 4):
      Total C4 Copies: 3
      Isotypes: 2 C4A, 1 C4B
      Sizes: 2 Long (with HERV-K), 1 Short (HERV-K deleted)

    CYP21 (Steroid 21-Hydroxylase):
      Total CYP21 Copies: 3
      Functional (CYP21A2): 1
      Pseudogene (CYP21A1P): 2
      Chimera Classification: None
      Conversion Events: None

    TNX (Tenascin X):
      Functional TNXB: 2
      CAH-X Classification: Class 3: CAH-X Chimeric Event (Unequal Crossover)
      Clinical Note: Unequal crossover with C4/CYP21A2 CN loss confirmed. High risk for EDS; evaluate CYP21A2 for concurrent CAH.
      Boundary Conversions:
        chr6:32041153-32045326 (alleles: +1)
```

## Validation Methodology

### Why HPRC Assemblies Are Unreliable Ground Truth for RCCX

While the Human Pangenome Reference Consortium (HPRC) provides unprecedented resolution of complex genomic regions, long-read assemblies are fundamentally unreliable as ground truth for the RCCX locus. The locus consists of ~30kb tandemly repeated modules (C4-CYP21-TNX) sharing >98% sequence homology. Even high-fidelity (HiFi) long reads, which typically average 15-20kb, are frequently insufficient to span multiple identical modules and uniquely anchor to flanking unique sequence.

When assembly algorithms construct de Bruijn or string graphs of this region, the lack of spanning reads causes the graph to loop back on itself. Assemblers "collapse" these homologous loops to resolve the path, systematically crushing multi-modular diploid genomes into artificially truncated haplotypes. Intermediate modules — which most frequently contain the pseudogenes CYP21A1P and TNXA — are systematically deleted from generated contigs, leading to massive underestimation of true copy number.

### Evidence of Assembly Collapse

The limitations of HPRC assemblies are evident when examining their internal structural consistency. In the RCCX tandem array, the biological definition of a structural module dictates a strict 1:1 ratio between C4 and CYP21 genes — a constraint enforced by the segdup caller via `cn_constraint_penalty()`.

However, specific HPRC samples reveal physically impossible genomic states. For example, sample HG00146 is annotated with 3 total modules and 3 copies of CYP21A2, yet only 2 total copies of C4 and 0 copies of CYP21A1P. This "orphan" CYP21 without a corresponding C4 violates the fundamental physical architecture of the locus. It is the bioinformatic signature of a graph resolution failure where the assembler artificially merged and dropped sequence blocks. Penalizing the segdup caller against such mathematically broken assemblies guarantees the accumulation of artificial false negatives.

### Validation Framework

Given the unreliability of assembly-based results, validation relies on biologically anchored methodologies.

**Orthogonal wet-lab validation** (MLPA, ddPCR) is the true gold standard for CAH diagnostics — probe hybridization cleanly separates CYP21A2, CYP21A1P, and chimeric fusions without read mapping artifacts. However, such assays are not available in the current validation scope.

The primary validation approach is **population genetics consistency**: the segdup caller's output distribution must recapitulate established epidemiological expectations (e.g., Hardy-Weinberg equilibrium, known carrier frequencies). An algorithm whose output matches population genetics is inherently more trustworthy than assembly results that systematically skew toward collapsed genotypes.

**Mendelian segregation in trios** provides supplementary evidence but is not definitive for RCCX. The two RCCX modules span up to ~100kb — larger than most long reads — and there is no long-range phasing across modules. Mendelian inheritance models can sometimes mathematically accommodate both a collapsed HPRC assembly result and the segdup caller's expanded prediction, so Mendelian concordance is supportive but cannot resolve the ambiguity.

### HPRC Population Statistics

To assess biological plausibility, we aggregated structural counts across the HPRC cohort (212 samples). Samples where the HPRC assembly exhibited discordant module and total C4 counts were filtered out to ensure a standardized baseline.

**Table 1: Overall RCCX Module Count Distribution**

| Module Count | HPRC Assembly (Samples) | Segdup Caller (Samples) |
| :---: | :---: | :---: |
| 2 | 101 | 3 |
| 3 | 72 | 29 |
| 4 | 32 | 149 |
| 5 | 6 | 26 |
| 6 | 1 | 4 |
| 7 | 0 | 1 |

**Table 2: Functional Gene (CYP21A2) vs. Pseudogene (CYP21A1P) Distribution**

| Copy Number | HPRC Assembly CYP21A2 | Segdup Caller CYP21A2 | HPRC Assembly CYP21A1P | Segdup Caller CYP21A1P |
| :---: | :---: | :---: | :---: | :---: |
| 0 | 0 | 0 | 104 | 2 |
| 1 | 0 | 2 | 71 | 39 |
| 2 | 207 | 194 | 30 | 146 |
| 3 | 5 | 15 | 6 | 20 |
| 4 | 0 | 1 | 1 | 5 |

**Table 3: Common HERV-K Structural Combinations (C4 Long/Short)**

*Minor structural states (<1% in both datasets) are omitted for brevity.*

| Structural State (Long / Short) | HPRC Assembly (%) | Segdup Caller (%) | Expected Biological Phenotype |
| :--- | :---: | :---: | :--- |
| 2 Long / 0 Short | 45.75% | 0.94% | Truncated / Collapsed Array |
| 3 Long / 0 Short | 30.66% | 3.77% | Truncated / Collapsed Array |
| 4 Long / 0 Short | 14.15% | 11.79% | "All Long" Canonical State (4 Modules) |
| 2 Long / 2 Short | 0.47% | 34.43% | Global Major Allele (4 Modules) |
| 3 Long / 1 Short | 0.94% | 22.64% | Heterozygous Short Allele (4 Modules) |
| 2 Long / 1 Short | 1.89% | 8.02% | Bimodular/Monomodular (3 Modules) |

### Biological Concordance Analysis

**Macro-Structural Resolution (Module Counts)**

As described above, assembly collapse systematically truncates multi-modular arrays. This artifact is visible in Table 1: the HPRC assembly assigns 4 modules to only 32 samples while assigning 2 or 3 modules to 173 samples. In contrast, the segdup caller predicts 4 modules as the dominant state (149 out of 212 samples, ~70%), closely matching the expected ~60-70% literature frequency for the Bimodular/Bimodular diploid state. Trimodular or expanded states (e.g., 5 modules: 25 samples) appear at frequencies expected in a diverse cohort.

**Gene-Level Resolution (Pseudogene Attrition)**

Table 2 shows that the functional gene CYP21A2 is highly concordant between the HPRC assembly and the segdup caller — both agree on CYP21A2 = 2 for the vast majority of samples (207 vs 194). The discrepancy lies in the pseudogene: the HPRC assembly annotates 104 samples with 0 copies of CYP21A1P, a state that is highly uncommon in healthy populations. This confirms that assembly collapse selectively deletes the interior pseudogene CYP21A1P while preserving the boundary-anchored functional CYP21A2. The segdup caller restores the expected 1:1 module pairing, identifying 2 copies of CYP21A1P in the majority of samples (146).

The segdup caller also identifies 2 samples with heterozygous CYP21A2 deletion (CYP21A2 = 1), consistent with the general population carrier frequency (~1.5%) for Congenital Adrenal Hyperplasia (CAH).

**HERV-K Structural Variant Segregation**

The resolution of the 6.4kb HERV-K endogenous retrovirus insertion (C4 Long vs C4 Short) provides a distinct marker for structural fidelity. Notably, C4 Long copy numbers are **100% concordant** between the HPRC assembly and the segdup caller — both datasets produce identical Long counts across all samples. The discrepancy lies entirely in C4 Short: the HPRC assembly frequently annotates arrays with zero C4 Short copies (representing over 90% of the dataset), which conflicts with established RCCX linkage disequilibrium patterns where the C4B gene is frequently Short. This indicates that assemblers faithfully resolve the Long (HERV-K present) alleles but systematically collapse the Short (HERV-K deleted) alleles, likely because the deleted region reduces the sequence diversity available for graph resolution.

The segdup caller's structural variant distribution aligns closely with known canonical human haplotypes. The segdup caller identifies the 2 Long / 2 Short state as the most common (34.43%), representing individuals heterozygous for the HERV-K insertion, followed by the 3 Long / 1 Short state (22.64%).

Direct long-read evidence corroborates the segdup caller's HERV-K calls. For example, sample HG03540 is annotated as having 0 C4 Short copies in the HPRC assembly, yet inspection of the PacBio HiFi reads for this sample reveals clear 6.4kb HERV-K deletion signatures — confirming that the assembler collapsed the Short allele despite the underlying reads containing the evidence. This demonstrates that the assembly's systematic bias toward "all Long" states is an artifact of graph resolution, not a biological reality.

### Limitations and Future Directions

While the segdup caller demonstrates consistency with population-level epidemiological data and internal structural constraints, this validation approach has inherent limitations:

1. **Phase Limitations**: Read-depth and unlinked HMM transitions cannot definitively phase multi-modular arrays across entire chromosomes. While absolute copy numbers and allele fractions can be quantified, assigning specific microconversions to specific *cis*-alleles remains inferred.
2. **Lack of Absolute Orthogonal Truth**: Statistical concordance is a proxy for accuracy, not a definitive measurement.

**Future Validation**: To establish absolute ground truth, future validation must incorporate orthogonal, targeted wet-lab structural assays. Techniques such as Multiplex Ligation-dependent Probe Amplification (MLPA) and Droplet Digital PCR (ddPCR) quantify absolute copy numbers via probe hybridization rather than read mapping or graph assembly, providing the necessary baseline to unambiguously validate the segdup caller's total module counts and chimeric breakpoints.

### Literature

1. Sekar A, Bialas AR, de Rivera H, et al. (2016) "Schizophrenia risk from complex variation of complement component 4" Nature 530(7589):177-183. https://doi.org/10.1038/nature16549
2. Chung EK, Yang Y, Rennebohm RM, et al. (2002) "Genetic sophistication of human complement components C4A and C4B and RP-C4-CYP21-TNX (RCCX) modules in congenital adrenal hyperplasia" The American Journal of Human Genetics 71(4):823-837. https://doi.org/10.1086/342775
3. Merke DP, Auchus RJ. (2020) "Congenital Adrenal Hyperplasia Due to 21-Hydroxylase Deficiency" New England Journal of Medicine 383(13):1248-1261. https://doi.org/10.1056/NEJMra1909786
4. Chen W, Lin-Su K, Hagman KD, et al. (2015) "TNXB mutations can cause Ehlers-Danlos syndrome in patients with congenital adrenal hyperplasia" The Journal of Clinical Investigation 125(3):1236-1245. https://doi.org/10.1172/JCI79194
5. Liao WW, Asber M, Didion JP, et al. (2023) "A draft human pangenome reference" *Nature* 617:312-324. https://doi.org/10.1038/s41586-023-05896-x
