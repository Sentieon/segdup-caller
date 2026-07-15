# RH (RHD / RHCE) — Rh blood group locus

## Overview

The Rh blood group is the most polymorphic and, after ABO, the most clinically
important blood group system. It is encoded by two paralogous genes on chromosome
1p36.11:

- RHD carries the D antigen. Presence or absence of a functional RHD gene is the basis
  of the RhD-positive / RhD-negative phenotype.
- RHCE carries the C/c and E/e antigens.

Clinical relevance: RhD matching for transfusion, and hemolytic disease of the fetus and
newborn — anti-D alloimmunization in an RhD-negative mother carrying an RhD-positive
fetus. RhCE variants and RHD/RHCE hybrids also cause alloimmunization and complicate
transfusion in chronically transfused patients such as those with sickle cell disease.

## What This Module Does

The RH module performs:

- RHD copy number and RhD zygosity: reports RHD copy number and predicted RhD status
  (0 = RhD-negative, 1 = hemizygous, 2 = homozygous RhD-positive). Whole-gene RHD deletion
  is the common cause of RhD-negative.
- RHCE copy number.
- Weak D and DEL typing: identifies reduced-D RHD alleles from their defining variants
  (curated from the ISBT RHD allele table), reporting the allele name, zygosity, and a
  clinical note — weak D types 1/2/3 are flagged as conventionally managed RhD-positive,
  while other weak-D types and DEL carry a caution. Genotype-level; no read phasing needed.
- RhCE E/e antigen typing: reports E-positive / e-positive from the c.676 determinant.
- Gene-conversion flagging: flags RHD-CE-D hybrid (partial-D) candidates. Preliminary —
  see Validation.
- Small variant calling across RHD and RHCE.

The RH module does not:

- Type RhCE C/c — its determinant lies in a part of RHCE that overlaps the RHD/RHCE
  conversion region, which the caller cannot resolve reliably. E/e is reported.
- Resolve exact partial-D alleles or the RHCE haplotype (Ce / cE / CE / ce) — these need
  read phasing.
- Detect RHD\*Ψ or other alleles defined by insertions/deletions.

Antigen predictions are genotype-based and should be confirmed serologically before
clinical use.

## Genomic structure (GRCh38)

RHD and RHCE arose by duplication and share ~97–98% nucleotide identity across the gene
body. They lie in opposite (inverted) orientation — a tail-to-tail arrangement with
their 3′ ends facing each other, separated by ~30 kb.

| Gene | GRCh38 region | Strand |
|---|---|---|
| RHD  | chr1:25,272,393–25,330,445 | + |
| RHCE | chr1:25,362,249–25,430,192 | − |

RHD is flanked by two ~9 kb near-identical segmental duplications, the "Rhesus boxes."
Unequal recombination between them deletes the entire RHD gene — the molecular basis of
the common RhD-negative haplotype.

## Allelic / structural landscape

- RHD whole-gene deletion → RhD-negative. The dominant RhD-negative mechanism in
  Europeans (deletion allele frequency ~40%; ~15% of Europeans are RhD-negative).
  Reported as RHD copy number (0 = D-negative, 1 = hemizygous, 2 = homozygous D-positive).
- Weak D — single amino-acid substitutions that reduce D antigen density (<1% of
  Europeans). Types 1, 2 and 3 are the common European types and are conventionally
  managed as RhD-positive; other types can form anti-D and are managed more cautiously.
- DEL — very low D expression, e.g. RHD\*1227A (c.1227G>A), common in East Asian
  RhD-negative typings.
- RHD/RHCE gene conversion (RHD-CE-D hybrids) → partial D. Segments of RHD are replaced
  by the corresponding RHCE segments, producing partial-D phenotypes (e.g. DIIIa, DVI).
- RHD\*Ψ (pseudogene) — a 37 bp duplication in exon 4 plus a downstream stop; a common
  RhD-negative background in individuals of African ancestry.

## Output Examples

Output from representative HPRC and GIAB samples. When gene-conversion events are
detected an additional `Gene conversion events:` section is appended; it is preliminary
(see Validation) and is omitted here for clarity.

### RhD-positive, conventional (HG002)

```yaml
RH:
  Copy numbers:
    RHD: 2
    RHCE: 2
  cn_interpretation: |-
    RH (Rh blood group) Analysis:
      RHD:  CN=2
      RHCE: CN=2
      Predicted RhD zygosity: RhD-positive, homozygous — two RHD copies
      RhD antigen markers: none of the curated weak-D/DEL alleles detected — conventional RHD.
      Predicted RhCE E/e antigen: E-e+  (no c.676 E-determinant variant called — E-negative only if the site is covered)
        NOTE: C/c is NOT typed — the RHCE exon-2 C determinant lies in the RHD->RHCE gene-conversion region, a caller blind spot. Exact RHCE haplotype and RHCE hybrids are not resolved here.
```

### RhD-positive, E-positive (HG00423)

```yaml
RH:
  Copy numbers:
    RHD: 2
    RHCE: 2
  cn_interpretation: |-
    RH (Rh blood group) Analysis:
      RHD:  CN=2
      RHCE: CN=2
      Predicted RhD zygosity: RhD-positive, homozygous — two RHD copies
      RhD antigen markers: none of the curated weak-D/DEL alleles detected — conventional RHD.
      Predicted RhCE E/e antigen: E+e+  (genotype Ee; c.676 Ala226Pro, phase-free)
        NOTE: C/c is NOT typed — the RHCE exon-2 C determinant lies in the RHD->RHCE gene-conversion region, a caller blind spot. Exact RHCE haplotype and RHCE hybrids are not resolved here.
```

### RhD-negative (HG00126)

```yaml
RH:
  Copy numbers:
    RHD: 0
    RHCE: 2
  cn_interpretation: |-
    RH (Rh blood group) Analysis:
      RHD:  CN=0
      RHCE: CN=2
      Predicted RhD zygosity: RhD-negative — no intact RHD (homozygous RHD deletion/absence)
      Predicted RhCE E/e antigen: E-e+  (no c.676 E-determinant variant called — E-negative only if the site is covered)
        NOTE: C/c is NOT typed — the RHCE exon-2 C determinant lies in the RHD->RHCE gene-conversion region, a caller blind spot. Exact RHCE haplotype and RHCE hybrids are not resolved here.
```

### Weak D carrier (HG01981)

```yaml
RH:
  Copy numbers:
    RHD: 2
    RHCE: 2
  cn_interpretation: |-
    RH (Rh blood group) Analysis:
      RHD:  CN=2
      RHCE: CN=2
      Predicted RhD zygosity: RhD-positive, homozygous — two RHD copies
      RhD antigen alleles detected:
        RHD*01W.3 (weak D type 3) — heterozygous
        => single allele heterozygous at CN2: D is likely expressed from the apparently-normal RHD in trans (carrier).
        NOTE: genotype-based prediction (require-all-variants); confirm serologically. Partial-D hybrids are flagged separately below.
      Predicted RhCE E/e antigen: E-e+  (no c.676 E-determinant variant called — E-negative only if the site is covered)
        NOTE: C/c is NOT typed — the RHCE exon-2 C determinant lies in the RHD->RHCE gene-conversion region, a caller blind spot. Exact RHCE haplotype and RHCE hybrids are not resolved here.
```

## Validation

Validated on a 226-sample HPRC cohort (GRCh38) using assembly-based truth, plus the GIAB
trio for small-variant accuracy.

Copy number is the primary output:

| metric | result |
|---|---|
| RHD copy number, exact | 222 / 223 (99.6%) |
| RHCE copy number, exact | 222 / 223 (99.6%) |
| RhD-negative (RHD CN0) detection | sensitivity and specificity 1.00 |

The cohort carries real RHD copy-number variation — 12 homozygous deletions, 59 hemizygous
samples, and one RHD triplication — so both deletion and duplication detection are genuinely
tested, not just specificity.

Small variants (GIAB trio, curated truth): with the population short-read model, RHD
SNV/indel calling is strong — F1 ≈ 0.95–0.97 at high precision (HG002, HG003). The result
is cross-platform: Ultima short reads reproduce the same RHD/RHCE copy number and RhD
zygosity and reach RHD F1 ≈ 0.92 (SNVs ≈ 0.93–0.95), with a modest indel-recall gap
consistent with the platform. RHCE calls the well-mappable fraction well (F1 ≈ 0.90) but is
incomplete over the most divergent part of the RHCE gene body; complete RHCE genotyping,
and with it C/c typing, is therefore limited.

Antigen typing: across the 226 samples the weak-D and DEL call rate matches population
expectation — 3 weak-D and 1 DEL carrier, with no evidence of over-calling — and the RhCE
E-positive frequency of 27% matches the expected population frequency. These are
consistency checks; the antigen calls have not been validated against serology.

Gene conversion is reported but preliminary. Because RHD and RHCE are ~97% identical, the
detector over-calls short normal tracts; it is tuned for specificity and is not
benchmarked against hybrid truth. Treat conversion output and partial-D flags as
candidates.

Cross-build: the RHD and RHCE reference sequences are identical between GRCh38 and GRCh37,
so the hg19 and b37 results reproduce the hg38 result.

Known gaps: RhCE C/c, partial-D resolution, RHD\*Ψ, complete RHCE small-variant genotyping
in the divergent core, and copy-number changes beyond RHD deletion.

## References

- Wagner FF, Flegel WA. *RHD gene deletion occurred in the Rhesus box.* Blood 2000.
- Flegel WA. *Molecular genetics and clinical applications for RH.* Transfus Apher Sci.
- Chou ST et al. (bloodTyper) and Lane WJ et al. — WGS-based RH genotyping.
