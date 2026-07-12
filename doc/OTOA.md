# OTOA (OTOA / OTOAP1) — DFNB22 nonsyndromic hearing loss locus

## Overview

DFNB22 is an autosomal-recessive, nonsyndromic sensorineural hearing loss caused by biallelic
loss of OTOA, which encodes otoancorin — a protein that anchors the inner-ear acellular gels
(the tectorial and otolithic membranes) to the sensory epithelium. OTOA copy-number loss is a
notable cause of recessive deafness: in a large hearing-loss cohort it was the second most
frequent single-gene cause (about 0.6%, 14 / 2262).

OTOA sits on chromosome 16p12.2 next to a nearly identical partial pseudogene, OTOAP1, about
0.8 Mb distal. The two share the last third of the gene as a ~68 kb, ~99%-identical direct
duplication. That arrangement produces two paralog-driven disease mechanisms that standard
short-read pipelines handle poorly:

- a recurrent ~500 kb deletion of the whole gene, formed by non-allelic homologous
  recombination between the OTOA and OTOAP1 repeats; and
- a gene-conversion allele, in which the OTOAP1 stop codon **p.Glu787\*** (OTOA c.2359G>T) is
  copied into OTOA.

Because OTOAP1 is so similar, ordinary pipelines both **miss the deletion** (OTOAP1 reads
mismap onto OTOA and refill the lost depth) and **mis-call the conversion stop** (the same
p.Glu787\* is a notorious short-read false positive that laboratories confirm with long-range
PCR). This module is paralog-aware and addresses both.

## What This Module Does

The OTOA module performs:

- Copy-number calling: OTOA and OTOAP1 copy number. OTOA copy number is anchored on the
  uniquely-mappable 5′ two-thirds of the gene, so it is measured directly and is not confounded
  by OTOAP1 copy-number variation.
- Recurrent-deletion typing: reports OTOA copy loss and summarizes a recessive-disease call —
  two pathogenic alleles → consistent with DFNB22 (bi-allelic), one → carrier, zero → negative.
- Conversion-variant typing: types the recurrent pathogenic conversion stop p.Glu787\*
  (c.2359G>T) by zygosity, and combines deletion and conversion alleles into the recessive
  call. Because this stop is also a common short-read artifact, an unreliable call is surfaced
  as low-confidence for orthogonal confirmation rather than reported as pathogenic.
- Small-variant calling: paralog-differentiated small variants across OTOA and OTOAP1.

The OTOA module does not:

- Type OTOA pathogenic alleles beyond the recurrent deletion and p.Glu787\* into a named
  disease call. Private point mutations still appear in the standard small-variant calls but
  are not summarized into the carrier/affected call.
- Resolve variants inside the ~15 kb PSV-desert in the middle of the 3′ duplication (5 coding
  exons), which is duplicated but has no reliable paralog-distinguishing sites — a dark region
  flagged for orthogonal confirmation.

## Genomic structure (GRCh38)

| Gene | GRCh38 region | Strand | Role |
|---|---|---|---|
| OTOA   | chr16:21,663,968–21,761,935 | + | functional gene (recipient) |
| OTOAP1 | chr16:22,546,964–22,576,682 | + | pseudogene (donor) |

The two paralogs are in the same (direct) orientation, ~0.8 Mb apart. Only OTOA's 3′ ~33 kb
(from ~chr16:21,729,043) is duplicated with OTOAP1 — a single collinear ~68 kb block at ~99%
identity, containing 94 curated paralog-distinguishing sites (PSVs). The 5′ ~65 kb of OTOA is
single-copy and uniquely mappable; the recurrent deletion removes the whole gene, and its loss
is read cleanly from that unique 5′ depth. Direct repeats recombine to a **deletion** (not an
inversion), so the event is copy-number-detectable.

## Pathogenic mechanisms

**Recurrent deletion (primary).** NAHR between the OTOA and OTOAP1 direct repeats deletes a
contiguous ~500 kb segment of 16p12.2 that removes the whole OTOA gene (and co-deletes
neighboring METTL9 / IGSF6). Biallelic deletion causes DFNB22; a single deletion is a carrier.
This copy-number change is the main clinically actionable signal.

**Gene-conversion stop p.Glu787\* (OTOA c.2359G>T).** The OTOAP1 allele carries a stop at this
codon; conversion copies it into OTOA. In GRCh38 the site is chr16:21,736,318 G>T (donor OTOAP1
chr16:22,552,441; ClinVar 218841). It is genuinely pathogenic when truly in OTOA, but is also a
frequent short-read artifact from OTOAP1 mismapping, so a positive call needs orthogonal
confirmation.

Compound heterozygosity across the two mechanisms (one deletion + one conversion allele) also
produces DFNB22.

## Output Examples

All examples are actual caller output. Copy number and a SUMMARY line are always reported; the
deletion line and the p.Glu787\* call appear when present. Ratio-based conversion tracts, when
shown, are supporting context and are non-diagnostic on their own. OTOAP1 copy number is
polymorphic in the population (commonly 4) and is not clinically actionable.

### Normal, no pathogenic allele (GIAB HG002)

HG002 shows a short conversion tract that is flagged non-diagnostic, so the SUMMARY stays
negative:

```yaml
OTOA:
  Copy numbers:
    OTOA: 2
    OTOAP1: 4
  cn_interpretation: |-
    OTOA (DFNB22 nonsyndromic hearing loss) Analysis:
      OTOA:   CN=2
      OTOAP1: CN=4

    Positive-direction conversion tract(s): 1 — ratio-based context, non-diagnostic alone (allelic dropout / OTOAP1 mismapping); pathogenic only where a tract spans a typed PSV (see above)
      Tract 1: chr16:21735751-21735889, converted_alleles=1

    SUMMARY: no recurrent OTOA deletion or conversion allele detected.
```

### Deletion carrier, one pathogenic allele (GIAB/1000G NA18948 — a real carrier)

```yaml
OTOA:
  Copy numbers:
    OTOA: 1
    OTOAP1: 4
  cn_interpretation: |-
    OTOA (DFNB22 nonsyndromic hearing loss) Analysis:
      OTOA:   CN=1
      OTOAP1: CN=4

    OTOA copy loss: 1 allele(s) deleted — recurrent DFNB22 ~500 kb NAHR deletion of the whole gene.

    SUMMARY: 1 pathogenic OTOA allele (deletion) — carrier (heterozygous).
```

### Conversion carrier, one pathogenic allele (simulated positive control)

```yaml
OTOA:
  Copy numbers:
    OTOA: 2
    OTOAP1: 4
  cn_interpretation: |-
    OTOA (DFNB22 nonsyndromic hearing loss) Analysis:
      OTOA:   CN=2
      OTOAP1: CN=4

    Recurrent OTOAP1->OTOA conversion variant(s) — tract-corroborated:
      mono-allelic OTOA_c.2359G>T (chr16:21736318 G>T; pseudogene-derived stop, conversion tract spans the PSV)

    Positive-direction conversion tract(s): 1 — ratio-based context, non-diagnostic alone (allelic dropout / OTOAP1 mismapping); pathogenic only where a tract spans a typed PSV (see above)
      Tract 1: chr16:21734548-21756274, converted_alleles=1

    SUMMARY: 1 pathogenic OTOA allele (conversion) — carrier (heterozygous).
```

### Affected, homozygous conversion (simulated positive control)

```yaml
OTOA:
  Copy numbers:
    OTOA: 2
    OTOAP1: 4
  cn_interpretation: |-
    OTOA (DFNB22 nonsyndromic hearing loss) Analysis:
      OTOA:   CN=2
      OTOAP1: CN=4

    Recurrent OTOAP1->OTOA conversion variant(s) — tract-corroborated:
      bi-allelic OTOA_c.2359G>T (chr16:21736318 G>T; pseudogene-derived stop, conversion tract spans the PSV)

    SUMMARY: 2 pathogenic OTOA allele(s) (deletion=0, conversion=2) — consistent with DFNB22 hearing loss (bi-allelic). Confirm; DFNB22 is autosomal recessive.
```

## Validation

### Copy number (population cohort)

Validated on 248 population samples (HPRC / 1000 Genomes assembly cohort, GRCh38). Because the
OTOA/OTOAP1 block is ~99% identical, the assemblies collapse there, so copy-number truth is
taken from genome-normalized read depth — OTOA from the uniquely-mappable 5′ region, OTOAP1
from the total 3′-block depth.

| metric | result |
|---|---|
| OTOA copy number (effective 5′, exact) | 248 / 248 (100%) |
| Recurrent-deletion carriers detected | 1 / 1 (NA18948; 0 false deletions) |
| OTOAP1 copy number vs independent block-depth | ~99% (semi-independent cross-check) |

The cohort contains one real OTOA deletion carrier (NA18948, half depth over the unique 5′,
confidently mapped) and one duplication (NA18976); both are called correctly, with no false
deletions anywhere else. OTOAP1 copy number is polymorphic (2–5, predominantly 4) and agrees
with an independent block-depth estimate on ~99% of samples; as with other pseudogenes, its
dosage carries no disease meaning and is reported only for completeness.

### Conversion specificity and sensitivity

No DFNB22-positive samples were available, so specificity is measured on the 248-sample
cohort (all expected negative for conversion) and sensitivity by simulation.

| metric | result |
|---|---|
| Conversion specificity (p.Glu787\*) | 0 / 248 false carriers |

Every cohort sample is correctly negative — including the samples where OTOAP1 mismapping
produces the p.Glu787\* short-read artifact that misleads standard pipelines. Sensitivity was
confirmed by injecting the conversion allele into reads at the OTOA locus:

| simulated genotype | caller call |
|---|---|
| carrier (heterozygous p.Glu787\*) | 1 pathogenic allele — carrier |
| affected (homozygous p.Glu787\*) | 2 alleles — consistent with DFNB22 (bi-allelic) |
| control (no conversion) | negative (ratio tract present but non-diagnostic) |

### Small-variant accuracy (GIAB, curated truth)

Benchmarked against NIST GIAB v4.2.1 within the OTOA gene ∩ high-confidence regions (v5.0q does
not cover the locus), HG001–HG005, population bundle. All five are OTOA CN=2 (diploid,
scoreable); OTOAP1 is a tetraploid pseudogene that GIAB excludes and is validated separately.

| region | result |
|---|---|
| OTOA full gene, overall F1 | ~0.91–0.96 (mean ~0.94) |
| OTOA 5′ unique region (recall) | complete — 0 false negatives across all 5 samples |

Every false negative falls in the 3′-duplicated region — the PSV-desert dark zone and, in the
two most divergent haplotypes (HG004/HG005), a dense cluster — where OTOA and OTOAP1 are not
short-read-separable. The uniquely-mappable 5′ region, which carries the clinically-relevant
point mutations and anchors the deletion call, is called completely. The population model
recovers a substantial fraction of the divergent-haplotype 3′ variants.

### Cross-build reproducibility (hg38 / hg19 / b37)

The OTOA and OTOAP1 reference sequences are byte-identical between GRCh38 and GRCh37, offset by
a constant +11,321 bp, so the GRCh37 (hg19 and b37) support is an exact coordinate lift of the
GRCh38 model (p.Glu787\* → 16:21,747,639). Runs on GIAB HG001–HG005 reproduce the GRCh38 result:
identical copy number (5/5) and concordant small-variant calls.

### Cross-platform (Ultima vs Illumina)

On GIAB HG001–HG005, Ultima and Illumina short reads give **identical** OTOA/OTOAP1 copy number
and identical deletion/conversion status (no false deletion or conversion on either). Small
variants show Ultima's expected, modest recall gap concentrated in homopolymer indels and the
near-paralogous 3′/dark region; in the uniquely-mappable 5′ clinical core Ultima SNV recall is
complete (0 false negatives, all 5 samples), and p.Glu787\* is correctly negative on both
platforms. The clinically actionable outputs are platform-invariant.

### Validation scope and known gaps

Validated (GRCh38 Illumina unless noted):

- OTOA copy number and recurrent-deletion detection — 100% exact on the 248-sample cohort, real
  carrier caught, 0 false deletions.
- The p.Glu787\* conversion allele (carrier / affected typing) — specificity 0/248 false
  positives; sensitivity confirmed on simulated positives.
- Small variants (SNV/indel) — full-gene F1 ~0.94 vs curated GIAB truth (HG001–HG005), with the
  uniquely-mappable 5′ region called completely.
- Cross-build (hg38 / hg19 / b37) and cross-platform (Ultima) reproducibility.

Not yet validated (deliberately out of v1 scope):

- OTOA pathogenic alleles beyond the recurrent deletion and p.Glu787\* — surfaced via the
  paralog-aware small-variant calls, but not typed into a named disease call.
- The ~15 kb PSV-desert dark region (5 coding exons of the 3′ duplication) — duplicated but
  paralog-indistinguishable; variants there are low-confidence and flagged for orthogonal
  confirmation.
- Conversion-tract detection — reported only as non-diagnostic supporting context, not
  benchmarked against conversion-tract truth; the pathogenic call rests on the p.Glu787\*
  genotype.
- Real DFNB22 patient material — conversion sensitivity is established by simulation; the
  deletion mechanism is confirmed on a real cohort carrier.

## References

- Zwaenepoel I, et al. *Otoancorin, an inner ear protein restricted to the interface between
  the apical surface of sensory epithelia and their overlying acellular gels, is defective in
  autosomal recessive deafness DFNB22.* PNAS 2002.
- Laurent S, et al. *OTOAP1 → OTOA gene conversion and the p.Glu787\* allele in DFNB22.* 2021
  (PMID 33492714).
- Kim BJ, et al. *OTOA p.Glu787\* as a short-read false positive requiring long-range
  confirmation.* 2025 (PMID 39943967).
- Sugimoto H, et al. *OTOA copy-number variation as a frequent cause of autosomal-recessive
  hearing loss* (Japanese cohort; PMID 31527525).
