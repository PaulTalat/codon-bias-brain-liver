# Codon Usage Bias: Prefrontal Cortex vs. Liver Gene Sets

## Question

Do genes highly expressed in the prefrontal cortex (e.g. NMDA receptor subunits, synaptic proteins) show a different codon usage bias profile than genes highly expressed in the liver?  [measured via Codon Adaptation Index (CAI) and Relative Synonymous Codon Usage (RSCU)]

## Method (Summary)

- **Gene sets:** 14 brain-associated genes (e.g. GRIN1, GRIN2A/B/C, GRIA1/2, SYN1, SYT1, SNAP25, CAMK2A, DLG4, SYP, NRGN, RBFOX3) and 14 liver-associated genes (e.g. ALB, APOB, CYP3A4, CYP2E1, SERPINA1, F2, HNF4A, TTR, APOA1, HAMP, CYP1A2, ASGR1, TAT, ARG1), each confirmed as tissue-enriched via GTEx/Human Protein Atlas before inclusion.
- **Sequence retrieval:** Coding sequences (CDS only, not full mRNA or genomic DNA) fetched programmatically via the NCBI Entrez API (Biopython), with an automated sanity check on every sequence (valid start codon, length divisible by 3, valid stop codon).
- **CAI reference weights:** Built from scratch using a set of 10 human ribosomal protein genes (RPL3, RPL4, RPL7, RPL8, RPL11, RPL13, RPS2, RPS3, RPS6, RPS9) as the highly-expressed reference set, following the standard Sharp & Li (1987) method — codon counts pooled across all 10 genes, then converted to relative adaptiveness (w) per codon.
- **Verification:** Before running the full pipeline, one gene (GRIN1) was independently checked by hand using a manual codon usage tool (CAIcal) to confirm the automated RSCU calculation was correct (matched exactly: ATC RSCU = 2.38 in both). The automated CAI value for GRIN1 (0.752) differs from the manual CAIcal value (0.853) because the manual tool uses a whole-genome reference set, while this pipeline deliberately uses a small, ribosomal-protein-only reference set — an expected and explainable difference, not an error.
- **Statistical testing:** Mann-Whitney U tests (non-parametric, appropriate for small sample sizes) were used to compare CAI distributions between groups, and to test individual codon RSCU divergences. A Bonferroni correction was applied when testing multiple codons simultaneously, to control for the increased chance of false positives from running several tests at once.

## Results

### 1. Overall codon optimization (CAI) does not differ significantly between groups

| Group | n | Mean CAI | Variance |
|---|---|---|---|
| Brain | 14 | 0.686 | 0.00230 |
| Liver | 14 | 0.672 | 0.00408 |

Mann-Whitney U = 105.0, **p = 0.765** — no significant difference. The spread of CAI values *within* each group (e.g. Brain ranges from 0.580 to 0.752) is substantially larger than the difference *between* group averages, indicating that individual gene identity is a stronger determinant of codon optimization than tissue category alone.

A check for a confounding relationship between gene length and CAI found only a weak correlation (r = -0.257 across all 28 genes; r = 0.082 when excluding the one notable length outlier, APOB), suggesting gene length is not meaningfully driving the CAI values observed.

### 2. Specific codon preferences do diverge between groups — but the strength of evidence varies by codon

The three codons with the largest raw average RSCU differences between groups were tested individually:

| Codon | Amino Acid | Brain avg RSCU | Liver avg RSCU | Raw p-value | Significant after Bonferroni correction (α = 0.0167)? |
|---|---|---|---|---|---|
| AGA | Arginine | 0.685 | 1.460 | 0.043 | No |
| CGG | Arginine | 1.617 | 0.914 | 0.147 | No |
| CTC | Leucine | 1.003 | 1.549 | 0.037 | No |

Notably, both AGA and CTC individually cleared the conventional p < 0.05 threshold — but neither survives correction for having tested three codons simultaneously. Applying a Bonferroni-adjusted threshold (0.05 ÷ 3 ≈ 0.0167) is the conservative, appropriate step here, since testing multiple codons increases the chance that at least one appears "significant" purely by chance.

## Honest Interpretation

- **No strong evidence** that overall codon adaptation (CAI) differs between genes highly expressed in prefrontal cortex versus liver, in this sample.
- **Suggestive, but not statistically confirmed**, evidence that brain and liver genes may prefer different synonymous codons for Arginine and Leucine specifically (notably a possible AGA ↔ CGG preference split for Arginine). These patterns are visually apparent and individually significant before correction, but should be treated as an exploratory observation warranting a larger sample, not a confirmed finding, given they do not survive correction for multiple comparisons.
- This project does not establish *why* any observed pattern exists — a plausible mechanism discussed in the literature is tissue-specific variation in tRNA pool abundance, but tRNA abundance was not measured here, and this remains a hypothesis for follow-up work, not a conclusion of this analysis.

## Limitations

- Small sample size (14 genes per group) limits statistical power, particularly for detecting subtle codon-level effects.
- Gene selection, while confirmed via GTEx/HPA tissue-expression data, was not exhaustive or randomly sampled — a different set of 14+14 genes could plausibly yield different results.
- CAI values are only meaningful relative to the specific reference set used (ribosomal proteins, in this case); they should not be directly compared to CAI values computed with a different reference set.
