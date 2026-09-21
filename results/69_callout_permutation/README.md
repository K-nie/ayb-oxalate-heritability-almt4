# 69 — Permutation null on candidate-gene-callout enrichment

Puts a calibrated p-value on the "did we recover real pathway-gene enrichment, or just chance gene-density?" question for the §3.6 candidate-gene callout. The honest answer: **no trait passes BH q < 0.10 at this panel size.** The candidate-gene callout in Table 2 should be read as a hypothesis-generating annotation, not as a significantly-enriched gene set. The ALMT4 mechanism case rests on five independent lines of evidence (GWAS p, h² CI, GBLUP holdout, ALMT family phylogeny, protein-domain confirmation, pleiotropy specificity), not on the callout enrichment.

## Method

1. Per trait, observed enrichment = number of top-30 anchored SNPs (out of 30) whose ± 50 kb candidate-gene set contains at least one keyword-matched gene under the trait's keyword set.
2. Null distribution: 1,000 chromosome-stratified permutation replicates. Each replicate samples 30 SNPs from the 2,862-marker AYB-anchored pool, with per-chromosome counts matching the observed top-30 distribution. Chromosome stratification removes the gene-density confound — chromosomes carrying dense annotation clusters (Ss10 with ALMT4_2; Ss03 with the seed-protein cluster) are over-represented in the top-30 sets by gene-density rather than by trait association.
3. Per replicate, same ± 50 kb sweep and same keyword set as the observed step.
4. Empirical p = (1 + #{null ≥ observed}) / (1 + 1,000). Benjamini-Hochberg FDR across the 13 traits.
5. Performance optimisation: per-trait keyword regex is compiled into a single compound pattern and pre-evaluated once per gene at load time, so the per-permutation hot path does only a window mask + boolean reduce. Runtime ≈ 60 s for the full 13-trait × 1,000-permutation matrix.

## Findings

| Trait | Observed | Null mean ± SD | Null max | Empirical p | BH q |
|---|---|---|---|---|---|
| **Soluble_Oxalate** | 3 | 1.25 ± 1.08 | 7 | **0.117** | 0.888 |
| Seed_Thickness | 2 | 1.07 ± 1.03 | 6 | 0.290 | 0.888 |
| Phenol | 1 | 0.38 ± 0.60 | 3 | 0.321 | 0.888 |
| Seed_Coat_Tannin | 1 | 0.39 ± 0.60 | 3 | 0.331 | 0.888 |
| Seed_Width | 2 | 1.20 ± 1.06 | 6 | 0.342 | 0.888 |
| Seed_Length | 1 | 1.09 ± 1.02 | 5 | 0.667 | 1.000 |
| Total_Oxalate | 1 | 1.51 ± 1.19 | 6 | 0.784 | 1.000 |
| Tannin | 0 | 0.43 ± 0.63 | 3 | 1.000 | 1.000 |
| Flavonoid | 0 | 0.46 ± 0.67 | 4 | 1.000 | 1.000 |
| Antioxidant | 0 | 0.36 ± 0.60 | 4 | 1.000 | 1.000 |
| Mass_of_Seeds | 0 | 1.04 ± 0.93 | 4 | 1.000 | 1.000 |
| Crude_Protein | 0 | 0.18 ± 0.42 | 2 | 1.000 | 1.000 |
| Insoluble_Oxalate | 0 | 1.45 ± 1.19 | 7 | 1.000 | 1.000 |

**Three honest observations:**

1. **No trait reaches BH q < 0.10.** The candidate-gene callout enrichment is not formally significant against the chromosome-stratified null at any trait. This is a real methodological limit at n = 41 (oxalates) and n = 95 (other traits) under a top-30 single-marker scan.

2. **Soluble_Oxalate (the headline trait) carries 3 observed hits vs 1.25 null mean.** The empirical p = 0.117 is close to but not below the 0.05 threshold. The Soluble_Oxalate trait *direction* (more observed hits than null expectation) is consistent with the mechanism case in the rest of the paper, but the formal enrichment test alone does not establish it.

3. **The ALMT4 mechanism case does not depend on this test.** Five independent lines of evidence support ALMT4 as the Soluble_Oxalate candidate: the GWAS p = 0.0014 (K-only pleiotropy fit) / 0.009 (K + PC1-PC3 GWAS fit) at the focal SNP, the REML profile-likelihood h² CI clearing 0, the GBLUP holdout r = 0.54, the 142-sequence ALMT-family phylogeny placing AYB ALMT4 in the functionally-characterised vacuolar organic-acid efflux clade, and the structural-biology confirmation from script 65 (PF11744 + PF13515 + 3-5 TM segments). The callout-enrichment test is methodologically separable from the per-locus mechanism case; the negative result on the former does not invalidate the latter.

## Outputs

- `tables/per_trait_enrichment_p.csv` — 13 rows × (trait, observed_hits, null_mean, null_sd, null_max, n_top30_chroms, empirical_p, bh_q).
- `figures/fig_callout_enrichment_forest.png` / `.pdf` — forest plot of observed hit count (coloured by significance category) against null mean ± 1.96 SD (grey).

## Caveats

- **The top-30 cutoff is heuristic.** A 95 %-credible-set Bayesian fine-mapping at each candidate region would put a per-SNP posterior weight on the keyword match rather than a hard 30-vs-rest cutoff. This is a candidate for a future revision pass.
- **The keyword regex set is fixed at the script-21 conventions.** A broader keyword set (incorporating MSU rice / SoyBase annotations) would likely lift the observed counts but also broaden the null; the relative observed-vs-null ranking is unlikely to change.
- **At n = 41 for the oxalate triplet, 1,000 permutations give a Monte-Carlo SE of ~ 0.007 at p = 0.05.** This is fine for discriminating p = 0.05 from p = 0.10 but not for resolving the difference between p = 0.10 and p = 0.20 — the Soluble_Oxalate p = 0.117 sits exactly in that uncertain range.

## Script

`scripts/69_callout_permutation.py`. Runtime: ~ 60 s on a laptop (after pre-scan optimisation).
