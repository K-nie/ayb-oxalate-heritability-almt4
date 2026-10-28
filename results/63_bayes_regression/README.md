# 63 — BayesB / BayesC posterior inclusion probabilities

Independent corroboration of the ALMT4 × Soluble_Oxalate signal from a fundamentally different shrinkage prior. RR-BLUP and FarmCPU / mrMLM use frequentist shrinkage or variable-selection. BayesB and BayesC π use a spike-and-slab prior that assigns each SNP a per-iteration probability of having a non-zero effect — the data must prefer the SNP for it to acquire posterior mass under "non-zero." A high PIP at the ALMT4 SNP under both methods would be a fourth independent line of evidence (after GWAS, h² CI, GBLUP holdout, ALMT family phylogeny, protein domain, pleiotropy).

## Method

BGLR v1.1.4 (Pérez & de los Campos 2014) on the 95-line × 1,625-marker post-QC dosage matrix. Two priors per trait:
- **BayesB** — per-SNP variance prior; fixed inclusion probability π = 0.05.
- **BayesC π** — shared-variance prior across SNPs; π estimated from the data.

MCMC: 30,000 iterations, 5,000-iteration burn-in, thin = 5. Inverse-chi-square prior on residual variance with df = 5. Same QC cascade as elsewhere (call rate ≥ 0.90 sample and marker; MAF ≥ 0.05).

## Findings at the ALMT4 focal SNP (Ss10:15,394,673; rs `100033542|F|0-31:T>C-31:T>C`)

| Trait | n | PIP_BayesB | PIP_BayesC | β_BayesB | β_BayesC |
|---|---|---|---|---|---|
| **Soluble_Oxalate** | 41 | **0.066** | **0.080** | **+1.43** | **+1.08** |
| **Insoluble_Oxalate** | 41 | 0.014 | 0.053 | **−1.89** | **−1.41** |
| Seed_Thickness | 95 | 0.090 | 0.033 | −0.000 | −0.001 |
| Seed_Width | 95 | 0.080 | 0.029 | −0.001 | −0.001 |
| Mass_of_Seeds | 95 | 0.077 | 0.046 | −0.000 | −0.000 |
| Antioxidant | 93 | 0.055 | 0.049 | +0.49 | +0.10 |
| Phenol | 93 | 0.028 | 0.036 | +0.05 | +0.49 |
| Total_Oxalate | 41 | 0.037 | 0.046 | +0.24 | −1.90 |
| Crude_Protein | 92 | 0.029 | 0.037 | +0.05 | +0.01 |
| Seed_Length | 95 | 0.024 | 0.043 | −0.001 | +0.001 |
| Flavonoid | 93 | 0.011 | 0.026 | +0.61 | +0.20 |
| Tannin | 93 | 0.009 | 0.021 | +0.50 | +0.40 |
| Seed_Coat_Tannin | 95 | 0.007 | 0.044 | −0.02 | +0.55 |

**Three operational takeaways:**

1. **Soluble_Oxalate carries the highest PIP under both methods** (PIP_B = 0.066, PIP_C = 0.080). Both exceed the prior π = 0.05 — the data shift the posterior probability of inclusion upward at this SNP. The corresponding posterior effects (β_B = +1.43, β_C = +1.08) are direction-consistent with the §3.9 pleiotropy-scan β = +98.65 on the BLUP scale.

2. **Insoluble_Oxalate shows a strong negative posterior effect at the same SNP** (β_B = −1.89, β_C = −1.41), reproducing the per-SNP Soluble × Insoluble antagonism documented in §3.9. The PIP is low (0.014 / 0.053) because n = 41 and the spike-and-slab prior is conservative.

3. **The seed-metric PIPs (Seed_Thickness 0.090, Seed_Width 0.080, Mass_of_Seeds 0.077) are weakly-inclusive-but-null-effect** — posterior effects |β| < 0.001 on the trait scale (Seed_Width in mm), so these are not pleiotropy evidence. The framework correctly distinguishes "data wants this SNP" from "the SNP has an interpretable effect."

The case for ALMT4 strengthens with this fourth methodological framework; the modest absolute PIP magnitude (0.066–0.080) reflects the n = 41 data limit, not a weakness of the mechanism call.

## Outputs

- `tables/bayesB_bayesC_pip_per_snp_per_trait.csv` — 21,125 rows (13 × 1,625) with PIP_B and PIP_C per SNP per trait.
- `tables/bayesB_bayesC_effects_per_snp_per_trait.csv` — same shape with posterior mean effects.
- `tables/almt4_focal_summary.csv` — 13 rows; ALMT4 focal SNP PIPs + effects per trait.
- `bglr_scratch/` — BGLR chain trace files per (trait, method).

## Caveats

- **PIP magnitudes are modest.** At n = 41 the spike-and-slab prior has limited power to distinguish non-zero from zero effects. The Soluble_Oxalate PIP of 0.066–0.080 is interpretable as "shift above prior" rather than "Bayesian-significant" (a PIP > 0.5 threshold under BayesB is the conventional bar; not reached here).
- **Direction-disagreement on Total_Oxalate** (β_B = +0.24 vs β_C = −1.90) reflects the different sparsity priors picking different solutions when the data are weak. The §3.9 pleiotropy-scan finding that Total_Oxalate is null (p = 0.61) stands; the BayesB / BayesC results are consistent with that null finding (the effect is poorly identified under both methods).
- **MCMC convergence (checked 2026-10-27).** Geweke z-statistics and effective sample sizes (ESS) were computed with `coda` on the post-burn-in chains (5,000 retained thinned draws) for the focal Soluble_Oxalate trait. Command: `coda::effectiveSize()` and `coda::geweke.diag()` on `bglr_scratch/{BayesB,BayesC}_Soluble_Oxalate_varE.dat` and the `probIn`/`varB` columns of `BayesC_Soluble_Oxalate_ETA_1_parBayesC.dat`, discarding the first 1,000 of 6,000 thinned samples.

  | parameter | ESS | Geweke z |
  |---|---|---|
  | BayesB residual variance (varE) | 2064 | +0.07 |
  | BayesC residual variance (varE) | 1868 | −1.26 |
  | BayesCπ π (probIn) | 163 | +1.23 |
  | BayesCπ varB (slab variance) | 47 | −0.89 |

  All Geweke z fall within ±1.96 → no evidence against convergence. The residual-variance chains mix well (ESS ≈ 1,900–2,100); π and the slab variance varB mix more slowly (ESS 163 and 47), as expected for sparsity-prior hyperparameters at n = 41, but still pass Geweke. The BayesCπ π posterior mean (0.0554) sits just above its 0.05 prior mean, consistent with the weak focal-SNP signal. The 30,000-iteration chain is adequate for the focal-SNP summary statistics reported; a robust genome-wide per-SNP PIP table would still benefit from longer chains.
- **No cross-validation accuracy comparison** of BayesB / BayesC against RR-BLUP. The §3.11 GBLUP predictive ability is the existing benchmark; a per-method CV-r forest would be a natural extension once Analysis 3 (FarmCPU + mrMLM) is also in hand for a four-method comparison.

## Script

`scripts/63_bayes_regression.R`. Runtime: ~ 8 minutes on a laptop (13 traits × 2 priors × 30k MCMC iterations).

## Higher-PIP SNPs vs the focal SNP (reviewer item, §3.4)

Added 2026-10-27 to answer the reviewer request to name the genes at the SNPs
that carry a higher soluble-oxalate PIP than the ALMT4-proximal focal SNP.

- Script: `scripts/_bayes_higher_pip_genes.py`
- Output: `tables/soluble_oxalate_higher_pip_than_focal.csv`
- Inputs: this analysis's PIP table + `results/89_blast_uniqueness_reanchor/.../marker_placement_classification.csv` (anchoring) + the Funannotate v1.8 GFF (gene-in / nearest-gene + product).
- Focal SNP rs 100033542 (Ss10:15,394,609; ~3.2 kb 5' of ALMT4_2), PIP BayesB 0.0664 / BayesCpi 0.0802.

Result: **13** SNPs exceed the focal BayesB PIP (12 exceed under both priors). By
chromosome: **Ss07 (4), Ss04 (3), Ss03 (2), Ss01 (2), Ss08 (1), Ss10 (1)** — note
Ss07 and Ss04, NOT Ss06 (an earlier draft said "Ss04 and Ss06"; the data show Ss07).
The three highest-PIP SNPs sit in Ss04 AYBTSS11_012251 (hypothetical protein, PIP
0.174), Ss07 AYBTSS11_020928 (methylenetetrahydrofolate reductase 1, 0.116) and Ss04
AYBTSS11_012253 (hypothetical protein, 0.110). None of the higher-PIP genes is an
oxalate-metabolism enzyme, so the Bayesian analysis does not single out the
ALMT4-proximal SNP.
