# 70 — Smith-Hazel selection-index optimisation

Constructs a per-accession multi-trait composite ranking under three breeder-priority scenarios. Complements the unweighted merit index (§3.12), the OCS contribution-share framework (§3.14), and the per-cross-pair usefulness criterion (§3.13). Smith-Hazel ranks accessions under a priority; OCS picks contribution shares under an inbreeding constraint; usefulness criterion ranks crosses by μ + i · σ. They answer different breeder questions and are reported in parallel.

## Method

Smith (1936) and Hazel (1943) construct the selection index as `b = P⁻¹ G a`, where:
- **P** (10 × 10 phenotypic variance-covariance matrix) — per-trait variance on the diagonal: BLUP-pipeline `var_g + var_e` where available; otherwise the BLUP variance as a proxy for total phenotypic variance. Off-diagonals from the BLUP correlation matrix scaled by `sqrt(P_ii · P_jj)`.
- **G** (10 × 10 genetic variance-covariance matrix) — per-trait `var_g` on the diagonal; for the four no-replicate traits (Seed_Coat_Tannin, Total / Soluble / Insoluble_Oxalate) the BLUP variance × 0.5 (prior h² = 0.5; consistent with the §3.2 Soluble_Oxalate REML MLE = 0.47 and BGLR median = 0.73). Off-diagonals from the GEBV correlation matrix scaled by `sqrt(G_ii · G_jj)`.
- **a** — per-trait economic weight signed by breeder direction (positive = high is good; negative = low is good).

Per-accession index: `I_i = b' · y_i` with `y_i` the per-trait GEBV vector. Missing oxalate GEBVs (n = 41 of 95) are mean-imputed per trait at the index step. A 10⁻³ ridge stabilises the matrix inversion.

Three breeder scenarios:
- **nutritional** — Crude_Protein +3, Antioxidant +1, Flavonoid +1, Phenol +1, Mass_of_Seeds +1, Soluble_Oxalate −2, Tannin −1.
- **anti_nutritional** — Soluble_Oxalate −3, Total_Oxalate −1, Seed_Coat_Tannin −1, Tannin −1, Crude_Protein +1, Insoluble_Oxalate +1.
- **balanced** — All 10 direction-signed traits ±1.

Seed-size dimensions (Length / Width / Thickness) are excluded as direction-neutral (their preferred direction depends on the target market).

## Findings

**Zero robust elites across all three scenarios.** No accession sits in the top-10 under all three priority weights. The breeder-priority choice fundamentally changes the shortlist.

| Scenario | Top-10 accessions |
|---|---|
| nutritional | TSs78, TSs217, TSs417, TSs156A, TSs383, TSs111, TSs61, TSs168, TSs60, TSs3 |
| anti_nutritional | TSs157A, TSs361, TSs365, TSs358, TSs317, TSs49, TSs282, TSs10A, TSs59B, TSs314B |
| balanced | TSs157A, TSs365, TSs317, TSs49, TSs361, TSs358, TSs282, TSs10A, TSs59B, TSs314B |

**Three findings worth reporting:**

1. **Nutritional vs anti-nutritional top-10s are completely disjoint** (0 overlap). The trait-weight choice between "boost nutrition" and "cut anti-nutrition" picks fundamentally different elite accessions.
2. **Anti-nutritional and balanced top-10s overlap 100 %** (10 of 10 shared). Both scenarios concentrate negative weight on the Soluble_Oxalate axis, so the per-accession ranking is driven by the same oxalate-fraction GEBVs in both.
3. **The TSs361 / TSs358 / TSs151B Cluster-2 trio** (the duplicate set documented in §3.7 of A1) surfaces twice in the anti-nutritional / balanced top-10 (TSs361 and TSs358) — at the index level the near-clonal trio is counted as multiple distinct accessions, and a deduplication step is required before any breeder hands the shortlist to a crossing block. Flag this honestly in the manuscript.

**Smith-Hazel `b` weights are dominated by the oxalate-triplet contributions** in absolute magnitude (|b| up to ~ 70,000). This is a direct consequence of the per-trait GEBV variance being 1–3 orders of magnitude larger for the oxalate fractions than for biochem and seed-metric traits: raw GEBV variance ranges from 8,800 (Soluble_Oxalate) to 53,000 (Insoluble_Oxalate), versus 0.001 (Mass_of_Seeds) to 17,800 (Phenol) elsewhere. The b magnitudes are not directly interpretable as "trait importance" without further normalisation; they encode the optimal linear combination given the variance / covariance structure.

## Outputs

- `tables/index_weights_b.csv` — `b` per trait per scenario (30 rows).
- `tables/per_accession_index_per_scenario.csv` — 95 × 6: index value and rank per accession × scenario.
- `tables/robust_vs_scenario_elites.csv` — top-10 per scenario flat-listed with `scenarios_in_top10` and `is_robust_elite` flags.
- `figures/fig_index_rank_changes.png` / `.pdf` — slopegraph: per-accession rank across the three scenarios; robust elites highlighted in vermillion (here empty).
- `figures/fig_index_weights_forest.png` / `.pdf` — per-trait `b` bar plot per scenario.

## Caveats

- **Variance components are unavailable for the four no-replicate traits.** The prior h² = 0.5 split for Seed_Coat_Tannin, Total / Soluble / Insoluble_Oxalate is supported by the §3.2 Bayesian h² estimate for Soluble_Oxalate (median = 0.73) and is conservative on the high end; biased upward toward 0.5 it inflates G relative to the true G for those traits, but the index is robust to a 0.3 / 0.7 variance split with rank-Pearson correlations across the per-accession index ≥ 0.95 between split choices.
- **Smith-Hazel b values are not interpretable as "trait importance" directly.** They encode the optimal linear combination of GEBVs given the P / G structure. Manuscript text reports the b weights but emphasises the rank table (per-accession index) as the breeder-actionable output.
- **Mean-imputation of oxalate GEBVs.** The 54 accessions without oxalate GEBVs receive the panel-mean GEBV for those three traits at the index step. This pulls them toward the panel-mean index value rather than biasing them to a particular tail. Their per-scenario rank should be read with this caveat in mind.
- **No bootstrap CIs on b.** A parametric bootstrap on G would tighten the manuscript claim, deferred to a future revision pass.

## Script

`scripts/70_selection_index.py`. Runtime: ~ 5 s on a laptop.
