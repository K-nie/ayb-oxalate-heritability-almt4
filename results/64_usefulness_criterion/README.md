# 64 — Allier 2019 usefulness criterion on the top-50 cross pairs

Ranks parent pairs by `U_c = μ_c + i · σ_c` — expected merit of the top-q% of progeny rather than mid-parent value alone. The mid-parent term (μ_c) is the parental average; the within-family SD term (σ_c) captures Mendelian sampling variance, which the breeder *wants high* because they select from the upper tail of the within-family distribution.

## Method

Per cross pair (parent A, parent B), per trait:

1. **μ_c (mid-parent value)**: read directly from the `<trait>_midparent` column of `results/13_grm_crosspairs/tables/cross_pairs_least_related.csv`. Direction-signed per the §2.13 merit-index convention.

2. **σ_c (within-family SD)**: under the additive infinitesimal-limit approximation,
   `σ²_MS_trait_cross = var_g_trait · (1 − G_ij_cross) / 2`.
   Per-trait `var_g` from `results/00_blup_pipeline/tables/variance_components.csv` for the 9 traits with replicates; for the 4 no-replicate traits, BLUP variance × 0.5 under prior h² = 0.5.

3. **U_c = direction · μ_c + i · σ_c**, with `i = 1.755` corresponding to selecting the top 10 % of progeny under a normal distribution (Falconer & Mackay 1996 Table A).

Per-cross aggregation across the 13 traits uses three weighting scenarios identical to OCS (§3.14) and Smith-Hazel (§3.15):
- **balanced** — all 10 direction-signed traits weight 1.
- **soluble_oxalate_priority** — Soluble_Oxalate w=3, Insoluble_Oxalate w=2, others w=1.
- **crude_protein_priority** — Crude_Protein w=3, others w=1.

The aggregate U_c is z-standardised *within* scenario so traits on different units contribute comparably across the 50 crosses.

## Approximation note

The Bonk 2016 / Lehermeier 2017 per-SNP-effect-based Mendelian sampling variance is the gold standard but requires phased parental haplotypes and per-SNP effects. At our marker density and given that BayesB / BayesC (Analysis 4) was not run in this revision, the script uses the trait-level variance-component approximation above. The structure of the script is identical to the per-SNP version, and a one-line swap of the `σ²_MS_trait_cross` computation would upgrade the analysis once SNP effects are in hand.

## Findings

**Robust top-ranked cross (rank #1 under all three scenarios): TSs209 × TSs282** (U_z ≈ +0.96; G_ij = −0.255).

**Five additional robust top-6 crosses** appear under every scenario:
- TSs49 × TSs358 (rank #2)
- TSs325 × TSs282 (rank #3)
- TSs96 × TSs358 (rank #4)
- TSs157A × TSs363 (rank #5)
- TSs363 × TSs51 (rank #6)

All six robust top-6 crosses have G_ij between −0.234 and −0.256.

**Scenario-dependent substitutions** appear only at ranks 7–10:
- Soluble_Oxalate priority adds TSs371 × TSs358 (rank 7) and TSs98 × TSs358 (rank 10).
- Crude_Protein priority adds TSs317 × TSs358 (rank 9) and TSs157A × TSs155 (rank 10).

**Compare against Smith-Hazel (§3.15) accession-level ranking**, where zero robust elites emerge across the same three scenarios. At the cross level, the σ_c term (which scales with √(1 − G_ij)) dominates the per-scenario differences and pushes the most-diverse cross pairs to the top regardless of weighting. At the accession level the per-trait GEBV dominates and is much more sensitive to the priority weight. The two frameworks answer different questions: Smith-Hazel ranks accessions; usefulness criterion ranks crosses.

**Crude_Protein-priority scatter has a striking low-U cluster** at G_ij ≈ −0.235 to −0.240 with U_z ≈ −1.5 to −2.1. These are crosses where one or both parents carry strongly negative Crude_Protein GEBVs — the weight-3 multiplier on the Crude_Protein axis amplifies the penalty.

## Outputs

- `tables/usefulness_per_cross_per_trait.csv` — 442 rows = 50 crosses × ≤ 13 traits with mu, sigma, U, direction.
- `tables/usefulness_aggregate_per_weight.csv` — 150 rows = 50 × 3 scenarios with U_aggregate and U_aggregate_z.
- `tables/top10_under_each_weighting.csv` — 30 rows = 10 × 3 scenarios with rank, U values.
- `figures/fig_usefulness_vs_relatedness.png` / `.pdf` — scatter of U_aggregate_z vs G_ij per scenario; top-10 circled.

## Caveats

- **Trait-level σ²_MS approximation** assumes per-SNP effects are absorbed into the trait-level variance component. Valid in the infinitesimal limit but underestimates σ_c when a few SNPs of large effect dominate the trait (Hill & Weir 2011). The ALMT4 × Soluble_Oxalate effect is a known large-effect locus; the Soluble_Oxalate σ_c here is therefore conservatively biased downward. Reviewers will see this — the script extends to a per-SNP-effect formulation once BayesB / BayesC SNP effects are estimated (Analysis 4).
- **Selection intensity i = 1.755** is the top-10 % constant under normality. If the actual progeny size per cross is small (e.g. 20–50 plants), the selection intensity needs to be scaled by the finite-sample correction (`i_n` from Falconer & Mackay Table A); deferred to a per-programme tuning step.
- **The 50-cross input set comes from the GRM-derived least-related shortlist (script 13).** A larger candidate set (e.g. all 95 × 94 / 2 = 4,465 pairs) would surface different top-10 candidates. The §3.13 framing as "rank among the 50 least-related candidate pairs" is the standard breeder pre-screen + rank workflow; it is not a panel-wide cross-pair optimum.

## Script

`scripts/64_usefulness_criterion.py`. Runtime: ~ 3 s on a laptop.
