# 24 — Profile-likelihood 95 % CIs on REML h² per trait

**Script:** `scripts/24_bootstrap_h2.py`

## What was done

Replaces sample-bootstrap (which fails with duplicate-row K) with profile-likelihood CIs on REML marker-based heritability — the canonical method for h² uncertainty in mixed-model genomic-prediction work (Lynch & Walsh 1998 ch. 27).

## Method

For each trait, REML log-likelihood evaluated on a fine grid of h² ∈ {0.001 … 0.999, step 0.005 with finer steps near boundaries}. MLE = argmin of −log L. 95 % profile-likelihood CI = the set of h² values whose −2 ∆log L is within χ²₁,₀.₉₅ = 3.841 of the minimum (Wilks' theorem; LRT at α = 0.05). Per-trait CI curves saved alongside the summary.

## Quick findings

- **Soluble_Oxalate is the only trait whose 95 % CI excludes 0**: MLE = 1.00, CI [0.14, 1.00].
- Antioxidant MLE = 0.41, CI [0.001, 0.999] — wide CI.
- Seed_Length MLE = 0.26, CI [0.004, 0.66] — narrowest non-Soluble CI.
- Tannin / Phenol / Flavonoid / Mass_of_Seeds / Seed_Coat_Tannin / Crude_Protein MLE ≈ 0, all CIs effectively span [0, 0.2–0.8] depending on n.
- This is the **honest** h² estimate set; point estimates from `03_gwas_mlm` and `15_gwas_new_traits` should be interpreted with these CIs attached.

## Caveats

- Profile-likelihood at the upper boundary (h² → 1) is conservative; the Soluble_Oxalate MLE at 1.00 should be read as "consistent with strong heritability, panel doesn't have enough information to localise above 0.14".
- Underlying EMMAX REML grid uses a single intercept; including PC1..PC5 as fixed covariates (M2 of script 03) would tighten some CIs but the qualitative ordering is preserved.
- Sample bootstrap was attempted and failed: resampling with replacement creates duplicate K rows, which forces REML to the boundary every time (CI always [1.0, 1.0]).

## Boundary likelihood-ratio test (added 2026-10-27, `scripts/24b_boundary_lrt.py`)

Reviewer 1 (R1.10) flagged the boundary MLE. The calibrated test of H0: h² = 0 for a
variance component on the boundary of the parameter space uses the ½χ²₀ + ½χ²₁ mixture
(Self & Liang, 1987) as the reference distribution. The LR statistic is read directly off
the saved profile: at the grid floor (h² ≈ 0) the column `neg2_delta_logL` = 2(logL_max − logL_0),
the LR statistic. `scripts/24b_boundary_lrt.py` post-processes `tables/h2_profile_curves.csv`
(no model refit) → `tables/h2_boundary_lrt.csv`.

| trait | LR (h²=0) | calibrated mixture p | BH q (13 traits) |
|---|---|---|---|
| Soluble_Oxalate | 6.02 | **0.0071** | 0.092 |
| Seed_Length | 4.00 | 0.023 | 0.148 |
| (all others) | ≤ 2.21 | ≥ 0.069 | ≥ 0.24 |

Soluble oxalate — the pre-specified focal trait — gives a calibrated nominal rejection of
zero (p = 0.0071) and the only substantial profile lower bound (0.14; seed length's two-sided
interval descends only to 0.004). Neither survives Benjamini–Hochberg correction across the
13 traits (q = 0.09, 0.15), as expected at n = 41; the result is reported as a focused,
pre-specified test, not a genome-wide-heritability discovery. Feeds Abstract, Results §3.2,
Figure 3, and response-letter R1.10.
