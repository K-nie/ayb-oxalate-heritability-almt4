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
