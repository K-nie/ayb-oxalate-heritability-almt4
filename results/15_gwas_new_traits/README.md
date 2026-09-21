# 15 — MLM GWAS + GBLUP for the 9 new traits

**Script:** `scripts/15_gwas_new_traits.py`

## What was done

EMMAX-style mixed-model GWAS and 5-fold cross-validated GBLUP for the 9 new traits: Seed_Length, Seed_Width, Seed_Thickness, Mass_of_Seeds, Seed_Coat_Tannin (n = 95) and Crude_Protein, Total_Oxalate, Soluble_Oxalate, Insoluble_Oxalate (n = 41 after QC of the 46-line subset).

## Method

Same EMMAX pipeline as `03_gwas_mlm.py` (VanRaden method-1 K, REML grid, spectral decomposition for fast GLS scan). For the protein/oxalate block, the dosage matrix and K are rebuilt on the trait-phenotyped lines only (n = 41), with markers re-filtered for MAF ≥ 0.05 on the subset (1,562 markers retained vs 1,625 panel-wide). GBLUP: kernel ridge with K, 5-fold CV × 50 reps + 100-perm null. Folds reduced to min(5, n/5) when n is small.

## Quick findings

- λ_GC range 0.85–1.14 across the 9 traits — well-controlled.
- **No marker reaches Bonferroni significance** for any trait.
- REML h² point estimates: Seed_Thickness 0.21 (highest seed trait); **Soluble_Oxalate 0.71 (highest of all 9 new traits)**; Total_Oxalate 0.17; most others < 0.10.
- GBLUP r: Total_Oxalate 0.11 (perm p = 0.18), Insoluble_Oxalate 0.08 (p = 0.30), most seed traits r < 0.05.
- Soluble_Oxalate single-trait GBLUP r is small but its multi-trait counterpart in `28_multi_trait_gblup/` gains modestly.

## Caveats

- Protein/oxalate n = 41 is well below any GWAS-power threshold. Frame as candidate-gene probe.
- h² point estimates here are imprecise; profile-likelihood CIs in `24_bootstrap_h2/` span 0 for everything except Soluble_Oxalate.
- The n = 41 subset uses 1,562 markers; the rest of the analyses use the full-panel 1,625.
