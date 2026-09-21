# 14 — 13-trait phenotype EDA

**Script:** `scripts/14_pheno_eda_13traits.py`

## What was done

Exploratory phenotype analysis across all 13 traits combined from three Excel sources: 4 wet-chemistry (n ≈ 105), 5 seed metrics (n = 105), 4 protein + oxalate (n = 46). Per-trait distributions, Shapiro-Wilk normality, 13 × 13 Pearson and Spearman correlation matrices with pairwise-complete observations, and focused two-trait scatter plots.

## Method

Phenotype loading via the unified `_pheno.py` loader (handles the three Excel files + applies the two flagged decimal-point corrections: TSs282 Seed_Thickness 63.747 → 6.3747 and TSs136 Seed_Coat_Tannin 7063.72 → 706.372). Per-trait Shapiro-Wilk W and p. Correlation matrices computed pairwise-complete (different n per cell, n_obs table preserved). Manual `imshow` + `ax.text` annotation grid replaces seaborn `heatmap(annot=True, square=True)` which silently drops lower-row annotations under `constrained_layout`.

## Quick findings

- Non-normal traits (W p < 0.05): Tannin, Phenol, Antioxidant, Mass_of_Seeds, Total_Oxalate, Insoluble_Oxalate.
- Normal-ish: Flavonoid, Seed_Length / Width / Thickness, Seed_Coat_Tannin, Crude_Protein, Soluble_Oxalate.
- **Bulk Tannin × Seed_Coat_Tannin Pearson r = 0.14** — barely correlated; the two tannin assays measure different phenotypes.
- Total_Oxalate × Insoluble_Oxalate r = **0.91** (insoluble dominates total).
- Soluble × Insoluble Oxalate r = **−0.38** (antagonistic; breeder-actionable).
- Seed-size cluster: Seed_Width × Seed_Thickness 0.74; Seed_Width × Mass 0.73.
- Antioxidant × Flavonoid 0.49; Antioxidant × Seed_Coat_Tannin 0.42.

## Caveats

- n per cell of the correlation matrix differs (4-biochem-only pairs use n = 103–105; protein-oxalate-involved pairs use n ≤ 46).
- Two decimal corrections were applied at the loader; the raw source spreadsheets are untouched.
- The 13-trait set has uneven coverage — Crude_Protein/Oxalate is a 46-line subset.
