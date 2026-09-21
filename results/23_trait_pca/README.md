# 23 — Multivariate phenotype PCA + biplot

**Script:** `scripts/23_trait_pca.py`

## What was done

Two phenotype PCAs: (A) full 13-trait PCA on lines with all 13 traits scored (n = 45); (B) 9-trait PCA on biochem + seed metrics only (n = 103). Both renderings include a biplot with trait-loading arrows.

## Method

Phenotype frame standardised (StandardScaler), PCA via scikit-learn 1.4. Biplot loadings overlaid as arrows from origin, scaled to fit the score scatter. Samples coloured by marker-PCA cluster from `01_qc_pca_power/`.

## Quick findings

- **9-trait PCA (n = 103)**: PC1 = 34.5 %, PC2 = 20.5 %, PC1+PC2 cumulative = 55.0 %, top-3 = 70.6 %.
  - PC1 dominated by seed-size traits: Seed_Width (+0.50), Mass_of_Seeds (+0.49), Seed_Thickness (+0.46) — a "seed size" axis.
- **13-trait PCA (n = 45)**: PC1 = 24.9 %, PC2 = 18.1 %, cumulative top-3 = 58.6 %.
- Marker-PCA clusters don't separate cleanly in trait PC space — confirms the cluster × trait result (`25_cluster_trait_tests/`): genetic structure is not strongly phenotypically differentiated.

## Caveats

- 13-trait PCA at n = 45 lines has marginal power; loadings should not be over-interpreted on individual trait scales.
- StandardScaler ignores trait-specific units (a fair choice for variance-explained PCA) but means the loading arrows aren't on comparable trait-magnitude scales.
- PC1 of the 9-trait PCA captures only the seed-size dimension; the biochem traits load mostly on PC2 and PC3.
