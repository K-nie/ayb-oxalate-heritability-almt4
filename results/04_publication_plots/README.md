# 04 — Publication-quality figure regeneration

**Script:** `scripts/04_publication_plots.py`

## What was done

Re-rendering of 14 panel figures at journal-standard 300 dpi using a unified style helper (`_plotstyle.py`): Wong/Okabe-Ito colour-blind-safe palette, no top/right spines, PDF + PNG outputs at every save.

## Method

All upstream tables (QC, PCA, power, AMOVA, GWAS) are loaded from `results/01_qc_pca_power/`, `results/02_amova/`, and `results/03_gwas_mlm/`. The script reproduces the headline figures of the paper with consistent fonts (Helvetica/Arial sans-serif), 300 dpi PNG plus PDF type-42 fonts (editable text in vector form).

## Quick findings

- 14 figures rendered: marker QC, sample QC, PCA scree, PCA scatter (3 colourings), phenotype distributions, trait correlation matrix, power heatmap, minimum-detectable-h², AMOVA null, QQ panel, Manhattan panel, GBLUP bars.
- **fig07 trait correlation heatmap** annotation rendering fixed: manual `ax.text` placement replaced seaborn `heatmap(annot=True, square=True, constrained_layout=True)` which was dropping annotations in lower rows.
- **Bonferroni labels** rewritten to show actual α values (e.g. α = 3.08 × 10⁻⁵) instead of "0.05 / 1672" formula notation.

## Caveats

- The AMOVA pie chart (`fig10_amova_pie`) was removed at user request — variance components are in the AMOVA table.
- 14 figures correspond to the n = 105 / 4-trait paper draft. The 9-new-trait expansion uses figures from `14_pheno_eda_13traits/`, `15_gwas_new_traits/`, and `21_candidate_genes_ayb/`.
