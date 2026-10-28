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

## Figure 1C — re-anchoring callout (script `04d_figure_1c_reanchor_callout.py`)

**Rebuilt 2026-09-22 for the revision (Reviewer 1 major point 3).** The submitted
Fig 1C plotted a single "89.3 % anchored (n = 2,862)" bar taken from the legacy
`-max_target_seqs 1` run (`refs/ayb_genome/ayb_marker_anchoring.csv`), which cannot
support a placement-uniqueness claim. The figure now separates two quantities on
two panels, each keyed to the correct source:

- **Top panel — placement confidence** (from script 89,
  `results/89_blast_uniqueness_reanchor/tables/marker_placement_classification.csv`,
  `-max_target_seqs 50`): stacked AYB bar of **2,603 strictly unique (81.2 %)** +
  **288 clear best locus (9.0 %)** = **2,891 confident (90.2 %)**, against the
  cowpea proxy **382 (11.9 %)**. Confident = unique locus, or a clear best locus
  with a best-vs-second bitscore gap ≥ 2 (identical rule to script 89). Lift =
  2,891 / 382 = **7.57 → 7.6×** (was 7.5× on the old 2,862 numerator).
- **Bottom panel + Table S6 — per-chromosome distribution** (from the working
  anchored set `ayb_marker_anchoring.csv`, `chr_ayb`): per-chromosome anchored
  counts **149 (Ss07) – 401 (Ss03), total 2,862**, with cowpea proxy on the same
  markers **6.0 % (Ss07) – 15.8 % (Ss09)**. This set is deliberately *not* switched
  to the script-89 confident set because every downstream position-based analysis
  (LD in `results/91`, marker density, selection scans; manuscript §3.2 and §3.3)
  used these `chr_ayb` positions. Keeping the per-chromosome view on `chr_ayb`
  keeps Fig 1C consistent with manuscript paras reporting the 149–401 range and the
  6.0–15.8 % cowpea range. The two panel totals (2,891 confident vs 2,862 anchored
  positions) therefore describe related but distinct quantities by design.

**Outputs:** `results/04_publication_plots/figures/fig01c_reanchor_callout.png/.pdf`,
`results/04_publication_plots/tables/Table_S6_reanchoring_per_marker_summary.csv`,
and the manuscript copies `manuscript/paper_A1/share/figures/main/A1_NEW_1_fig01c_reanchor_callout.*`
and `share/tables/supp/TableS6_Table_S6_reanchoring_per_marker_summary.csv`.

**Accompanying manuscript tracked-change edits (2026-09-22):** §3.2 header and §3.3
body "7.5-fold" → "7.6-fold"; Methods "bitscore gap below 1" → "below 2" (the
reported 200 ambiguous / 288 clear-best counts require the ≥ 2-bit threshold, so
"below 1" was an internal inconsistency).

**Reproduce:** `cd scripts && /Users/black_einstein/miniconda3/envs/yeast-viz/bin/python 04d_figure_1c_reanchor_callout.py`
(Python 3.9.23 `yeast-viz`; pandas 2.3.1, matplotlib 3.8.4).
