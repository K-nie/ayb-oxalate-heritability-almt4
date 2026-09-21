#!/usr/bin/env python3
"""
Assemble the canonical Paper 1 figure set into
`results/04_publication_plots/figures/`.

Reads from the per-analysis directories (the authoritative copies) and
writes paper-numbered duplicates so the manuscript can reference a single
`fig01..fig20` series. Source paths are explicit, so any later swap of
a figure version is a one-line edit here.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
RES  = ROOT / "results"
DST  = RES / "04_publication_plots" / "figures"
DST.mkdir(parents=True, exist_ok=True)

# (paper-figure number, short slug, source basename without extension,
#  source per-analysis dir).
# Listed in the draft-v0 figure table order.
FIGS = [
    (1,  "marker_qc",                "fig01_marker_qc",                          "04_publication_plots"),
    (2,  "sample_qc",                "fig02_sample_qc",                          "04_publication_plots"),
    (3,  "pca_scree",                "fig03_pca_scree",                          "04_publication_plots"),
    (4,  "pca_clusters",             "fig04_pca_clusters",                       "04_publication_plots"),
    (5,  "pca_phenotype",            "fig05_pca_phenotype",                      "04_publication_plots"),
    # admixture
    (6,  "admixture_cv",             "fig22_admixture_cv_error",                 "09_admixture"),
    (7,  "admixture_Q",              "fig24_admixture_Q_bars_all_K",             "09_admixture"),
    # ld decay
    (8,  "ld_decay",                 "fig27_ld_decay",                           "11_ld_decay"),
    (9,  "ld_decay_per_chrom",       "fig28_ld_decay_per_chrom",                 "11_ld_decay"),
    # phylogeny tanglegrams
    (10, "tangle_iqtree_raxml",      "fig69_tanglegram_iqtree_vs_raxml",         "19_phylogenetic_trees"),
    (11, "tangle_iqtree_upgma",      "fig70_tanglegram_iqtree_vs_upgma",         "19_phylogenetic_trees"),
    (12, "tangle_raxml_upgma",       "fig71_tanglegram_raxml_vs_upgma",          "19_phylogenetic_trees"),
    # selection scans
    (13, "pcadapt_manhattan",        "fig64_pcadapt_manhattan",                  "27_pcadapt_outliers"),
    (14, "pcadapt_qq",               "fig65_pcadapt_qq",                         "27_pcadapt_outliers"),
    (15, "dapc_density",             "fig67_dapc_density",                       "29_dapc"),
    (16, "dapc_manhattan",           "fig68_dapc_marker_manhattan",              "29_dapc"),
    # GRM
    (17, "grm_heatmap",              "fig31_grm_heatmap",                        "13_grm_crosspairs"),
    (18, "grm_offdiag",              "fig32_grm_offdiag_distribution",           "13_grm_crosspairs"),
    # phenotype eda
    (19, "pheno_distributions",      "fig34_phenotype_distributions_13",         "14_pheno_eda_13traits"),
    (20, "corr_pearson",             "fig35_correlation_matrix_pearson",         "14_pheno_eda_13traits"),
    (21, "corr_spearman",            "fig36_correlation_matrix_spearman",        "14_pheno_eda_13traits"),
    # trait PCA
    (22, "trait_pca_all13",          "fig59_trait_pca_biplot_all13",             "23_trait_pca"),
    (23, "trait_pca_n9",             "fig60_trait_pca_biplot_n9",                "23_trait_pca"),
    # cluster trait tests
    (24, "cluster_trait_effects",    "fig62_cluster_trait_effects",              "25_cluster_trait_tests"),
    # GWAS panels
    (25, "qq_biochem",               "fig12_qq_M1_K",                            "04_publication_plots"),
    (26, "manhattan_biochem",        "fig13_manhattan_M1_K",                     "04_publication_plots"),
    (27, "qq_seed",                  "fig39_qq_panels_seed",                     "15_gwas_new_traits"),
    (28, "qq_protox",                "fig40_qq_panels_protox",                   "15_gwas_new_traits"),
    (29, "manhattan_seed",           "fig41_manhattan_panels_seed",              "15_gwas_new_traits"),
    (30, "manhattan_protox",         "fig42_manhattan_panels_protox",            "15_gwas_new_traits"),
    # heritability
    (31, "h2_intervals",             "fig61_h2_profile_intervals",               "24_bootstrap_h2"),
    (32, "h2_curves",                "fig62_h2_profile_curves",                  "24_bootstrap_h2"),
    # GBLUP
    (33, "gblup_biochem",            "fig14_gblup_prediction",                   "04_publication_plots"),
    (34, "gblup_all9",               "fig43_gblup_all_9",                        "15_gwas_new_traits"),
    (35, "mt_vs_st_gblup",           "fig66_mt_vs_st_gblup",                     "28_multi_trait_gblup"),
    # AYB candidate-gene Manhattans (top 3 by trait interest)
    (36, "manhattan_ayb_Soluble_Oxalate",   "fig46_manhattan_ayb_Soluble_Oxalate",   "21_candidate_genes_ayb"),
    (37, "manhattan_ayb_Insoluble_Oxalate", "fig46_manhattan_ayb_Insoluble_Oxalate", "21_candidate_genes_ayb"),
    (38, "manhattan_ayb_Crude_Protein",     "fig46_manhattan_ayb_Crude_Protein",     "21_candidate_genes_ayb"),
    # locus zoom
    (39, "zoom_ALMT4",               "fig63_zoom_ALMT4_locus",                   "26_mate_locus_zoom"),
    (40, "zoom_MATE_DETOX",          "fig63_zoom_MATE_DETOX_locus",              "26_mate_locus_zoom"),
    (41, "ldheatmap_ALMT4",          "fig63b_ldheatmap_ALMT4_locus",             "26_mate_locus_zoom"),
    (42, "ldheatmap_MATE",           "fig63b_ldheatmap_MATE_DETOX_locus",        "26_mate_locus_zoom"),
    # ALMT family tree (Paper 1 — supports the Soluble_Oxalate / ALMT4 story)
    (43, "almt_family_tree",         "fig72_almt_family_tree",                   "33_almt_phylogeny"),
    # GO enrichment (topGO combined panel)
    (44, "go_topgo_combined",        "fig75_go_topgo_combined",                  "38_go_enrichment"),
    # F_IS + F_ROH reconciliation
    (45, "fis_distribution",         "fig56_fis_distribution",                   "22_fis_inbreeding"),
    (46, "f_per_accession",          "fig57_f_per_accession",                    "22_fis_inbreeding"),
    (47, "roh_length",               "fig73_roh_length_distribution",            "34_roh"),
    (48, "f_roh_per_sample",         "fig74_f_roh_per_sample",                   "34_roh"),
    # GEBV + merit index
    (49, "merit_index_top20",        "fig48_merit_index_top20",                  "18_gebv"),
    (50, "merit_radar_top10",        "fig49_merit_index_radar_top10",            "18_gebv"),
    # core collection
    (51, "core_stratified",          "fig20_stratified_core_pca",                "06_breeding"),
    (52, "core_comparison",          "fig21_core_comparison_pca",                "06_breeding"),
]


def copy_pair(src_base: Path, dst_base: Path):
    """Copy both .png and .pdf if either exists; skip if src == dst."""
    for ext in (".png", ".pdf"):
        src = src_base.with_suffix(ext)
        dst = dst_base.with_suffix(ext)
        if not src.exists():
            continue
        if src.resolve() == dst.resolve():
            continue
        shutil.copy2(src, dst)


print(f"[curate] target = {DST}")
manifest_rows = []
missing = []
for num, slug, src_basename, src_dir in FIGS:
    src_base = RES / src_dir / "figures" / src_basename
    src_png = src_base.with_suffix(".png")
    if not src_png.exists():
        missing.append(f"  fig{num:02d}_{slug}: {src_png}")
        continue
    dst_base = DST / f"fig{num:02d}_{slug}"
    copy_pair(src_base, dst_base)
    manifest_rows.append((num, slug, src_basename, src_dir,
                          str(dst_base.with_suffix('.png').relative_to(ROOT))))

manifest_path = DST.parent / "figure_manifest.csv"
with open(manifest_path, "w") as fh:
    fh.write("fig_num,slug,source_basename,source_dir,destination_path\n")
    for row in manifest_rows:
        fh.write(",".join(str(v) for v in row) + "\n")

print(f"[curate] copied {len(manifest_rows)} figure pairs")
print(f"[curate] manifest -> {manifest_path}")
if missing:
    print("[!] missing sources:")
    for m in missing:
        print(m)
