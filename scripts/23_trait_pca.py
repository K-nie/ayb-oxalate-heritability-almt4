#!/usr/bin/env python3
"""
Multivariate PCA + biplot on the 13-trait AYB phenotype matrix.

Trait PCA is the standard multivariate breeder visualisation — it shows
which accessions are multi-trait elite, which traits drive each axis, and
how the genetic-cluster partition (PCA on markers) maps onto trait space.

Two PCAs are run:
  A. Full 13-trait PCA on lines that have ALL traits scored (typically the
     ~41 lines in the protein/oxalate subset).
  B. 9-trait PCA dropping the protein/oxalate block, so all 95 lines with
     biochem + seed metrics contribute.

For each, render a biplot with trait loadings overlaid as arrows, and
colour the accession scatter by PCA-on-markers cluster from script 01.

Outputs (results/23_trait_pca/)
-------------------------------
tables/
    trait_pca_loadings_all13.csv
    trait_pca_loadings_n9.csv
    trait_pca_scores_n9.csv
figures/
    fig59_trait_pca_biplot_all13.png/.pdf
    fig60_trait_pca_biplot_n9.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from _plotstyle import apply, WONG, CLUSTER_PAL
from _figstyle import adjust_labels, label_top_n, publishable_axes
from _pheno import (load_phenotypes, TRAITS_ALL, TRAITS_BIOCHEM,
                    TRAITS_SEED, TRAITS_PROT_OX)
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "23_trait_pca"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# Breeding-relevant accessions to label on the n9 scores panel (merit-index /
# OCS elites, the top usefulness cross partners, and the near-clone trio).
# TSs157A (top merit-index accession) is emphasised. (dx, dy in points, ha).
HIGHLIGHT_OFFSETS = {
    "TSs217": (24, 14, "left"),
    "TSs291": (-16, -20, "right"),
    "TSs157A": (-30, 30, "right"), "TSs282": (22, 20, "left"),
    "TSs38": (-34, -16, "right"), "TSs49": (4, -32, "center"),
    "TSs22A": (12, -18, "left"), "TSs151B": (12, 16, "left"),
    "TSs209": (14, 6, "left"), "TSs358": (12, 10, "left"),
    "TSs361": (14, -12, "left"),
}
# Two merit leaders are marked on the n9 scores panel:
#   TSs217  -> star    : leader of the restricted (primary) merit index
#   TSs157A -> diamond : leader of the all-trait merit index (kept marked
#                        because Reviewer 2 point 26 asked whether it appears in the PCA)
STAR_ACC, DIAMOND_ACC = "TSs217", "TSs157A"


def _conf_ellipse(ax, x, y, color):
    """95% confidence ellipse (chi-square, 2 d.f.) for a 2-D point cloud."""
    from matplotlib.patches import Ellipse
    if len(x) < 3:
        return
    xy = np.column_stack([x, y])
    cov = np.cov(xy, rowvar=False)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    w, h = 2.0 * np.sqrt(5.991 * vals)       # chi2_0.95(2 d.f.) = 5.991
    ax.add_patch(Ellipse(xy.mean(axis=0), w, h, angle=ang,
                         facecolor=color, alpha=0.10, edgecolor=color,
                         lw=1.0, zorder=0))


def biplot(scores, loadings, traits, samples, clusters, title, fname,
           scale_arrow=4.0, highlight=None, group_style=None, trait_labels=None):
    """Two-panel biplot: scores scatter + loading arrows overlay.

    group_style: optional {cluster_code: (colour, legend_prefix, draw_ellipse)}.
      Lets the n9 panel colour the ungenotyped accessions in a neutral colour as a
      separate "Not genotyped" group rather than folding them into a marker cluster.
    trait_labels: optional {raw_name: plain_label} for the loadings arrows.
    """
    def _style(k):
        if group_style and k in group_style:
            return group_style[k]
        return (CLUSTER_PAL[k], f"Cluster {k + 1}", True)

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    # left panel: PC1 vs PC2 scores
    ax = axes[0]
    for k in sorted(set(clusters)):
        m = np.array(clusters) == k
        color, prefix, draw_ell = _style(k)
        ax.scatter(scores[m, 0], scores[m, 1], s=44,
                   color=color, edgecolor="white", linewidth=0.5,
                   label=f"{prefix} (n = {m.sum()})", alpha=0.85)
        if draw_ell:
            _conf_ellipse(ax, scores[m, 0], scores[m, 1], color)
    if highlight:
        # Label the named breeding-relevant accessions (leader lines + offsets).
        samp = np.array(samples)
        for name in highlight:
            idx = np.where(samp == name)[0]
            if len(idx) == 0:
                continue
            j = idx[0]
            is_star = (name == STAR_ACC)
            is_diamond = (name == DIAMOND_ACC)
            focal = is_star or is_diamond
            if is_star:
                ax.scatter([scores[j, 0]], [scores[j, 1]], marker="*",
                           s=260, color="black", zorder=6)
            elif is_diamond:
                ax.scatter([scores[j, 0]], [scores[j, 1]], marker="D",
                           s=90, facecolor="none", edgecolor="black",
                           linewidths=1.5, zorder=6)
            dx, dy, ha = HIGHLIGHT_OFFSETS.get(name, (10, 10, "left"))
            ax.annotate(name, (scores[j, 0], scores[j, 1]),
                        xytext=(dx, dy), textcoords="offset points",
                        fontsize=9 if focal else 7.3,
                        fontweight="bold" if focal else "normal",
                        color="black", zorder=7, ha=ha, va="center",
                        arrowprops=dict(arrowstyle="-", color="0.5",
                                        lw=0.5, shrinkA=0, shrinkB=2))
    else:
        # fall back to labelling the most outlying accessions (all13 panel).
        score_pc = np.abs(scores[:, 0]) + np.abs(scores[:, 1])
        label_top_n(ax, scores[:, 0], scores[:, 1], samples,
                    n=10, score=score_pc, fontsize=7.0, color="#222222",
                    force_text=(1.0, 1.4), expand=(1.30, 1.55))
    # Key the two merit-leader symbols alongside the cluster colours.
    if highlight and (STAR_ACC in samples or DIAMOND_ACC in samples):
        ax.scatter([], [], marker="*", s=180, color="black",
                   label=f"{STAR_ACC} (restricted-index leader)")
        ax.scatter([], [], marker="D", s=55, facecolor="none",
                   edgecolor="black", linewidths=1.5,
                   label=f"{DIAMOND_ACC} (all-trait-index leader)")
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    # No baked title -- caption carries it (no-title rule). The legend keys the
    # marker-PCA cluster colours and the two merit-leader symbols.
    ax.legend(loc="upper right", framealpha=0.9, edgecolor="none", fontsize=7.5)
    publishable_axes(ax, grid="both", grid_alpha=0.20)

    # right panel: biplot with loadings
    ax = axes[1]
    for k in sorted(set(clusters)):
        m = np.array(clusters) == k
        ax.scatter(scores[m, 0], scores[m, 1], s=22,
                   color=_style(k)[0], edgecolor="white", linewidth=0.4,
                   alpha=0.45)
    # loadings arrows -- collect trait labels then push them apart with
    # adjustText so the 13 arrowheads do not stack labels on top of one
    # another and do not occlude the scatter cloud at the origin.
    trait_texts, tips = [], {}
    for i, t in enumerate(traits):
        x = loadings[i, 0] * scale_arrow
        y = loadings[i, 1] * scale_arrow
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=WONG["vermillion"],
                                    lw=1.0, alpha=0.85))
        lab = (trait_labels or {}).get(t, t)
        txt = ax.text(x * 1.12, y * 1.12, lab, fontsize=8,
                      color=WONG["vermillion"], fontstyle="italic",
                      ha="center", va="center", zorder=10)
        trait_texts.append(txt)
        tips[t] = (x, y, txt)
    adjust_labels(trait_texts, ax=ax,
                  only_move={"text": "xy", "static": "xy"},
                  force_text=(1.6, 2.0), expand=(1.5, 1.7))
    # Pin the near-horizontal Seed_Thickness label clear of the crowded arrow
    # origin (its tiny PC2 component otherwise leaves the label sitting at 0,0).
    if "Seed_Thickness" in tips:
        x, y, txt = tips["Seed_Thickness"]
        lx, ly = x * 1.12 + 0.15, min(y * 1.12, -0.55)
        txt.set_position((lx, ly)); txt.set_ha("left"); txt.set_va("top")
        # short grey leader line from the arrowhead to the offset label
        ax.plot([x, lx], [y, ly], color="0.55", lw=0.6, zorder=9,
                solid_capstyle="round")
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    # No baked title -- caption carries it (no-title rule).
    publishable_axes(ax, grid="both", grid_alpha=0.20)

    save(fig, fname)
    print(f"[fig] {fname}")


# ---------------------------------------------------------------------------
# Load phenotype + marker-PCA cluster assignment
# ---------------------------------------------------------------------------
print("[load] phenotype + marker-PCA clusters...")
pheno = load_phenotypes()
pca_marker = pd.read_csv(PCA_CSV).set_index("sample")


# ---------------------------------------------------------------------------
# (A) Full 13-trait PCA on lines with ALL traits scored
# ---------------------------------------------------------------------------
print("\n[A] full 13-trait PCA...")
sub13 = pheno.dropna(how="any")
print(f"  lines with all 13 traits: {len(sub13)}")
if len(sub13) >= 5:
    scaler = StandardScaler()
    Xs = scaler.fit_transform(sub13.values)
    pca13 = PCA(n_components=min(len(sub13) - 1, 13))
    scores13 = pca13.fit_transform(Xs)
    loadings13 = pca13.components_.T  # 13 x ncomp
    var13 = pca13.explained_variance_ratio_

    print(f"  variance explained: PC1 = {var13[0]:.1%}, "
          f"PC2 = {var13[1]:.1%}, cumulative top-3 = {sum(var13[:3]):.1%}")

    loadings_df = pd.DataFrame(loadings13[:, :5],
                               index=sub13.columns,
                               columns=[f"PC{i+1}" for i in range(5)])
    loadings_df["variance_explained"] = (
        loadings_df.index.map(lambda t: np.nan))
    loadings_df.to_csv(TAB / "trait_pca_loadings_all13.csv")

    var_row = pd.DataFrame([var13], columns=[f"PC{i+1}" for i in range(len(var13))],
                           index=["var_explained"])
    var_row.to_csv(TAB / "trait_pca_loadings_all13.csv", mode="a")

    clu13 = [int(pca_marker.loc[s, "cluster"]) if s in pca_marker.index else 0
             for s in sub13.index]
    biplot(scores13[:, :2], loadings13[:, :2], list(sub13.columns),
           list(sub13.index), clu13,
           f"Trait PCA, all 13 traits (n = {len(sub13)}; PC1 + PC2 = {sum(var13[:2]):.0%} variance)",
           "fig59_trait_pca_biplot_all13", scale_arrow=2.0)
else:
    print("  too few lines with all 13 traits, skipping panel A")


# ---------------------------------------------------------------------------
# (B) 9-trait PCA dropping protein/oxalate subset
# ---------------------------------------------------------------------------
print("\n[B] 9-trait PCA (biochem + seed metrics)...")
keep_traits = TRAITS_BIOCHEM + TRAITS_SEED  # 9 traits
sub9 = pheno[keep_traits].dropna(how="any")
print(f"  lines with all 9 biochem + seed traits: {len(sub9)}")

scaler = StandardScaler()
Xs = scaler.fit_transform(sub9.values)
pca9 = PCA(n_components=min(len(sub9) - 1, 9))
scores9 = pca9.fit_transform(Xs)
loadings9 = pca9.components_.T
var9 = pca9.explained_variance_ratio_
print(f"  variance explained: PC1 = {var9[0]:.1%}, "
      f"PC2 = {var9[1]:.1%}, cumulative top-3 = {sum(var9[:3]):.1%}")

loadings_df = pd.DataFrame(loadings9[:, :5], index=sub9.columns,
                           columns=[f"PC{i+1}" for i in range(5)])
loadings_df.to_csv(TAB / "trait_pca_loadings_n9.csv")
scores_df = pd.DataFrame(scores9[:, :5], index=sub9.index,
                         columns=[f"PC{i+1}" for i in range(5)])
scores_df["marker_cluster"] = [
    int(pca_marker.loc[s, "cluster"]) + 1 if s in pca_marker.index else 0
    for s in sub9.index]
scores_df.to_csv(TAB / "trait_pca_scores_n9.csv")

print("  top |loading| per PC:")
for i in range(min(3, loadings9.shape[1])):
    abs_load = np.abs(loadings9[:, i])
    order = np.argsort(-abs_load)
    top_str = ", ".join(f"{sub9.columns[j]} ({loadings9[j,i]:+.2f})"
                         for j in order[:3])
    print(f"    PC{i+1} ({var9[i]:.1%}): {top_str}")

# Marker-cluster code per accession: 0 = not genotyped (no marker cluster),
# 1 = minor marker cluster, 2 = main marker cluster. The ungenotyped accessions
# are kept as their own neutral-coloured group, NOT folded into a marker cluster.
clu9 = [int(pca_marker.loc[s, "cluster"]) + 1 if s in pca_marker.index else 0
        for s in sub9.index]
GROUP_STYLE_N9 = {
    0: ("#B0B0B0", "Not genotyped", False),   # the 10 ungenotyped accessions
    1: (WONG["blue"], "Cluster 1", True),       # minor marker cluster (n = 11)
    2: (WONG["vermillion"], "Cluster 2", True),  # main marker cluster (n = 84)
}
NICE_TRAITS = {
    "Seed_Length": "Seed length", "Seed_Width": "Seed width",
    "Seed_Thickness": "Seed thickness", "Mass_of_Seeds": "Seed mass",
    "Seed_Coat_Tannin": "Seed coat tannin", "Crude_Protein": "Crude protein",
    "Total_Oxalate": "Total oxalate", "Soluble_Oxalate": "Soluble oxalate",
    "Insoluble_Oxalate": "Insoluble oxalate", "Tannin": "Tannin",
    "Phenol": "Phenol", "Flavonoid": "Flavonoid", "Antioxidant": "Antioxidant",
}
# Breeding-relevant accessions to label on the n9 panel: the merit-index /
# OCS elites (TSs217 restricted leader, TSs157A all-trait leader), the top
# usefulness cross partners, the near-clone trio, and the lone PC1 outlier TSs291.
HIGHLIGHT_N9 = ["TSs217", "TSs157A", "TSs282", "TSs38", "TSs22A", "TSs49",
                "TSs209", "TSs151B", "TSs358", "TSs361", "TSs291"]
biplot(scores9[:, :2], loadings9[:, :2], list(sub9.columns),
       list(sub9.index), clu9,
       f"Trait PCA, 9 traits (n = {len(sub9)}; PC1 + PC2 = {sum(var9[:2]):.0%} variance)",
       "fig60_trait_pca_biplot_n9", scale_arrow=3.0, highlight=HIGHLIGHT_N9,
       group_style=GROUP_STYLE_N9, trait_labels=NICE_TRAITS)

print(f"\nOutputs in: {OUT}")
