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


def biplot(scores, loadings, traits, samples, clusters, title, fname,
           scale_arrow=4.0):
    """Two-panel biplot: scores scatter + loading arrows overlay."""
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), constrained_layout=True)
    # left panel: PC1 vs PC2 scores
    ax = axes[0]
    for k in sorted(set(clusters)):
        m = np.array(clusters) == k
        ax.scatter(scores[m, 0], scores[m, 1], s=44,
                   color=CLUSTER_PAL[k], edgecolor="white", linewidth=0.5,
                   label=f"Cluster {k + 1} (n = {m.sum()})", alpha=0.85)
    # Label only the 12 most outlying accessions by |PC1| + |PC2| so the
    # interior cloud is readable.
    score_pc = np.abs(scores[:, 0]) + np.abs(scores[:, 1])
    label_top_n(ax, scores[:, 0], scores[:, 1], samples,
                n=10, score=score_pc, fontsize=7.0, color="#222222",
                force_text=(1.0, 1.4), expand=(1.30, 1.55))
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA scores -- coloured by marker-PCA cluster")
    ax.legend(loc="upper right", framealpha=0.9, edgecolor="none")
    publishable_axes(ax, grid="both", grid_alpha=0.20)

    # right panel: biplot with loadings
    ax = axes[1]
    for k in sorted(set(clusters)):
        m = np.array(clusters) == k
        ax.scatter(scores[m, 0], scores[m, 1], s=22,
                   color=CLUSTER_PAL[k], edgecolor="white", linewidth=0.4,
                   alpha=0.45)
    # loadings arrows -- collect trait labels then push them apart with
    # adjustText so the 13 arrowheads do not stack labels on top of one
    # another and do not occlude the scatter cloud at the origin.
    trait_texts = []
    for i, t in enumerate(traits):
        x = loadings[i, 0] * scale_arrow
        y = loadings[i, 1] * scale_arrow
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=WONG["vermillion"],
                                    lw=1.0, alpha=0.85))
        trait_texts.append(ax.text(x * 1.08, y * 1.08, t,
                                    fontsize=8, color=WONG["vermillion"],
                                    fontstyle="italic",
                                    ha="center", va="center", zorder=10))
    adjust_labels(trait_texts, ax=ax,
                  only_move={"text": "xy", "static": "xy"},
                  force_text=(1.2, 1.6), expand=(1.35, 1.55))
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"{title} -- biplot")
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

clu9 = [int(pca_marker.loc[s, "cluster"]) if s in pca_marker.index else 0
        for s in sub9.index]
biplot(scores9[:, :2], loadings9[:, :2], list(sub9.columns),
       list(sub9.index), clu9,
       f"Trait PCA, 9 traits (n = {len(sub9)}; PC1 + PC2 = {sum(var9[:2]):.0%} variance)",
       "fig60_trait_pca_biplot_n9", scale_arrow=3.0)

print(f"\nOutputs in: {OUT}")
