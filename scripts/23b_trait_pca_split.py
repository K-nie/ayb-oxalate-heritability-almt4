#!/usr/bin/env python3
"""
Split the 13-trait AYB trait-PCA biplot into two independent figures:

  A. PCA scores plot -- accession scatter coloured by marker-PCA cluster,
     with top-outlier accession labels.
  B. Loadings biplot -- accession scatter (dimmer) overlaid with trait
     loading arrows and italicised trait labels.

Renders both PDF and PNG. Same data + PCA as script 23; only the plotting
layout changes so each panel is a standalone plate.

Outputs
-------
results/23_trait_pca/figures/
    fig59a_trait_pca_scores_all13.png/.pdf
    fig59b_trait_pca_loadings_all13.png/.pdf

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
from _pheno import load_phenotypes
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "23_trait_pca"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def scores_plot(scores, samples, clusters, title, fname):
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for k in sorted(set(clusters)):
        m = np.array(clusters) == k
        ax.scatter(scores[m, 0], scores[m, 1], s=48,
                   color=CLUSTER_PAL[k], edgecolor="white", linewidth=0.6,
                   label=f"Cluster {k + 1} (n = {m.sum()})", alpha=0.90)
    score_pc = np.abs(scores[:, 0]) + np.abs(scores[:, 1])
    label_top_n(ax, scores[:, 0], scores[:, 1], samples,
                n=12, score=score_pc, fontsize=7.5, color="#222222",
                force_text=(1.0, 1.4), expand=(1.30, 1.55))
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title)
    ax.legend(loc="upper right", framealpha=0.9, edgecolor="none")
    publishable_axes(ax, grid="both", grid_alpha=0.20)
    save(fig, fname)
    print(f"[fig] {fname}")


def loadings_plot(scores, loadings, traits, clusters, title, fname,
                  scale_arrow=2.0):
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    for k in sorted(set(clusters)):
        m = np.array(clusters) == k
        ax.scatter(scores[m, 0], scores[m, 1], s=22,
                   color=CLUSTER_PAL[k], edgecolor="white", linewidth=0.4,
                   alpha=0.40)
    trait_texts = []
    for i, t in enumerate(traits):
        x = loadings[i, 0] * scale_arrow
        y = loadings[i, 1] * scale_arrow
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=WONG["vermillion"],
                                    lw=1.1, alpha=0.90))
        trait_texts.append(ax.text(x * 1.08, y * 1.08, t,
                                    fontsize=8.5, color=WONG["vermillion"],
                                    fontstyle="italic",
                                    ha="center", va="center", zorder=10))
    adjust_labels(trait_texts, ax=ax,
                  only_move={"text": "xy", "static": "xy"},
                  force_text=(1.2, 1.6), expand=(1.35, 1.55))
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title)
    publishable_axes(ax, grid="both", grid_alpha=0.20)
    save(fig, fname)
    print(f"[fig] {fname}")


print("[load] phenotype + marker-PCA clusters...")
pheno = load_phenotypes()
pca_marker = pd.read_csv(PCA_CSV).set_index("sample")

sub13 = pheno.dropna(how="any")
print(f"[A] lines with all 13 traits: {len(sub13)}")
if len(sub13) < 5:
    raise SystemExit("too few complete lines")

scaler = StandardScaler()
Xs = scaler.fit_transform(sub13.values)
pca13 = PCA(n_components=min(len(sub13) - 1, 13))
scores13 = pca13.fit_transform(Xs)
loadings13 = pca13.components_.T
var13 = pca13.explained_variance_ratio_

print(f"  variance explained: PC1 = {var13[0]:.1%}, PC2 = {var13[1]:.1%}, "
      f"cumulative top-3 = {sum(var13[:3]):.1%}")

clu13 = [int(pca_marker.loc[s, "cluster"]) if s in pca_marker.index else 0
         for s in sub13.index]

n = len(sub13)
pc12_pct = f"{sum(var13[:2]):.0%}"
scores_plot(scores13[:, :2], list(sub13.index), clu13,
            f"Trait PCA scores (n = {n}, all 13 traits, PC1 + PC2 = {pc12_pct})",
            "fig59a_trait_pca_scores_all13")
loadings_plot(scores13[:, :2], loadings13[:, :2], list(sub13.columns), clu13,
              f"Trait PCA loadings biplot (n = {n}, all 13 traits)",
              "fig59b_trait_pca_loadings_all13",
              scale_arrow=2.0)
print("[done]")
