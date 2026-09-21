#!/usr/bin/env python3
"""
Extended phenotype EDA across all 13 AYB traits.

Inputs come from `_pheno.load_phenotypes()` which already harmonises the three
source spreadsheets and applies the two flagged decimal-point corrections
(TSs282 Seed_Thickness and TSs136 Seed_Coat_Tannin).

Three trait groups are handled separately so figures and tables stay
interpretable:
  * BIOCHEM  — Tannin, Phenol, Flavonoid, Antioxidant         (n approx 105)
  * SEED     — Length, Width, Thickness, Mass, Seed_Coat_Tannin (n = 105)
  * PROT/OX  — Crude_Protein, Total / Soluble / Insoluble Oxalate (n = 46)

Outputs (results/14_pheno_eda_13traits/)
---------------------------------------
tables/
    trait_summary_n_mean_sd.csv     per-trait n, mean, sd, min, max + Shapiro p
    pearson_13x13.csv               Pearson r (pairwise complete obs)
    spearman_13x13.csv              Spearman rho (pairwise complete obs)
    pair_n_obs.csv                  n of pairwise complete obs per cell
figures/
    fig34_phenotype_distributions_13.png/.pdf   small-multiple histograms
    fig35_correlation_matrix_pearson.png/.pdf   13x13 Pearson heatmap
    fig36_correlation_matrix_spearman.png/.pdf  13x13 Spearman heatmap
    fig37_seed_coat_vs_bulk_tannin.png/.pdf     scatter: do they agree?
    fig38_protein_vs_oxalate.png/.pdf           Crude Protein vs each oxalate

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import shapiro, pearsonr, spearmanr

from _plotstyle import apply, WONG, TRAIT_PAL
from _pheno import (load_phenotypes, trait_n_table,
                    TRAITS_BIOCHEM, TRAITS_SEED, TRAITS_PROT_OX, TRAITS_ALL)
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUT = ROOT / "results" / "14_pheno_eda_13traits"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# extend the trait palette to cover the new traits, keeping the four original
# colours where they were
EXTRA = {
    "Seed_Length":       "#1f78b4",
    "Seed_Width":        "#33a02c",
    "Seed_Thickness":    "#ff7f00",
    "Mass_of_Seeds":     "#6a3d9a",
    "Seed_Coat_Tannin":  "#b15928",
    "Crude_Protein":     "#a6cee3",
    "Total_Oxalate":     "#b2df8a",
    "Soluble_Oxalate":   "#fb9a99",
    "Insoluble_Oxalate": "#fdbf6f",
}
PALETTE = {**TRAIT_PAL, **EXTRA}


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
df = load_phenotypes()
print(f"[load] {df.shape[0]} samples x {df.shape[1]} traits")
summary = trait_n_table(df)
# Shapiro-Wilk per trait
sw_p = []
for t in TRAITS_ALL:
    s = df[t].dropna()
    if 3 <= len(s) <= 5000:
        sw_p.append(shapiro(s).pvalue)
    else:
        sw_p.append(np.nan)
summary["shapiro_p"] = sw_p
summary["shapiro_normal_at_0.05"] = summary["shapiro_p"] > 0.05
summary.to_csv(TAB / "trait_summary_n_mean_sd.csv", index=False)
print(summary.to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 34. Phenotype distributions, three rows by trait group
# ---------------------------------------------------------------------------
print("[fig] phenotype distributions across all 13 traits...")
groups = [
    ("Biochemistry (n approx 105)", TRAITS_BIOCHEM),
    ("Seed metrics (n = 105)", TRAITS_SEED),
    ("Protein and oxalate (n = 46)", TRAITS_PROT_OX),
]
max_per_row = max(len(g[1]) for g in groups)
fig, axes = plt.subplots(len(groups), max_per_row,
                         figsize=(2.4 * max_per_row, 2.4 * len(groups)),
                         constrained_layout=True)
for r, (title, traits) in enumerate(groups):
    for c in range(max_per_row):
        ax = axes[r, c]
        if c >= len(traits):
            ax.set_visible(False)
            continue
        t = traits[c]
        s = df[t].dropna()
        ax.hist(s, bins=18, color=PALETTE.get(t, WONG["blue"]),
                edgecolor="white", linewidth=0.4)
        sw_row = summary[summary["trait"] == t].iloc[0]
        sig = "*" if sw_row["shapiro_p"] < 0.05 else ""
        ax.set_title(f"{t}\nn = {int(sw_row['n'])}, "
                     f"SW p = {sw_row['shapiro_p']:.2g}{sig}",
                     fontsize=9)
        ax.set_xlabel(t)
        ax.set_ylabel("count")
        ax.grid(True, axis="y")
    axes[r, 0].annotate(title, xy=(-0.32, 0.5), xycoords="axes fraction",
                        rotation=90, ha="center", va="center",
                        fontsize=10, fontweight="bold")
fig.suptitle("Phenotype distributions across the three data sources",
             fontsize=12, y=1.02)
save(fig, "fig34_phenotype_distributions_13")


# ---------------------------------------------------------------------------
# 13 x 13 correlation matrices (pairwise complete obs)
# ---------------------------------------------------------------------------
print("[corr] computing 13x13 correlation matrices...")
n_obs = pd.DataFrame(np.zeros((len(TRAITS_ALL), len(TRAITS_ALL)), dtype=int),
                     index=TRAITS_ALL, columns=TRAITS_ALL)
P = pd.DataFrame(np.nan, index=TRAITS_ALL, columns=TRAITS_ALL)
S = pd.DataFrame(np.nan, index=TRAITS_ALL, columns=TRAITS_ALL)
for i, a in enumerate(TRAITS_ALL):
    for j, b in enumerate(TRAITS_ALL):
        s_ab = df[[a, b]].dropna()
        n_obs.loc[a, b] = len(s_ab)
        if len(s_ab) >= 3:
            if a == b:
                P.loc[a, b] = 1.0
                S.loc[a, b] = 1.0
            else:
                P.loc[a, b] = pearsonr(s_ab[a], s_ab[b]).statistic
                S.loc[a, b] = spearmanr(s_ab[a], s_ab[b]).correlation
P.to_csv(TAB / "pearson_13x13.csv")
S.to_csv(TAB / "spearman_13x13.csv")
n_obs.to_csv(TAB / "pair_n_obs.csv")


def heatmap(mat, n_mat, title, fname):
    fig, ax = plt.subplots(figsize=(10.5, 9), constrained_layout=True)
    cmap = plt.get_cmap("RdBu_r")
    norm = plt.Normalize(vmin=-1, vmax=1)
    arr = mat.values.astype(float)
    im = ax.imshow(arr, cmap=cmap, norm=norm, aspect="equal")
    n = arr.shape[0]
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(mat.columns, rotation=40, ha="right", fontsize=9)
    ax.set_yticklabels(mat.index, fontsize=9)
    ax.set_xticks(np.arange(n + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(n):
        for j in range(n):
            r = arr[i, j]
            if np.isnan(r):
                ax.text(j, i, "NA", ha="center", va="center",
                        color="black", fontsize=7)
                continue
            txt = "white" if abs(r) > 0.55 else "black"
            ax.text(j, i, f"{r:.2f}", ha="center", va="center",
                    color=txt, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.7, label="correlation")
    ax.set_title(title, fontsize=12)
    save(fig, fname)


heatmap(P, n_obs, "Pearson r — pairwise complete obs (n shown in audit table)",
        "fig35_correlation_matrix_pearson")
print("[fig] fig35_correlation_matrix_pearson")
heatmap(S, n_obs, "Spearman rho — pairwise complete obs",
        "fig36_correlation_matrix_spearman")
print("[fig] fig36_correlation_matrix_spearman")


# ---------------------------------------------------------------------------
# Fig 37. Seed-coat vs bulk Tannin
# ---------------------------------------------------------------------------
both = df[["Tannin", "Seed_Coat_Tannin"]].dropna()
r_p = pearsonr(both["Tannin"], both["Seed_Coat_Tannin"])
r_s = spearmanr(both["Tannin"], both["Seed_Coat_Tannin"])
fig, ax = plt.subplots(figsize=(5.6, 5.0), constrained_layout=True)
ax.scatter(both["Tannin"], both["Seed_Coat_Tannin"], s=30,
           color=WONG["blue"], edgecolor="black", linewidth=0.4, alpha=0.85)
ax.set_xlabel("bulk Tannin (Wet Chemistry)")
ax.set_ylabel("Seed Coat Tannin (Seed Metrics)")
ax.set_title(f"Do the two tannin assays agree?\n"
             f"Pearson r = {r_p.statistic:.2f} (p = {r_p.pvalue:.2g}), "
             f"Spearman rho = {r_s.correlation:.2f}, "
             f"n = {len(both)}")
ax.grid(True)
save(fig, "fig37_seed_coat_vs_bulk_tannin")
print(f"[fig] fig37_seed_coat_vs_bulk_tannin "
      f"(r = {r_p.statistic:.2f}, n = {len(both)})")


# ---------------------------------------------------------------------------
# Fig 38. Crude Protein vs each oxalate fraction
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8), constrained_layout=True,
                         sharey=True)
for ax, t in zip(axes, ["Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"]):
    s = df[["Crude_Protein", t]].dropna()
    if len(s) < 3:
        ax.set_visible(False)
        continue
    rp = pearsonr(s["Crude_Protein"], s[t])
    ax.scatter(s["Crude_Protein"], s[t], s=30,
               color=PALETTE.get(t), edgecolor="black", linewidth=0.4,
               alpha=0.85)
    ax.set_xlabel("Crude Protein (%)")
    ax.set_ylabel(t)
    ax.set_title(f"{t}\nr = {rp.statistic:.2f} (p = {rp.pvalue:.2g}), n = {len(s)}",
                 fontsize=10)
    ax.grid(True)
fig.suptitle("Crude protein vs oxalate fractions (n = 46 subset)",
             y=1.04, fontsize=12)
save(fig, "fig38_protein_vs_oxalate")
print("[fig] fig38_protein_vs_oxalate")

print(f"\nOutputs in: {OUT}")
