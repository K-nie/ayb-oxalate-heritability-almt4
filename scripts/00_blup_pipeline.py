#!/usr/bin/env python3
"""
BLUP pipeline driver for the AYB DArTseq panel.

Fits a per-trait REML linear mixed model `value ~ 1 | sample` on the
replicate-level data from `Updated_Wet_Chemistry_Data.xlsx`, with
moisture-content (MC) as a fixed covariate when modelling
Crude_Protein. Saves:

    results/00_blup_pipeline/tables/
        phenotype_blups.csv          one row per sample, 13 trait columns
        variance_components.csv      var_g, var_e, n_rep_median, H2_Holland
        replicate_data.csv           the long-format replicate table used
        typo_corrections_log.csv     audit trail for the two replicate-level
                                     typos (TSs49 / TSs282 Seed_Thickness)
    figures/
        fig00_h2_holland_per_trait.png/.pdf

This is the new canonical phenotype source for the project. All
downstream scripts that call `from _pheno import load_phenotypes`
pick up these BLUPs automatically because the default behaviour of
`load_phenotypes()` is now `use_blups=True`.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
from _pheno import (load_blups, load_replicates,
                     TYPO_CORRECTIONS, REPLICATE_TYPO_CORRECTIONS,
                     TRAITS_ALL, TRAITS_REPLICATED, TRAITS_MEANS_ONLY)
apply()
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUT  = ROOT / "results" / "00_blup_pipeline"
TAB  = OUT / "tables"
FIG  = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fit BLUPs + variance components (with Rep-effect LRT sensitivity +
# Cullis H2 + per-sample reliability for Garrick de-regression)
# ---------------------------------------------------------------------------
print("[blup] fitting REML LMM per trait (with Rep-effect LRT sensitivity)...")
blups, vc, reliab = load_blups(use_mc_covariate_for_protein=True,
                                 apply_corrections=True,
                                 return_variance_components=True,
                                 return_reliability=True,
                                 run_rep_lrt=True)
blups.to_csv(TAB / "phenotype_blups.csv")
vc.to_csv(TAB / "variance_components.csv", index=False)
reliab.to_csv(TAB / "reliability_per_sample.csv", index=False)
print(f"[blup] wrote {blups.shape[0]} samples x {blups.shape[1]} traits "
       f"to phenotype_blups.csv")
print(f"[blup] wrote {reliab.shape[0]} (sample, trait) reliability rows "
       f"for Garrick de-regression")
print("\n[blup] variance-component table:")
print(vc.to_string(index=False))

# Long-format replicate table (with typo corrections applied)
reps = load_replicates(apply_corrections=True)
reps.to_csv(TAB / "replicate_data.csv", index=False)
print(f"[blup] wrote {reps.shape[0]} replicate-level rows to replicate_data.csv")

# Typo corrections audit trail
typo_rows = []
for c in TYPO_CORRECTIONS:
    typo_rows.append({"level": "means", **c})
for c in REPLICATE_TYPO_CORRECTIONS:
    typo_rows.append({"level": "replicate", **c})
pd.DataFrame(typo_rows).to_csv(TAB / "typo_corrections_log.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 00. Holland vs Cullis H2 per trait (side-by-side bars per trait)
# ---------------------------------------------------------------------------
print("[fig] H2 bar plot (Holland vs Cullis)")
vc_plot = vc.copy().sort_values("H2_cullis", ascending=False,
                                  na_position="last")

fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
xpos = np.arange(len(vc_plot))
w = 0.4
ax.bar(xpos - w/2, vc_plot["H2_holland"].fillna(0), width=w,
        color=WONG["blue"], edgecolor="black", linewidth=0.4,
        label="Holland (n_rep_median)")
ax.bar(xpos + w/2, vc_plot["H2_cullis"].fillna(0), width=w,
        color=WONG["vermillion"], edgecolor="black", linewidth=0.4,
        label="Cullis (avg PEV)")
for i, t in enumerate(vc_plot["trait"]):
    h_holl = vc_plot["H2_holland"].iloc[i]
    h_cull = vc_plot["H2_cullis"].iloc[i]
    if np.isnan(h_holl):
        ax.text(i, 0.02, "no reps", ha="center", va="bottom",
                 rotation=90, fontsize=7, color="#4d4d4d")
    else:
        ax.text(i - w/2, h_holl + 0.01, f"{h_holl:.2f}",
                 ha="center", va="bottom", fontsize=6.5)
        ax.text(i + w/2, h_cull + 0.01, f"{h_cull:.2f}",
                 ha="center", va="bottom", fontsize=6.5)
ax.set_xticks(xpos)
ax.set_xticklabels(vc_plot["trait"], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("broad-sense $H^2$")
ax.set_ylim(0, 1.05)
ax.set_title("Per-trait broad-sense heritability: Holland (median n_rep) "
              "vs Cullis (avg PEV) line-mean $H^2$")
ax.legend(loc="upper right", frameon=True, framealpha=0.9)
ax.grid(True, axis="y", alpha=0.4)
save(fig, "fig00_h2_holland_per_trait")
print(f"[done] outputs in {OUT}")
