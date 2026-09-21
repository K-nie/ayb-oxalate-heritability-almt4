"""
A2 Figure 3C: FarmCPU + mrMLM + BayesB/BayesC concordance with single-marker
MLM at the ALMT4 focal locus. Forest panel showing per-trait MLM beta vs
FarmCPU effect vs BayesB PIP-weighted effect, with consistent colour for
the headline trait (Soluble_Oxalate, the only Holm-significant signal).

Reads:  results/61_farmcpu_mrmlm/tables/three_method_concordance.csv
Writes: results/61_farmcpu_mrmlm/figures/fig_farmcpu_mrmlm_concordance.png/.pdf
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from _plotstyle import apply, WONG
apply()

PROJ = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB  = PROJ / "results/61_farmcpu_mrmlm/tables/three_method_concordance.csv"
FIG  = PROJ / "results/61_farmcpu_mrmlm/figures"
FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(TAB)
df["trait_label"] = df["trait"].str.replace("_", " ")
# Highlight Soluble_Oxalate (the Holm-significant MLM signal)
df["highlight"] = df["trait"] == "Soluble_Oxalate"
# Order: headline first, then by |mlm_beta| descending
df = df.sort_values(
    by=["highlight", "mlm_beta"],
    ascending=[False, False],
).reset_index(drop=True)

fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.5),
                          constrained_layout=True, sharey=True)
ypos = np.arange(len(df))[::-1]
hl_color = WONG["vermillion"]; bg_color = "#777777"
cols = [hl_color if h else bg_color for h in df["highlight"]]

# Panel A: single-marker MLM beta (Wald) — the reference signal
ax = axes[0]
ax.scatter(df["mlm_beta"], ypos, color=cols, s=100, edgecolor="black",
           linewidth=0.6, zorder=3)
for y, p in zip(ypos, df["mlm_wald_p_holm"]):
    p_label = f"Holm p = {p:.3f}" if p < 0.10 else f"Holm p = {p:.2f}"
    ax.text(ax.get_xlim()[1] * 0.98, y, p_label,
            ha="right", va="center", fontsize=7.5, color="#444")
ax.axvline(0, color="black", linewidth=0.5)
ax.set_yticks(ypos)
ax.set_yticklabels(df["trait_label"], fontsize=9)
ax.set_xlabel("Single-marker MLM β at ALMT4 SNP")
ax.set_title("(A) MLM Wald (reference)", fontsize=10)
ax.grid(True, axis="x", alpha=0.3)

# Panel B: FarmCPU effect at the same SNP (multi-locus model)
ax = axes[1]
ax.scatter(df["farmcpu_effect"], ypos, color=cols, s=100,
           edgecolor="black", linewidth=0.6, zorder=3)
for y, p in zip(ypos, df["farmcpu_p"]):
    p_label = f"p = {p:.3f}" if p < 0.05 else f"p = {p:.2f}"
    ax.text(ax.get_xlim()[1] * 0.98, y, p_label,
            ha="right", va="center", fontsize=7.5, color="#444")
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel("FarmCPU effect at ALMT4 SNP")
ax.set_title("(B) FarmCPU (multi-locus)", fontsize=10)
ax.grid(True, axis="x", alpha=0.3)

# Panel C: BayesB / BayesC PIP — posterior inclusion probability
ax = axes[2]
ax.scatter(df["pip_bayesB_focal"], ypos - 0.14, color=cols, s=80,
           edgecolor="black", linewidth=0.6, marker="o", zorder=3,
           label="BayesB PIP")
ax.scatter(df["pip_bayesC_focal"], ypos + 0.14, color=cols, s=80,
           edgecolor="black", linewidth=0.6, marker="s", zorder=3,
           label="BayesCπ PIP")
ax.axvline(0.05, color="grey", linestyle="--", linewidth=0.6,
           label="PIP = 0.05 (weak)")
ax.axvline(0.10, color="grey", linestyle=":", linewidth=0.6,
           label="PIP = 0.10 (moderate)")
ax.set_xlabel("Posterior inclusion probability")
ax.set_title("(C) Bayesian shrinkage (BayesB / BayesCπ)", fontsize=10)
ax.set_xlim(0, max(0.15, df[["pip_bayesB_focal",
                              "pip_bayesC_focal"]].max().max() * 1.15))
ax.grid(True, axis="x", alpha=0.3)
ax.legend(loc="lower right", fontsize=8, framealpha=0.9)

fig.suptitle(
    "Four-method concordance at the ALMT4_2 focal SNP (Ss10:15,394,673)\n"
    "Soluble_Oxalate (vermillion) is the only Holm-significant single-"
    "marker signal; FarmCPU corroborates with p < 0.02, "
    "BayesB / BayesCπ PIPs are weak but non-zero",
    fontsize=11)
out = FIG / "fig_farmcpu_mrmlm_concordance"
fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"[fig] {out}.png / .pdf")
