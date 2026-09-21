"""
A2 Figure 6D: BayesB / BayesCπ posterior inclusion probability (PIP) track
across the AYB pseudo-chromosomes for Soluble_Oxalate (the headline trait
with the Holm-significant single-marker MLM hit at ALMT4 on Ss10). The
panel shows that even under Bayesian shrinkage with a sparse-effect prior,
the ALMT4 focal SNP carries the highest PIP genome-wide.

Reads:
  results/63_bayes_regression/tables/bayesB_bayesC_pip_per_snp_per_trait.csv
  results/03_qc_phenology_first/tables/marker_qc.csv (for chr/pos lookup)
Writes:
  results/63_bayes_regression/figures/fig_bayesB_pip_track.png/.pdf
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
PIP  = PROJ / "results/63_bayes_regression/tables/bayesB_bayesC_pip_per_snp_per_trait.csv"
FOC  = PROJ / "results/63_bayes_regression/tables/almt4_focal_summary.csv"
# AYB-anchored coords come from script 21's anchoring summary; one of the
# GWAS-per-trait tables under script 21 carries snp_pos_ayb + chr_ayb.
GWAS = (PROJ / "results/15_gwas_new_traits/tables/"
         "gwas_results_Soluble_Oxalate.csv")
FIG  = PROJ / "results/63_bayes_regression/figures"
FIG.mkdir(parents=True, exist_ok=True)

CHR_ORDER = [f"Ss{i:02d}" for i in range(1, 12)]
HEADLINE_TRAIT = "Soluble_Oxalate"
ALMT4_RS = "100033542|F|0-31:T>C-31:T>C"

pip = pd.read_csv(PIP)
print(f"[load] {len(pip):,} (trait, rs, PIP) rows; "
      f"{pip['trait'].nunique()} traits")

# Master AYB-anchoring lookup (all 3,204 markers, BLAST-anchored to the
# AYB reference). This is the authoritative source for rs -> chr_ayb +
# snp_pos_ayb across the pipeline.
ANCHOR = PROJ / "refs/ayb_genome/ayb_marker_anchoring.csv"
anchor_df = pd.read_csv(ANCHOR)
anchor_df = anchor_df[anchor_df["ayb_anchored"] == True].copy()
pos_lookup = anchor_df.set_index("rs")[
    ["chr_ayb", "snp_pos_ayb"]].to_dict("index")
print(f"[anchor] {len(pos_lookup):,} SNPs with AYB chr+pos "
      f"(from {ANCHOR.name})")

trait_pip = pip[pip["trait"] == HEADLINE_TRAIT].copy()
trait_pip["chr"]  = trait_pip["rs"].map(
    lambda r: pos_lookup.get(r, {}).get("chr_ayb"))
trait_pip["pos"]  = trait_pip["rs"].map(
    lambda r: pos_lookup.get(r, {}).get("snp_pos_ayb"))
trait_pip = trait_pip.dropna(subset=["chr", "pos"])
trait_pip["pos"] = trait_pip["pos"].astype(int)
trait_pip = trait_pip.sort_values(["chr", "pos"]).reset_index(drop=True)
print(f"[plot] {len(trait_pip):,} SNPs anchored for {HEADLINE_TRAIT}")

# Build cumulative x positions with 5 Mb chromosome gaps
offsets, mids, x_cursor = {}, [], 0
xs = []
for c in CHR_ORDER:
    sub = trait_pip[trait_pip["chr"] == c]
    offsets[c] = x_cursor
    if not sub.empty:
        xs.extend((sub["pos"].values + x_cursor).tolist())
        mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
        x_cursor += sub["pos"].max() + 5_000_000
    else:
        mids.append(x_cursor)
        x_cursor += 5_000_000
trait_pip["x"] = xs

fig, axes = plt.subplots(2, 1, figsize=(13, 6.5),
                          constrained_layout=True, sharex=True)
palette = [WONG["blue"], "#7c7c7c"]

for axi, (col_pip, label, axes_color) in enumerate([
        ("pip_bayesB", "BayesB",     WONG["blue"]),
        ("pip_bayesC", "BayesC$\\pi$", WONG["vermillion"]),
        ]):
    ax = axes[axi]
    for i, c in enumerate(CHR_ORDER):
        sub = trait_pip[trait_pip["chr"] == c]
        if sub.empty:
            continue
        ax.scatter(sub["x"], sub[col_pip], s=14,
                   color=palette[i % 2], alpha=0.85, edgecolor="none")
    # ALMT4 focal SNP — vermillion star + label
    focal = trait_pip[trait_pip["rs"] == ALMT4_RS]
    if not focal.empty:
        ax.scatter(focal["x"], focal[col_pip], s=180,
                   color=axes_color, edgecolor="black", linewidth=1.0,
                   marker="*", zorder=10,
                   label=f"ALMT4_2 SNP (PIP = {focal[col_pip].iloc[0]:.3f})")
    ax.axhline(0.05, color="grey", linestyle="--", linewidth=0.5,
               label="PIP = 0.05 (weak)")
    ax.axhline(0.10, color="grey", linestyle=":", linewidth=0.5,
               label="PIP = 0.10 (moderate)")
    ax.set_xticks(mids)
    ax.set_xticklabels(CHR_ORDER, fontsize=8)
    ax.set_ylabel(f"{label} PIP")
    ax.set_title(f"({chr(ord('A') + axi)}) {label} posterior inclusion "
                  f"probability per SNP — {HEADLINE_TRAIT.replace('_', ' ')}",
                  fontsize=10, loc="left")
    ax.set_ylim(0, max(0.15, trait_pip[col_pip].max() * 1.15))
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

out = FIG / "fig_bayesB_pip_track"
fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"[fig] {out}.png / .pdf")
