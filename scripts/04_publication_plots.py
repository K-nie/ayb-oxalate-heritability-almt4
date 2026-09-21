#!/usr/bin/env python3
"""
Re-render all figures from the AYB analysis pipeline at publication quality.

Reads the CSV outputs of scripts 01-03 plus the source HapMap + wet-chemistry
sheet, applies the shared journal style, and writes new high-resolution PNG
(and matching PDF) files into results/04_publication_plots/figures/.

New panels added beyond the draft pass
--------------------------------------
- Phenotype distribution panel (Tannin / Phenol / Flavonoid / Antioxidant)
  with Shapiro-Wilk p-values and KDE overlay.
- Trait correlation heatmap (Pearson + Spearman, paired).
- QQ panels with 95 % concentration band (Beta order statistics).
- Manhattan with cleaned chromosome labels (parens stripped), explicit
  Bonferroni and FDR thresholds, alternate-colour two-tone strips.
- GBLUP prediction figure with p-value annotations placed cleanly.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from _plotstyle import apply, WONG, TRAIT_PAL, CLUSTER_PAL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PHENO_XLSX = ROOT / "Updated_Wet_Chemistry_Data.xlsx"
R1 = ROOT / "results" / "01_qc_pca_power"
R2 = ROOT / "results" / "02_amova"
R3 = ROOT / "results" / "03_gwas_mlm"
OUT = ROOT / "results" / "04_publication_plots"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

TRAITS = ["Tannin", "Phenol", "Flavonoid", "Antioxidant"]


def save(fig, name):
    """Save figure as both PNG and PDF in the publication folder."""
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load tables and source files
# ---------------------------------------------------------------------------
marker_qc = pd.read_csv(R1 / "tables/marker_qc.csv")
sample_qc = pd.read_csv(R1 / "tables/sample_qc.csv")
pca_var = pd.read_csv(R1 / "tables/pca_variance.csv")
pca_coords = pd.read_csv(R1 / "tables/pca_coords.csv")
power_df = pd.read_csv(R1 / "tables/power_curve.csv")
mdh2 = pd.read_csv(R1 / "tables/min_detectable_h2.csv")
amova_null = pd.read_csv(R2 / "tables/amova_null_distribution.csv")
trait_sum = pd.read_csv(R3 / "tables/trait_model_summary.csv")
pred = pd.read_csv(R3 / "tables/gblup_prediction_accuracy.csv")

# AMOVA observed Phi_ST from summary file
amova_summary = (R2 / "tables/amova_summary.txt").read_text()
phi_obs = float([ln for ln in amova_summary.splitlines()
                 if ln.startswith("Phi_ST (observed):")][0].split(":")[1])
p_amova = float([ln for ln in amova_summary.splitlines()
                 if ln.startswith("Empirical p-value (Phi_ST):")][0].split(":")[1])
pct_among = float([ln for ln in amova_summary.splitlines()
                   if ln.startswith("% variance among groups:")][0].split(":")[1])
pct_within = 100 - pct_among

pheno = pd.read_excel(PHENO_XLSX, sheet_name="Means")
pheno = pheno.rename(columns={"Genotypes": "sample", "Flavinoid": "Flavonoid"})

print("[load] all tables in.")


# ---------------------------------------------------------------------------
# Fig 1. Marker QC panels
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7, 5.5), constrained_layout=True)
for ax, (col, title, cutoff, x) in zip(
    axes.ravel(),
    [("maf",      "Minor allele frequency",        0.05, "MAF"),
     ("call_rate", "Marker call rate",             0.90, "Call rate"),
     ("het",      "Observed heterozygosity",       None, "H$_{obs}$"),
     ("pic",      "Polymorphism information content", None, "PIC")],
):
    data = marker_qc[col].dropna().values
    ax.hist(data, bins=45, color=WONG["blue"], edgecolor="white", linewidth=0.4)
    if cutoff is not None:
        ax.axvline(cutoff, color=WONG["vermillion"], ls="--", lw=1.2,
                   label=f"cutoff = {cutoff:g}")
        ax.legend(loc="upper right")
    ax.set_xlabel(x)
    ax.set_ylabel("SNP count")
    ax.set_title(title)
fig.suptitle(f"DArTseq marker QC — {len(marker_qc):,} SNPs across {len(sample_qc)} AYB lines",
             y=1.03, fontsize=12)
save(fig, "fig01_marker_qc")
print("[fig] fig01_marker_qc")


# ---------------------------------------------------------------------------
# Fig 2. Sample QC
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4), constrained_layout=True)
ax = axes[0]
ax.hist(sample_qc["call_rate"], bins=22, color=WONG["vermillion"],
        edgecolor="white", linewidth=0.4)
ax.axvline(0.90, color="black", ls="--", lw=1.2, label="cutoff = 0.90")
ax.set_xlabel("Per-sample call rate")
ax.set_ylabel("Sample count")
ax.set_title("Sample call rate")
ax.legend(loc="upper left")

ax = axes[1]
ax.scatter(sample_qc["call_rate"], sample_qc["observed_het"],
           s=22, color=WONG["vermillion"], edgecolor="black", linewidth=0.4,
           alpha=0.85)
ax.set_xlabel("Per-sample call rate")
ax.set_ylabel("Observed heterozygosity")
ax.set_title("Call rate vs heterozygosity")
fig.suptitle(f"AYB sample QC — {len(sample_qc)} lines", y=1.05, fontsize=12)
save(fig, "fig02_sample_qc")
print("[fig] fig02_sample_qc")


# ---------------------------------------------------------------------------
# Fig 3. PCA scree
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.2, 3.6), constrained_layout=True)
xs = np.arange(len(pca_var))
ax.bar(xs, pca_var["variance_explained"] * 100, color=WONG["blue"],
       edgecolor="white", linewidth=0.5, width=0.7)
ax2 = ax.twinx()
ax2.plot(xs, pca_var["cumulative"] * 100, color=WONG["vermillion"],
         marker="o", lw=1.8, ms=4)
ax2.set_ylabel("Cumulative variance (%)", color=WONG["vermillion"])
ax2.tick_params(axis="y", colors=WONG["vermillion"])
ax2.spines["right"].set_visible(True)
ax2.spines["right"].set_color(WONG["vermillion"])
ax.set_xticks(xs)
ax.set_xticklabels(pca_var["PC"])
ax.set_ylabel("Variance per PC (%)")
ax.set_title(f"PCA scree — {len(marker_qc[marker_qc['maf'].ge(0.05) & marker_qc['call_rate'].ge(0.90)]):,} QC-filtered SNPs")
save(fig, "fig03_pca_scree")
print("[fig] fig03_pca_scree")


# ---------------------------------------------------------------------------
# Fig 4. PCA PC1 vs PC2 colored by k-means cluster
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.6, 5.0), constrained_layout=True)
for k in sorted(pca_coords["cluster"].unique()):
    m = pca_coords["cluster"] == k
    ax.scatter(pca_coords.loc[m, "PC1"], pca_coords.loc[m, "PC2"],
               s=44, color=CLUSTER_PAL[k], edgecolor="black",
               linewidth=0.5, alpha=0.9,
               label=f"Cluster {k + 1} (n = {m.sum()})")
v1 = pca_var.iloc[0]["variance_explained"] * 100
v2 = pca_var.iloc[1]["variance_explained"] * 100
ax.set_xlabel(f"PC1 ({v1:.1f} %)")
ax.set_ylabel(f"PC2 ({v2:.1f} %)")
ax.set_title("Genetic-PCA structure of the 105-line AYB panel")
ax.legend(loc="upper right")
ax.grid(True)
save(fig, "fig04_pca_clusters")
print("[fig] fig04_pca_clusters")


# ---------------------------------------------------------------------------
# Fig 5. PCA colored by phenotype (4 panels)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.8), constrained_layout=True)
for ax, trait in zip(axes.ravel(), TRAITS):
    sc = ax.scatter(pca_coords["PC1"], pca_coords["PC2"],
                    c=pca_coords[trait], s=40, cmap="viridis",
                    edgecolor="black", linewidth=0.4)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label(trait, rotation=270, labelpad=14)
    cb.ax.tick_params(labelsize=8)
    ax.set_xlabel(f"PC1 ({v1:.1f} %)")
    ax.set_ylabel(f"PC2 ({v2:.1f} %)")
    ax.set_title(trait)
fig.suptitle("PCA scores coloured by seed-biochemistry trait values",
             y=1.02, fontsize=12)
save(fig, "fig05_pca_phenotype")
print("[fig] fig05_pca_phenotype")


# ---------------------------------------------------------------------------
# Fig 6. Phenotype distributions (NEW)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.8), constrained_layout=True)
for ax, trait in zip(axes.ravel(), TRAITS):
    v = pheno[trait].dropna().values
    sh = stats.shapiro(v)
    ax.hist(v, bins=22, color=TRAIT_PAL[trait], edgecolor="white",
            linewidth=0.4, alpha=0.85, density=True)
    xs = np.linspace(v.min(), v.max(), 200)
    kde = stats.gaussian_kde(v)
    ax.plot(xs, kde(xs), color="black", lw=1.5)
    ax.set_xlabel(trait)
    ax.set_ylabel("Density")
    ax.set_title(f"{trait}   n = {len(v)}  Shapiro W = {sh.statistic:.3f}, p = {sh.pvalue:.2g}")
fig.suptitle("Seed-biochemistry trait distributions across 105 AYB accessions",
             y=1.02, fontsize=12)
save(fig, "fig06_phenotype_distributions")
print("[fig] fig06_phenotype_distributions")


# ---------------------------------------------------------------------------
# Fig 7. Phenotype correlation heatmap (Pearson + Spearman)
# Manual imshow + ax.text annotations — seaborn.heatmap with square=True drops
# annotations on lower rows under constrained_layout, so we render the grid
# ourselves to guarantee every cell shows its r value.
# ---------------------------------------------------------------------------
pearson = pheno[TRAITS].corr(method="pearson")
spearman = pheno[TRAITS].corr(method="spearman")
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
cmap = plt.get_cmap("RdBu_r")
norm = plt.Normalize(vmin=-1, vmax=1)
for ax, mat, lab in zip(axes, [pearson, spearman], ["Pearson", "Spearman"]):
    M_ = mat.values
    im = ax.imshow(M_, cmap=cmap, norm=norm, aspect="equal")
    n = M_.shape[0]
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(mat.columns, rotation=30, ha="right")
    ax.set_yticklabels(mat.index, rotation=0)
    # white minor-grid lines between cells
    ax.set_xticks(np.arange(n + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # annotate every cell — pick text colour by |r| so it stays legible
    for i in range(n):
        for j in range(n):
            r = M_[i, j]
            txt_color = "white" if abs(r) > 0.55 else "black"
            ax.text(j, i, f"{r:.2f}", ha="center", va="center",
                    color=txt_color, fontsize=10)
    ax.set_title(f"{lab} correlation")
    fig.colorbar(im, ax=ax, shrink=0.7, label="r")
fig.suptitle("Pairwise trait correlations across the 105-line AYB panel",
             y=1.04, fontsize=12)
save(fig, "fig07_trait_correlation")
print("[fig] fig07_trait_correlation")


# ---------------------------------------------------------------------------
# Fig 8. GWAS power heatmap
# ---------------------------------------------------------------------------
piv = power_df.pivot(index="h2_qtl", columns="maf", values="power")
fig, ax = plt.subplots(figsize=(6.4, 4.6), constrained_layout=True)
sns.heatmap(piv, annot=True, fmt=".2f", cmap="rocket_r", vmin=0, vmax=1,
            ax=ax, cbar_kws={"label": "Statistical power", "shrink": 0.8},
            linewidths=0.4, linecolor="white",
            annot_kws={"fontsize": 9})
ax.invert_yaxis()
ax.set_xlabel("Minor allele frequency (MAF)")
ax.set_ylabel("QTL variance fraction $h^2_{qtl}$")
ax.set_title(f"GWAS power — n = 105, Bonferroni $\\alpha$ = {0.05/1672:.2e}")
save(fig, "fig08_power_heatmap")
print("[fig] fig08_power_heatmap")


# ---------------------------------------------------------------------------
# Fig 9. Minimum-detectable QTL effect curve
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.4, 3.6), constrained_layout=True)
ax.plot(mdh2["maf"], mdh2["h2_for_power_0.80"], marker="o",
        color=WONG["blue"], lw=2, ms=7)
for _, r in mdh2.iterrows():
    ax.annotate(f"{r['h2_for_power_0.80']:.2f}",
                xy=(r["maf"], r["h2_for_power_0.80"]),
                xytext=(0, 8), textcoords="offset points",
                ha="center", fontsize=8)
ax.set_xlabel("Minor allele frequency (MAF)")
ax.set_ylabel("Required $h^2_{qtl}$ for 80 % power")
ax.set_title("Minimum-detectable QTL effect at n = 105")
ax.grid(True)
save(fig, "fig09_min_detectable_h2")
print("[fig] fig09_min_detectable_h2")


# Fig 10. AMOVA pie chart removed at user request 2026-05-16.
# Variance components are reported in the AMOVA table (results/02_amova/
# tables/amova_table.csv); fig11 shows the permutation null with Phi_ST_obs.


# ---------------------------------------------------------------------------
# Fig 11. AMOVA permutation null
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
ax.hist(amova_null["phi_st_null"], bins=45, color=WONG["skyblue"],
        edgecolor="white", linewidth=0.4)
ax.axvline(phi_obs, color=WONG["vermillion"], lw=2.0,
           label=f"Observed $\\Phi_{{ST}}$ = {phi_obs:.3f}")
ax.set_xlabel("$\\Phi_{ST}$ under permuted group labels")
ax.set_ylabel("Count")
ax.set_title(f"AMOVA permutation null distribution (B = {len(amova_null)})")
ax.legend(loc="upper right")
ax.grid(True)
save(fig, "fig11_amova_null")
print("[fig] fig11_amova_null")


# ---------------------------------------------------------------------------
# Fig 12. QQ panels with 95 % CI band
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.8), constrained_layout=True)
for ax, trait in zip(axes.ravel(), TRAITS):
    df = pd.read_csv(R3 / f"tables/gwas_{trait}_M1_K.csv")
    p = np.sort(df["p"].values)
    m = len(p)
    k = np.arange(1, m + 1)
    exp = -np.log10((k - 0.5) / m)
    obs = -np.log10(p)
    # 95% CI via Beta order statistics: k-th of m uniforms ~ Beta(k, m-k+1)
    lo = -np.log10(stats.beta.ppf(0.975, k, m - k + 1))
    hi = -np.log10(stats.beta.ppf(0.025, k, m - k + 1))
    ax.fill_between(exp, lo, hi, color=WONG["skyblue"], alpha=0.25,
                    label="95 % concentration band")
    lam_row = trait_sum[(trait_sum["trait"] == trait) & (trait_sum["model"] == "M1_K")].iloc[0]
    ax.scatter(exp, obs, s=10, color=WONG["blue"], edgecolor="none")
    lim = max(exp.max(), obs.max(), hi.max()) + 0.2
    ax.plot([0, lim], [0, lim], color="black", ls="--", lw=0.9)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Expected $-\\log_{10}(p)$")
    ax.set_ylabel("Observed $-\\log_{10}(p)$")
    ax.set_title(f"{trait}   $\\lambda_{{GC}}$ = {lam_row['lambda_GC']:.2f}")
fig.suptitle("QQ — MLM (K kinship) per trait", y=1.02, fontsize=12)
save(fig, "fig12_qq_M1_K")
print("[fig] fig12_qq_M1_K")


# ---------------------------------------------------------------------------
# Fig 13. Manhattan, anchored markers (clean labels)
# ---------------------------------------------------------------------------
import re
def clean_chrom(c):
    if pd.isna(c):
        return None
    s = str(c)
    s = re.sub(r"\s*\(.*\)\s*", "", s)  # strip "(old3)" parts
    return s

fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.5), constrained_layout=True)
for ax, trait in zip(axes.ravel(), TRAITS):
    df = pd.read_csv(R3 / f"tables/gwas_{trait}_M1_K.csv")
    df["chrom_c"] = df["chrom"].map(clean_chrom)
    df = df.dropna(subset=["chrom_c"]).copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    if df.empty:
        ax.set_title(f"{trait} — no anchored markers")
        continue
    df = df.sort_values(["chrom_c", "pos"])
    chroms = sorted(df["chrom_c"].unique())
    palette = [WONG["blue"], WONG["skyblue"]]
    xs = []
    cum = 0
    xticks, xlabels = [], []
    for i, ch in enumerate(chroms):
        sub = df[df["chrom_c"] == ch]
        x = sub["pos"].values - sub["pos"].values.min() + cum
        ax.scatter(x, -np.log10(sub["p"].values), s=12,
                   color=palette[i % 2], edgecolor="none")
        xticks.append(cum + (sub["pos"].max() - sub["pos"].min()) / 2)
        xlabels.append(ch)
        cum += (sub["pos"].max() - sub["pos"].min()) + 5e6
    m_eff = trait_sum.loc[(trait_sum["trait"] == trait) & (trait_sum["model"] == "M1_K"),
                          "n_hits_FDR_0.10"]
    alpha_bonf = 0.05 / 1672
    bonf = -np.log10(alpha_bonf)
    fdr10 = -np.log10(0.10 / 1672 * 100)  # rough FDR display line
    ax.axhline(bonf, color=WONG["vermillion"], ls="--", lw=1,
               label=f"Bonferroni $\\alpha$ = {alpha_bonf:.2e}")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("$-\\log_{10}(p)$")
    ax.set_title(f"{trait}   {len(df)} cowpea-anchored SNPs")
    ax.legend(loc="upper right", fontsize=8)
fig.suptitle("Manhattan — MLM (K), cowpea-anchored marker subset",
             y=1.02, fontsize=12)
save(fig, "fig13_manhattan_M1_K")
print("[fig] fig13_manhattan_M1_K")


# ---------------------------------------------------------------------------
# Fig 14. GBLUP prediction accuracy — per-rep paired-violin distribution
# (observed vs permuted null) per trait. Reads the long-form per-rep table
# saved by script 03_gwas_mlm.py.
# ---------------------------------------------------------------------------
per_rep_path = R3 / "tables/gblup_per_rep_long.csv"
if per_rep_path.exists():
    per_rep = pd.read_csv(per_rep_path)
else:
    per_rep = None
    print("[warn] gblup_per_rep_long.csv not found -- "
          "fig14 will fall back to bar+errbar; re-run script 03 to enable violins")

fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
xpos = np.arange(len(TRAITS))
if per_rep is not None:
    for i, (_, row) in enumerate(pred.iterrows()):
        trait = row["trait"]
        obs_vals  = per_rep.loc[(per_rep["trait"] == trait) &
                                 (per_rep["kind"] == "observed"),
                                 "r"].values
        perm_vals = per_rep.loc[(per_rep["trait"] == trait) &
                                 (per_rep["kind"] == "permuted_null"),
                                 "r"].values
        if len(obs_vals) > 1:
            v1 = ax.violinplot([obs_vals], positions=[i - 0.20],
                                widths=0.34, showmeans=True,
                                showmedians=False, showextrema=False)
            for body in v1["bodies"]:
                body.set_facecolor(WONG["blue"]); body.set_edgecolor("black")
                body.set_alpha(0.75); body.set_linewidth(0.6)
            v1["cmeans"].set_color("black"); v1["cmeans"].set_linewidth(1.0)
        if len(perm_vals) > 1:
            v2 = ax.violinplot([perm_vals], positions=[i + 0.20],
                                widths=0.34, showmeans=True,
                                showmedians=False, showextrema=False)
            for body in v2["bodies"]:
                body.set_facecolor("#c9c9c9"); body.set_edgecolor("black")
                body.set_alpha(0.75); body.set_linewidth(0.6)
            v2["cmeans"].set_color("black"); v2["cmeans"].set_linewidth(1.0)
else:
    w = 0.36
    ax.bar(xpos - w / 2, pred["mean_r"], yerr=pred["sd_r"], width=w,
           color=WONG["blue"], edgecolor="black", linewidth=0.6, capsize=4)
    ax.bar(xpos + w / 2, pred["perm_mean"], yerr=pred["perm_sd"], width=w,
           color="#c9c9c9", edgecolor="black", linewidth=0.6, capsize=4)
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(xpos)
ax.set_xticklabels(TRAITS)
ax.set_ylabel("Predictive ability (Pearson $r$)")
ax.set_title("GBLUP genomic prediction across four seed-biochemistry traits "
             "(per-rep distribution)")
ax.bar([np.nan], [np.nan], color=WONG["blue"], edgecolor="black",
       label="Observed (5-fold CV x 50 reps)")
ax.bar([np.nan], [np.nan], color="#c9c9c9", edgecolor="black",
       label="Permuted null")
ax.legend(loc="upper left")
y_top = (pred["mean_r"] + pred["sd_r"]).max() + 0.05
for i, r in pred.iterrows():
    sig = "*" if r["perm_p"] < 0.10 else ""
    ax.text(xpos[i], y_top, f"p = {r['perm_p']:.2f}{sig}",
            ha="center", fontsize=9,
            color=WONG["vermillion"] if r["perm_p"] < 0.10 else "black")
ax.set_ylim(top=y_top + 0.08)
ax.grid(True, axis="y")
save(fig, "fig14_gblup_prediction")
print("[fig] fig14_gblup_prediction")


print("\nAll figures written to:")
print(f"  {FIG}")
