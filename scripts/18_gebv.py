#!/usr/bin/env python3
"""
Per-accession genomic-estimated breeding values (GEBVs / EBVs) across all 13
AYB traits, with a multi-trait merit index for breeder ranking.

For each trait:
  1. Restrict to phenotyped genotyped lines (n approx 95 for seed + biochem,
     41 for protein + oxalate).
  2. Fit GBLUP at the full panel (no CV split), using kernel-ridge regression
     with K (VanRaden method 1) as the precomputed kernel; alpha tuned by
     leave-one-out within-panel Pearson r over {0.5, 1, 2, 5, 10, 25}.
  3. Predict each line's BV from its own row of K (the GBLUP-implied
     additive breeding value).
  4. Standardise within trait (z-score) and convert to a favourable-direction
     score (positive = breeder-desirable) using the trait-direction table
     from script 17.
  5. Aggregate the per-trait standardised scores into one multi-trait merit
     index (equal weights across the favourable-direction signed traits),
     and a separate "neutral" composite excluding the size-direction traits.

Outputs (results/18_gebv/)
--------------------------
tables/
    gebv_per_accession.csv             one row per accession, one column per trait + index
    merit_index_top20.csv              top-20 ranked elites by merit
    merit_index_top20_radar.csv        normalised values used for radar
figures/
    fig47_gebv_distribution_per_trait.png/.pdf
    fig48_merit_index_top20.png/.pdf   horizontal bar chart of top-20 merit
    fig49_merit_index_radar_top10.png/.pdf  radar plot of top-10 multi-trait

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.kernel_ridge import KernelRidge

from _plotstyle import apply, WONG
from _pheno import (load_phenotypes, TRAITS_BIOCHEM, TRAITS_SEED,
                    TRAITS_PROT_OX, TRAITS_ALL)
apply()
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
DIR_DIR_CSV = ROOT / "results" / "17_trait_extremes_13" / "tables" / "trait_directions.csv"
OUT = ROOT / "results" / "18_gebv"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
ALPHAS = [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Build filtered dosage + K (VanRaden method 1)
# ---------------------------------------------------------------------------
print("[load] HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)
ref_alt = hm["alleles"].str.split("/", expand=True)
ref_alt.columns = ["ref", "alt"]
dosage = np.full((len(hm), len(sample_cols)), np.nan)
for i in range(len(hm)):
    ref, alt = ref_alt.iloc[i]
    arr = calls.iloc[i].values
    dosage[i, arr == ref + ref] = 0.0
    dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
    dosage[i, arr == alt + alt] = 2.0
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)

call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values
M_all = dosage.loc[keep_m, keep_s].T
M_all = M_all.fillna(M_all.mean(axis=0)).values
samples_all = dosage.columns[keep_s].tolist()
print(f"[load] post-QC: {M_all.shape[0]} samples x {M_all.shape[1]} markers")

p_all = M_all.mean(axis=0) / 2.0
Z_all = M_all - 2.0 * p_all[None, :]
denom_all = 2.0 * (p_all * (1.0 - p_all)).sum()
K_all = (Z_all @ Z_all.T) / denom_all


# ---------------------------------------------------------------------------
# Load phenotypes + trait directions
# ---------------------------------------------------------------------------
pheno = load_phenotypes()
print(f"[pheno] {pheno.shape[0]} samples loaded for 13 traits")

if DIR_DIR_CSV.exists():
    direction_df = pd.read_csv(DIR_DIR_CSV)
    DIRECTION = dict(zip(direction_df["trait"], direction_df["favourable_extreme"]))
else:
    print("[dir] trait_directions.csv missing, falling back to script-17 defaults")
    DIRECTION = {
        "Tannin": "low", "Phenol": "high", "Flavonoid": "high",
        "Antioxidant": "high", "Seed_Length": "size", "Seed_Width": "size",
        "Seed_Thickness": "size", "Mass_of_Seeds": "high",
        "Seed_Coat_Tannin": "low", "Crude_Protein": "high",
        "Total_Oxalate": "low", "Soluble_Oxalate": "low",
        "Insoluble_Oxalate": "high",
    }


def make_K(M):
    p = M.mean(axis=0) / 2.0
    p = np.clip(p, 1e-6, 1 - 1e-6)
    Z = M - 2.0 * p[None, :]
    return (Z @ Z.T) / (2.0 * (p * (1.0 - p)).sum())


def fisher_z_ci(r, n, alpha=0.05):
    """Two-sided Fisher-z 95% CI for a Pearson correlation."""
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return float("nan"), float("nan")
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - z_crit * se)), float(np.tanh(z + z_crit * se))


def fit_gblup_full(y, K):
    """LOO-tuned alpha by minimum RMSE, then refit on full data and
    return per-line BV.

    BV_i = sum_j K_ij * c_j where c = (K + alpha I)^{-1} y, the closed-form
    kernel-ridge prediction at training points. Alpha is selected by
    minimum leave-one-out RMSE — the proper objective for kernel ridge,
    not Pearson r (which is scale-invariant and can pick an alpha whose
    prediction has correct shape but biased magnitude).
    """
    n = len(y)
    best_alpha, best_rmse = None, np.inf
    best_score = -np.inf
    for alpha in ALPHAS:
        preds = np.zeros(n)
        for i in range(n):
            tr = np.arange(n) != i
            kr = KernelRidge(kernel="precomputed", alpha=alpha)
            kr.fit(K[np.ix_(tr, tr)], y[tr])
            preds[i] = kr.predict(K[i:i+1, tr])[0]
        rmse = float(np.sqrt(np.mean((preds - y) ** 2)))
        if rmse < best_rmse:
            best_rmse = rmse; best_alpha = alpha
            if np.std(preds) > 0:
                best_score = float(stats.pearsonr(preds, y)[0])
    kr = KernelRidge(kernel="precomputed", alpha=best_alpha or 5.0)
    kr.fit(K, y)
    bv = kr.predict(K)
    r_lo, r_hi = fisher_z_ci(best_score, n)
    return bv, float(best_alpha or 5.0), float(best_score), r_lo, r_hi


# ---------------------------------------------------------------------------
# Per-trait GEBV
# ---------------------------------------------------------------------------
results = {}
print("\n[gebv] fitting full-panel GBLUP per trait...")
for trait in TRAITS_ALL:
    y_full = pheno[trait].reindex(samples_all).values.astype(float)
    keep_y = ~np.isnan(y_full)
    samples_used = [samples_all[i] for i, k in enumerate(keep_y) if k]
    if len(samples_used) < 20:
        print(f"  {trait:>18}: only {len(samples_used)} samples, skipping")
        continue
    # subset dosage + K to the used samples; recompute K so kinship reflects
    # actual covariance among the trait-phenotyped lines
    dose_sub = dosage.loc[keep_m, samples_used].T.values.astype(float)
    col_mean = np.nanmean(dose_sub, axis=0)
    inds = np.where(np.isnan(dose_sub))
    dose_sub[inds] = np.take(col_mean, inds[1])
    K_sub = make_K(dose_sub)
    y = y_full[keep_y]
    bv, alpha, loo_r, loo_lo, loo_hi = fit_gblup_full(y, K_sub)
    bv_series = pd.Series(bv, index=samples_used, name=trait)
    results[trait] = {
        "bv": bv_series, "alpha": alpha, "loo_r": loo_r,
        "loo_r_ci_low": loo_lo, "loo_r_ci_high": loo_hi,
        "n": len(samples_used),
    }
    print(f"  {trait:>18}: n = {len(samples_used)}, alpha = {alpha:.2f}, "
          f"LOO r = {loo_r:+.3f}  Fisher-z 95% CI [{loo_lo:+.3f}, {loo_hi:+.3f}]")


# ---------------------------------------------------------------------------
# Assemble GEBV table + favourable-direction signed score
# ---------------------------------------------------------------------------
loo_table = pd.DataFrame([{
    "trait": t,
    "n": r["n"],
    "alpha": r["alpha"],
    "loo_r": r["loo_r"],
    "loo_r_ci95_low": r["loo_r_ci_low"],
    "loo_r_ci95_high": r["loo_r_ci_high"],
    "lower_excludes_zero": (r["loo_r_ci_low"] > 0
                            if not np.isnan(r["loo_r_ci_low"]) else False),
} for t, r in results.items()])
loo_table.to_csv(TAB / "gebv_loo_accuracy.csv", index=False)
print("\n[loo] per-trait LOO Fisher-z CIs:")
print(loo_table.to_string(index=False))

gebv = pd.DataFrame(index=samples_all)
gebv.index.name = "sample"
for t, r in results.items():
    gebv[t] = r["bv"]
gebv = gebv.round(4)

# Standardise + sign: positive score means breeder-desirable
favourable_traits = [t for t in TRAITS_ALL if DIRECTION.get(t) in ("high", "low")]
z = gebv[favourable_traits].copy()
for t in favourable_traits:
    s = z[t]
    if s.std() == 0 or s.notna().sum() < 5: continue
    zt = (s - s.mean()) / s.std()
    if DIRECTION[t] == "low":
        zt = -zt
    z[t] = zt
gebv["merit_index"] = z.mean(axis=1, skipna=True)
gebv["merit_index_n_traits"] = z.notna().sum(axis=1)
gebv = gebv.sort_values("merit_index", ascending=False)
gebv.to_csv(TAB / "gebv_per_accession.csv")

top20 = gebv.head(20).copy()
top20.to_csv(TAB / "merit_index_top20.csv")
print("\n[merit] top-20 elite accessions:")
print(top20[["merit_index", "merit_index_n_traits"] +
            favourable_traits].round(3).to_string())


# ---------------------------------------------------------------------------
# Fig 47 GEBV distributions per trait
# ---------------------------------------------------------------------------
n_cols = 4
n_rows = int(np.ceil(len(TRAITS_ALL) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.2 * n_cols, 2.4 * n_rows),
                         constrained_layout=True)
PALETTE = {
    "Tannin": "#D55E00", "Phenol": "#0072B2", "Flavonoid": "#009E73",
    "Antioxidant": "#CC79A7", "Seed_Length": "#1f78b4",
    "Seed_Width": "#33a02c", "Seed_Thickness": "#ff7f00",
    "Mass_of_Seeds": "#6a3d9a", "Seed_Coat_Tannin": "#b15928",
    "Crude_Protein": "#7570b3", "Total_Oxalate": "#1b9e77",
    "Soluble_Oxalate": "#d95f02", "Insoluble_Oxalate": "#e7298a",
}
for ax, t in zip(axes.ravel(), TRAITS_ALL):
    if t not in results:
        ax.set_visible(False); continue
    bv = results[t]["bv"]
    ax.hist(bv, bins=15, color=PALETTE.get(t, WONG["blue"]),
            edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set_xlabel(f"{t} GEBV")
    ax.set_title(f"{t}\nn = {results[t]['n']}, LOO r = {results[t]['loo_r']:.2f}",
                 fontsize=9)
    ax.set_ylabel("count")
    ax.grid(True, axis="y")
for ax in axes.ravel()[len(TRAITS_ALL):]:
    ax.set_visible(False)
fig.suptitle("Per-accession GBLUP genomic-estimated breeding values (GEBVs)",
             y=1.01, fontsize=12)
save(fig, "fig47_gebv_distribution_per_trait")
print("[fig] fig47_gebv_distribution_per_trait")


# ---------------------------------------------------------------------------
# Fig 48 top-20 merit-index bars
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 8), constrained_layout=True)
top_show = top20.iloc[::-1]
y = np.arange(len(top_show))
colors = [WONG["blue"] if v > 0 else WONG["vermillion"]
          for v in top_show["merit_index"]]
ax.barh(y, top_show["merit_index"], color=colors, edgecolor="black",
        linewidth=0.4)
ax.set_yticks(y)
ax.set_yticklabels(top_show.index, fontsize=8)
ax.set_xlabel("Multi-trait merit index (avg favourable-direction z-score GEBV)")
ax.axvline(0, color="black", linewidth=0.6)
ax.set_title("Top-20 elite accessions by multi-trait merit index")
ax.grid(True, axis="x")
save(fig, "fig48_merit_index_top20")
print("[fig] fig48_merit_index_top20")


# ---------------------------------------------------------------------------
# Fig 49 radar plot of top-10 across favourable-direction traits
# ---------------------------------------------------------------------------
top10 = gebv.head(10)
# z-score-and-sign matrix, NaN-filled with 0 for radar
zmat = z.loc[top10.index, favourable_traits].copy().fillna(0).values
# normalise across all panel samples (0-1 range per trait) for radar polygons
zmin = z[favourable_traits].min().values
zmax = z[favourable_traits].max().values
rng = np.where(zmax > zmin, zmax - zmin, 1.0)
zmat_norm = (zmat - zmin) / rng
zmat_norm = np.clip(zmat_norm, 0.0, 1.0)

theta = np.linspace(0, 2*np.pi, len(favourable_traits), endpoint=False)
theta_closed = np.concatenate([theta, theta[:1]])

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection="polar"),
                       constrained_layout=True)
COLORS_RADAR = plt.cm.tab10.colors
for i in range(len(top10)):
    vals = np.concatenate([zmat_norm[i], zmat_norm[i:i+1, 0].flatten()])
    ax.plot(theta_closed, vals, linewidth=1.5, color=COLORS_RADAR[i % 10],
            label=top10.index[i])
    ax.fill(theta_closed, vals, alpha=0.08, color=COLORS_RADAR[i % 10])
ax.set_xticks(theta)
ax.set_xticklabels(favourable_traits, fontsize=8)
ax.set_yticks([0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels(["", "", "", ""])
ax.set_ylim(0, 1.05)
ax.set_title("Top-10 elites — favourable-direction GEBV profile (normalised)",
             pad=18, fontsize=11)
ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.05), fontsize=8,
          frameon=False)
save(fig, "fig49_merit_index_radar_top10")
print("[fig] fig49_merit_index_radar_top10")

# save the radar-input matrix for caption / supplement
radar_df = pd.DataFrame(zmat_norm, index=top10.index, columns=favourable_traits)
radar_df.to_csv(TAB / "merit_index_top20_radar.csv")

print(f"\nOutputs in: {OUT}")
