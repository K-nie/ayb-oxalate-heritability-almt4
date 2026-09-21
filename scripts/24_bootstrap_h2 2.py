#!/usr/bin/env python3
"""
Profile-likelihood 95 % CIs on REML marker-based heritability per AYB trait.

At n = 41-95, single-point REML h^2 estimates are imprecise; reporting a
single decimal masks the fact that the CI typically spans most of [0, 1] for
traits with low marker-trait signal. Per-trait CIs communicate this honestly.

Method
------
For each trait:
  1. Use the EMMAX pipeline as in script 03 (REML h^2 from spectral
     decomposition of K, K = VanRaden method 1).
  2. Compute the REML log-likelihood across a fine h^2 grid (0.001..0.999,
     step 0.005). The 95 % profile-likelihood CI is the set of h^2 values
     whose -2 log L is within chi^2_{1, 0.95} = 3.841 of the minimum
     (equivalent to a likelihood-ratio test at alpha = 0.05).
  3. Report MLE + grid-derived 95 % CI. This is the canonical CI for h^2
     in REML (Lynch & Walsh 1998 chapter 27; Pritchard & Donnelly 2001).

Sample-bootstrap doesn't work here because resampling with replacement
creates duplicate rows in K, which forces REML to the boundary every time.

Outputs (results/24_bootstrap_h2/)
----------------------------------
tables/
    h2_profile_summary.csv          point, lo, hi per trait
    h2_profile_curves.csv           full likelihood grid (long-form)
figures/
    fig61_h2_profile_intervals.png/.pdf
    fig62_h2_profile_curves.png/.pdf  per-trait LRT-band curve

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from _plotstyle import apply, WONG
from _pheno import load_phenotypes, TRAITS_ALL
apply()
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results" / "24_bootstrap_h2"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
H2_GRID = np.concatenate([
    np.linspace(0.001, 0.10, 30),
    np.linspace(0.11, 0.90, 80),
    np.linspace(0.91, 0.999, 30),
])
LRT_THRESHOLD = 3.841 / 2.0  # chi^2_{1, 0.95} -- threshold for 2 * Delta-log-L


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage
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
dose_filt = dosage.loc[keep_m, keep_s]
samples_all = list(dose_filt.columns)
print(f"[load] {len(samples_all)} samples x {dose_filt.shape[0]} markers")


# ---------------------------------------------------------------------------
# Phenotypes
# ---------------------------------------------------------------------------
pheno = load_phenotypes()


# ---------------------------------------------------------------------------
# REML h^2 estimator (same as script 03)
# ---------------------------------------------------------------------------
def make_K(M):
    p = M.mean(axis=0) / 2.0
    p = np.clip(p, 1e-6, 1 - 1e-6)
    Z = M - 2.0 * p[None, :]
    return (Z @ Z.T) / (2.0 * (p * (1.0 - p)).sum())


def reml_nll(h2, y_t, X_t, d):
    n_ = len(y_t); p_ = X_t.shape[1]
    v = h2 * d + (1.0 - h2)
    Vinv = 1.0 / v
    XtVinvX = (X_t.T * Vinv) @ X_t
    sign, logdet = np.linalg.slogdet(XtVinvX)
    if sign <= 0: return 1e12
    XtVinvy = (X_t.T * Vinv) @ y_t
    beta = np.linalg.solve(XtVinvX, XtVinvy)
    resid = y_t - X_t @ beta
    rss = float((resid * resid * Vinv).sum())
    sigma2 = rss / (n_ - p_)
    return 0.5 * (np.log(v).sum() + (n_ - p_) * np.log(sigma2) + logdet)


def h2_reml(y, K):
    eig_vals, U = np.linalg.eigh(K)
    eig_vals = np.maximum(eig_vals, 1e-8)
    y_t = U.T @ y
    X = np.ones((len(y), 1))
    X_t = U.T @ X
    res = minimize_scalar(reml_nll, args=(y_t, X_t, eig_vals),
                          bounds=(1e-4, 1 - 1e-4), method="bounded",
                          options={"xatol": 1e-4})
    return float(res.x)


def h2_profile_ci(y, K, h2_grid=H2_GRID, lrt_thresh=LRT_THRESHOLD):
    """Profile-likelihood 95 % CI on h^2 with boundary diagnostics.

    Returns a dict with:
      h2_mle           grid MLE
      ci_low, ci_high  95 % profile-likelihood interval (closed at the
                       grid edges, not extrapolated past them)
      on_lower_bound   True if the MLE sits at the bottom of the grid
                       (interpret CI as a one-sided upper bound)
      on_upper_bound   True if the MLE sits at the top of the grid
                       (interpret CI as a one-sided lower bound)
      ci_kind          one of {"two_sided", "one_sided_lower",
                       "one_sided_upper", "uninformative"}
      delta            full -2 (log L - max log L) curve

    Boundary behaviour. At small n / low marker count, the REML
    likelihood is often flat at one or both ends of [0, 1]. A grid
    MLE that lands at the edge is not "h^2 ~ 1" or "h^2 ~ 0" — it is
    "the data cannot distinguish that boundary from any nearby value."
    The CI kind flag is what the Methods section should report.
    """
    eig_vals, U = np.linalg.eigh(K)
    eig_vals = np.maximum(eig_vals, 1e-8)
    y_t = U.T @ y
    X = np.ones((len(y), 1))
    X_t = U.T @ X
    nlls = np.array([reml_nll(h2, y_t, X_t, eig_vals) for h2 in h2_grid])
    min_idx = int(np.argmin(nlls))
    h2_mle = float(h2_grid[min_idx])
    delta = 2.0 * (nlls - nlls[min_idx])
    inside = delta < 3.841
    lo_idx = min_idx
    while lo_idx > 0 and inside[lo_idx - 1]:
        lo_idx -= 1
    hi_idx = min_idx
    while hi_idx < len(h2_grid) - 1 and inside[hi_idx + 1]:
        hi_idx += 1

    on_lower = (min_idx == 0) or (h2_mle <= h2_grid[1] and inside[0])
    on_upper = (min_idx == len(h2_grid) - 1) or (
        h2_mle >= h2_grid[-2] and inside[-1])
    if on_lower and on_upper:
        ci_kind = "uninformative"
    elif on_upper:
        ci_kind = "one_sided_lower"
    elif on_lower:
        ci_kind = "one_sided_upper"
    else:
        ci_kind = "two_sided"
    return {
        "h2_mle": h2_mle,
        "ci_low": float(h2_grid[lo_idx]),
        "ci_high": float(h2_grid[hi_idx]),
        "on_lower_bound": bool(on_lower),
        "on_upper_bound": bool(on_upper),
        "ci_kind": ci_kind,
        "delta": delta,
    }


# ---------------------------------------------------------------------------
# Per-trait MLE + profile-likelihood CI
# ---------------------------------------------------------------------------
print(f"\n[h2] computing profile-likelihood CIs per trait...")
rows = []
curves = []
for trait in TRAITS_ALL:
    y_full = pheno[trait].reindex(samples_all).values.astype(float)
    keep_y = ~np.isnan(y_full)
    samples_used = [samples_all[i] for i, k in enumerate(keep_y) if k]
    n = len(samples_used)
    if n < 20:
        print(f"  {trait:>18}: only {n} samples, skipping"); continue

    dose_sub = dose_filt[samples_used].T.values.astype(float)
    col_mean = np.nanmean(dose_sub, axis=0)
    inds = np.where(np.isnan(dose_sub))
    dose_sub[inds] = np.take(col_mean, inds[1])
    maf_sub = np.minimum(dose_sub.mean(axis=0) / 2.0,
                         1.0 - dose_sub.mean(axis=0) / 2.0)
    keep_mark = maf_sub >= MIN_MAF
    M_sub = dose_sub[:, keep_mark]
    K_sub = make_K(M_sub)
    y = y_full[keep_y]

    out = h2_profile_ci(y, K_sub)
    rows.append({"trait": trait, "n_samples": n,
                 "n_markers_subset": int(keep_mark.sum()),
                 "h2_mle": out["h2_mle"],
                 "h2_ci95_low": out["ci_low"],
                 "h2_ci95_high": out["ci_high"],
                 "ci_width": out["ci_high"] - out["ci_low"],
                 "ci_kind": out["ci_kind"],
                 "on_lower_bound": out["on_lower_bound"],
                 "on_upper_bound": out["on_upper_bound"]})
    for h2_val, d in zip(H2_GRID, out["delta"]):
        curves.append({"trait": trait, "h2": float(h2_val),
                       "neg2_delta_logL": float(d)})
    tag = {"two_sided": "two-sided",
            "one_sided_lower": "one-sided lower bound (CI upper edge is the grid limit)",
            "one_sided_upper": "one-sided upper bound (CI lower edge is the grid limit)",
            "uninformative": "uninformative — CI spans the grid"}.get(out["ci_kind"], "")
    print(f"  {trait:>18}: n = {n:>3d}, MLE h^2 = {out['h2_mle']:.3f}, "
          f"95 % CI [{out['ci_low']:.3f}, {out['ci_high']:.3f}]  "
          f"[{tag}]")

df = pd.DataFrame(rows)
df.to_csv(TAB / "h2_profile_summary.csv", index=False)
pd.DataFrame(curves).to_csv(TAB / "h2_profile_curves.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 61. Forest plot of h^2 MLE with profile-likelihood 95 % CI
# ---------------------------------------------------------------------------
PALETTE = {
    "Tannin": "#D55E00", "Phenol": "#0072B2", "Flavonoid": "#009E73",
    "Antioxidant": "#CC79A7", "Seed_Length": "#1f78b4",
    "Seed_Width": "#33a02c", "Seed_Thickness": "#ff7f00",
    "Mass_of_Seeds": "#6a3d9a", "Seed_Coat_Tannin": "#b15928",
    "Crude_Protein": "#7570b3", "Total_Oxalate": "#1b9e77",
    "Soluble_Oxalate": "#d95f02", "Insoluble_Oxalate": "#e7298a",
}
fig, ax = plt.subplots(figsize=(8.5, 6), constrained_layout=True)
df_show = df.sort_values("h2_mle", ascending=True).reset_index(drop=True)
y = np.arange(len(df_show))
for i, row in df_show.iterrows():
    color = PALETTE.get(row["trait"], WONG["blue"])
    # ensure non-negative xerr (clip to MLE in case of grid edge cases)
    xerr_lo = max(0.0, row["h2_mle"] - row["h2_ci95_low"])
    xerr_hi = max(0.0, row["h2_ci95_high"] - row["h2_mle"])
    ax.errorbar(row["h2_mle"], i,
                xerr=[[xerr_lo], [xerr_hi]],
                fmt="o", color=color, ecolor=color,
                markersize=8, capsize=4, elinewidth=1.2)
    ax.text(1.04, i, f"n = {row['n_samples']}", va="center", fontsize=8)
ax.set_yticks(y)
ax.set_yticklabels(df_show["trait"], fontsize=9)
ax.set_xlabel(r"REML marker-based heritability h$^2$")
ax.axvline(0, color="black", linewidth=0.5)
ax.axvline(1, color="black", linewidth=0.5)
ax.set_xlim(-0.05, 1.15)
ax.set_title("REML h$^2$ MLE + 95 % profile-likelihood CI")
ax.grid(True, axis="x")
save(fig, "fig61_h2_profile_intervals")
print("[fig] fig61_h2_profile_intervals")

# Fig 62. Profile-likelihood curves per trait (4 x 4 small multiples)
all_curves = pd.DataFrame(curves)
n_cols = 4
n_rows = int(np.ceil(len(TRAITS_ALL) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.0 * n_cols, 2.2 * n_rows),
                         constrained_layout=True)
for ax, trait in zip(axes.ravel(), TRAITS_ALL):
    sub = all_curves[all_curves["trait"] == trait]
    if sub.empty:
        ax.set_visible(False); continue
    color = PALETTE.get(trait, WONG["blue"])
    ax.plot(sub["h2"], sub["neg2_delta_logL"], linewidth=1.5, color=color)
    ax.axhline(3.841, color="black", linewidth=0.6, linestyle="--",
               label=r"$\chi^2_{1, 0.95}$ = 3.841")
    row = df[df["trait"] == trait].iloc[0]
    ax.axvspan(row["h2_ci95_low"], row["h2_ci95_high"],
               color=color, alpha=0.18)
    ax.axvline(row["h2_mle"], color=color, linewidth=1.2)
    ax.set_xlabel(r"h$^2$")
    ax.set_ylabel(r"$-2\,\Delta \log L$")
    ax.set_title(f"{trait}\nMLE = {row['h2_mle']:.2f}, "
                 f"CI [{row['h2_ci95_low']:.2f}, {row['h2_ci95_high']:.2f}]",
                 fontsize=8.5)
    ax.set_ylim(0, max(8, sub["neg2_delta_logL"].max() + 0.5))
    ax.grid(True)
for ax in axes.ravel()[len(TRAITS_ALL):]:
    ax.set_visible(False)
fig.suptitle("REML profile-likelihood curves for h$^2$ (per trait)",
             fontsize=12)
save(fig, "fig62_h2_profile_curves")
print("[fig] fig62_h2_profile_curves")

print(f"\nOutputs in: {OUT}")
