#!/usr/bin/env python3
"""
Pairwise LD (r^2) decay across the cowpea-anchored AYB SNPs.

1. Restrict to filtered markers with a v1.0 chromosomal position (the only
   ones for which physical distance is meaningful).
2. Compute pairwise r^2 within each chromosome for all SNP pairs <= MAX_DIST_BP.
3. Bin pairs by physical distance, compute mean r^2 and 95 % CI in each bin.
4. Fit the Hill & Weir (1988) expectation under recombination/drift:
       E[r^2] = (10 + C) / ((2 + C)(11 + C))   ,  C = 4 Ne c
   where c is recombination distance in Morgans. Treating c = beta * d for an
   unknown rho = 4 Ne, fit beta by non-linear least squares against the
   distance-binned means. The half-decay distance (where E[r^2] = max/2)
   summarises decay rate.
5. Output the per-bin curve + the model fit, plus an Ne estimate assuming
   1 cM per Mb in cowpea (Lonardi 2019, ~1.06 cM/Mb panel-wide).

Outputs (results/11_ld_decay/)
------------------------------
tables/
    ld_pairs.csv                         all per-pair r^2 by chrom + distance
    ld_decay_bins.csv                    distance bin, mean / sd / n
    ld_summary.csv                       half-decay BP, rho, Ne (1 cM/Mb)
figures/
    fig27_ld_decay.png/.pdf              binned r^2 + Hill-Weir fit
    fig28_ld_decay_per_chrom.png/.pdf    11 small multiples

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results" / "11_ld_decay"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
MAX_DIST_BP = 5_000_000
N_BINS = 40
CM_PER_MB = 1.06   # Lonardi et al. 2019 cowpea panel-wide average


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered, anchored dosage
# ---------------------------------------------------------------------------
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
dos_f = dosage.loc[keep_m, keep_s]

chr_norm = hm.loc[keep_m, "chrom"].astype(str).str.extract(
    r"(Vu\d+)", expand=False).values
pos = pd.to_numeric(hm.loc[keep_m, "pos"], errors="coerce").values
anc_mask = pd.notna(chr_norm) & np.isfinite(pos)
dos_anc = dos_f.iloc[anc_mask].copy()
dos_anc = dos_anc.fillna(dos_anc.mean(axis=1))    # mean-impute for r^2 calc
chr_anc = chr_norm[anc_mask]
pos_anc = pos[anc_mask].astype(int)
rs_anc = dos_anc.index.tolist()
print(f"[load] anchored markers: {len(rs_anc)}")


# ---------------------------------------------------------------------------
# Pairwise r^2 within each chromosome, up to MAX_DIST_BP
# ---------------------------------------------------------------------------
print(f"[ld] computing pairwise r^2 within chromosomes "
      f"(<= {MAX_DIST_BP / 1e6:.1f} Mb)...")
records = []
for c in sorted(set(chr_anc), key=lambda s: int(re.search(r"\d+", s).group())):
    idx = np.where(chr_anc == c)[0]
    if len(idx) < 2:
        continue
    p = pos_anc[idx]
    order = np.argsort(p)
    idx = idx[order]; p = p[order]
    # standardise dosages for correlation: r = z_i^T z_j / (n-1) (for two SNPs)
    X = dos_anc.values[idx]            # markers x samples
    Xz = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True,
                                                     ddof=1)
    Xz = np.nan_to_num(Xz, nan=0.0)
    L, N = Xz.shape
    for i in range(L - 1):
        # only walk forward and stop when distance exceeds MAX_DIST_BP
        for j in range(i + 1, L):
            d = p[j] - p[i]
            if d > MAX_DIST_BP:
                break
            r = (Xz[i] * Xz[j]).sum() / (N - 1)
            records.append((c, rs_anc[idx[i]], rs_anc[idx[j]],
                            int(p[i]), int(p[j]), int(d), r * r))
    print(f"  {c}: {len(idx)} SNPs")
pairs = pd.DataFrame(records, columns=["chr", "rs_i", "rs_j",
                                        "pos_i", "pos_j",
                                        "dist_bp", "r2"])
pairs.to_csv(TAB / "ld_pairs.csv", index=False)
print(f"[ld] {len(pairs):,} SNP pairs computed")


# ---------------------------------------------------------------------------
# Bin by distance, mean r^2 with 95 % CI per bin
# ---------------------------------------------------------------------------
edges = np.linspace(0, MAX_DIST_BP, N_BINS + 1)
mids = 0.5 * (edges[:-1] + edges[1:])
bin_idx = np.digitize(pairs["dist_bp"], edges, right=False) - 1
pairs["bin"] = bin_idx
bins = (pairs.groupby("bin")["r2"]
              .agg(["mean", "std", "count"])
              .reindex(range(N_BINS))
              .fillna({"mean": np.nan, "std": np.nan, "count": 0}))
bins["mid_bp"] = mids
bins["se"] = bins["std"] / np.sqrt(bins["count"].replace(0, np.nan))
bins.to_csv(TAB / "ld_decay_bins.csv")


# ---------------------------------------------------------------------------
# Hill & Weir (1988) fit
# ---------------------------------------------------------------------------
# Reparameterise C = rho * d_morgans = rho * (cm_per_mb / 100) * d_bp
# Fit rho by least squares on the bin means.
cm_per_bp = CM_PER_MB / 1e6
morgans_per_bp = cm_per_bp / 100.0


def hill_weir(d_bp, rho):
    C = rho * morgans_per_bp * d_bp
    return (10.0 + C) / ((2.0 + C) * (11.0 + C))


valid = bins["mean"].notna()
xfit = bins.loc[valid, "mid_bp"].values
yfit = bins.loc[valid, "mean"].values
try:
    popt, _ = curve_fit(hill_weir, xfit, yfit, p0=[100.0],
                        bounds=(1e-3, 1e6), maxfev=10000)
    rho_hat = float(popt[0])
    print(f"[fit] Hill-Weir rho = 4 Ne (in Morgans^-1) = {rho_hat:.2f}")
    Ne_hat = rho_hat / 4.0
    # E[r^2] at d=0 is 10/22 = 0.4545 (under H&W); decay-to-half = where
    # E[r^2] = max/2. Solve numerically over a dense grid.
    d_grid = np.linspace(1, MAX_DIST_BP, 200_000)
    y_grid = hill_weir(d_grid, rho_hat)
    half_target = y_grid[0] / 2.0
    half_idx = np.argmin(np.abs(y_grid - half_target))
    half_bp = float(d_grid[half_idx])
    print(f"[fit] half-decay distance = {half_bp:.0f} bp "
          f"({half_bp / 1e3:.1f} kb)")
    # Distance at which r^2 < 0.1
    r2_01 = d_grid[np.argmax(y_grid < 0.1)] if (y_grid < 0.1).any() else np.nan
    print(f"[fit] distance at which r^2 < 0.1 = {r2_01:.0f} bp")
except Exception as e:
    print(f"[fit] failed: {e}")
    rho_hat = Ne_hat = half_bp = r2_01 = np.nan


pd.DataFrame([{
    "n_pairs": len(pairs),
    "max_dist_bp": MAX_DIST_BP,
    "rho_hat": rho_hat,
    "Ne_estimate": Ne_hat,
    "half_decay_bp": half_bp,
    "dist_r2_below_0.1_bp": r2_01,
    "cm_per_mb_assumed": CM_PER_MB,
}]).to_csv(TAB / "ld_summary.csv", index=False)


# ---------------------------------------------------------------------------
# Fig 27. Global decay + fit
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.0, 4.4), constrained_layout=True)
ax.scatter(pairs["dist_bp"] / 1e3, pairs["r2"], s=4, alpha=0.06,
           color=WONG["skyblue"], edgecolor="none", rasterized=True,
           label="per-pair r$^2$")
ax.errorbar(bins["mid_bp"] / 1e3, bins["mean"],
            yerr=1.96 * bins["se"],
            fmt="o", color=WONG["blue"], markersize=4,
            elinewidth=0.8, capsize=2, label="bin mean +/- 95 % CI")
if np.isfinite(rho_hat):
    d_plot = np.linspace(1, MAX_DIST_BP, 1000)
    ax.plot(d_plot / 1e3, hill_weir(d_plot, rho_hat),
            color=WONG["vermillion"], linewidth=1.6,
            label=f"Hill-Weir fit, $\\hat\\rho$ = {rho_hat:.1f}; "
                  f"half-decay = {half_bp / 1e3:.0f} kb")
ax.set_xlabel("physical distance (kb)")
ax.set_ylabel(r"r$^2$")
ax.set_title(f"LD decay across cowpea-anchored AYB SNPs "
             f"({len(pairs):,} pairs)")
ax.set_xlim(0, MAX_DIST_BP / 1e3)
ax.set_ylim(0, 1.02)
ax.legend(loc="upper right", fontsize=8)
ax.grid(True)
save(fig, "fig27_ld_decay")
print("[fig] fig27_ld_decay")


# ---------------------------------------------------------------------------
# Fig 28. Per-chromosome small multiples
# ---------------------------------------------------------------------------
chrs = sorted(pairs["chr"].unique(),
              key=lambda s: int(re.search(r"\d+", s).group()))
ncol = 4
nrow = int(np.ceil(len(chrs) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.4 * nrow),
                         constrained_layout=True, sharex=True, sharey=True)
for ax, c in zip(axes.ravel(), chrs):
    sub = pairs[pairs["chr"] == c]
    ax.scatter(sub["dist_bp"] / 1e3, sub["r2"], s=6, alpha=0.15,
               color=WONG["skyblue"], edgecolor="none", rasterized=True)
    if len(sub) > 50:
        edges_c = np.linspace(0, MAX_DIST_BP, 21)
        mids_c = 0.5 * (edges_c[:-1] + edges_c[1:])
        bin_c = np.digitize(sub["dist_bp"], edges_c) - 1
        b = sub.groupby(bin_c)["r2"].mean()
        b = b.reindex(range(20))
        ax.plot(mids_c / 1e3, b.values, "o-", color=WONG["blue"],
                markersize=3, linewidth=1)
    ax.set_title(c)
    ax.grid(True)
for ax in axes.ravel()[len(chrs):]:
    ax.set_visible(False)
for ax in axes[-1, :]:
    ax.set_xlabel("distance (kb)")
for ax in axes[:, 0]:
    ax.set_ylabel("r$^2$")
fig.suptitle("LD decay per cowpea pseudo-chromosome",
             y=1.01, fontsize=12)
save(fig, "fig28_ld_decay_per_chrom")
print("[fig] fig28_ld_decay_per_chrom")

print(f"\nOutputs in: {OUT}")
