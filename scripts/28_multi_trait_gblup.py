#!/usr/bin/env python3
"""
Multi-trait (multivariate) GBLUP for the AYB panel.

Single-trait GBLUP at n = 41 (protein + oxalate subset) has near-zero
predictive ability. Multi-trait GBLUP can borrow strength across correlated
traits — particularly relevant here because:
  Total_Oxalate x Insoluble_Oxalate    r = 0.91   (insoluble dominates total)
  Soluble_Oxalate x Insoluble_Oxalate  r = -0.38  (antagonistic; breeder-actionable)

Method (Calus & Veerkamp 2011, Jia & Jannink 2012):
  y_t = mu_t + g_t + e_t with g ~ MVN(0, G_genetic x K) and e ~ MVN(0, R x I).
We use a kernelised reduced-rank multi-output ridge with the precomputed
kinship K as the joint covariance basis. Implementation: vectorised
multi-output kernel-ridge with a single trait-by-trait alpha grid; alpha is
chosen by leave-one-out (Pearson r averaged across traits).

For each trait pair / triple of interest, we compare:
  ST = single-trait GBLUP (existing in script 15)
  MT = multi-trait GBLUP across the chosen trait set

Outputs (results/28_multi_trait_gblup/)
---------------------------------------
tables/
    mt_gblup_predictions.csv     per-trait obs vs MT-predicted vs ST-predicted
    mt_gblup_summary.csv         delta-r per trait (MT - ST)
figures/
    fig66_mt_vs_st_gblup.png/.pdf

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
from sklearn.model_selection import KFold

from _plotstyle import apply, WONG
from _pheno import load_phenotypes, TRAITS_PROT_OX
apply()
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results" / "28_multi_trait_gblup"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_FOLDS = 5
N_CV_REPS = 50
N_BOOTSTRAP = 500         # sample-level paired bootstrap for delta-r CI
N_BOOT_INNER_REPS = 5     # CV reps inside each bootstrap (keep small for runtime)
ALPHAS = [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]
RNG = np.random.default_rng(20250517)

# Trait groups for joint analysis. Strongest gain expected for the protein /
# oxalate block because of the strong inter-trait correlations.
TRAIT_GROUPS = {
    "oxalate_triplet": ["Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"],
    "protein_oxalate_full": ["Crude_Protein", "Total_Oxalate",
                              "Soluble_Oxalate", "Insoluble_Oxalate"],
}


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def make_K(M):
    p = M.mean(axis=0) / 2.0
    p = np.clip(p, 1e-6, 1 - 1e-6)
    Z = M - 2.0 * p[None, :]
    return (Z @ Z.T) / (2.0 * (p * (1.0 - p)).sum())


def generate_fold_splits(n, n_folds, n_reps, seed):
    """Pre-generate deterministic KFold partitions so ST and MT can be
    compared on exactly the same fold structure within each rep."""
    rs = np.random.default_rng(seed)
    splits = []
    for r in range(n_reps):
        kf = KFold(n_splits=n_folds, shuffle=True,
                   random_state=int(rs.integers(0, 2**31 - 1)))
        splits.append(list(kf.split(np.arange(n))))
    return splits


def single_trait_gblup_one_rep(y, K, fold_split):
    """One CV-rep ST GBLUP. Alpha selected by minimum-RMSE on the
    held-out fold (proper objective for kernel ridge; Pearson r is
    scale-invariant and can pick an alpha whose prediction is on the
    wrong scale).
    """
    preds = np.full(len(y), np.nan)
    for tr, te in fold_split:
        best_a, best_rmse = None, np.inf
        for a in ALPHAS:
            kr = KernelRidge(kernel="precomputed", alpha=a)
            kr.fit(K[np.ix_(tr, tr)], y[tr])
            yhat = kr.predict(K[np.ix_(te, tr)])
            rmse = float(np.sqrt(np.mean((yhat - y[te]) ** 2)))
            if rmse < best_rmse:
                best_rmse = rmse; best_a = a
        kr = KernelRidge(kernel="precomputed", alpha=best_a or 5.0)
        kr.fit(K[np.ix_(tr, tr)], y[tr])
        preds[te] = kr.predict(K[np.ix_(te, tr)])
    ok = ~np.isnan(preds)
    if ok.sum() > 5 and np.std(preds[ok]) > 0:
        return float(stats.pearsonr(preds[ok], y[ok])[0])
    return float("nan")


def single_trait_gblup_cv(y, K, n_folds=5, n_reps=50, seed=0):
    splits = generate_fold_splits(len(y), n_folds, n_reps, seed)
    return np.array([single_trait_gblup_one_rep(y, K, s) for s in splits])


def multi_trait_gblup_one_rep(Y, K, fold_split):
    """One CV-rep MT GBLUP. Alpha shared across traits, selected by
    minimum mean-RMSE across traits on the held-out fold (proper
    objective; correlation-averaging is scale-invariant and biases the
    chosen alpha).
    """
    n, T = Y.shape
    preds = np.full((n, T), np.nan)
    for tr, te in fold_split:
        best_a, best_rmse = None, np.inf
        for a in ALPHAS:
            kr = KernelRidge(kernel="precomputed", alpha=a)
            kr.fit(K[np.ix_(tr, tr)], Y[tr])
            yhat = kr.predict(K[np.ix_(te, tr)])
            # Per-trait RMSE on the held-out fold, normalised by per-trait
            # SD so traits on different scales contribute comparably
            rmse_sum = 0.0; n_ok = 0
            for t in range(T):
                sd = np.std(Y[tr, t])
                if sd > 0:
                    rmse_sum += np.sqrt(np.mean(
                        ((yhat[:, t] - Y[te, t]) / sd) ** 2))
                    n_ok += 1
            rmse_mean = rmse_sum / n_ok if n_ok > 0 else np.inf
            if rmse_mean < best_rmse:
                best_rmse = rmse_mean; best_a = a
        kr = KernelRidge(kernel="precomputed", alpha=best_a or 5.0)
        kr.fit(K[np.ix_(tr, tr)], Y[tr])
        preds[te] = kr.predict(K[np.ix_(te, tr)])
    out = np.full(T, np.nan)
    for t in range(T):
        ok = ~np.isnan(preds[:, t])
        if ok.sum() > 5 and np.std(preds[ok, t]) > 0:
            out[t] = stats.pearsonr(preds[ok, t], Y[ok, t])[0]
    return out


def multi_trait_gblup_cv(Y, K, n_folds=5, n_reps=50, seed=0):
    """Multi-output kernel-ridge regression.

    Y is n x T (samples x traits). The same alpha applies across traits;
    sklearn KernelRidge with a vector target solves T independent ridges
    sharing the same K, which is the closed-form analytical multi-trait
    GBLUP solution when the genetic covariance G is taken to be the
    identity scaled by alpha (equal heritability across traits). For
    unequal heritability this is a defensible approximation; the real
    win at small n is the shared regularisation, which lets shrinkage
    borrow strength across traits.
    """
    splits = generate_fold_splits(Y.shape[0], n_folds, n_reps, seed)
    return np.vstack([multi_trait_gblup_one_rep(Y, K, s) for s in splits])


def paired_st_mt_cv(Y, K, traits, n_folds=5, n_reps=50, seed=0):
    """Paired ST vs MT CV — both passes run on the same fold partitions
    so per-rep delta-r is a proper within-pair difference, not a
    difference of independent samples. Returns (st_acc, mt_acc, delta_acc)
    each of shape (n_reps, T).
    """
    n, T = Y.shape
    splits = generate_fold_splits(n, n_folds, n_reps, seed)
    st_mat = np.full((n_reps, T), np.nan)
    mt_mat = np.full((n_reps, T), np.nan)
    for r, fold_split in enumerate(splits):
        # ST: one independent ridge per trait on this fold split
        for j in range(T):
            st_mat[r, j] = single_trait_gblup_one_rep(Y[:, j], K, fold_split)
        # MT: joint ridge across all T traits on the same fold split
        mt_mat[r] = multi_trait_gblup_one_rep(Y, K, fold_split)
    delta_mat = mt_mat - st_mat
    return st_mat, mt_mat, delta_mat


def paired_bootstrap_delta_r(Y, K, traits, n_folds=5, n_reps=10,
                              n_boot=1000, seed=0):
    """Sample-level paired bootstrap for the mean delta-r per trait.

    Resamples genotypes (with replacement) at the sample level — this
    captures the dominant source of uncertainty at small n. For each
    bootstrap resample the kinship and the phenotype matrix are rebuilt
    on the bootstrap sample, both ST and MT are CV-fit on a small number
    of reps (n_reps small to keep runtime reasonable), and the mean
    delta-r per trait is recorded. The bootstrap distribution gives a
    2.5/97.5-percentile CI on the mean delta-r.
    """
    rs = np.random.default_rng(seed + 7919)
    n, T = Y.shape
    boots = np.full((n_boot, T), np.nan)
    for b in range(n_boot):
        sel = rs.integers(0, n, size=n)
        Y_b = Y[sel]
        K_b = K[np.ix_(sel, sel)]
        # add a tiny ridge to K_b to ensure invertibility — bootstrap
        # samples duplicate rows so K_b is singular without it
        K_b = K_b + 1e-6 * np.eye(n)
        _, _, delta = paired_st_mt_cv(Y_b, K_b, traits,
                                       n_folds=n_folds,
                                       n_reps=n_reps,
                                       seed=seed + b)
        boots[b] = np.nanmean(delta, axis=0)
    return boots


# ---------------------------------------------------------------------------
# Load + filter
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
samples_all = dose_filt.columns.tolist()
pheno = load_phenotypes()


# ---------------------------------------------------------------------------
# Run per-group MT vs ST comparison
# ---------------------------------------------------------------------------
all_rows = []
all_per_rep_rows = []
for group_name, traits in TRAIT_GROUPS.items():
    print(f"\n=== {group_name}: traits = {traits} ===")
    sub = pheno[traits].dropna(how="any")
    sub_samples = [s for s in sub.index if s in samples_all]
    sub = sub.loc[sub_samples]
    n = len(sub_samples)
    if n < 20:
        print(f"  too few samples ({n}); skipping group"); continue
    # rebuild K on this trait-subset panel
    dose_sub = dose_filt[sub_samples].T.values.astype(float)
    col_mean = np.nanmean(dose_sub, axis=0)
    inds = np.where(np.isnan(dose_sub))
    dose_sub[inds] = np.take(col_mean, inds[1])
    maf_sub = np.minimum(dose_sub.mean(axis=0) / 2.0,
                         1.0 - dose_sub.mean(axis=0) / 2.0)
    M_sub = dose_sub[:, maf_sub >= MIN_MAF]
    K_sub = make_K(M_sub)
    Y = sub.values
    # Paired ST vs MT CV — both passes on the same fold partitions.
    print(f"  n = {n}; paired ST vs MT CV (n_reps={N_CV_REPS}, same folds)...")
    st_mat, mt_mat, delta_mat = paired_st_mt_cv(
        Y, K_sub, traits, n_folds=N_FOLDS, n_reps=N_CV_REPS, seed=42)

    # Sample-level paired bootstrap on the mean delta-r.
    print(f"  paired bootstrap on delta-r (n_boot={N_BOOTSTRAP}, "
          f"n_reps_inner=10)...")
    boot_delta = paired_bootstrap_delta_r(
        Y, K_sub, traits, n_folds=N_FOLDS, n_reps=N_BOOT_INNER_REPS,
        n_boot=N_BOOTSTRAP, seed=42)

    # Persist per-rep raw r values for downstream violin rendering
    for j, t in enumerate(traits):
        for r_idx in range(st_mat.shape[0]):
            all_per_rep_rows.append({
                "group": group_name, "trait": t, "n": n,
                "rep": r_idx,
                "ST_r": float(st_mat[r_idx, j]),
                "MT_r": float(mt_mat[r_idx, j]),
                "delta_r": float(delta_mat[r_idx, j]),
            })

    for j, t in enumerate(traits):
        st = st_mat[:, j]; mt = mt_mat[:, j]; delta = delta_mat[:, j]
        st_ok = st[~np.isnan(st)]; mt_ok = mt[~np.isnan(mt)]
        delta_ok = delta[~np.isnan(delta)]
        # Paired t-test on delta-r across CV reps
        if len(delta_ok) > 2 and np.std(delta_ok) > 0:
            t_stat, p_paired = stats.ttest_1samp(delta_ok, 0.0)
        else:
            t_stat, p_paired = float("nan"), float("nan")
        # Bootstrap CI on the mean delta-r
        bj = boot_delta[:, j]
        bj = bj[~np.isnan(bj)]
        if len(bj) > 10:
            boot_lo, boot_hi = float(np.percentile(bj, 2.5)), float(np.percentile(bj, 97.5))
        else:
            boot_lo, boot_hi = float("nan"), float("nan")
        verdict = ("benefit"
                    if boot_lo > 0
                    else "penalty" if boot_hi < 0
                    else "no detectable benefit")
        print(f"    {t:>18}: ST r={st_ok.mean():+.3f}, "
              f"MT r={mt_ok.mean():+.3f}, delta r={delta_ok.mean():+.3f}  "
              f"95 % bootstrap CI [{boot_lo:+.3f}, {boot_hi:+.3f}]  "
              f"paired-t p={p_paired:.3f}  -> {verdict}")
        all_rows.append({"group": group_name, "trait": t, "n": n,
                         "ST_mean_r": float(st_ok.mean()),
                         "ST_sd_r": float(st_ok.std()),
                         "MT_mean_r": float(mt_ok.mean()),
                         "MT_sd_r": float(mt_ok.std()),
                         "delta_r": float(delta_ok.mean()),
                         "delta_r_sd_across_reps": float(delta_ok.std()),
                         "delta_r_boot_ci95_low": boot_lo,
                         "delta_r_boot_ci95_high": boot_hi,
                         "paired_t_p": float(p_paired),
                         "verdict": verdict,
                         "MT_better": bool(delta_ok.mean() > 0)})

df = pd.DataFrame(all_rows)
df.to_csv(TAB / "mt_gblup_summary.csv", index=False)
per_rep_df = pd.DataFrame(all_per_rep_rows)
per_rep_df.to_csv(TAB / "mt_gblup_per_rep_long.csv", index=False)
print(f"[per-rep] saved {len(per_rep_df)} rows to mt_gblup_per_rep_long.csv")
print("\n[summary]")
print(df.round(3).to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 66. Side-by-side ST vs MT predictive ability per trait
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.6), constrained_layout=True,
                         gridspec_kw={"width_ratios": [1.0, 0.9]})

# Panel A: ST vs MT mean r per trait, with per-rep SD as the error bar
ax = axes[0]
labels = [f"{r['trait'].replace('_', ' ')}\n({r['group']}, n={r['n']})"
          for _, r in df.iterrows()]
x = np.arange(len(df))
# Per-rep distributions as paired violins (ST = grey, MT = blue) — the bars
# previously hid the within-trait CV variance; violins expose it.
for i, (_, row) in enumerate(df.iterrows()):
    sub = per_rep_df[(per_rep_df["group"] == row["group"]) &
                     (per_rep_df["trait"] == row["trait"])]
    st_vals = sub["ST_r"].dropna().values
    mt_vals = sub["MT_r"].dropna().values
    if len(st_vals) > 1:
        v1 = ax.violinplot([st_vals], positions=[i - 0.20],
                           widths=0.34, showmeans=True, showmedians=False,
                           showextrema=False)
        for body in v1["bodies"]:
            body.set_facecolor("#7c7c7c"); body.set_edgecolor("black")
            body.set_alpha(0.7); body.set_linewidth(0.6)
        v1["cmeans"].set_color("black"); v1["cmeans"].set_linewidth(1.0)
    if len(mt_vals) > 1:
        v2 = ax.violinplot([mt_vals], positions=[i + 0.20],
                           widths=0.34, showmeans=True, showmedians=False,
                           showextrema=False)
        for body in v2["bodies"]:
            body.set_facecolor(WONG["blue"]); body.set_edgecolor("black")
            body.set_alpha(0.75); body.set_linewidth(0.6)
        v2["cmeans"].set_color("black"); v2["cmeans"].set_linewidth(1.0)
# Legend proxies for the violin fills
from matplotlib.patches import Patch
ax.bar([np.nan], [np.nan], color="#7c7c7c", edgecolor="black",
       label="Single-trait GBLUP (per-rep distribution)")
ax.bar([np.nan], [np.nan], color=WONG["blue"], edgecolor="black",
       label="Multi-trait GBLUP (per-rep distribution)")
ax.axhline(0, color="black", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8.0)
plt.setp(ax.get_xticklabels(), rotation_mode="anchor")
ax.set_ylabel("Predictive ability (Pearson r)")
ax.set_title("(A) ST vs MT predictive ability per trait")
ax.legend(framealpha=0.9, edgecolor="none")
ax.grid(True, axis="y", alpha=0.25)

# Panel B: paired delta-r with bootstrap 95 % CI
ax = axes[1]
y = np.arange(len(df))
for i, (_, row) in enumerate(df.iterrows()):
    lo = row["delta_r_boot_ci95_low"]
    hi = row["delta_r_boot_ci95_high"]
    color = (WONG["vermillion"]
              if (not np.isnan(lo) and lo > 0)
              else WONG["blue"]
              if (not np.isnan(hi) and hi < 0)
              else "#7c7c7c")
    ax.errorbar(row["delta_r"], i,
                 xerr=[[max(row["delta_r"] - lo, 0)],
                        [max(hi - row["delta_r"], 0)]],
                 fmt="o", color=color, ecolor=color,
                 markersize=7, capsize=4, elinewidth=1.4)
    ax.text(hi + 0.02 if hi > row["delta_r"] else row["delta_r"] + 0.02,
             i, row["verdict"], va="center", fontsize=7, color=color)
ax.axvline(0, color="black", linewidth=0.6)
ax.set_yticks(y)
ax.set_yticklabels([f"{r['trait'].replace('_', ' ')}\n{r['group']}"
                     for _, r in df.iterrows()],
                    fontsize=8.5)
ax.set_xlabel(r"$\Delta r$ (MT $-$ ST)  with 95 % paired-bootstrap CI")
ax.set_title("(B) Paired-bootstrap CI on the MT minus ST predictive-ability gain")
ax.grid(True, axis="x", alpha=0.25)

save(fig, "fig66_mt_vs_st_gblup")
print("[fig] fig66_mt_vs_st_gblup")

print(f"\nOutputs in: {OUT}")
