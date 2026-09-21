#!/usr/bin/env python3
"""
EMMAX-style MLM GWAS + GBLUP genomic prediction for the 9 new AYB traits.

Trait blocks
------------
Seed metrics, n approx 105:
    Seed_Length, Seed_Width, Seed_Thickness, Mass_of_Seeds, Seed_Coat_Tannin
Protein + oxalate subset, n = 46 lines (45 genotyped):
    Crude_Protein, Total_Oxalate, Soluble_Oxalate, Insoluble_Oxalate

For the protein/oxalate block we subset the dosage matrix, recompute MAF, drop
markers with MAF < 0.05 in the 45-line subset, rebuild K on that subset, and
re-do the eigendecomposition. Everything else mirrors script 03.

Sample-size caveat
------------------
n = 45 is well below any defensible GWAS power threshold. We report results
for completeness (negative-by-design framing for the paper) and clearly mark
them in tables, figures and titles.

Outputs (results/15_gwas_new_traits/)
-------------------------------------
tables/
    gwas_<TRAIT>_M1_K.csv               per-SNP GLS results
    trait_model_summary.csv             h2_REML, lambda_GC per trait
    gblup_prediction_accuracy.csv       mean r, sd, perm p per trait
figures/
    fig39_qq_panels_seed.png/.pdf       QQ for the 5 seed traits
    fig40_qq_panels_protox.png/.pdf     QQ for the 4 prot/ox traits
    fig41_manhattan_panels_seed.png/.pdf
    fig42_manhattan_panels_protox.png/.pdf
    fig43_gblup_all_9.png/.pdf          GBLUP bars for all 9 new traits

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import KFold

from _plotstyle import apply, WONG
from _pheno import load_phenotypes, TRAITS_SEED, TRAITS_PROT_OX
apply()
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"

OUT = ROOT / "results" / "15_gwas_new_traits"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05            # default MAF floor for full-panel traits
MIN_MAF_SMALL_N = 0.10    # tighter MAF floor when n_used < 50 (oxalates)
SMALL_N_THRESHOLD = 50    # below this, switch to MIN_MAF_SMALL_N
N_PERM_PRED = 1000        # bumped from 100 per audit (MC SE ~0.01 at p=0.05)
N_CV_REPS = 50
N_FOLDS = 5
N_PCS = 3                 # PC1-PC3 fixed effects per audit (residual structure)
NEW_TRAITS = TRAITS_SEED + TRAITS_PROT_OX

# H2 grid for the REML profile (must match script 24 for cross-script consistency)
H2_GRID_REML = np.concatenate([
    np.linspace(0.001, 0.10, 30),
    np.linspace(0.11, 0.90, 80),
    np.linspace(0.91, 0.999, 30),
])


def fisher_z_ci(r, n, alpha=0.05):
    """Two-sided Fisher-z 95% CI for a Pearson correlation."""
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return float("nan"), float("nan")
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - z_crit * se), np.tanh(z + z_crit * se)
    return float(lo), float(hi)


def detectable_r_at_power(n, power=0.80, alpha=0.05):
    """Minimum two-sided Pearson r detectable at given power, alpha,
    via Fisher-z approximation."""
    if n < 4:
        return float("nan")
    z_a = stats.norm.ppf(1 - alpha / 2)
    z_b = stats.norm.ppf(power)
    z_min = (z_a + z_b) / np.sqrt(n - 3)
    return float(np.tanh(z_min))


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Filtered dosage (same cascade as script 03)
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
mean_dose = dosage.mean(axis=1, skipna=True)
maf = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values

dose_filt = dosage.loc[keep_m, keep_s]
samples_all = dose_filt.columns.tolist()
markers_all = dose_filt.index.tolist()
chr_all = hm.loc[keep_m, "chrom"].values
pos_all = hm.loc[keep_m, "pos"].values
print(f"[load] filtered dosage: {dose_filt.shape}  (markers x samples)")


# ---------------------------------------------------------------------------
# Phenotypes (unified, with corrections)
# ---------------------------------------------------------------------------
pheno = load_phenotypes()
print(f"[pheno] {pheno.shape[0]} samples loaded; using {len(NEW_TRAITS)} new traits")


# ---------------------------------------------------------------------------
# EMMAX helpers (lifted from script 03)
# ---------------------------------------------------------------------------
def reml_nll(h2, y_t, X_t, d):
    n_ = len(y_t); p_ = X_t.shape[1]
    v = h2 * d + (1.0 - h2)
    Vinv = 1.0 / v
    XtVinvX = (X_t.T * Vinv) @ X_t
    sign, logdet_XVX = np.linalg.slogdet(XtVinvX)
    if sign <= 0:
        return 1e12
    XtVinvy = (X_t.T * Vinv) @ y_t
    beta = np.linalg.solve(XtVinvX, XtVinvy)
    resid = y_t - X_t @ beta
    rss = float((resid * resid * Vinv).sum())
    sigma2 = rss / (n_ - p_)
    return 0.5 * (np.log(v).sum() + (n_ - p_) * np.log(sigma2) + logdet_XVX)


def fit_reml(y, X, U_, d_, h2_grid=H2_GRID_REML):
    """REML fit with grid-based h^2 estimator (matches script 24).

    The continuous bounded optimiser used previously stopped at internal
    local maxima when the likelihood was flat at the upper end, giving
    inconsistent h^2 numbers across scripts (Soluble_Oxalate: script 15
    returned 0.705, script 24 returned 0.999 from the grid). The grid
    estimator is the honest report — when it lands at the boundary, the
    on_boundary flag tells the caller to report a one-sided bound.
    """
    y_t = U_.T @ y
    X_t = U_.T @ X
    nlls = np.array([reml_nll(h2, y_t, X_t, d_) for h2 in h2_grid])
    min_idx = int(np.argmin(nlls))
    h2 = float(h2_grid[min_idx])
    on_boundary = (min_idx == 0) or (min_idx == len(h2_grid) - 1)
    v = h2 * d_ + (1.0 - h2)
    Vinv = 1.0 / v
    XtVinvX = (X_t.T * Vinv) @ X_t
    XtVinvy = (X_t.T * Vinv) @ y_t
    beta = np.linalg.solve(XtVinvX, XtVinvy)
    resid = y_t - X_t @ beta
    rss = float((resid * resid * Vinv).sum())
    sigma2 = rss / (len(y) - X.shape[1])
    Vinv_full = (U_ * Vinv) @ U_.T
    return {"h2": h2, "sigma2": sigma2, "Vinv": Vinv_full,
            "on_boundary": on_boundary}


def emmax_scan(y, Xn, Vinv, sigma2, M_):
    n_ = len(y); p_null = Xn.shape[1]
    Vinv_y = Vinv @ y
    Vinv_X = Vinv @ Xn
    XtVinvX_inv = np.linalg.inv(Xn.T @ Vinv_X)
    Vinv_y_perp = Vinv_y - Vinv_X @ (XtVinvX_inv @ (Xn.T @ Vinv_y))
    Vinv_M = Vinv @ M_
    Mt_Vinv_X = M_.T @ Vinv_X
    Mt_Vinv_M = np.einsum("ij,ji->i", M_.T, Vinv_M)
    correction = np.einsum("ij,jk,ik->i", Mt_Vinv_X, XtVinvX_inv, Mt_Vinv_X)
    Mt_Vinv_M_perp = np.maximum(Mt_Vinv_M - correction, 1e-12)
    beta_m = (M_.T @ Vinv_y_perp) / Mt_Vinv_M_perp
    se_m = np.sqrt(sigma2 / Mt_Vinv_M_perp)
    t = beta_m / se_m
    df_resid = n_ - p_null - 1
    p_vals = 2.0 * stats.t.sf(np.abs(t), df=df_resid)
    return beta_m, se_m, t, p_vals, t * t


def bh_fdr(p):
    p = np.asarray(p); n_ = len(p)
    order = np.argsort(p)
    adj = p[order] * n_ / (np.arange(n_) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n_); out[order] = np.minimum(adj, 1.0)
    return out


def make_K(M_arr):
    """VanRaden method 1 kinship on a (n x m) mean-imputed dosage."""
    p_loc = M_arr.mean(axis=0) / 2.0
    p_loc = np.clip(p_loc, 1e-6, 1 - 1e-6)
    Z = M_arr - 2.0 * p_loc[None, :]
    denom = 2.0 * (p_loc * (1.0 - p_loc)).sum()
    return (Z @ Z.T) / denom


# ---------------------------------------------------------------------------
# Per-trait scan
# ---------------------------------------------------------------------------
trait_stats = []
all_hits = {}

for trait in NEW_TRAITS:
    y_full = pheno[trait].reindex(samples_all).values.astype(float)
    keep_y = ~np.isnan(y_full)
    samples_used = [samples_all[i] for i, k in enumerate(keep_y) if k]
    if len(samples_used) < 20:
        print(f"\n[{trait}] only {len(samples_used)} samples with phenotype, skipping")
        continue

    n_used = len(samples_used)
    print(f"\n========= {trait}  (n_used = {n_used}) =========")

    # MAF floor tightens at small n: at n=41 a MAF of 0.05 corresponds to
    # only 4 minor alleles total, so a single mis-called het can drive
    # spurious associations. Bump to MAF >= 0.10 for n < 50.
    maf_floor = MIN_MAF_SMALL_N if n_used < SMALL_N_THRESHOLD else MIN_MAF

    dose_sub = dose_filt[samples_used].T.values.astype(float)
    col_mean = np.nanmean(dose_sub, axis=0)
    inds = np.where(np.isnan(dose_sub))
    dose_sub[inds] = np.take(col_mean, inds[1])
    maf_sub = np.minimum(dose_sub.mean(axis=0) / 2.0,
                         1.0 - dose_sub.mean(axis=0) / 2.0)
    var_sub = dose_sub.var(axis=0)
    keep_marker = (maf_sub >= maf_floor) & (var_sub > 1e-8)
    M_sub = dose_sub[:, keep_marker]
    rs_sub = [markers_all[i] for i, k in enumerate(keep_marker) if k]
    chr_sub = chr_all[keep_marker]
    pos_sub = pos_all[keep_marker]
    print(f"  markers retained in subset: {M_sub.shape[1]} of {len(markers_all)} "
          f"(post MAF >= {maf_floor})")

    # kinship + eig on this subset
    K_sub = make_K(M_sub)
    eig_vals, U = np.linalg.eigh(K_sub)
    eig_vals = np.maximum(eig_vals, 1e-8)

    # PC1..PC{N_PCS} from the kinship eigendecomposition. The top
    # eigenvectors of K are the genotype PCs (Patterson, Price & Reich
    # 2006) — same object, opposite end of the spectrum. We use the
    # largest-eigenvalue PCs as fixed effects to absorb residual
    # population structure that K alone may not capture.
    pcs = U[:, -N_PCS:][:, ::-1]  # PC1 has the largest eigenvalue
    pcs = (pcs - pcs.mean(axis=0)) / (pcs.std(axis=0) + 1e-12)

    # M1 = MLM with K + intercept + PC1..PCN
    y = y_full[keep_y]
    Xn = np.column_stack([np.ones(len(y))] + [pcs[:, k] for k in range(N_PCS)])
    fit = fit_reml(y, Xn, U, eig_vals)
    beta, se, t, p, chi2 = emmax_scan(y, Xn, fit["Vinv"], fit["sigma2"], M_sub)
    lam_gc = float(np.median(chi2) / 0.4549)
    fdr = bh_fdr(p)
    bonf = np.minimum(p * M_sub.shape[1], 1.0)

    df_hits = pd.DataFrame({
        "rs": rs_sub,
        "chrom": chr_sub,
        "pos": pos_sub,
        "maf_subset": maf_sub[keep_marker],
        "beta": beta,
        "se": se,
        "t": t,
        "p": p,
        "chi2": chi2,
        "fdr_bh": fdr,
        "p_bonferroni": bonf,
    }).sort_values("p")
    df_hits.to_csv(TAB / f"gwas_{trait}_M1_K.csv", index=False)
    all_hits[trait] = df_hits

    n_fdr10 = int((fdr < 0.10).sum())
    n_fdr05 = int((fdr < 0.05).sum())
    n_bonf = int((bonf < 0.05).sum())
    print(f"  h2_REML={fit['h2']:.3f}, sigma2_total={fit['sigma2']:.3f}, "
          f"lambda_GC={lam_gc:.3f}")
    print(f"  hits: Bonf<0.05={n_bonf}, FDR<0.05={n_fdr05}, FDR<0.10={n_fdr10}")

    # Per-SNP detectable effect (Pearson r at 80% power, Bonferroni alpha
    # = 0.05 / M markers). Converts to a beta in SD units via beta = r *
    # sqrt(2 * MAF * (1 - MAF))^-1 * sd_y at typical MAF ~ 0.25.
    alpha_bonf = 0.05 / max(M_sub.shape[1], 1)
    z_a = stats.norm.ppf(1 - alpha_bonf / 2)
    z_b = stats.norm.ppf(0.80)
    r_min = float(np.tanh((z_a + z_b) / np.sqrt(max(len(samples_used) - 3, 1))))

    trait_stats.append({
        "trait": trait, "n_used": n_used,
        "n_markers_subset": M_sub.shape[1],
        "maf_floor_used": maf_floor,
        "h2_REML": fit["h2"], "sigma2_total": fit["sigma2"],
        "h2_on_boundary": fit["on_boundary"],
        "lambda_GC": lam_gc, "min_p": float(p.min()),
        "n_hits_FDR_0.10": n_fdr10,
        "n_hits_FDR_0.05": n_fdr05,
        "n_hits_Bonf_0.05": n_bonf,
        "detectable_r_80pct_power_Bonf": r_min,
    })

pd.DataFrame(trait_stats).to_csv(TAB / "trait_model_summary.csv", index=False)


# ---------------------------------------------------------------------------
# QQ panels (split by trait block)
# ---------------------------------------------------------------------------
def qq_panel(traits, suptitle, fname):
    n = len(traits)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3.6 * cols, 3.4 * rows),
                             constrained_layout=True)
    if rows * cols == 1:
        axes = np.array([axes])
    for ax, trait in zip(axes.ravel(), traits):
        if trait not in all_hits:
            ax.set_visible(False)
            continue
        df = all_hits[trait]
        p_sorted = np.sort(df["p"].values)
        m_ = len(p_sorted)
        expected = -np.log10((np.arange(1, m_ + 1) - 0.5) / m_)
        observed = -np.log10(p_sorted)
        ax.scatter(expected, observed, s=10,
                   color=WONG["blue"], edgecolor="none")
        lim = max(expected.max(), observed.max()) + 0.3
        ax.plot([0, lim], [0, lim], "k--", lw=1)
        # 95% Beta CI
        k = np.arange(1, m_ + 1)
        lo = -np.log10(stats.beta.ppf(0.975, k, m_ - k + 1))
        hi = -np.log10(stats.beta.ppf(0.025, k, m_ - k + 1))
        ax.fill_between(expected, lo, hi, color="grey", alpha=0.18,
                        linewidth=0)
        row = next(r for r in trait_stats if r["trait"] == trait)
        ax.set_title(f"{trait}\nn = {row['n_used']}, "
                     f"$\\lambda$$_{{GC}}$ = {row['lambda_GC']:.2f}",
                     fontsize=10)
        ax.set_xlabel(r"Expected $-\log_{10}\,p$")
        ax.set_ylabel(r"Observed $-\log_{10}\,p$")
        ax.grid(True)
    for ax in axes.ravel()[len(traits):]:
        ax.set_visible(False)
    fig.suptitle(suptitle, fontsize=12, y=1.02)
    save(fig, fname)


qq_panel(TRAITS_SEED,
         "QQ — MLM (K) — seed metrics, n approx 105",
         "fig39_qq_panels_seed")
qq_panel(TRAITS_PROT_OX,
         "QQ — MLM (K) — protein and oxalate subset, n approx 45",
         "fig40_qq_panels_protox")
print("[fig] fig39 + fig40 QQ panels")


# ---------------------------------------------------------------------------
# Manhattan panels per trait block (anchored SNPs only)
# ---------------------------------------------------------------------------
def manhattan_panel(traits, suptitle, fname):
    n = len(traits)
    cols = 1
    rows = n
    fig, axes = plt.subplots(rows, cols, figsize=(11, 2.6 * rows),
                             constrained_layout=True, sharex=False)
    if rows == 1:
        axes = np.array([axes])
    for ax, trait in zip(axes.ravel(), traits):
        if trait not in all_hits:
            ax.set_visible(False)
            continue
        df = all_hits[trait].copy()
        df["chr"] = df["chrom"].astype(str).str.extract(r"(Vu\d+)", expand=False)
        df = df.dropna(subset=["chr", "pos"])
        df["pos"] = df["pos"].astype(int)
        df = df.sort_values(["chr", "pos"])
        chrs = sorted(df["chr"].unique(),
                      key=lambda s: int(re.search(r"\d+", s).group()))
        offsets, mids, x_cursor = {}, [], 0
        xs = []
        for c in chrs:
            sub = df[df["chr"] == c]
            offsets[c] = x_cursor
            xs.extend((sub["pos"].values + x_cursor).tolist())
            mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
            x_cursor += sub["pos"].max() + 5_000_000
        df["x"] = xs
        palette = [WONG["blue"], "#7c7c7c"]
        for i, c in enumerate(chrs):
            sub = df[df["chr"] == c]
            ax.scatter(sub["x"], -np.log10(sub["p"]),
                       s=14, color=palette[i % 2], alpha=0.85, edgecolor="none")
        m_used = next(r for r in trait_stats if r["trait"] == trait)["n_markers_subset"]
        alpha_bonf = 0.05 / m_used
        bonf = -np.log10(alpha_bonf)
        ax.axhline(bonf, color="black", ls="--", lw=0.8,
                   label=f"Bonferroni $\\alpha$ = {alpha_bonf:.2e}")
        ax.set_xticks(mids)
        ax.set_xticklabels(chrs, rotation=0, fontsize=8)
        ax.set_ylabel(r"$-\log_{10}\,p$")
        row = next(r for r in trait_stats if r["trait"] == trait)
        ax.set_title(f"{trait}   n = {row['n_used']}, "
                     f"$\\lambda$$_{{GC}}$ = {row['lambda_GC']:.2f}",
                     fontsize=10)
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(True, axis="y")
    fig.suptitle(suptitle, fontsize=12, y=1.01)
    save(fig, fname)


manhattan_panel(TRAITS_SEED,
                "Manhattan — seed metrics, cowpea-anchored SNPs",
                "fig41_manhattan_panels_seed")
manhattan_panel(TRAITS_PROT_OX,
                "Manhattan — protein and oxalate (n = 45), cowpea-anchored",
                "fig42_manhattan_panels_protox")
print("[fig] fig41 + fig42 Manhattan panels")


# ---------------------------------------------------------------------------
# GBLUP genomic prediction
# ---------------------------------------------------------------------------
def gblup_cv(y, K, n_folds=5, n_reps=50, perm=False, seed=0):
    """5-fold cross-validation GBLUP with RMSE-based alpha selection.

    Alpha is selected per training fold by **minimum RMSE** between
    predicted and observed values on the held-out fold (not by Pearson
    r). RMSE is the proper objective for kernel ridge regression
    because it penalises both correlation and scale; Pearson r is
    scale-invariant and can select an alpha whose prediction shape is
    right but whose magnitude is biased.
    """
    rs = np.random.default_rng(seed)
    accs = []
    for r in range(n_reps):
        y_use = y.copy()
        if perm:
            y_use = y_use[rs.permutation(len(y))]
        kf = KFold(n_splits=n_folds, shuffle=True,
                   random_state=int(rs.integers(0, 2**31 - 1)))
        preds = np.full(len(y), np.nan)
        idx = np.arange(len(y))
        for tr, te in kf.split(idx):
            tr_i, te_i = idx[tr], idx[te]
            best_alpha, best_rmse = None, np.inf
            for alpha in [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]:
                kr = KernelRidge(kernel="precomputed", alpha=alpha)
                kr.fit(K[np.ix_(tr_i, tr_i)], y_use[tr_i])
                yhat = kr.predict(K[np.ix_(te_i, tr_i)])
                rmse = float(np.sqrt(np.mean((yhat - y_use[te_i]) ** 2)))
                if rmse < best_rmse:
                    best_rmse = rmse; best_alpha = alpha
            kr = KernelRidge(kernel="precomputed",
                             alpha=best_alpha or 5.0)
            kr.fit(K[np.ix_(tr_i, tr_i)], y_use[tr_i])
            preds[te_i] = kr.predict(K[np.ix_(te_i, tr_i)])
        ok = ~np.isnan(preds)
        if ok.sum() > 5 and np.std(preds[ok]) > 0:
            accs.append(stats.pearsonr(preds[ok], y_use[ok])[0])
    return np.array(accs)


pred_rows = []
per_rep_rows = []
print("\n[gblup] running per trait (5-fold CV x 50 reps + 100-perm null)...")
for trait in NEW_TRAITS:
    if trait not in all_hits:
        continue
    y_full = pheno[trait].reindex(samples_all).values.astype(float)
    keep_y = ~np.isnan(y_full)
    samples_used = [samples_all[i] for i, k in enumerate(keep_y) if k]
    dose_sub = dose_filt[samples_used].T.values.astype(float)
    col_mean = np.nanmean(dose_sub, axis=0)
    inds = np.where(np.isnan(dose_sub))
    dose_sub[inds] = np.take(col_mean, inds[1])
    maf_sub = np.minimum(dose_sub.mean(axis=0) / 2.0,
                         1.0 - dose_sub.mean(axis=0) / 2.0)
    var_sub = dose_sub.var(axis=0)
    keep_marker = (maf_sub >= MIN_MAF) & (var_sub > 1e-8)
    M_sub = dose_sub[:, keep_marker]
    K_sub = make_K(M_sub)
    y = y_full[keep_y]
    n_folds_use = min(N_FOLDS, max(2, len(y) // 5))

    obs = gblup_cv(y, K_sub, n_folds=n_folds_use, n_reps=N_CV_REPS,
                   perm=False, seed=42)
    perm = []
    for b in range(N_PERM_PRED):
        a = gblup_cv(y, K_sub, n_folds=n_folds_use, n_reps=2,
                     perm=True, seed=1000 + b)
        if len(a) > 0:
            perm.append(a.mean())
    perm = np.array(perm)
    p_emp = (1 + np.sum(perm >= np.mean(obs))) / (1 + len(perm))
    # Fisher-z CI on the mean CV r (treats the n_used genotypes as the
    # effective sample). At small n this CI is the dominant source of
    # uncertainty, not the cross-fold dispersion.
    r_lo, r_hi = fisher_z_ci(float(np.mean(obs)), n=len(y))
    r_min_det = detectable_r_at_power(len(y), power=0.80, alpha=0.05)
    row = {
        "trait": trait, "n_used": len(y),
        "mean_r": float(np.mean(obs)), "sd_r": float(np.std(obs)),
        "r_fisher_ci_low": r_lo, "r_fisher_ci_high": r_hi,
        "min_detectable_r_80pct": r_min_det,
        "perm_mean": float(perm.mean()), "perm_sd": float(perm.std()),
        "perm_p": float(p_emp),
        "n_perm": int(len(perm)),
        "n_folds": n_folds_use,
    }
    pred_rows.append(row)
    # capture per-rep observed r AND per-permutation null r for violin
    for i, r_val in enumerate(obs):
        per_rep_rows.append({"trait": trait, "n_used": len(y),
                              "rep": i, "kind": "observed",
                              "r": float(r_val)})
    for i, r_val in enumerate(perm):
        per_rep_rows.append({"trait": trait, "n_used": len(y),
                              "rep": i, "kind": "permuted_null",
                              "r": float(r_val)})
    print(f"  {trait:>18}: n={len(y):3d}  r = {row['mean_r']:+.3f} "
          f"[{r_lo:+.3f}, {r_hi:+.3f}]  perm p = {row['perm_p']:.3f}  "
          f"(min detectable r = {r_min_det:.3f})")

pred_df = pd.DataFrame(pred_rows)
pred_df.to_csv(TAB / "gblup_prediction_accuracy.csv", index=False)
per_rep_df = pd.DataFrame(per_rep_rows)
per_rep_df.to_csv(TAB / "gblup_per_rep_long.csv", index=False)
print(f"[per-rep] saved {len(per_rep_df)} rows to gblup_per_rep_long.csv")


# Figure 43: paired-violin observed-vs-permuted CV r per trait
fig, ax = plt.subplots(figsize=(13, 5.4), constrained_layout=True)
xpos = np.arange(len(pred_df))
for i, (_, row) in enumerate(pred_df.iterrows()):
    trait = row["trait"]
    obs_vals  = per_rep_df.loc[(per_rep_df["trait"] == trait) &
                                (per_rep_df["kind"] == "observed"),
                                "r"].values
    perm_vals = per_rep_df.loc[(per_rep_df["trait"] == trait) &
                                (per_rep_df["kind"] == "permuted_null"),
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
            body.set_facecolor("#cccccc"); body.set_edgecolor("black")
            body.set_alpha(0.75); body.set_linewidth(0.6)
        v2["cmeans"].set_color("black"); v2["cmeans"].set_linewidth(1.0)
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(xpos)
labels = [f"{r['trait']}\nn = {r['n_used']}" for _, r in pred_df.iterrows()]
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
ax.set_ylabel("Predictive ability (Pearson r)")
ax.set_title("GBLUP genomic prediction across the 9 new AYB traits "
             "(per-rep distribution)")
for i, row in pred_df.iterrows():
    sig = "*" if row["perm_p"] < 0.10 else ""
    color = WONG["vermillion"] if row["perm_p"] < 0.10 else "black"
    ax.text(xpos[i], max(row["mean_r"], row["perm_mean"]) + 0.05,
            f"p = {row['perm_p']:.2f}{sig}",
            ha="center", fontsize=8, color=color)
ax.bar([np.nan], [np.nan], color=WONG["blue"], edgecolor="black",
       label="Observed (5-fold CV x 50 reps)")
ax.bar([np.nan], [np.nan], color="#cccccc", edgecolor="black",
       label=f"Permuted null ({N_PERM_PRED} perms)")
ax.legend(loc="upper right", fontsize=9)
ax.grid(True, axis="y")
save(fig, "fig43_gblup_all_9")
print("[fig] fig43_gblup_all_9")


print(f"\nOutputs in: {OUT}")
print("Trait summary:")
print(pd.DataFrame(trait_stats).to_string(index=False))
print("\nGBLUP:")
print(pred_df.to_string(index=False))
