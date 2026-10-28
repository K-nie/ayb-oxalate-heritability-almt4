#!/usr/bin/env python3
"""
Mixed-model GWAS + genomic prediction for AYB seed-biochemistry traits.

Methods
-------
Kinship matrix K
    VanRaden (2008) method 1:
        Z = M - 2P   (centred dosage; P = column of allele frequencies)
        K = Z Z^T / (2 sum_l p_l (1 - p_l))

Mixed model (per trait y)
    y = Xb + Z u + e,  u ~ N(0, sigma_g^2 K),  e ~ N(0, sigma_e^2 I)
    Spectral decomposition K = U D U^T turns this into n independent
    Gaussian observations with variance (h^2 D_i + (1 - h^2)) * sigma^2.
    REML estimate of h^2 by grid search (EMMA, Kang et al. 2008).

Per-SNP test (EMMAX, Kang et al. 2010)
    Variance components from the null model are reused at every SNP, so the
    per-SNP fit reduces to a generalised least-squares regression with the
    pre-computed V^{-1}. Test statistic for marker effect is t (df = n - p),
    converted to chi^2 with 1 df for genomic-control lambda.

Two model variants
    M1 = MLM with K only.
    M2 = MLM with K + PC1..PC3 as fixed-effect covariates (Q + K).

Multiple-testing
    Bonferroni (alpha / m_eff) and Benjamini-Hochberg FDR.
    Genomic control lambda_GC = median(chi^2_obs) / 0.4549.

Genomic prediction
    GBLUP via kernel ridge regression with the kinship K as the kernel
    (mathematically equivalent to RR-BLUP on the centred dosage matrix).
    5-fold cross-validation x 50 random splits. Predictive ability =
    Pearson r(predicted, observed). Empirical null = same procedure on
    100 random label permutations; one-sided p-value.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
PHENO_XLSX = ROOT / "Updated_Wet_Chemistry_Data.xlsx"
PCA_TAB = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"

OUT = ROOT / "results" / "03_gwas_mlm"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

sns.set_context("talk")
sns.set_style("whitegrid")

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
N_PERM_PRED = 100
N_CV_REPS = 50
N_FOLDS = 5
RNG = np.random.default_rng(0)

TRAITS = ["Tannin", "Phenol", "Flavonoid", "Antioxidant"]


# ---------------------------------------------------------------------------
# Load + filter (mirrors script 01)
# ---------------------------------------------------------------------------
print("[load] reading HapMap...")
hm = pd.read_csv(HAPMAP, low_memory=False)
META_COLS = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
             "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META_COLS]
calls = hm[sample_cols].astype(str)


def dosage_row(row, alleles):
    ref, alt = alleles.split("/")
    out = np.full(len(row), np.nan)
    arr = row.values
    out[arr == ref + ref] = 0.0
    out[(arr == ref + alt) | (arr == alt + ref)] = 1.0
    out[arr == alt + alt] = 2.0
    return out


print("[load] converting calls to dosage...")
dosage = np.vstack([dosage_row(calls.iloc[i], hm["alleles"].iloc[i])
                    for i in range(len(hm))])
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)

call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
mean_dose = dosage.mean(axis=1, skipna=True)
maf = np.minimum(mean_dose / 2.0, 1.0 - mean_dose / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]

keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values

dose_filt = dosage.loc[keep_m, keep_s]
M = dose_filt.T.values.astype(float)
samples = dose_filt.columns.tolist()
markers = dose_filt.index.tolist()
n, m = M.shape
# mean-impute
col_mean = np.nanmean(M, axis=0)
inds = np.where(np.isnan(M))
M[inds] = np.take(col_mean, inds[1])
print(f"[load] n samples = {n}, m markers = {m}")

# allele freq per marker after imputation
p = M.mean(axis=0) / 2.0
p = np.clip(p, 1e-6, 1 - 1e-6)
maf_v = np.minimum(p, 1 - p)

# kinship K (VanRaden method 1)
Z = M - 2.0 * p[None, :]
denom = 2.0 * (p * (1.0 - p)).sum()
K = (Z @ Z.T) / denom
print(f"[kinship] K shape={K.shape}, trace/n={np.trace(K)/n:.3f}")


# ---------------------------------------------------------------------------
# Phenotypes
# ---------------------------------------------------------------------------
pheno = pd.read_excel(PHENO_XLSX, sheet_name="Means")
pheno = pheno.rename(columns={"Genotypes": "sample", "Flavinoid": "Flavonoid"})
pheno = pheno.set_index("sample").reindex(samples)
missing = pheno.isna().all(axis=1).sum()
print(f"[pheno] aligned phenotypes for {n - missing} of {n} samples "
      f"(genotype-only with no phenotype: {missing})")


# ---------------------------------------------------------------------------
# PCs for Q + K variant
# ---------------------------------------------------------------------------
pca_df = pd.read_csv(PCA_TAB).set_index("sample").loc[samples]
PCs = pca_df[[f"PC{i+1}" for i in range(3)]].values  # n x 3 (PC1-PC3, matching §3.3 and script 15)


# ---------------------------------------------------------------------------
# Spectral decomposition for fast EMMA
# ---------------------------------------------------------------------------
print("[emma] eigendecomposition of K...")
eig_vals, U = np.linalg.eigh(K)
# numerical floor
eig_vals = np.maximum(eig_vals, 1e-8)


def reml_nll(h2, y_t, X_t, d):
    """Negative restricted log-likelihood given h^2 in (0,1)."""
    n_ = len(y_t)
    p_ = X_t.shape[1]
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
    nll = 0.5 * (np.log(v).sum() + (n_ - p_) * np.log(sigma2) + logdet_XVX)
    return nll


def fit_reml(y, X, U_, d_):
    """Find REML-optimal h^2 and return variance components and V^{-1}."""
    y_t = U_.T @ y
    X_t = U_.T @ X
    res = minimize_scalar(reml_nll, args=(y_t, X_t, d_),
                          bounds=(1e-4, 1 - 1e-4), method="bounded",
                          options={"xatol": 1e-4})
    h2_hat = float(res.x)
    v = h2_hat * d_ + (1.0 - h2_hat)
    Vinv = 1.0 / v
    XtVinvX = (X_t.T * Vinv) @ X_t
    XtVinvy = (X_t.T * Vinv) @ y_t
    beta = np.linalg.solve(XtVinvX, XtVinvy)
    resid = y_t - X_t @ beta
    rss = float((resid * resid * Vinv).sum())
    sigma2 = rss / (len(y) - X.shape[1])
    Vinv_full = (U_ * Vinv) @ U_.T
    return {"h2": h2_hat, "sigma2": sigma2, "Vinv": Vinv_full,
            "beta_null": beta, "rss_null": rss, "v_diag": v}


def emmax_scan(y, X_null, Vinv, sigma2, M_):
    """Per-SNP GLS scan with pre-computed V^{-1} from the null model."""
    n_ = len(y)
    p_null = X_null.shape[1]
    Vinv_y = Vinv @ y
    Vinv_X = Vinv @ X_null
    XtVinvX_inv = np.linalg.inv(X_null.T @ Vinv_X)
    # Project out null covariates from y and from markers:
    # y_perp = (I - X (X' Vinv X)^{-1} X' Vinv) y    in V-norm sense
    proj = Vinv_X @ XtVinvX_inv @ X_null.T   # n x n  (acting via Vinv on the right)
    # For each marker m_j (n-vector):
    # beta_j = m' Vinv y_perp / (m' Vinv m_perp)
    # where m_perp is m residualised against X under Vinv.
    # Implementation: precompute Vinv_y_perp = Vinv y - Vinv X (X' Vinv X)^{-1} X' Vinv y
    Vinv_y_perp = Vinv_y - Vinv_X @ (XtVinvX_inv @ (X_null.T @ Vinv_y))
    # Per-marker:
    # Mt_Vinv_M_diag = sum over i,j: M_ij * Vinv * M_ij ; here per column
    # We need m_j' Vinv m_j_perp -- equivalently:
    Vinv_M = Vinv @ M_                # n x m_
    Mt_Vinv_X = M_.T @ Vinv_X         # m_ x p_null
    # cross-product correction
    Mt_Vinv_M = np.einsum("ij,ji->i", M_.T, Vinv_M)   # m_
    correction = np.einsum("ij,jk,ik->i", Mt_Vinv_X, XtVinvX_inv, Mt_Vinv_X)
    Mt_Vinv_M_perp = Mt_Vinv_M - correction
    Mt_Vinv_M_perp = np.maximum(Mt_Vinv_M_perp, 1e-12)
    beta_m = (M_.T @ Vinv_y_perp) / Mt_Vinv_M_perp
    se_m = np.sqrt(sigma2 / Mt_Vinv_M_perp)
    t = beta_m / se_m
    df_resid = n_ - p_null - 1
    p_vals = 2.0 * stats.t.sf(np.abs(t), df=df_resid)
    chi2 = t * t
    return beta_m, se_m, t, p_vals, chi2


def bh_fdr(p):
    p = np.asarray(p)
    n_ = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n_ / (np.arange(n_) + 1)
    # enforce monotonicity
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n_)
    out[order] = np.minimum(adj, 1.0)
    return out


# ---------------------------------------------------------------------------
# GWAS per trait, two model variants
# ---------------------------------------------------------------------------
all_hits = {}
trait_stats = []

intercept = np.ones((n, 1))

for trait in TRAITS:
    y = pheno[trait].values.astype(float)
    keep_y = ~np.isnan(y)
    print(f"\n=========== {trait} ({keep_y.sum()} non-missing) ===========")

    # build X for two variants
    Xnull_M1 = intercept
    Xnull_M2 = np.hstack([intercept, PCs])

    for label, Xnull in [("M1_K", Xnull_M1), ("M2_KQ", Xnull_M2)]:
        # subset to non-missing (rare for this dataset but handled cleanly)
        y_s = y[keep_y]
        Xn_s = Xnull[keep_y, :]
        M_s = M[keep_y, :]
        # rebuild eig for sub-sample (needed only if dropping rows; we keep U on the full set)
        # For 105 samples with no missing phenotype values this branch is identical to the full set.
        if keep_y.sum() < n:
            K_sub = K[np.ix_(keep_y, keep_y)]
            d_use, U_use = np.linalg.eigh(K_sub)
            d_use = np.maximum(d_use, 1e-8)
        else:
            d_use, U_use = eig_vals, U
        fit = fit_reml(y_s, Xn_s, U_use, d_use)
        beta_m, se_m, t_m, p_m, chi2_m = emmax_scan(y_s, Xn_s, fit["Vinv"], fit["sigma2"], M_s)
        lam_gc = float(np.median(chi2_m) / 0.4549)
        fdr = bh_fdr(p_m)
        bonf = np.minimum(p_m * m, 1.0)

        df_hits = pd.DataFrame({
            "rs": markers,
            "chrom": hm.set_index("rs#").loc[markers, "chrom"].values,
            "pos": hm.set_index("rs#").loc[markers, "pos"].values,
            "maf": maf_v,
            "beta": beta_m,
            "se": se_m,
            "t": t_m,
            "p": p_m,
            "chi2": chi2_m,
            "fdr_bh": fdr,
            "p_bonferroni": bonf,
        }).sort_values("p")

        df_hits.to_csv(TAB / f"gwas_{trait}_{label}.csv", index=False)
        all_hits[(trait, label)] = df_hits

        n_fdr10 = int((fdr < 0.10).sum())
        n_fdr05 = int((fdr < 0.05).sum())
        n_bonf  = int((bonf < 0.05).sum())
        print(f"  {label}: h2_null={fit['h2']:.3f}, sigma2_e_total={fit['sigma2']:.2f}, "
              f"lambda_GC={lam_gc:.3f}")
        print(f"          hits: Bonf<0.05={n_bonf}, FDR<0.05={n_fdr05}, FDR<0.10={n_fdr10}")

        trait_stats.append({
            "trait": trait, "model": label,
            "h2_REML": fit["h2"], "sigma2_total": fit["sigma2"],
            "lambda_GC": lam_gc, "min_p": float(p_m.min()),
            "n_hits_FDR_0.10": n_fdr10,
            "n_hits_FDR_0.05": n_fdr05,
            "n_hits_Bonf_0.05": n_bonf,
        })

pd.DataFrame(trait_stats).to_csv(TAB / "trait_model_summary.csv", index=False)


# ---------------------------------------------------------------------------
# Per-trait QQ and Manhattan plots (primary model = M1_K)
# ---------------------------------------------------------------------------
print("\n[fig] QQ panels...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, trait in zip(axes.ravel(), TRAITS):
    df = all_hits[(trait, "M1_K")]
    p_sorted = np.sort(df["p"].values)
    expected = -np.log10((np.arange(1, len(p_sorted) + 1) - 0.5) / len(p_sorted))
    observed = -np.log10(p_sorted)
    ax.scatter(expected, observed, s=14, color="#3b6fb6", edgecolor="none")
    lim = max(expected.max(), observed.max()) + 0.3
    ax.plot([0, lim], [0, lim], "k--", lw=1)
    lam = trait_stats[[i for i, r in enumerate(trait_stats) if r["trait"] == trait and r["model"] == "M1_K"][0]]["lambda_GC"]
    ax.set_title(f"{trait}   lambda_GC = {lam:.2f}")
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
fig.suptitle("QQ — MLM (K) per trait", y=1.02)
fig.tight_layout()
fig.savefig(FIG / "qq_panels_M1_K.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "qq_panels_M1_K.pdf", bbox_inches="tight")
plt.close(fig)

print("[fig] Manhattan panels (anchored markers only)...")
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
for ax, trait in zip(axes.ravel(), TRAITS):
    df = all_hits[(trait, "M1_K")].dropna(subset=["chrom", "pos"]).copy()
    if df.empty:
        ax.set_title(f"{trait} — no anchored markers")
        continue
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["chrom_str"] = df["chrom"].astype(str)
    chroms = sorted(df["chrom_str"].unique())
    cum = 0
    xticks, xlabels, xs = [], [], []
    palette = sns.color_palette("Paired", n_colors=max(2, len(chroms)))
    for i, ch in enumerate(chroms):
        sub = df[df["chrom_str"] == ch].sort_values("pos")
        x = sub["pos"].values + cum
        ax.scatter(x, -np.log10(sub["p"].values), s=10, color=palette[i % len(palette)])
        xticks.append(cum + sub["pos"].max() / 2)
        xlabels.append(ch)
        cum += sub["pos"].max() + 1e6
    alpha_bonf = 0.05 / m
    bonf = -np.log10(alpha_bonf)
    ax.axhline(bonf, color="red", ls="--", lw=1,
               label=f"Bonferroni $\\alpha$ = {alpha_bonf:.2e}")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("-log10(p)")
    ax.set_title(f"{trait} — {len(df)} anchored SNPs")
fig.suptitle("Manhattan — MLM (K), cowpea-anchored subset", y=1.01)
fig.tight_layout()
fig.savefig(FIG / "manhattan_panels_M1_K.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "manhattan_panels_M1_K.pdf", bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------------------
# Genomic prediction (GBLUP via kernel ridge with K as kernel)
# ---------------------------------------------------------------------------
# Choose alpha (lambda) for KernelRidge per trait from REML: alpha = (1 - h^2) / h^2
def gblup_cv_accuracy(y_obs, K_mat, n_folds=5, n_reps=50, perm=False, seed=0):
    """5-fold CV x n_reps GBLUP accuracies. Different seeds give independent
    CV partitions (and, when perm=True, independent label shuffles).

    Inner alpha tuning chooses among {0.5, 1, 2, 5, 10, 25} by fold-internal
    Pearson r. Returns one accuracy per CV rep (CV-aggregated across folds).
    """
    rs = np.random.default_rng(seed)
    accs = []
    for r in range(n_reps):
        y_use = y_obs.copy()
        if perm:
            y_use = y_use[rs.permutation(len(y_obs))]
        kf = KFold(n_splits=n_folds, shuffle=True,
                   random_state=int(rs.integers(0, 2**31 - 1)))
        preds = np.full(len(y_obs), np.nan)
        idx = np.arange(len(y_obs))
        for tr, te in kf.split(idx):
            tr_i = idx[tr]; te_i = idx[te]
            best_alpha, best_score = None, -np.inf
            for alpha in [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]:
                kr = KernelRidge(kernel="precomputed", alpha=alpha)
                kr.fit(K_mat[np.ix_(tr_i, tr_i)], y_use[tr_i])
                yhat = kr.predict(K_mat[np.ix_(te_i, tr_i)])
                if np.std(yhat) > 0:
                    sc = stats.pearsonr(yhat, y_use[te_i])[0]
                    if sc > best_score:
                        best_score = sc
                        best_alpha = alpha
            kr = KernelRidge(kernel="precomputed", alpha=best_alpha or 5.0)
            kr.fit(K_mat[np.ix_(tr_i, tr_i)], y_use[tr_i])
            preds[te_i] = kr.predict(K_mat[np.ix_(te_i, tr_i)])
        ok = ~np.isnan(preds)
        if ok.sum() > 5 and np.std(preds[ok]) > 0:
            accs.append(stats.pearsonr(preds[ok], y_use[ok])[0])
    return np.array(accs)


pred_rows = []
per_rep_rows = []
print("\n[gblup] running 5-fold CV x 50 reps per trait + 100-perm null...")
for trait in TRAITS:
    y = pheno[trait].values.astype(float)
    keep_y = ~np.isnan(y)
    y_s = y[keep_y]
    K_s = K[np.ix_(keep_y, keep_y)]
    obs_acc = gblup_cv_accuracy(y_s, K_s, n_folds=N_FOLDS, n_reps=N_CV_REPS,
                                perm=False, seed=42)
    perm_acc = []
    for b in range(N_PERM_PRED):
        a = gblup_cv_accuracy(y_s, K_s, n_folds=N_FOLDS, n_reps=2,
                              perm=True, seed=1000 + b)
        if len(a) > 0:
            perm_acc.append(a.mean())
    perm_acc = np.array(perm_acc)
    p_emp = (1 + np.sum(perm_acc >= np.mean(obs_acc))) / (1 + len(perm_acc))
    row = {
        "trait": trait,
        "mean_r": float(np.mean(obs_acc)),
        "sd_r": float(np.std(obs_acc)),
        "n_reps": len(obs_acc),
        "perm_mean": float(perm_acc.mean()),
        "perm_sd": float(perm_acc.std()),
        "perm_p": float(p_emp),
    }
    pred_rows.append(row)
    for i, r_val in enumerate(obs_acc):
        per_rep_rows.append({"trait": trait, "rep": i,
                              "kind": "observed", "r": float(r_val)})
    for i, r_val in enumerate(perm_acc):
        per_rep_rows.append({"trait": trait, "rep": i,
                              "kind": "permuted_null", "r": float(r_val)})
    print(f"  {trait}: mean r = {row['mean_r']:.3f} +/- {row['sd_r']:.3f}, "
          f"perm p = {row['perm_p']:.3f}")
pred_df = pd.DataFrame(pred_rows)
pred_df.to_csv(TAB / "gblup_prediction_accuracy.csv", index=False)
per_rep_df = pd.DataFrame(per_rep_rows)
per_rep_df.to_csv(TAB / "gblup_per_rep_long.csv", index=False)
print(f"[per-rep] saved {len(per_rep_df)} rows to gblup_per_rep_long.csv")


# prediction figure (paired violins: observed vs permuted-null per trait)
fig, ax = plt.subplots(figsize=(10, 5.5))
xpos = np.arange(len(TRAITS))
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
            body.set_facecolor("#3b6fb6"); body.set_edgecolor("black")
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
ax.set_xticklabels(TRAITS, rotation=20, ha="right")
ax.set_ylabel("Predictive ability (Pearson r)")
ax.set_title(f"GBLUP genomic prediction "
             f"(5-fold CV x {N_CV_REPS} reps, per-rep distribution)")
for i, row in pred_df.iterrows():
    ax.text(xpos[i], max(row["mean_r"], row["perm_mean"]) + 0.04,
            f"p={row['perm_p']:.2g}",
            ha="center", fontsize=10)
ax.bar([np.nan], [np.nan], color="#3b6fb6", edgecolor="black",
       label="Observed")
ax.bar([np.nan], [np.nan], color="#cccccc", edgecolor="black",
       label="Permuted null")
ax.legend(loc="best")
fig.tight_layout()
fig.savefig(FIG / "gblup_prediction_accuracy.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "gblup_prediction_accuracy.pdf", bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------------------
# Summary print
# ---------------------------------------------------------------------------
print("\n=========== TRAIT SUMMARY ===========")
print(pd.DataFrame(trait_stats).to_string(index=False))
print("\n=========== GBLUP PREDICTION ===========")
print(pred_df.to_string(index=False))
print(f"\nOutputs in: {OUT}")
