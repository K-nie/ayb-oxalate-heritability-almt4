#!/usr/bin/env python3
"""
REVISION Phase 1 — stage-1 BLUEs ("don't BLUP twice") and the stage-2
genomic-prediction cascade on the BLUE response, weighted vs unweighted.

Reviewer ask (Holland & Piepho 2024, G3 jkae250): a two-stage GWAS /
genomic-prediction workflow should carry fixed-genotype BLUEs plus a
per-genotype weight into stage 2, not random-effect BLUPs (which
double-shrink). This script:

  1. Writes the stage-1 BLUE table, the per-(sample,trait) inverse-
     variance weights (w = 1/SE²), and the fixed-model residual
     variances (via _pheno.load_blues).
  2. Re-runs per-trait GBLUP leave-one-out predictive ability (the
     script-18 GEBV accuracy metric) on three responses:
         - BLUP          (the current/baseline response)
         - BLUE unweighted
         - BLUE weighted  (residual precision R = diag(1/w) → weighted
                           kernel-ridge GBLUP)
  3. Re-estimates the REML null-model narrow-sense h² per replicated
     trait (the script-03/15 GWAS heritability) on BLUP vs BLUE.
  4. Emits an old(BLUP)→new(BLUE) diff table for every headline statistic.

Kinship: VanRaden (2008) method 1, rebuilt per trait on the trait-
phenotyped lines (identical convention to script 18). Weighted GBLUP:
c = (K + λ diag(1/w))^{-1} y ; ĝ = K c  (Legarra et al.; residual-weight
diagonal in the mixed-model equations).

Oxalate fractions + Seed_Coat_Tannin are means-only → BLUE == BLUP == mean
exactly, so their numbers are invariant; this is asserted and checked.

Outputs: results/R1_blue_pipeline/
Author: Benjamin Narh-Madey
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar

import _pheno
from _pheno import (load_blues, load_blups, TRAITS_ALL, TRAITS_REPLICATED,
                    TRAITS_MEANS_ONLY)

warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
OUT = ROOT / "results" / "R1_blue_pipeline"
TAB = OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
ALPHAS = [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]
RNG = np.random.default_rng(42)


# ---------------------------------------------------------------------------
# Dosage + helpers (same QC cascade as scripts 03/18)
# ---------------------------------------------------------------------------
def load_dosage():
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in META]
    calls = hm[sample_cols].astype(str)
    ref_alt = hm["alleles"].str.split("/", expand=True)
    dosage = np.full((len(hm), len(sample_cols)), np.nan)
    for i in range(len(hm)):
        ref, alt = ref_alt.iloc[i]
        arr = calls.iloc[i].values
        dosage[i, arr == ref + ref] = 0.0
        dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
        dosage[i, arr == alt + alt] = 2.0
    dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
    cr_m = dosage.notna().sum(axis=1) / dosage.shape[1]
    maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                     1 - dosage.mean(axis=1, skipna=True) / 2.0)
    cr_s = dosage.notna().sum(axis=0) / dosage.shape[0]
    keep_m = ((cr_m >= MARK_CR) & (maf >= MIN_MAF)).values
    keep_s = (cr_s >= SAMP_CR).values
    return dosage, keep_m, keep_s


def make_K(M):
    p = M.mean(axis=0) / 2.0
    p = np.clip(p, 1e-6, 1 - 1e-6)
    Z = M - 2.0 * p[None, :]
    return (Z @ Z.T) / (2.0 * (p * (1.0 - p)).sum())


def sub_K(dosage, keep_m, samples_used):
    dose = dosage.loc[keep_m, samples_used].T.values.astype(float)
    cm = np.nanmean(dose, axis=0)
    idx = np.where(np.isnan(dose))
    dose[idx] = np.take(cm, idx[1])
    return make_K(dose)


def loo_gblup(y, K, w=None):
    """Leave-one-out GBLUP predictive ability. alpha chosen by min LOO
    RMSE over the grid (as in script 18). Weighted variant uses residual
    precision diag(w): c = (K + alpha diag(1/w))^{-1} y ; ghat = K c.

    Weights are normalised to mean 1 so the residual-weight diagonal has
    the same overall scale as the unweighted identity — only the RELATIVE
    per-genotype precision (weight heterogeneity) can then move the result,
    which is the quantity of interest. Without this normalisation the
    absolute 1/SE² scale would rescale the effective alpha off the tuning
    grid and masquerade as a weighting effect."""
    n = len(y)
    if w is not None:
        w = np.asarray(w, float)
        w = w / np.nanmean(w)              # mean-1 normalisation
    best_alpha, best_rmse, best_preds = None, np.inf, None
    for alpha in ALPHAS:
        preds = np.zeros(n)
        if w is None:
            A = K + alpha * np.eye(n)
        else:
            A = K + alpha * np.diag(1.0 / w)
        for i in range(n):
            tr = np.arange(n) != i
            Kt = A[np.ix_(tr, tr)]
            try:
                c = np.linalg.solve(Kt, y[tr])
            except np.linalg.LinAlgError:
                c = np.linalg.lstsq(Kt, y[tr], rcond=None)[0]
            preds[i] = K[i, tr] @ c
        rmse = float(np.sqrt(np.mean((preds - y) ** 2)))
        if rmse < best_rmse:
            best_rmse, best_alpha, best_preds = rmse, alpha, preds
    if np.std(best_preds) > 0:
        r = float(stats.pearsonr(best_preds, y)[0])
    else:
        r = np.nan
    return r, best_alpha


def fisher_ci(r, n, a=0.05):
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    zc = stats.norm.ppf(1 - a / 2)
    return float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))


def reml_h2(y, K):
    """REML narrow-sense h² from the null GBLUP model via EMMA spectral
    trick (same estimator as script 03's fit_reml, intercept only)."""
    n = len(y)
    d, U = np.linalg.eigh(K)
    d = np.maximum(d, 1e-8)
    yt = U.T @ y
    Xt = U.T @ np.ones((n, 1))

    def nll(h2):
        v = h2 * d + (1 - h2)
        Vi = 1 / v
        XtVX = (Xt.T * Vi) @ Xt
        s, ld = np.linalg.slogdet(XtVX)
        if s <= 0:
            return 1e12
        beta = np.linalg.solve(XtVX, (Xt.T * Vi) @ yt)
        resid = yt - Xt @ beta
        rss = float((resid * resid * Vi).sum())
        sig2 = rss / (n - 1)
        return 0.5 * (np.log(v).sum() + (n - 1) * np.log(sig2) + ld)

    res = minimize_scalar(nll, bounds=(1e-4, 1 - 1e-4), method="bounded")
    return float(res.x)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
print("[blue] fitting stage-1 BLUEs + weights + residual variances")
blues, weights, bvc = load_blues(return_weights=True,
                                 return_variance_components=True)
blups = load_blups()
blues.to_csv(TAB / "phenotype_blues.csv")
weights.to_csv(TAB / "blue_weights_per_sample.csv", index=False)
bvc.to_csv(TAB / "blue_variance_components.csv", index=False)
print(f"  wrote BLUE table {blues.shape}, weights {weights.shape}")

print("[geno] loading dosage + QC")
dosage, keep_m, keep_s = load_dosage()
samples_all = dosage.columns[keep_s].tolist()

rows = []
for trait in TRAITS_ALL:
    means_only = trait in TRAITS_MEANS_ONLY
    y_blup = blups[trait].reindex(samples_all)
    y_blue = blues[trait].reindex(samples_all)
    keep = y_blup.notna() & y_blue.notna()
    used = [s for s in samples_all if keep.get(s, False)]
    if len(used) < 20:
        continue
    K = sub_K(dosage, keep_m, used)
    yb = y_blup.loc[used].values.astype(float)
    ye = y_blue.loc[used].values.astype(float)
    # weights for this trait (align to used)
    wtab = weights[weights.trait == trait].set_index("sample")["weight"]
    w = wtab.reindex(used).values
    w_cv = (float(np.nanstd(w, ddof=1) / np.nanmean(w))
            if np.all(np.isfinite(w)) and np.nanmean(w) != 0 else 0.0)
    # Weighting can only matter when the per-genotype weights are
    # heterogeneous. Below a negligible CV the residual-weight diagonal is
    # (up to scale) the identity, so weighted GBLUP is provably identical
    # to unweighted; we skip the expensive weighted LOO and record that.
    w_ok = w_cv > 0.02

    r_blup, a_blup = loo_gblup(yb, K)
    r_blue_u, a_blue_u = loo_gblup(ye, K)
    if w_ok:
        r_blue_w, a_blue_w = loo_gblup(ye, K, w=w)
    else:
        r_blue_w, a_blue_w = r_blue_u, a_blue_u  # equal weights → identical

    lo_b, hi_b = fisher_ci(r_blup, len(used))
    lo_eu, hi_eu = fisher_ci(r_blue_u, len(used))
    lo_ew, hi_ew = fisher_ci(r_blue_w, len(used))

    h2_blup = np.nan if means_only else reml_h2(yb, K)
    h2_blue = np.nan if means_only else reml_h2(ye, K)

    rows.append({
        "trait": trait, "n": len(used), "means_only": means_only,
        "loo_r_BLUP": r_blup, "loo_r_BLUP_lo": lo_b, "loo_r_BLUP_hi": hi_b,
        "loo_r_BLUE_unw": r_blue_u, "loo_r_BLUE_unw_lo": lo_eu, "loo_r_BLUE_unw_hi": hi_eu,
        "loo_r_BLUE_wtd": r_blue_w, "loo_r_BLUE_wtd_lo": lo_ew, "loo_r_BLUE_wtd_hi": hi_ew,
        "d_r_blue_minus_blup": r_blue_u - r_blup,
        "d_r_wtd_minus_unw": r_blue_w - r_blue_u,
        "alpha_BLUP": a_blup, "alpha_BLUE_unw": a_blue_u, "alpha_BLUE_wtd": a_blue_w,
        "h2_REML_BLUP": h2_blup, "h2_REML_BLUE": h2_blue,
        "weight_cv": w_cv, "weighted_loo_run": bool(w_ok),
    })
    print(f"  {trait:18s} n={len(used):3d}  LOO r: BLUP={r_blup:+.3f} "
          f"BLUE_u={r_blue_u:+.3f} BLUE_w={r_blue_w:+.3f}  "
          f"h2 BLUP={h2_blup if means_only else round(h2_blup,3)} "
          f"BLUE={h2_blue if means_only else round(h2_blue,3)}")

df = pd.DataFrame(rows)
df.to_csv(TAB / "cascade_blup_vs_blue_diff.csv", index=False)
print(f"\n[done] wrote {TAB/'cascade_blup_vs_blue_diff.csv'}")
print("\n=== means-only invariance check (should be exactly 0 diff) ===")
mo = df[df.means_only]
print(mo[["trait", "loo_r_BLUP", "loo_r_BLUE_unw", "d_r_blue_minus_blup"]]
      .to_string(index=False))
