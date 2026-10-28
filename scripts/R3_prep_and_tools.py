#!/usr/bin/env python3
"""
REVISION Phase 3 (b,c,d) + data export for 3a — breeding-tool fixes on the
BLUE response.

Reviewer asks:
 (b) Usefulness criterion must include accuracy: UC = mu + i * h * sigma_G
     (Schnell & Utz 1975), where h = per-trait prediction accuracy (CV/LOO
     r, clamped at 0). The current script-64 UC = mu + i*sigma_G omits h.
 (c) Smith-Hazel P must be built from PHENOTYPIC (co)variances, not a
     BLUP-correlation proxy, and economic weights expressed per genetic SD
     so unstandardised oxalate traits don't dominate by scale. Re-derive
     b = P^-1 G a and the rankings.
 (d) Crude_Protein has CV r < 0 (worse than random). Produce merit-index +
     OCS variants that EXCLUDE negative-CV traits, alongside the full ones,
     so Benjamin can choose.

Also exports the trait-subset data + kinship that R3a_mtm_sommer.R consumes
for the proper multivariate GBLUP.

Inputs : BLUEs (results/R1_blue_pipeline), per-trait LOO accuracy
         (results/R1_blue_pipeline cascade), GRM (results/13), var
         components (results/00), cross pairs (results/13), GEBVs
         (results/18, re-signed merit from BLUEs here).
Output : results/R3_breeding_tools/
Author : Benjamin Narh-Madey
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge

import _pheno
from _pheno import TRAITS_ALL, TRAITS_REPLICATED, TRAITS_MEANS_ONLY
warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
GRM_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"
CROSS_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "cross_pairs_least_related.csv"
BLUP_VC = ROOT / "results" / "00_blup_pipeline" / "tables" / "variance_components.csv"
CASCADE = ROOT / "results" / "R1_blue_pipeline" / "tables" / "cascade_blup_vs_blue_diff.csv"
OUT = ROOT / "results" / "R3_breeding_tools"
TAB = OUT / "tables"; FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

ALPHAS = [0.5, 1.0, 2.0, 5.0, 10.0, 25.0]
SELECTION_INTENSITY = 1.755     # top 10% under normality (Falconer & Mackay)

# Direction-signed traits (+ want high, - want low); seed-size dims neutral.
TRAIT_DIR = {
    "Tannin": -1, "Phenol": +1, "Flavonoid": +1, "Antioxidant": +1,
    "Mass_of_Seeds": +1, "Seed_Coat_Tannin": -1, "Crude_Protein": +1,
    "Total_Oxalate": -1, "Soluble_Oxalate": -1, "Insoluble_Oxalate": +1,
}
INDEX_TRAITS = list(TRAIT_DIR.keys())

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
blues = _pheno.load_blues()                      # BLUE response
grm = pd.read_csv(GRM_CSV, index_col=0)
vc_blup = pd.read_csv(BLUP_VC).set_index("trait")
cascade = pd.read_csv(CASCADE).set_index("trait")

# Per-trait prediction accuracy h = LOO r (BLUE, unweighted), clamped >= 0.
acc = cascade["loo_r_BLUE_unw"].to_dict()
h_acc = {t: max(float(acc.get(t, np.nan)), 0.0) if pd.notna(acc.get(t, np.nan)) else 0.0
         for t in TRAITS_ALL}
neg_cv_traits = sorted([t for t in TRAITS_ALL
                        if pd.notna(acc.get(t, np.nan)) and acc[t] < 0])
print(f"[acc] negative-CV traits (LOO r<0, to exclude in 3d): {neg_cv_traits}")

# Per-trait genetic variance: var_g from the BLUP REML where present;
# else 0.5 * phenotypic (BLUE) variance under prior h2=0.5 (script-70 rule).
var_g = {}
for t in TRAITS_ALL:
    if t in vc_blup.index and pd.notna(vc_blup.loc[t, "var_g"]):
        var_g[t] = float(vc_blup.loc[t, "var_g"])
    else:
        var_g[t] = 0.5 * float(blues[t].var(ddof=1))

# ---------------------------------------------------------------------------
# Full-panel GEBVs on the BLUE response (re-fit GBLUP per trait, as script 18)
# ---------------------------------------------------------------------------
def load_dosage():
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sc = [c for c in hm.columns if c not in META]
    calls = hm[sc].astype(str)
    ra = hm["alleles"].str.split("/", expand=True)
    dose = np.full((len(hm), len(sc)), np.nan)
    for i in range(len(hm)):
        ref, alt = ra.iloc[i]
        arr = calls.iloc[i].values
        dose[i, arr == ref + ref] = 0.0
        dose[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
        dose[i, arr == alt + alt] = 2.0
    dose = pd.DataFrame(dose, index=hm["rs#"].values, columns=sc)
    cr_m = dose.notna().sum(axis=1) / dose.shape[1]
    maf = np.minimum(dose.mean(axis=1, skipna=True) / 2, 1 - dose.mean(axis=1, skipna=True) / 2)
    cr_s = dose.notna().sum(axis=0) / dose.shape[0]
    keep_m = ((cr_m >= 0.90) & (maf >= 0.05)).values
    keep_s = (cr_s >= 0.90).values
    return dose, keep_m, keep_s


def make_K(M):
    p = M.mean(axis=0) / 2.0
    p = np.clip(p, 1e-6, 1 - 1e-6)
    Z = M - 2.0 * p[None, :]
    return (Z @ Z.T) / (2.0 * (p * (1.0 - p)).sum())


dosage, keep_m, keep_s = load_dosage()
samples_all = dosage.columns[keep_s].tolist()

gebv = pd.DataFrame(index=samples_all)
for t in TRAITS_ALL:
    y = blues[t].reindex(samples_all)
    used = [s for s in samples_all if pd.notna(y.get(s))]
    if len(used) < 20:
        continue
    dose = dosage.loc[keep_m, used].T.values.astype(float)
    cm = np.nanmean(dose, axis=0)
    idx = np.where(np.isnan(dose)); dose[idx] = np.take(cm, idx[1])
    K = make_K(dose)
    yv = y.loc[used].values.astype(float)
    # LOO-RMSE-tuned alpha (as script 18)
    best_a, best_rmse = 5.0, np.inf
    for a in ALPHAS:
        pr = np.zeros(len(used))
        A = K + a * np.eye(len(used))
        for i in range(len(used)):
            tr = np.arange(len(used)) != i
            c = np.linalg.solve(A[np.ix_(tr, tr)], yv[tr])
            pr[i] = K[i, tr] @ c
        rmse = np.sqrt(np.mean((pr - yv) ** 2))
        if rmse < best_rmse:
            best_rmse, best_a = rmse, a
    kr = KernelRidge(kernel="precomputed", alpha=best_a); kr.fit(K, yv)
    gebv.loc[used, t] = kr.predict(K)
gebv.index.name = "sample"
gebv.to_csv(TAB / "gebv_per_accession_BLUE.csv")

# ===========================================================================
# (b) Usefulness criterion WITH accuracy: UC = mu + i * h * sigma_G
# ===========================================================================
crosses = pd.read_csv(CROSS_CSV)
uc_rows = []
for _, row in crosses.iterrows():
    cid = f"{row['parent_A']}_x_{row['parent_B']}"
    gij = float(row["G_ij"])
    for t in TRAITS_ALL:
        mid = row.get(f"{t}_midparent")
        if pd.isna(mid):
            continue
        d = TRAIT_DIR.get(t, 0)
        if d == 0:
            continue
        sigma_ms = np.sqrt(max(var_g[t] * (1.0 - gij) / 2.0, 0.0))
        h = h_acc[t]
        uc_noacc = d * float(mid) + SELECTION_INTENSITY * sigma_ms            # old (script 64)
        uc_acc = d * float(mid) + SELECTION_INTENSITY * h * sigma_ms          # Schnell & Utz 1975
        uc_rows.append({"cross_id": cid, "parent_A": row["parent_A"],
                        "parent_B": row["parent_B"], "G_ij": gij, "trait": t,
                        "mu_mid": float(mid), "sigma_MS": sigma_ms, "accuracy_h": h,
                        "UC_no_accuracy_old": uc_noacc, "UC_with_accuracy_new": uc_acc,
                        "direction": d})
uc_df = pd.DataFrame(uc_rows)
uc_df.to_csv(TAB / "usefulness_with_accuracy_per_cross_per_trait.csv", index=False)
# per-cross balanced aggregate (standardised within scenario), old vs new
agg = (uc_df.groupby("cross_id")
       .agg(G_ij=("G_ij", "first"),
            UC_old=("UC_no_accuracy_old", "mean"),
            UC_new=("UC_with_accuracy_new", "mean")).reset_index())
for col in ["UC_old", "UC_new"]:
    agg[col + "_z"] = (agg[col] - agg[col].mean()) / agg[col].std(ddof=1)
agg["rank_old"] = agg["UC_old_z"].rank(ascending=False).astype(int)
agg["rank_new"] = agg["UC_new_z"].rank(ascending=False).astype(int)
agg["rank_shift"] = (agg["rank_old"] - agg["rank_new"]).abs()
agg = agg.sort_values("rank_new")
agg.to_csv(TAB / "usefulness_with_accuracy_aggregate.csv", index=False)
spear = agg[["UC_old_z", "UC_new_z"]].corr(method="spearman").iloc[0, 1]
print(f"[b] UC old-vs-new (balanced) Spearman rank corr = {spear:.3f}; "
      f"max rank shift = {int(agg['rank_shift'].max())} of {len(agg)} crosses")

# ===========================================================================
# (c) Smith-Hazel with PHENOTYPIC P + economic weights per genetic SD
# ===========================================================================
# P = phenotypic variance-covariance among genotype means (BLUEs), pairwise-
#     complete; the honest phenotypic (co)variance, not a BLUP-corr proxy.
# G = genetic variance-covariance: diag = var_g; off-diag = genetic corr
#     (approx by GEBV corr) * sqrt(var_g_i var_g_j).
# economic weights per GENETIC SD: a_t = direction / sigma_g_t so one index
# weight unit = one genetic SD, removing the raw-scale dominance of oxalate.
#
# NOTE: including all three oxalate fractions (Total = Soluble + Insoluble
# EXACTLY) makes P and G singular → b = P^-1 G a explodes (|b| ~ 1e2-1e3).
# We therefore run TWO index definitions: the FULL 10-trait set (reported
# only as a cautionary ill-conditioned diagnostic) and the NON-REDUNDANT
# set that drops the linearly-dependent Total_Oxalate (the stable primary).
def smith_hazel(index_traits, label):
    P_emp = blues[index_traits].cov(min_periods=10)
    g_corr = gebv[index_traits].corr(method="pearson")
    gvar = np.array([var_g[t] for t in index_traits])
    sigma_g = np.sqrt(gvar)
    G_mat = g_corr.values * np.sqrt(np.outer(gvar, gvar))
    P_mat = np.nan_to_num(P_emp.values.copy(), nan=0.0)
    ridge = 1e-6 * np.trace(P_mat) / len(index_traits)
    P_mat = P_mat + ridge * np.eye(len(index_traits))
    G_mat = G_mat + ridge * np.eye(len(index_traits))
    scen = {
        "balanced_per_geneticSD": {t: 1.0 for t in index_traits},
        "nutritional": {"Crude_Protein": 3, "Antioxidant": 1, "Flavonoid": 1,
                        "Phenol": 1, "Mass_of_Seeds": 1, "Soluble_Oxalate": 2, "Tannin": 1},
        "anti_nutritional": {"Soluble_Oxalate": 3, "Total_Oxalate": 1,
                             "Seed_Coat_Tannin": 1, "Tannin": 1, "Crude_Protein": 1,
                             "Insoluble_Oxalate": 1},
    }
    rows, idx_vals = [], {}
    for sc, w in scen.items():
        a = np.array([TRAIT_DIR[t] * w.get(t, 0) / sigma_g[i]
                      for i, t in enumerate(index_traits)])
        b = np.linalg.solve(P_mat, G_mat @ a)
        for t, bv in zip(index_traits, b):
            rows.append({"index_set": label, "scenario": sc, "trait": t,
                         "b": float(bv), "a_per_geneticSD": float(a[index_traits.index(t)])})
        y = gebv[index_traits].copy().fillna(gebv[index_traits].mean())
        idx_vals[sc] = pd.Series(y.values @ b, index=gebv.index)
    return (pd.DataFrame(rows), pd.DataFrame(idx_vals),
            float(np.linalg.cond(P_mat)), float(np.linalg.cond(G_mat)))

INDEX_NONREDUNDANT = [t for t in INDEX_TRAITS if t != "Total_Oxalate"]
sh_w_full, sh_i_full, condP_f, condG_f = smith_hazel(INDEX_TRAITS, "full_10trait")
sh_w_nr, sh_i_nr, condP_nr, condG_nr = smith_hazel(INDEX_NONREDUNDANT, "nonredundant_drop_Total")
sh_weights = pd.concat([sh_w_nr, sh_w_full], ignore_index=True)
sh_weights.to_csv(TAB / "smith_hazel_weights_phenotypicP.csv", index=False)
sh_idx = sh_i_nr.copy()
for sc in sh_idx.columns:
    sh_idx[f"rank_{sc}"] = sh_idx[sc].rank(ascending=False).astype(int)
sh_idx.index.name = "sample"
sh_idx.to_csv(TAB / "smith_hazel_index_per_accession.csv")
diag = pd.DataFrame([
    {"index_set": "full_10trait", "P_condition_number": condP_f, "G_condition_number": condG_f,
     "note": "ILL-CONDITIONED: Total=Soluble+Insoluble collinear; b weights explode"},
    {"index_set": "nonredundant_drop_Total", "P_condition_number": condP_nr,
     "G_condition_number": condG_nr, "note": "stable primary index (collinearity removed)"},
])
diag.to_csv(TAB / "smith_hazel_PG_diagnostics.csv", index=False)
print(f"[c] Smith-Hazel full-index P(cond)={condP_f:.3g} (ill-conditioned); "
      f"non-redundant P(cond)={condP_nr:.1f} (stable primary)")

# ===========================================================================
# (d) Merit index + OCS EXCLUDING negative-CV traits
# ===========================================================================
def composite_merit(gebvs, traits_use):
    z = pd.DataFrame(index=gebvs.index)
    for t in traits_use:
        s = gebvs[t]
        if s.notna().sum() < 5 or s.std() == 0:
            continue
        zt = (s - s.mean()) / s.std()
        if TRAIT_DIR[t] == -1:
            zt = -zt
        z[t] = zt
    return z.mean(axis=1, skipna=True), z.notna().sum(axis=1)

all_signed = [t for t in INDEX_TRAITS]
pos_cv = [t for t in INDEX_TRAITS if t not in neg_cv_traits]
merit_full, n_full = composite_merit(gebv, all_signed)
merit_poscv, n_pos = composite_merit(gebv, pos_cv)
merit = pd.DataFrame({"merit_full_allsigned": merit_full,
                      "merit_full_n": n_full,
                      "merit_excl_negCV": merit_poscv,
                      "merit_excl_negCV_n": n_pos})
merit["rank_full"] = merit["merit_full_allsigned"].rank(ascending=False)
merit["rank_excl_negCV"] = merit["merit_excl_negCV"].rank(ascending=False)
merit = merit.sort_values("merit_excl_negCV", ascending=False)
merit.index.name = "sample"
merit.to_csv(TAB / "merit_index_full_vs_exclNegCV.csv")
sp_merit = merit[["merit_full_allsigned", "merit_excl_negCV"]].corr(method="spearman").iloc[0, 1]
# top-10 overlap
top_full = set(merit.sort_values("merit_full_allsigned", ascending=False).head(10).index)
top_excl = set(merit.sort_values("merit_excl_negCV", ascending=False).head(10).index)
print(f"[d] merit full-vs-exclNegCV Spearman={sp_merit:.3f}; "
      f"top-10 overlap = {len(top_full & top_excl)}/10; excluded traits={neg_cv_traits}")

# OCS excluding negative-CV traits (simple Lagrangian QP at one moderate lambda)
from scipy.optimize import minimize
common = sorted(set(gebv.index) & set(grm.index))
G_ocs = grm.loc[common, common].values + 1e-4 * np.eye(len(common))
def build_ebv(traits_use):
    comp = np.zeros(len(common)); wsum = np.zeros(len(common))
    for t in traits_use:
        s = gebv.loc[common, t]
        if s.notna().sum() < 5 or s.std() == 0:
            continue
        z = ((s - s.mean()) / s.std()).fillna(0).values
        m = s.notna().values.astype(float)
        comp += TRAIT_DIR[t] * z * m; wsum += m
    out = np.zeros(len(common)); ok = wsum > 0
    out[ok] = comp[ok] / wsum[ok]
    return out
def solve_ocs(ebv, lam=1.0):
    n = len(ebv)
    c0 = np.full(n, 1 / n)
    res = minimize(lambda c: -(ebv @ c - 0.5 * lam * c @ G_ocs @ c), c0,
                   jac=lambda c: -(ebv - lam * G_ocs @ c), method="SLSQP",
                   bounds=[(0, 1)] * n,
                   constraints=[{"type": "eq", "fun": lambda c: c.sum() - 1,
                                 "jac": lambda c: np.ones(n)}],
                   options={"maxiter": 500, "ftol": 1e-9})
    c = np.maximum(res.x, 0); c = c / c.sum()
    return c
ebv_full = build_ebv(all_signed); ebv_excl = build_ebv(pos_cv)
c_full = solve_ocs(ebv_full); c_excl = solve_ocs(ebv_excl)
ocs = pd.DataFrame({"sample": common,
                    "OCS_contribution_full": c_full,
                    "OCS_contribution_exclNegCV": c_excl}).sort_values(
    "OCS_contribution_exclNegCV", ascending=False)
ocs.to_csv(TAB / "ocs_full_vs_exclNegCV.csv", index=False)
print(f"[d] OCS (lambda=1): full n_active={int((c_full>1e-3).sum())}, "
      f"exclNegCV n_active={int((c_excl>1e-3).sum())}")

# ===========================================================================
# Export trait-subset data + K for R3a_mtm_sommer.R
# ===========================================================================
for tag, traits in [("oxalate_triplet", ["Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"]),
                    ("oxalate_pair", ["Soluble_Oxalate", "Insoluble_Oxalate"]),
                    ("protein_oxalate", ["Crude_Protein", "Soluble_Oxalate", "Insoluble_Oxalate"])]:
    sub = blues[traits].dropna(how="any")
    used = [s for s in sub.index if s in samples_all]
    sub = sub.loc[used]
    dose = dosage.loc[keep_m, used].T.values.astype(float)
    cm = np.nanmean(dose, axis=0); idx = np.where(np.isnan(dose)); dose[idx] = np.take(cm, idx[1])
    K = make_K(dose)
    sub.reset_index().rename(columns={"index": "sample"}).to_csv(
        TAB / f"mtm_data_{tag}.csv", index=False)
    pd.DataFrame(K, index=used, columns=used).to_csv(TAB / f"mtm_K_{tag}.csv")
    print(f"[export] mtm_data_{tag}.csv  n={len(used)} traits={traits}")

# record the negative-CV flag table
pd.DataFrame([{"trait": t, "loo_r_BLUE": acc.get(t, np.nan),
               "accuracy_h_clamped": h_acc[t],
               "negative_CV": t in neg_cv_traits} for t in TRAITS_ALL]
             ).to_csv(TAB / "per_trait_accuracy_and_negCV_flag.csv", index=False)
print("\n[done] Phase-3 b/c/d complete; run R3a_mtm_sommer.R and R3e_clone_blocked_cv.py")
