#!/usr/bin/env python3
"""
REVISION Phase 2 — honest tests at the focal candidate SNP
Ss10:15,394,673 (rs 100033542|F|0-31:T>C-31:T>C; ~3.2 kb upstream of the
ALMT4_2 homolog AYBTSS11_029726).

Reviewer asks (L436-477, L152-161):
 (a) genotype-class counts, a per-class dot/box plot of Soluble_Oxalate,
     and the variance explained (R²) by the SNP.
 (b) an EXACT likelihood test to replace the EMMAX approximation — run in
     R/sommer (R2_focal_snp_sommer.R); here we export the tidy data and
     the kinship subset it consumes, and record the EMMAX baselines.
 (c) because Insoluble = Total − Soluble EXACTLY (max dev 5.9e-12), the
     "partitions without changing total" result is an arithmetic identity
     (β_total ≡ β_soluble + β_insoluble). We replace it with two honest
     tests: SNP vs (i) the Soluble/Total ratio, (ii) Soluble with Total as
     a covariate. OLS cross-checks here; kinship-controlled mixed-model
     versions in the sommer script.

Outputs: results/R2_focal_snp_honest/
Author: Benjamin Narh-Madey
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf

import _pheno
from _plotstyle import apply, WONG
apply()
warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
GRM_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"
OUT = ROOT / "results" / "R2_focal_snp_honest"
TAB = OUT / "tables"; FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

FOCAL_RS = "100033542|F|0-31:T>C-31:T>C"
FOCAL_POS = 15_394_673
# EMMAX baselines to compare the exact test against (from script 60 / §3.7)
EMMAX_P_KONLY = 0.0013865108686350964     # results/60 almt4_per_trait_effects.csv
EMMAX_BETA_KONLY = 98.65151553777642
EMMAX_P_KPC = 0.009                        # §3.7 GWAS, K + PC1-3


def focal_dose():
    hm = pd.read_csv(HAPMAP, low_memory=False)
    META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sc = [c for c in hm.columns if c not in META]
    row = hm.loc[hm["rs#"] == FOCAL_RS].iloc[0]
    ref, alt = row["alleles"].split("/")
    calls = row[sc].astype(str).values
    d = np.full(len(sc), np.nan)
    d[calls == ref + ref] = 0.0
    d[(calls == ref + alt) | (calls == alt + ref)] = 1.0
    d[calls == alt + alt] = 2.0
    return pd.Series(d, index=sc, name="dose"), row["alleles"]


# ---------------------------------------------------------------------------
dose, alleles = focal_dose()
means = _pheno.load_means()
grm = pd.read_csv(GRM_CSV, index_col=0)

sol = means["Soluble_Oxalate"]; tot = means["Total_Oxalate"]; ins = means["Insoluble_Oxalate"]

# Honest-test analysis set: oxalate lines with a called genotype AND a GRM
# row (matches the script-60 n=41 intersection used for the mixed model).
common = sorted(set(sol.dropna().index) & set(dose.dropna().index) & set(grm.index))
dat = pd.DataFrame({
    "sample": common,
    "dose": dose.loc[common].values,
    "Soluble_Oxalate": sol.loc[common].values,
    "Total_Oxalate": tot.loc[common].values,
    "Insoluble_Oxalate": ins.loc[common].values,
})
dat["ratio_sol_tot"] = dat["Soluble_Oxalate"] / dat["Total_Oxalate"]
dat.to_csv(TAB / "focal_oxalate_data.csv", index=False)
grm.loc[common, common].to_csv(TAB / "focal_grm_subset.csv")
n = len(dat)
print(f"[focal] rs={FOCAL_RS}  alleles={alleles}  n(analysis set)={n}")

# ---- (a) genotype-class counts + R² ---------------------------------------
classes = dat["dose"].value_counts().reindex([0.0, 1.0, 2.0]).fillna(0).astype(int)
maf = min(dat["dose"].mean() / 2, 1 - dat["dose"].mean() / 2)
cls_rows = []
for c in [0.0, 1.0, 2.0]:
    sub = dat.loc[dat["dose"] == c, "Soluble_Oxalate"]
    cls_rows.append({"genotype_class": int(c),
                     "label": {0: "TT (hom ref)", 1: "TC (het)", 2: "CC (hom alt)"}[int(c)],
                     "n": int(classes[c]),
                     "SolOx_mean": float(sub.mean()) if len(sub) else np.nan,
                     "SolOx_sd": float(sub.std(ddof=1)) if len(sub) > 1 else np.nan})
cls_df = pd.DataFrame(cls_rows)

# R² = variance explained by additive dose (simple linear regression)
lm = smf.ols("Soluble_Oxalate ~ dose", data=dat).fit()
r2 = float(lm.rsquared)
r2_adj = float(lm.rsquared_adj)
beta_ols = float(lm.params["dose"]); se_ols = float(lm.bse["dose"]); p_ols = float(lm.pvalues["dose"])

print("\n[a] genotype-class counts (analysis set):")
print(cls_df.to_string(index=False))
print(f"    MAF(set)={maf:.3f}   additive R²={r2:.3f} (adj {r2_adj:.3f})   "
      f"OLS β={beta_ols:+.2f} SE={se_ols:.2f} p={p_ols:.4f}")

# ---- (c) honest OLS cross-checks: ratio + Soluble|Total -------------------
lm_ratio = smf.ols("ratio_sol_tot ~ dose", data=dat).fit()
lm_cov = smf.ols("Soluble_Oxalate ~ dose + Total_Oxalate", data=dat).fit()
honest_rows = [
    {"test": "Soluble ~ dose (marginal, OLS)", "beta": beta_ols, "se": se_ols,
     "p": p_ols, "r2": r2, "note": "same quantity as EMMAX but no kinship control"},
    {"test": "Soluble/Total ratio ~ dose (OLS)",
     "beta": float(lm_ratio.params["dose"]), "se": float(lm_ratio.bse["dose"]),
     "p": float(lm_ratio.pvalues["dose"]), "r2": float(lm_ratio.rsquared),
     "note": "honest fraction test (dimensionless)"},
    {"test": "Soluble ~ dose + Total (OLS, Total as covariate)",
     "beta": float(lm_cov.params["dose"]), "se": float(lm_cov.bse["dose"]),
     "p": float(lm_cov.pvalues["dose"]), "r2": float(lm_cov.rsquared),
     "note": "honest partitioning test; Total held fixed"},
]
honest_df = pd.DataFrame(honest_rows)
cls_df.to_csv(TAB / "genotype_class_counts.csv", index=False)
honest_df.to_csv(TAB / "honest_tests_ols.csv", index=False)
print("\n[c] honest OLS tests (ratio / covariate-adjusted):")
print(honest_df.to_string(index=False))

# identity demonstration row
ident = pd.DataFrame([{
    "max_abs_dev_Ins_minus_(Tot-Sol)": float((ins - (tot - sol)).abs().max()),
    "EMMAX_beta_total": -43.07746708504585,
    "EMMAX_beta_soluble": 98.65151553777642,
    "EMMAX_beta_insoluble": -142.54903278185876,
    "beta_sol_plus_beta_ins": 98.65151553777642 + (-142.54903278185876),
    "note": "beta_total ~= beta_soluble + beta_insoluble because Total=Sol+Ins; "
            "the insoluble 'opposite effect' is not independent evidence."}])
ident.to_csv(TAB / "arithmetic_identity_demo.csv", index=False)

# ---- (a) figure: per-class stripplot + box of Soluble_Oxalate -------------
fig, ax = plt.subplots(figsize=(5.2, 4.3), constrained_layout=True)
pal = {0: WONG["blue"], 1: WONG["orange"], 2: WONG["vermillion"]}
for c in [0.0, 1.0, 2.0]:
    sub = dat.loc[dat["dose"] == c, "Soluble_Oxalate"].values
    x = np.full(len(sub), c) + np.random.default_rng(int(c) + 1).normal(0, 0.05, len(sub))
    ax.scatter(x, sub, s=36, color=pal[int(c)], edgecolor="black",
               linewidth=0.4, zorder=3, alpha=0.9)
    if len(sub) >= 1:
        ax.hlines(sub.mean(), c - 0.2, c + 0.2, color="black", lw=1.6, zorder=4)
bx = [dat.loc[dat["dose"] == c, "Soluble_Oxalate"].values for c in [0.0, 1.0, 2.0]]
ax.boxplot(bx, positions=[0, 1, 2], widths=0.5, showfliers=False,
           boxprops=dict(color="#888888"), medianprops=dict(color="#888888"),
           whiskerprops=dict(color="#888888"), capprops=dict(color="#888888"), zorder=1)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([f"TT\nn={int(classes[0.0])}", f"TC\nn={int(classes[1.0])}",
                    f"CC\nn={int(classes[2.0])}"])
ax.set_xlabel(f"Genotype at Ss10:{FOCAL_POS:,} (alt-allele dose)")
ax.set_ylabel("Soluble oxalate (mg/100 g)")
ax.text(0.02, 0.98, f"additive $R^2$ = {r2:.2f}\nMAF = {maf:.3f},  n = {n}",
        transform=ax.transAxes, va="top", ha="left", fontsize=9)
fig.savefig(FIG / "focal_snp_genotype_class_boxplot.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "focal_snp_genotype_class_boxplot.pdf", bbox_inches="tight")
plt.close(fig)
print(f"\n[fig] focal_snp_genotype_class_boxplot  -> {FIG}")
print(f"[done] Phase-2 python; run R2_focal_snp_sommer.R for the exact REML test")
