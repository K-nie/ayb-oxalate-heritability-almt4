#!/usr/bin/env python3
"""
ALMT4 SNP pleiotropy scan across the 13 panel traits.

The Soluble_Oxalate top SNP at Ss10:15,394,673 (rs
100033542|F|0-31:T>C-31:T>C) lies 3,174 bp upstream of AYB ALMT4_2
(AYBTSS11_029726). The candidate-gene call rests on (i) p = 0.009 at
n = 41 for Soluble_Oxalate, (ii) the 142-sequence ALMT phylogeny, (iii)
the protein-domain confirmation from script 65. What the candidate-gene
call does NOT directly establish is whether ALMT4 acts on Soluble_Oxalate
specifically or whether the effect spills onto correlated traits
(Total_Oxalate, Insoluble_Oxalate, Crude_Protein, the seed-size cluster).
A pleiotropy scan settles this directly.

Method
------
For each of the 13 BLUP-derived trait series, fit an EMMAX-style mixed
model

    y = mu + beta * SNP + u + e,   u ~ N(0, sigma_g^2 K),
                                    e ~ N(0, sigma_e^2 I)

with K = VanRaden method-1 kinship from script 13. REML variance
components estimated by grid search over h^2 = sigma_g^2 /
(sigma_g^2 + sigma_e^2) on a 50-point grid in [0.001, 0.999]. Wald
p-value from beta / SE(beta) under chi^2_1.

Multiple-testing correction: Holm-Bonferroni across the 13 traits, with
the raw p-values reported alongside the adjusted ones so the reader can
see both layers.

Outputs (results/60_almt4_pleiotropy_scan/)
-------------------------------------------
tables/
    almt4_per_trait_effects.csv      beta / SE / Wald p / Holm-p / n per trait
figures/
    fig_almt4_pleiotropy_forest.png  forest plot of per-trait beta with 95% CI

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
BLUP_CSV = ROOT / "results" / "00_blup_pipeline" / "tables" / "phenotype_blups.csv"
GRM_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"

OUT = ROOT / "results" / "60_almt4_pleiotropy_scan"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# Focal SNP from the manuscript: Ss10:15,394,673; rs identifier confirmed
# against results/21_candidate_genes_ayb/tables/snp_top30_anchored_Soluble_Oxalate.csv
FOCAL_RS = "100033542|F|0-31:T>C-31:T>C"
FOCAL_CHROM = "Ss10"
FOCAL_POS = 15_394_673

TRAITS = [
    "Tannin", "Phenol", "Flavonoid", "Antioxidant",
    "Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
    "Seed_Coat_Tannin", "Crude_Protein",
    "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate",
]
BREEDER_PRIORITY_TRAITS = {"Soluble_Oxalate", "Insoluble_Oxalate",
                            "Total_Oxalate", "Crude_Protein"}

SAMP_CR = 0.90
H2_GRID = np.linspace(0.001, 0.999, 50)


# ---------------------------------------------------------------------------
# Load focal SNP dosage from HapMap
# ---------------------------------------------------------------------------
def load_focal_snp_dosage() -> pd.Series:
    """Return per-sample 0/1/2 dosage at the focal ALMT4 SNP."""
    print(f"[load] reading HapMap for focal SNP {FOCAL_RS}")
    hm = pd.read_csv(HAPMAP, low_memory=False)
    meta_cols = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                 "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
    sample_cols = [c for c in hm.columns if c not in meta_cols]
    focal_row = hm.loc[hm["rs#"] == FOCAL_RS]
    if len(focal_row) != 1:
        raise RuntimeError(
            f"Expected exactly 1 row for rs={FOCAL_RS}, got {len(focal_row)}")
    focal = focal_row.iloc[0]
    alleles = focal["alleles"]
    ref, alt = alleles.split("/")
    het1, het2 = ref + alt, alt + ref
    hom_ref = ref + ref
    hom_alt = alt + alt
    calls = focal[sample_cols].astype(str).values
    dosage = np.full(len(sample_cols), np.nan)
    dosage[calls == hom_ref] = 0.0
    dosage[(calls == het1) | (calls == het2)] = 1.0
    dosage[calls == hom_alt] = 2.0
    series = pd.Series(dosage, index=sample_cols, name="ALMT4_SNP_dose")
    n_called = int(series.notna().sum())
    p_alt = float(series.mean(skipna=True)) / 2.0
    print(f"  alleles = {alleles}, n_called = {n_called}/{len(series)}, alt-allele freq = {p_alt:.3f}")
    return series


# ---------------------------------------------------------------------------
# EMMAX-style fit at one SNP for one trait
# ---------------------------------------------------------------------------
def fit_emmax_one_snp(y: np.ndarray, X: np.ndarray, K: np.ndarray,
                      h2_grid: np.ndarray) -> dict:
    """Fit y = X b + u + e with u ~ N(0, sigma_g^2 K), e ~ N(0, sigma_e^2 I)
    under REML grid search on h^2 = sigma_g^2 / (sigma_g^2 + sigma_e^2).

    Returns dict with effect-size estimate, standard error, Wald p-value
    on the SNP column (assumed to be the last column of X), h^2 MLE."""
    n, p = X.shape
    eigvals, U = np.linalg.eigh(K)
    # Eigenvalues can dip slightly negative under numerical noise; clip
    # at zero so the weight terms below stay positive.
    eigvals = np.maximum(eigvals, 0.0)

    y_rot = U.T @ y
    X_rot = U.T @ X

    best_ll = -np.inf
    best_h2 = h2_grid[0]
    best_beta = None
    best_se = None
    best_sigma2 = None

    for h2 in h2_grid:
        w = h2 * eigvals + (1.0 - h2)
        # All weights stay positive because eigvals >= 0 and h2 in (0, 1).
        log_det_w = np.sum(np.log(w))
        Wy = y_rot / w
        WX = X_rot / w[:, None]
        XtWX = X_rot.T @ WX
        XtWy = X_rot.T @ Wy
        try:
            beta_hat = np.linalg.solve(XtWX, XtWy)
        except np.linalg.LinAlgError:
            continue
        resid = y_rot - X_rot @ beta_hat
        rss = float(np.sum(resid ** 2 / w))
        # REML log-likelihood up to constant.
        # ll = -0.5 (n-p) log(sigma2_hat) - 0.5 log|XtWX| - 0.5 log|V|
        # with V = sigma^2 (h2 K + (1-h2) I); sigma^2 is concentrated out.
        sigma2_hat = rss / max(n - p, 1)
        sign, log_det_XtWX = np.linalg.slogdet(XtWX)
        if sign <= 0:
            continue
        ll = (-0.5 * (n - p) * np.log(sigma2_hat)
              - 0.5 * log_det_w
              - 0.5 * log_det_XtWX)
        if ll > best_ll:
            best_ll = ll
            best_h2 = h2
            best_beta = beta_hat.copy()
            best_se = np.sqrt(np.diag(sigma2_hat * np.linalg.inv(XtWX)))
            best_sigma2 = sigma2_hat

    beta_snp = float(best_beta[-1])
    se_snp = float(best_se[-1])
    z = beta_snp / se_snp if se_snp > 0 else 0.0
    wald_p = float(stats.chi2.sf(z ** 2, df=1))
    return {
        "beta": beta_snp,
        "se": se_snp,
        "z": z,
        "wald_p": wald_p,
        "h2_mle": float(best_h2),
        "sigma2_hat": float(best_sigma2),
        "ll_max": float(best_ll),
    }


# ---------------------------------------------------------------------------
# Holm-Bonferroni
# ---------------------------------------------------------------------------
def holm_bonferroni(pvals: np.ndarray) -> np.ndarray:
    """Family-wise-error-rate adjusted p-values under Holm-Bonferroni
    across the supplied p-value vector."""
    m = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    adj_sorted = np.zeros(m)
    running = 0.0
    for k, p in enumerate(sorted_p):
        candidate = min(1.0, (m - k) * p)
        running = max(running, candidate)
        adj_sorted[k] = running
    adj = np.zeros(m)
    adj[order] = adj_sorted
    return adj


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def plot_forest(results: pd.DataFrame, out_path: Path) -> None:
    """Per-trait beta forest with 95 % CI. Each trait is on its own row;
    beta is on the BLUP scale (units differ across traits) so we also
    render a standardised-beta column (beta / trait_SD) so the user can
    compare effect sizes across traits on a common scale."""
    results = results.copy()
    results = results.sort_values("wald_p", ascending=False).reset_index(drop=True)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 5.5),
                                       sharey=True,
                                       gridspec_kw={"width_ratios": [1.0, 1.0]})

    y = np.arange(len(results))
    colours = [WONG["vermillion"] if t in BREEDER_PRIORITY_TRAITS
                else WONG["blue"]
                for t in results["trait"]]

    ci_lo = results["beta"] - 1.96 * results["se"]
    ci_hi = results["beta"] + 1.96 * results["se"]

    # Panel A: beta on the BLUP / native scale.
    for i, (lo, hi, beta, c) in enumerate(zip(ci_lo, ci_hi,
                                                results["beta"], colours)):
        ax_a.plot([lo, hi], [i, i], color=c, linewidth=1.6)
        ax_a.scatter([beta], [i], color=c, s=55, edgecolor="black",
                      linewidth=0.6, zorder=3)
    ax_a.axvline(0, color="grey", linewidth=0.6)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(results["trait"])
    ax_a.set_xlabel("Estimated beta on BLUP scale (per ALMT4 alt-allele dose)")
    ax_a.set_title("Per-trait ALMT4 SNP effect (BLUP scale)")
    ax_a.invert_yaxis()

    # Panel B: standardised beta = beta / trait_SD.
    std_beta = results["beta"] / results["trait_sd"]
    std_se = results["se"] / results["trait_sd"]
    std_lo = std_beta - 1.96 * std_se
    std_hi = std_beta + 1.96 * std_se
    for i, (lo, hi, sb, c) in enumerate(zip(std_lo, std_hi, std_beta, colours)):
        ax_b.plot([lo, hi], [i, i], color=c, linewidth=1.6)
        ax_b.scatter([sb], [i], color=c, s=55, edgecolor="black",
                      linewidth=0.6, zorder=3)
    ax_b.axvline(0, color="grey", linewidth=0.6)
    ax_b.set_xlabel("Standardised beta (in trait-SD units per alt-allele)")
    ax_b.set_title("Per-trait standardised ALMT4 SNP effect")

    # Panel-spanning annotation for the Soluble_Oxalate row.
    fig.suptitle(
        f"ALMT4 SNP (Ss10:{FOCAL_POS:,}; {FOCAL_RS}) pleiotropy across 13 traits\n"
        "Forest plot ordered by Wald p (most significant at bottom); "
        "breeder-priority traits in vermillion",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"[load] focal SNP {FOCAL_RS} at {FOCAL_CHROM}:{FOCAL_POS:,}")
    snp_dose = load_focal_snp_dosage()

    print(f"[load] BLUP table from {BLUP_CSV.name}")
    blups = pd.read_csv(BLUP_CSV).set_index("sample")
    print(f"  BLUPs for {blups.shape[0]} accessions x {blups.shape[1]} traits")

    print(f"[load] GRM from {GRM_CSV.name}")
    grm = pd.read_csv(GRM_CSV, index_col=0)
    print(f"  GRM {grm.shape[0]} x {grm.shape[1]}")

    results = []
    for trait in TRAITS:
        # Per-trait sample intersection: those with a BLUP, a focal-SNP dose,
        # and a GRM row.
        per_trait = blups[trait].dropna()
        common = per_trait.index.intersection(snp_dose.dropna().index)
        common = common.intersection(grm.index)
        common = sorted(common)
        n = len(common)
        if n < 20:
            print(f"  [skip] {trait}: n = {n} < 20; insufficient for the fit")
            continue

        y = per_trait.loc[common].values.astype(float)
        x_snp = snp_dose.loc[common].values.astype(float)
        K_sub = grm.loc[common, common].values

        X = np.column_stack([np.ones(n), x_snp])

        # Centre the SNP column so the intercept absorbs the mean dose;
        # this makes beta directly interpretable as "effect per unit dose
        # above the panel mean", and stabilises the EMMAX fit at low MAF.
        x_centered = x_snp - np.mean(x_snp)
        X = np.column_stack([np.ones(n), x_centered])

        fit = fit_emmax_one_snp(y, X, K_sub, H2_GRID)
        results.append({
            "trait": trait,
            "n": n,
            "trait_mean": float(np.mean(y)),
            "trait_sd": float(np.std(y, ddof=1)),
            "beta": fit["beta"],
            "se": fit["se"],
            "z": fit["z"],
            "wald_p": fit["wald_p"],
            "h2_mle": fit["h2_mle"],
        })

    results_df = pd.DataFrame(results)
    results_df["wald_p_holm"] = holm_bonferroni(results_df["wald_p"].values)
    results_df = results_df.sort_values("wald_p").reset_index(drop=True)

    out_csv = TAB / "almt4_per_trait_effects.csv"
    results_df.to_csv(out_csv, index=False)
    print(f"\n[result] wrote {out_csv}")
    print(results_df[["trait", "n", "beta", "se", "wald_p", "wald_p_holm", "h2_mle"]]
          .to_string(index=False))

    print("\n[plot] forest")
    plot_forest(results_df, FIG / "fig_almt4_pleiotropy_forest")

    print("\n=== manuscript-text summary ===")
    print(f"  focal SNP                                    : {FOCAL_RS}")
    print(f"  position                                     : {FOCAL_CHROM}:{FOCAL_POS:,}")
    sol = results_df.loc[results_df["trait"] == "Soluble_Oxalate"]
    if not sol.empty:
        r = sol.iloc[0]
        print(f"  Soluble_Oxalate                              : "
              f"beta = {r['beta']:+.3f}, p = {r['wald_p']:.4f}, Holm = {r['wald_p_holm']:.4f} (n = {int(r['n'])})")
    holm_sig = results_df.loc[results_df["wald_p_holm"] < 0.10, "trait"].tolist()
    print(f"  Holm-adjusted p < 0.10 traits                : {holm_sig if holm_sig else 'none'}")
    raw_sig = results_df.loc[results_df["wald_p"] < 0.10, "trait"].tolist()
    print(f"  raw p < 0.10 traits                          : {raw_sig if raw_sig else 'none'}")


if __name__ == "__main__":
    main()
