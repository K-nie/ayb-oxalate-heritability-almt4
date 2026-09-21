#!/usr/bin/env python3
"""
Optimal Contribution Selection (Meuwissen 1997; Woolliams et al. 2015).

Picks per-accession contribution shares c in [0, 1] summing to 1 that
maximise expected multi-trait genetic gain subject to a constraint on the
per-generation inbreeding increase. Distinct from the merit-index z-score
composite (§3.12 unweighted ranking) and from the Allier 2019 usefulness
criterion (Analysis 5, cross-pair ranking by mu + i*sigma): OCS is the
parental-contribution layer that a working breeder consumes at a crossing
block.

Formulation
-----------
    maximise   c' EBV - lambda/2 * c' G c
    subject to sum(c) = 1, c_i >= 0

where G is the VanRaden method-1 GRM (script 13) and EBV is a per-
accession direction-signed merit z-score under one of three breeder
priority scenarios. The Lagrangian-relaxation form yields the Pareto
front in (achieved gain, achieved inbreeding) by sweeping lambda. We
sweep lambda on a log-spaced grid and report:

  - Pareto front trace per scenario
  - top-15 contributing accessions at the moderate operating point
    (lambda picked to give roughly the median achieved inbreeding)

Outputs (results/68_ocs/)
-------------------------
tables/
    ocs_contributions_per_priority.csv   c per accession x priority x lambda
    pareto_front.csv                       per-lambda achieved gain / inbreeding
    top_contributions_per_priority.csv    top-15 accessions per scenario at moderate lambda
figures/
    fig_ocs_pareto.png                     gain vs inbreeding per scenario
    fig_ocs_top_contributions.png         top-15 per scenario bar plot

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
GEBV_CSV = ROOT / "results" / "18_gebv" / "tables" / "gebv_per_accession.csv"
GRM_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "grm_vanraden.csv"
OUT = ROOT / "results" / "68_ocs"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# Direction-signed traits: + when breeder wants high, - when wants low.
# Direction-neutral traits (seed-size dimensions) are excluded from the
# composite because their preferred direction depends on the target
# market.
TRAIT_DIRECTIONS = {
    "Tannin": -1, "Phenol": +1, "Flavonoid": +1, "Antioxidant": +1,
    "Mass_of_Seeds": +1, "Seed_Coat_Tannin": -1,
    "Crude_Protein": +1,
    "Total_Oxalate": -1, "Soluble_Oxalate": -1, "Insoluble_Oxalate": +1,
}

# Three breeder priority scenarios; weights multiply the direction-signed
# per-trait z-score before averaging into the composite EBV.
SCENARIOS = {
    "balanced": {t: 1.0 for t in TRAIT_DIRECTIONS},
    "soluble_oxalate_priority": {
        **{t: 1.0 for t in TRAIT_DIRECTIONS},
        "Soluble_Oxalate": 3.0,
        "Insoluble_Oxalate": 2.0,
    },
    "crude_protein_priority": {
        **{t: 1.0 for t in TRAIT_DIRECTIONS},
        "Crude_Protein": 3.0,
    },
}

# Lambda grid (log-spaced). Lower lambda = more aggressive gain, higher
# inbreeding; higher lambda = more diversity-preserving, lower gain.
LAMBDA_GRID = np.logspace(-2, 2, 20)

# Ridge added to G to ensure positive-definiteness for the QP. At our
# panel size some eigenvalues fall near zero.
G_RIDGE = 1e-4


# ---------------------------------------------------------------------------
# Build the priority-weighted composite EBV per accession
# ---------------------------------------------------------------------------
def build_composite_ebv(gebvs: pd.DataFrame, weights: dict[str, float]) -> np.ndarray:
    """Direction-signed z-score composite per accession across the
    direction-signed traits present in `weights`.

    For each trait t with direction d in {+1, -1} and weight w,
    contribution to EBV_i is `w * d * z(GEBV_t)`. The composite EBV is
    the weighted-average of available signed z-scores per accession.
    Accessions missing some traits (e.g. oxalate triplet n = 41 of 95)
    fall back to averaging across only the available traits."""
    samples = gebvs.index.tolist()
    n = len(samples)
    composite = np.zeros(n)
    weight_sum = np.zeros(n)

    for trait, direction in TRAIT_DIRECTIONS.items():
        if trait not in weights:
            continue
        w = weights[trait]
        series = gebvs[trait]
        valid = series.notna()
        if not valid.any():
            continue
        z = (series - series.mean(skipna=True)) / series.std(skipna=True, ddof=1)
        contribution = direction * w * z.fillna(0.0).values
        mask = valid.values.astype(float)
        composite = composite + contribution * mask
        weight_sum = weight_sum + abs(w) * mask

    safe = weight_sum > 0
    out = np.zeros(n)
    out[safe] = composite[safe] / weight_sum[safe]
    return out


# ---------------------------------------------------------------------------
# QP solver via scipy SLSQP
# ---------------------------------------------------------------------------
def solve_ocs(ebv: np.ndarray, G: np.ndarray, lam: float,
              c_init: np.ndarray = None) -> dict:
    """Maximise c' ebv - lam/2 * c' G c subject to sum(c) = 1, c >= 0.

    SLSQP via scipy.optimize.minimize; objective passed as the negative
    of the maximand."""
    n = len(ebv)
    if c_init is None:
        c_init = np.full(n, 1.0 / n)

    def neg_obj(c):
        return -(ebv @ c - 0.5 * lam * c @ G @ c)

    def neg_obj_grad(c):
        return -(ebv - lam * G @ c)

    constraints = [{"type": "eq", "fun": lambda c: np.sum(c) - 1.0,
                     "jac": lambda c: np.ones_like(c)}]
    bounds = [(0.0, 1.0)] * n

    result = minimize(neg_obj, c_init, jac=neg_obj_grad,
                       method="SLSQP", bounds=bounds,
                       constraints=constraints,
                       options={"maxiter": 500, "ftol": 1e-9})
    c = result.x
    # Numerical clean-up.
    c = np.maximum(c, 0.0)
    if c.sum() > 0:
        c = c / c.sum()
    gain = float(ebv @ c)
    inbreeding = float(0.5 * c @ G @ c)
    return {
        "c": c,
        "gain": gain,
        "inbreeding": inbreeding,
        "converged": bool(result.success),
        "n_active": int(np.sum(c > 1e-4)),
    }


# ---------------------------------------------------------------------------
# Run all scenarios across the lambda grid
# ---------------------------------------------------------------------------
def run_all(gebvs: pd.DataFrame, G_ridged: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (pareto_df with per-lambda summary, contributions_df with
    per-accession c at every (scenario, lambda)) pair."""
    pareto_rows = []
    contrib_rows = []
    samples = gebvs.index.tolist()

    for scenario, weights in SCENARIOS.items():
        ebv = build_composite_ebv(gebvs, weights)
        c_init = np.full(len(ebv), 1.0 / len(ebv))
        for lam in LAMBDA_GRID:
            result = solve_ocs(ebv, G_ridged, lam, c_init=c_init)
            pareto_rows.append({
                "scenario": scenario,
                "lambda": float(lam),
                "gain": result["gain"],
                "inbreeding": result["inbreeding"],
                "n_active_accessions": result["n_active"],
                "converged": result["converged"],
            })
            for s, c_i in zip(samples, result["c"]):
                contrib_rows.append({
                    "scenario": scenario,
                    "lambda": float(lam),
                    "sample": s,
                    "contribution": float(c_i),
                })
            # Warm-start the next lambda.
            c_init = result["c"]
    pareto_df = pd.DataFrame(pareto_rows)
    contrib_df = pd.DataFrame(contrib_rows)
    return pareto_df, contrib_df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
SCENARIO_COLOUR = {
    "balanced": WONG["blue"],
    "soluble_oxalate_priority": WONG["vermillion"],
    "crude_protein_priority": WONG["green"],
}
SCENARIO_LABEL = {
    "balanced": "Balanced (all traits weight 1)",
    "soluble_oxalate_priority": "Soluble_Oxalate priority (S_Ox w=3, I_Ox w=2)",
    "crude_protein_priority": "Crude_Protein priority (CP w=3)",
}


def plot_pareto(pareto_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for scenario, sub in pareto_df.groupby("scenario"):
        sub = sub.sort_values("inbreeding")
        ax.plot(sub["inbreeding"], sub["gain"],
                "-o", color=SCENARIO_COLOUR[scenario],
                markersize=6, linewidth=1.5,
                label=SCENARIO_LABEL[scenario])
    ax.set_xlabel("Achieved inbreeding (c' G c / 2)")
    ax.set_ylabel("Achieved genetic gain (composite EBV)")
    ax.set_title("OCS Pareto front: genetic gain vs inbreeding\n"
                  "Lambda swept on log scale from 0.01 (gain-aggressive) to 100 (diversity-preserving)")
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def plot_top_contributions(contrib_df: pd.DataFrame,
                            moderate_lambda: dict[str, float],
                            out_path: Path) -> None:
    fig, axes = plt.subplots(1, len(SCENARIOS), figsize=(15, 5.5),
                              sharex=False, sharey=True)
    for ax, scenario in zip(axes, SCENARIOS.keys()):
        lam = moderate_lambda[scenario]
        sub = contrib_df.loc[(contrib_df["scenario"] == scenario)
                              & (np.isclose(contrib_df["lambda"], lam))]
        top = sub.nlargest(15, "contribution").sort_values("contribution",
                                                            ascending=True)
        ax.barh(top["sample"], top["contribution"],
                color=SCENARIO_COLOUR[scenario],
                edgecolor="black", linewidth=0.4)
        ax.set_title(SCENARIO_LABEL[scenario], fontsize=9)
        ax.set_xlabel(f"OCS contribution at lambda = {lam:.3g}")
    axes[0].set_ylabel("Accession")
    fig.suptitle("Top-15 OCS contributions per breeder priority scenario\n"
                  "(at moderate-lambda operating point)",
                  fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"[load] GEBV table from {GEBV_CSV.name}")
    gebvs_full = pd.read_csv(GEBV_CSV)
    gebvs = gebvs_full.set_index("sample")[list(TRAIT_DIRECTIONS.keys())]
    print(f"  {gebvs.shape[0]} accessions x {gebvs.shape[1]} traits")

    print(f"[load] GRM from {GRM_CSV.name}")
    grm = pd.read_csv(GRM_CSV, index_col=0)
    print(f"  GRM {grm.shape[0]} x {grm.shape[1]}")

    # Restrict to samples present in both.
    common = sorted(set(gebvs.index) & set(grm.index))
    gebvs = gebvs.loc[common]
    G = grm.loc[common, common].values
    print(f"  intersected to {len(common)} accessions")

    # Ridge for QP positive-definiteness.
    G_ridged = G + G_RIDGE * np.eye(len(common))

    print(f"[solve] OCS sweep across {len(LAMBDA_GRID)} lambda values x "
          f"{len(SCENARIOS)} scenarios")
    pareto_df, contrib_df = run_all(gebvs, G_ridged)

    pareto_path = TAB / "pareto_front.csv"
    pareto_df.to_csv(pareto_path, index=False)
    contrib_path = TAB / "ocs_contributions_per_priority.csv"
    contrib_df.to_csv(contrib_path, index=False)
    print(f"  wrote {pareto_path}")
    print(f"  wrote {contrib_path}")

    # Moderate-lambda operating point per scenario: pick lambda whose
    # achieved inbreeding sits at the median of the swept range.
    moderate_lambda = {}
    for scenario in SCENARIOS:
        sub = pareto_df.loc[pareto_df["scenario"] == scenario]
        median_inbreeding = sub["inbreeding"].median()
        closest = sub.iloc[(sub["inbreeding"] - median_inbreeding).abs().argmin()]
        moderate_lambda[scenario] = float(closest["lambda"])
    print(f"\n[moderate-lambda picks per scenario]")
    for scen, lam in moderate_lambda.items():
        print(f"  {scen:35s} -> lambda = {lam:.4g}")

    # Top-15 contributions at moderate lambda.
    top_rows = []
    for scenario, lam in moderate_lambda.items():
        sub = contrib_df.loc[(contrib_df["scenario"] == scenario)
                              & (np.isclose(contrib_df["lambda"], lam))]
        top = sub.nlargest(15, "contribution")
        for _, row in top.iterrows():
            top_rows.append({
                "scenario": scenario,
                "lambda_moderate": lam,
                "sample": row["sample"],
                "contribution": row["contribution"],
            })
    top_df = pd.DataFrame(top_rows)
    top_csv = TAB / "top_contributions_per_priority.csv"
    top_df.to_csv(top_csv, index=False)
    print(f"  wrote {top_csv}")

    print("\n[plot] Pareto front")
    plot_pareto(pareto_df, FIG / "fig_ocs_pareto")
    print("[plot] top contributions per scenario")
    plot_top_contributions(contrib_df, moderate_lambda,
                             FIG / "fig_ocs_top_contributions")

    print("\n=== manuscript-text summary ===")
    for scenario in SCENARIOS:
        sub = pareto_df.loc[pareto_df["scenario"] == scenario]
        print(f"  {scenario}")
        print(f"    gain range  : [{sub['gain'].min():.3f}, {sub['gain'].max():.3f}]")
        print(f"    inbr range  : [{sub['inbreeding'].min():.4f}, {sub['inbreeding'].max():.4f}]")
        print(f"    moderate-lambda n_active : "
              f"{int(sub.loc[sub['lambda'] == moderate_lambda[scenario], 'n_active_accessions'].values[0])}")


if __name__ == "__main__":
    main()
