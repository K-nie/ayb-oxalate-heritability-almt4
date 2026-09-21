#!/usr/bin/env python3
"""
Allier-Lehermeier-Bonk usefulness criterion on the top-50 cross pairs.

Ranks parent pairs not by mid-parent value alone but by the expected
merit of the top-q% of their progeny, mu_c + i * sigma_c, where mu_c is
the cross mid-parent value, sigma_c is the within-family Mendelian
sampling SD, and i is the standardised selection intensity.

At our marker density, the per-SNP-effect-based Bonk 2016 / Lehermeier
2017 formulation is approximated with a kinship-derived Mendelian
sampling variance:

    sigma^2_MS_trait_cross = var_g_trait * (1 - G_ij_cross) / 2

This is the additive infinitesimal-limit approximation (Hill & Weir 2011;
Falconer & Mackay 1996); valid when per-SNP effects are absorbed into
the trait-level variance component. The exact Bonk 2016 formula
requires phased parental haplotypes and per-SNP effects that the kernel
GBLUP pipeline does not directly produce. Once BayesB / BayesC SNP
effects are in hand (Analysis 4), the script extends naturally by
swapping in the marker-resolution Mendelian sampling variance.

Selection intensity: i = 1.755 corresponds to selecting the top 10%
of progeny under a normal distribution (Falconer & Mackay 1996 Table A).

Cross-merit aggregation across the 13 traits uses three breeder weights:
balanced; Soluble_Oxalate priority (weight 3 on Sol_Ox-low,
weight 2 on Insol_Ox-high, weight 1 elsewhere); Crude_Protein priority
(weight 3 on CP-high, weight 1 elsewhere). Direction signs follow the
§2.13 merit-index convention from the manuscript.

Outputs (results/64_usefulness_criterion/)
------------------------------------------
tables/
    usefulness_per_cross_per_trait.csv   50 crosses x 13 traits x (mu, sigma, U)
    usefulness_aggregate_per_weight.csv    per-cross aggregate U per scenario
    top10_under_each_weighting.csv         top-10 crosses per scenario
figures/
    fig_usefulness_vs_relatedness.png      U scatter vs G_ij with top-10 marked

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
CROSS_CSV = ROOT / "results" / "13_grm_crosspairs" / "tables" / "cross_pairs_least_related.csv"
BLUP_CSV = ROOT / "results" / "00_blup_pipeline" / "tables" / "phenotype_blups.csv"
VARCOMP_CSV = ROOT / "results" / "00_blup_pipeline" / "tables" / "variance_components.csv"
OUT = ROOT / "results" / "64_usefulness_criterion"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

TRAITS = [
    "Tannin", "Phenol", "Flavonoid", "Antioxidant",
    "Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
    "Seed_Coat_Tannin", "Crude_Protein",
    "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate",
]
TRAIT_DIRECTIONS = {
    "Tannin": -1, "Phenol": +1, "Flavonoid": +1, "Antioxidant": +1,
    "Seed_Length": 0, "Seed_Width": 0, "Seed_Thickness": 0,
    "Mass_of_Seeds": +1, "Seed_Coat_Tannin": -1,
    "Crude_Protein": +1,
    "Total_Oxalate": -1, "Soluble_Oxalate": -1, "Insoluble_Oxalate": +1,
}

SELECTION_INTENSITY = 1.755   # top 10% under normality (Falconer & Mackay)

WEIGHTING_SCENARIOS = {
    "balanced": {t: 1.0 for t in TRAITS if TRAIT_DIRECTIONS[t] != 0},
    "soluble_oxalate_priority": {
        **{t: 1.0 for t in TRAITS if TRAIT_DIRECTIONS[t] != 0},
        "Soluble_Oxalate": 3.0,
        "Insoluble_Oxalate": 2.0,
    },
    "crude_protein_priority": {
        **{t: 1.0 for t in TRAITS if TRAIT_DIRECTIONS[t] != 0},
        "Crude_Protein": 3.0,
    },
}


def load_var_g_per_trait() -> dict[str, float]:
    """Per-trait genetic variance. For traits without replicates (the
    four no-replicate phenotypes), fall back to BLUP-variance * 0.5
    under prior h^2 = 0.5 -- same convention as script 70."""
    vc = pd.read_csv(VARCOMP_CSV).set_index("trait")
    blups = pd.read_csv(BLUP_CSV).set_index("sample")
    out = {}
    for t in TRAITS:
        if t in vc.index and pd.notna(vc.loc[t, "var_g"]):
            out[t] = float(vc.loc[t, "var_g"])
        else:
            out[t] = 0.5 * float(blups[t].var(ddof=1))
    return out


def main() -> None:
    print(f"[load] {CROSS_CSV.name}")
    crosses = pd.read_csv(CROSS_CSV)
    print(f"  {len(crosses)} cross pairs")

    var_g = load_var_g_per_trait()
    print(f"[var] var_g per trait:")
    for t in TRAITS:
        print(f"    {t:20s} = {var_g[t]:.3f}")

    # Per-trait per-cross mu, sigma, U.
    per_pair_rows = []
    aggregate_rows = []
    for idx, row in crosses.iterrows():
        cross_id = f"{row['parent_A']}_x_{row['parent_B']}"
        g_ij = float(row["G_ij"])

        # Mendelian sampling variance per trait under the infinitesimal
        # additive approximation.
        per_trait = {}
        for t in TRAITS:
            mid_col = f"{t}_midparent"
            mu = row.get(mid_col)
            if pd.isna(mu):
                # Mid-parent not computable when one parent lacks the
                # trait BLUP (common for the n=46 oxalate triplet on the
                # 95-line panel).
                continue
            mu = float(mu)
            sigma2 = max(var_g[t] * (1.0 - g_ij) / 2.0, 0.0)
            sigma = float(np.sqrt(sigma2))
            d = TRAIT_DIRECTIONS[t]
            U = d * mu + SELECTION_INTENSITY * sigma if d != 0 else np.nan
            per_pair_rows.append({
                "cross_id": cross_id,
                "parent_A": row["parent_A"],
                "parent_B": row["parent_B"],
                "G_ij": g_ij,
                "trait": t,
                "mu": mu,
                "sigma": sigma,
                "U": U,
                "direction": d,
            })
            per_trait[t] = {"mu": mu, "sigma": sigma, "U": U, "d": d}

        # Aggregate per weighting scenario: weighted standardised U on
        # the cross-pool z-score scale so traits on different units
        # contribute comparably. The standardisation is computed at the
        # end, but for now record per-cross unweighted sums.
        for scenario, weights in WEIGHTING_SCENARIOS.items():
            num = 0.0
            denom = 0.0
            for t, w in weights.items():
                if t not in per_trait:
                    continue
                u = per_trait[t]["U"]
                if u is None or pd.isna(u):
                    continue
                num += w * u
                denom += abs(w)
            U_agg = num / denom if denom > 0 else np.nan
            aggregate_rows.append({
                "cross_id": cross_id,
                "parent_A": row["parent_A"],
                "parent_B": row["parent_B"],
                "G_ij": g_ij,
                "scenario": scenario,
                "U_aggregate": U_agg,
            })

    per_pair_df = pd.DataFrame(per_pair_rows)
    aggregate_df = pd.DataFrame(aggregate_rows)

    # Standardise the aggregate U within each scenario so the ranking is
    # not dominated by a single trait scale.
    for scenario in WEIGHTING_SCENARIOS:
        mask = aggregate_df["scenario"] == scenario
        vals = aggregate_df.loc[mask, "U_aggregate"]
        z = (vals - vals.mean()) / vals.std(ddof=1)
        aggregate_df.loc[mask, "U_aggregate_z"] = z

    per_pair_df.to_csv(TAB / "usefulness_per_cross_per_trait.csv", index=False)
    aggregate_df.to_csv(TAB / "usefulness_aggregate_per_weight.csv", index=False)

    # Top-10 per scenario by U_aggregate_z.
    top_rows = []
    for scenario in WEIGHTING_SCENARIOS:
        sub = aggregate_df.loc[aggregate_df["scenario"] == scenario]
        top = sub.nlargest(10, "U_aggregate_z").copy()
        top["rank"] = np.arange(1, len(top) + 1)
        for _, r in top.iterrows():
            top_rows.append({
                "scenario": scenario,
                "rank": int(r["rank"]),
                "cross_id": r["cross_id"],
                "parent_A": r["parent_A"],
                "parent_B": r["parent_B"],
                "G_ij": r["G_ij"],
                "U_aggregate": r["U_aggregate"],
                "U_aggregate_z": r["U_aggregate_z"],
            })
    top_df = pd.DataFrame(top_rows)
    top_df.to_csv(TAB / "top10_under_each_weighting.csv", index=False)

    # Scatter: U vs G_ij per scenario.
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    colour_map = {
        "balanced": WONG["blue"],
        "soluble_oxalate_priority": WONG["vermillion"],
        "crude_protein_priority": WONG["green"],
    }
    label_map = {
        "balanced": "Balanced (all weight 1)",
        "soluble_oxalate_priority": "Soluble_Oxalate priority (w=3, Insol_Ox w=2)",
        "crude_protein_priority": "Crude_Protein priority (w=3)",
    }
    for scenario in WEIGHTING_SCENARIOS:
        sub = aggregate_df.loc[aggregate_df["scenario"] == scenario]
        ax.scatter(sub["G_ij"], sub["U_aggregate_z"],
                    color=colour_map[scenario], s=42,
                    edgecolor="white", linewidth=0.4, alpha=0.75,
                    label=label_map[scenario])
        top10 = sub.nlargest(10, "U_aggregate_z")
        ax.scatter(top10["G_ij"], top10["U_aggregate_z"],
                    facecolor="none", edgecolor=colour_map[scenario],
                    s=140, linewidth=1.5, zorder=4)
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("Parental G_ij (relatedness; lower = more diverse)")
    ax.set_ylabel("Aggregate usefulness criterion (standardised z within scenario)")
    ax.set_title("Allier 2019 usefulness criterion across the top-50 least-related cross pairs\n"
                  "Circled points = top-10 under each scenario")
    ax.legend(loc="lower left", fontsize=9, framealpha=0.95,
              edgecolor="none")
    fig.tight_layout()
    fig.savefig(FIG / "fig_usefulness_vs_relatedness.png", dpi=300,
                bbox_inches="tight")
    fig.savefig(FIG / "fig_usefulness_vs_relatedness.pdf",
                bbox_inches="tight")
    plt.close(fig)

    print(f"\n[result] wrote tables to {TAB}")
    print(f"  per_pair_df    : {per_pair_df.shape}")
    print(f"  aggregate_df   : {aggregate_df.shape}")
    print(f"  top_df         : {top_df.shape}")

    print("\n=== manuscript-text summary ===")
    for scenario in WEIGHTING_SCENARIOS:
        top10 = top_df.loc[top_df["scenario"] == scenario].sort_values("rank")
        print(f"\n  top-10 under {scenario}:")
        for _, r in top10.iterrows():
            print(f"    #{int(r['rank']):2d}  G_ij={r['G_ij']:+.3f}  U_z={r['U_aggregate_z']:+.3f}  "
                  f"{r['parent_A']} x {r['parent_B']}")


if __name__ == "__main__":
    main()
