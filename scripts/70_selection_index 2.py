#!/usr/bin/env python3
"""
Smith-Hazel / Pesek-Baker selection index optimisation.

Constructs a per-accession composite EBV under three breeder-priority
scenarios using the per-trait genetic variance, the genetic correlation
matrix, and per-trait economic / nutritional weights. The output is a
per-accession ranking under each scenario; the robust elites (top-10
under all three scenarios) and the scenario-dependent elites (top-10
under only one or two) discriminate between accessions whose superiority
is broad versus weight-specific.

Smith (1936) and Hazel (1943) construct the selection index as

    b = P^{-1} G a

where P is the phenotypic variance-covariance matrix, G is the genetic
variance-covariance matrix, and a is the vector of economic weights.
The per-accession index is

    I_i = b' y_i

with y_i the vector of trait BLUPs (or GEBVs) for accession i.

Pesek-Baker (1969) replaces a with desired-genetic-gain weights -- the
practical alternative when economic weights are not well-defined for
orphan crops.

Three scenarios encode the AYB-relevant breeder priorities:
  (a) nutritional        : Crude_Protein +3, Antioxidant +1, Flavonoid +1,
                            Phenol +1, Mass_of_Seeds +1, Soluble_Oxalate -2,
                            Tannin -1.
  (b) anti_nutritional   : Soluble_Oxalate -3, Total_Oxalate -1,
                            Seed_Coat_Tannin -1, Tannin -1, Crude_Protein +1,
                            Insoluble_Oxalate +1.
  (c) balanced           : all 10 direction-signed traits +/- 1.

Outputs (results/70_selection_index/)
-------------------------------------
tables/
    index_weights_b.csv                     b per trait per scenario
    per_accession_index_per_scenario.csv     95 x 3 index values + ranks
    robust_vs_scenario_elites.csv           top-10 per scenario with overlap
figures/
    fig_index_rank_changes.png              slopegraph of rank across 3 scenarios
    fig_index_weights_forest.png             per-trait b weight per scenario

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
GEBV_CSV = ROOT / "results" / "18_gebv" / "tables" / "gebv_per_accession.csv"
BLUP_CSV = ROOT / "results" / "00_blup_pipeline" / "tables" / "phenotype_blups.csv"
VARCOMP_CSV = ROOT / "results" / "00_blup_pipeline" / "tables" / "variance_components.csv"
OUT = ROOT / "results" / "70_selection_index"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# Traits used in the index. Seed-size dimensions excluded as direction-
# neutral (their preferred direction depends on the target market).
INDEX_TRAITS = [
    "Tannin", "Phenol", "Flavonoid", "Antioxidant",
    "Mass_of_Seeds", "Seed_Coat_Tannin",
    "Crude_Protein",
    "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate",
]
TRAIT_DIRECTIONS = {
    "Tannin": -1, "Phenol": +1, "Flavonoid": +1, "Antioxidant": +1,
    "Mass_of_Seeds": +1, "Seed_Coat_Tannin": -1,
    "Crude_Protein": +1,
    "Total_Oxalate": -1, "Soluble_Oxalate": -1, "Insoluble_Oxalate": +1,
}

# Three economic-weight scenarios. Weights are direction-signed magnitudes
# (the trait direction sign is applied separately in the EBV vector). The
# direction-signed final weight `a_t` enters Smith-Hazel as `a_t = sign * |w|`.
SCENARIOS = {
    "nutritional": {
        "Crude_Protein": 3, "Antioxidant": 1, "Flavonoid": 1,
        "Phenol": 1, "Mass_of_Seeds": 1,
        "Soluble_Oxalate": 2, "Tannin": 1,
    },
    "anti_nutritional": {
        "Soluble_Oxalate": 3, "Total_Oxalate": 1,
        "Seed_Coat_Tannin": 1, "Tannin": 1,
        "Crude_Protein": 1, "Insoluble_Oxalate": 1,
    },
    "balanced": {t: 1 for t in INDEX_TRAITS},
}


# ---------------------------------------------------------------------------
# Estimate per-trait genetic and phenotypic variances + covariance matrices
# ---------------------------------------------------------------------------
def load_variance_components() -> pd.DataFrame:
    vc = pd.read_csv(VARCOMP_CSV)
    print(f"[load] variance components for {len(vc)} traits")
    return vc


def build_PG_matrices(blups: pd.DataFrame, gebvs: pd.DataFrame,
                       vc: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Construct P (phenotypic) and G (genetic) variance-covariance
    matrices on the INDEX_TRAITS subset.

    P -- estimated from the BLUP-stage Cullis-line-mean phenotypic
         variance per trait on the diagonal, with off-diagonals from the
         observed BLUP correlation scaled by sqrt(P_ii * P_jj).
    G -- estimated from the per-trait sigma_g^2 reported in
         variance_components.csv on the diagonal, with off-diagonals from
         the GEBV correlation scaled by sqrt(G_ii * G_jj). GEBV
         correlation is a pragmatic proxy for the genetic correlation
         when a full multi-trait BLUP solver is not in hand.
    """
    vc = vc.set_index("trait")
    p_blups = blups[INDEX_TRAITS]
    p_gebv = gebvs[INDEX_TRAITS]

    # Phenotypic correlation matrix from BLUPs.
    p_corr = p_blups.corr(method="pearson").values
    g_corr = p_gebv.corr(method="pearson").values

    # Per-trait genetic and residual variance from the BLUP pipeline.
    # Traits without replicates (Seed_Coat_Tannin, Total/Soluble/Insoluble
    # Oxalate) appear in variance_components.csv as NaN rows; fall through
    # to the BLUP-variance-based fallback in that case.
    g_var = np.zeros(len(INDEX_TRAITS))
    e_var = np.zeros(len(INDEX_TRAITS))
    for i, t in enumerate(INDEX_TRAITS):
        has_vc = (t in vc.index
                   and pd.notna(vc.loc[t, "var_g"])
                   and pd.notna(vc.loc[t, "var_e"]))
        if has_vc:
            g_var[i] = float(vc.loc[t, "var_g"])
            e_var[i] = float(vc.loc[t, "var_e"])
        else:
            # Fallback: BLUP variance as a proxy for total phenotypic
            # variance, split 50/50 between genetic and residual under a
            # prior h^2 = 0.5. The Soluble_Oxalate h^2 estimate (REML MLE
            # = 0.47, BGLR median = 0.73) supports the 0.5 fallback.
            total_var = float(p_blups[t].var(ddof=1))
            g_var[i] = 0.5 * total_var
            e_var[i] = 0.5 * total_var

    p_var = g_var + e_var

    # Build G and P from correlations * sqrt(var_i * var_j).
    P = p_corr * np.sqrt(np.outer(p_var, p_var))
    G = g_corr * np.sqrt(np.outer(g_var, g_var))

    # Ridge for invertibility at small n.
    ridge = 1e-3 * np.eye(len(INDEX_TRAITS))
    P = P + ridge
    G = G + ridge
    return P, G


# ---------------------------------------------------------------------------
# Smith-Hazel solver
# ---------------------------------------------------------------------------
def solve_smith_hazel(P: np.ndarray, G: np.ndarray,
                        a: np.ndarray) -> np.ndarray:
    """b = P^{-1} G a. Solved via np.linalg.solve."""
    return np.linalg.solve(P, G @ a)


def signed_weights(scenario_weights: dict[str, int]) -> np.ndarray:
    """Construct the per-trait `a` vector from a scenario weight dict.
    Direction is taken from TRAIT_DIRECTIONS; magnitude from the weight
    value. Traits absent from the scenario receive weight zero."""
    a = np.zeros(len(INDEX_TRAITS))
    for i, t in enumerate(INDEX_TRAITS):
        w = scenario_weights.get(t, 0)
        a[i] = TRAIT_DIRECTIONS[t] * w
    return a


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def plot_index_rank_changes(ranks: pd.DataFrame, robust: list[str],
                              out_path: Path) -> None:
    """Slopegraph: per-accession rank across the three scenarios.

    Three layers:
      * background = all accessions in light grey (high enough contrast to
        actually appear in print) so the panel is full of context, not blank
      * scenario-elite = top-10 in any single scenario, drawn in that
        scenario's Wong colour
      * robust elite = top-10 in all three scenarios, drawn in vermillion
        with bold per-accession labels on the right edge

    When `robust` is empty (common given the panel's small overlap between
    nutritional and anti-nutritional priorities), we still surface the
    scenario-specific elites so the figure is interpretable.
    """
    from _figstyle import publishable_axes, adjust_labels

    scenarios = list(SCENARIOS.keys())
    scenario_colour = {
        s: c for s, c in zip(scenarios,
                              [WONG["blue"], WONG["vermillion"],
                               WONG["green"]])
    }

    # Per-scenario top-10 elites (excluding the intersection)
    scenario_elites: dict[str, set[str]] = {}
    for s in scenarios:
        top10 = ranks.sort_values(f"rank_{s}").head(10).index.tolist()
        scenario_elites[s] = set(top10) - set(robust)

    fig, ax = plt.subplots(figsize=(11.5, 10.5))
    x = np.arange(len(scenarios))

    # Layer 1: every accession as a thin mid-grey line (always visible)
    for sample, row in ranks.iterrows():
        y = [row[f"rank_{s}"] for s in scenarios]
        if sample in robust:
            continue
        any_scenario = next((s for s in scenarios
                              if sample in scenario_elites[s]), None)
        if any_scenario is not None:
            continue
        ax.plot(x, y, color="#B0B0B0", linewidth=0.55, alpha=0.55, zorder=1)

    # Layer 2: per-scenario top-10 lines coloured by their winning scenario
    for s in scenarios:
        for sample in scenario_elites[s]:
            row = ranks.loc[sample]
            y = [row[f"rank_{sc}"] for sc in scenarios]
            ax.plot(x, y, color=scenario_colour[s], linewidth=1.4,
                     alpha=0.85, zorder=3)

    # Layer 3: robust elites (top-10 in ALL three scenarios) in vermillion
    robust_label_texts = []
    for sample in robust:
        row = ranks.loc[sample]
        y = [row[f"rank_{sc}"] for sc in scenarios]
        ax.plot(x, y, color=WONG["vermillion"], linewidth=2.0,
                 alpha=1.0, zorder=5)
        robust_label_texts.append(
            ax.text(len(scenarios) - 1 + 0.07, y[-1], sample,
                     fontsize=8.5, va="center", ha="left",
                     color=WONG["vermillion"], fontweight="bold", zorder=6))
    adjust_labels(robust_label_texts, ax=ax,
                   only_move={"text": "y", "static": "y"},
                   force_text=(0.0, 0.9))

    # Per-scenario elite labels on the right edge, in their scenario colour.
    # Multiple scenarios may share elites -- collect each accession once,
    # picking the colour of the scenario where it ranks highest.
    edge_labels: dict[str, tuple[float, str]] = {}
    for s in scenarios:
        for sample in scenario_elites[s]:
            row = ranks.loc[sample]
            y_end = row[f"rank_{scenarios[-1]}"]
            if sample not in edge_labels or \
                    row[f"rank_{s}"] < ranks.loc[sample,
                                                  f"rank_{edge_labels[sample][1]}"]:
                edge_labels[sample] = (y_end, s)

    # Lay edge labels out as a vertical ladder split into TOP and BOTTOM
    # half-ladders so dense top-rank elites and dense bottom-rank elites
    # never collide. Each label is connected to its data tip by a thin
    # leader line in the scenario colour.
    labels_sorted = sorted(edge_labels.items(),
                            key=lambda kv: kv[1][0])
    n_labels = len(labels_sorted)
    n_total  = max(len(ranks), 1)
    # Anchor row height to font size, not data: 7.5pt font + 1.4pt leading
    # = ~9pt per row; over a 95-rank y-axis on a 10-in canvas the y unit
    # ~ 0.105 inch ~ 7.6pt, so 1 rank ~ 1 line. Bump to 3.0 ranks per row
    # to add safety margin.
    row_height_in_ranks = max(3.0, n_total * 0.032)
    mid = n_labels // 2
    for i, (sample, (y_end, s)) in enumerate(labels_sorted):
        if i < mid:
            # top half-ladder: anchor at topmost label y, descend
            anchor_y = labels_sorted[0][1][0]
            ladder_y = anchor_y + i * row_height_in_ranks
        else:
            # bottom half-ladder: anchor at bottommost label y, ascend
            anchor_y = labels_sorted[-1][1][0]
            ladder_y = anchor_y - (n_labels - 1 - i) * row_height_in_ranks
        x_label = len(scenarios) - 1 + 0.15
        ax.text(x_label, ladder_y, sample,
                 fontsize=7.5, va="center", ha="left",
                 color=scenario_colour[s], zorder=6)
        # leader line from end-of-data tip to label anchor
        ax.plot([len(scenarios) - 1, x_label - 0.01],
                [y_end, ladder_y],
                color=scenario_colour[s], linewidth=0.5,
                alpha=0.55, zorder=4)

    ax.set_xticks(np.arange(len(scenarios)))
    ax.set_xticklabels([s.replace("_", " ").title() for s in scenarios],
                         fontsize=10.5)
    ax.set_ylabel("Per-scenario rank (1 = top elite)")
    ax.invert_yaxis()

    if robust:
        subtitle = (f"Robust elites (top-10 in all three scenarios; n = "
                     f"{len(robust)}) in vermillion -- scenario-specific "
                     f"elites coloured by winning scenario")
    else:
        subtitle = ("No accession is top-10 under all three scenarios -- "
                     "scenario-specific elites coloured by winning scenario")
    # Use suptitle so the multi-line title sits in the figure margin and
    # cannot land on top of the topmost ladder label.
    fig.suptitle("Per-accession rank under three Smith-Hazel "
                  "selection-index scenarios\n" + subtitle,
                  fontsize=10.5, y=0.995)
    # Top y-limit at -6 (visually above rank 1) gives the topmost label
    # ~6 ranks of head-room beneath the suptitle.
    ax.set_ylim(top=-6)
    ax.set_xlim(-0.25, len(scenarios) - 1 + 0.85)
    publishable_axes(ax, grid="y", grid_alpha=0.20)

    # Scenario-colour legend placed in the figure margin below the x-axis
    # tick labels so it never overlaps the ladder labels at the top-left
    # of the data area.
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], color=scenario_colour[s], lw=2.0,
                label=f"{s.replace('_', ' ').title()} top-10")
        for s in scenarios
    ]
    if robust:
        legend_handles.insert(0, Line2D([0], [0], color=WONG["vermillion"],
                                          lw=2.5, label="Robust (all-3) elites"))
    # Reserve bottom 10% of the figure for the legend; place axis above it.
    fig.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.10)
    fig.legend(handles=legend_handles,
                loc="lower center",
                ncol=len(legend_handles),
                bbox_to_anchor=(0.5, 0.005),
                fontsize=9.5, framealpha=0.9, edgecolor="none")

    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_weights_forest(weights_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.5))
    scenarios = list(SCENARIOS.keys())
    width = 0.27
    x = np.arange(len(INDEX_TRAITS))
    for k, (scenario, colour) in enumerate(zip(scenarios,
                                                  [WONG["blue"], WONG["vermillion"],
                                                   WONG["green"]])):
        b = [weights_df.loc[(weights_df["scenario"] == scenario)
                              & (weights_df["trait"] == t), "b"].values[0]
              for t in INDEX_TRAITS]
        ax.bar(x + (k - 1) * width, b, width=width,
                color=colour, edgecolor="black", linewidth=0.4,
                label=scenario.replace("_", " ").title())
    ax.axhline(0, color="grey", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(INDEX_TRAITS, rotation=30, ha="right")
    ax.set_ylabel("Smith-Hazel index weight b")
    ax.set_title("Smith-Hazel selection-index weights per trait per scenario")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"[load] GEBVs and BLUPs")
    gebvs = pd.read_csv(GEBV_CSV).set_index("sample")
    blups = pd.read_csv(BLUP_CSV).set_index("sample")
    common = sorted(set(gebvs.index) & set(blups.index))
    gebvs = gebvs.loc[common]
    blups = blups.loc[common]
    print(f"  {len(common)} accessions x {len(INDEX_TRAITS)} index traits")

    vc = load_variance_components()
    P, G = build_PG_matrices(blups, gebvs, vc)
    print(f"  P diag mean = {np.mean(np.diag(P)):.3g}; G diag mean = {np.mean(np.diag(G)):.3g}")

    # Build b per scenario.
    weights_rows = []
    index_rows = {s: pd.Series(np.zeros(len(common)), index=common)
                   for s in SCENARIOS}
    for scenario, w_dict in SCENARIOS.items():
        a = signed_weights(w_dict)
        b = solve_smith_hazel(P, G, a)
        for t, bv in zip(INDEX_TRAITS, b):
            weights_rows.append({"scenario": scenario, "trait": t, "b": float(bv)})

        # Per-accession index value. The GEBVs already encode breeding
        # values on the trait scale; the index combines them via b.
        y = gebvs[INDEX_TRAITS].copy()
        # Mean-impute missing GEBVs per trait so accessions with the
        # oxalate triplet missing (n = 41 of 95 panel) still receive an
        # index value built from their available trait GEBVs. Reviewers
        # will see the imputation choice; flag this in the manuscript.
        y = y.fillna(y.mean())
        I = y.values @ b
        index_rows[scenario] = pd.Series(I, index=common)

    weights_df = pd.DataFrame(weights_rows)
    weights_csv = TAB / "index_weights_b.csv"
    weights_df.to_csv(weights_csv, index=False)
    print(f"[result] wrote {weights_csv}")

    per_accession = pd.DataFrame(index_rows)
    per_accession.columns = [f"index_{s}" for s in SCENARIOS]
    # Defensive: report any accessions whose index value came out NaN and
    # drop them before ranking. NaN can occur if the variance-component
    # estimation produced NaN for a trait and propagated through b.
    for s in SCENARIOS:
        n_nan = per_accession[f"index_{s}"].isna().sum()
        if n_nan:
            print(f"  [warn] scenario {s}: {n_nan} accessions with NaN index value; dropping")
    finite = per_accession.dropna(how="any")
    print(f"  per-accession index built for {len(finite)} of {len(per_accession)} accessions")
    per_accession = finite
    for s in SCENARIOS:
        per_accession[f"rank_{s}"] = per_accession[f"index_{s}"].rank(
            method="min", ascending=False).astype(int)
    per_accession_csv = TAB / "per_accession_index_per_scenario.csv"
    per_accession.to_csv(per_accession_csv)
    print(f"[result] wrote {per_accession_csv}")

    # Robust vs scenario-dependent elites.
    top10_per_scenario = {
        s: set(per_accession.sort_values(f"index_{s}", ascending=False).head(10).index)
        for s in SCENARIOS
    }
    robust_elites = set.intersection(*top10_per_scenario.values())
    elite_rows = []
    for s in SCENARIOS:
        for sample in top10_per_scenario[s]:
            scenarios_in = sum(sample in top10_per_scenario[ss] for ss in SCENARIOS)
            elite_rows.append({
                "sample": sample,
                "in_scenario": s,
                "scenarios_in_top10": scenarios_in,
                "is_robust_elite": sample in robust_elites,
                f"index": per_accession.loc[sample, f"index_{s}"],
                f"rank": per_accession.loc[sample, f"rank_{s}"],
            })
    elite_df = pd.DataFrame(elite_rows)
    elite_csv = TAB / "robust_vs_scenario_elites.csv"
    elite_df.to_csv(elite_csv, index=False)
    print(f"[result] wrote {elite_csv}")

    print(f"\n[plot] index rank slopegraph")
    plot_index_rank_changes(per_accession, list(robust_elites),
                             FIG / "fig_index_rank_changes")
    print(f"[plot] weights forest")
    plot_weights_forest(weights_df, FIG / "fig_index_weights_forest")

    print("\n=== manuscript-text summary ===")
    print(f"  robust elites (top-10 under all 3 scenarios): {sorted(robust_elites)}")
    for s in SCENARIOS:
        print(f"  top-10 under {s} : "
              f"{sorted(top10_per_scenario[s], key=lambda x: per_accession.loc[x, f'rank_{s}'])}")


if __name__ == "__main__":
    main()
