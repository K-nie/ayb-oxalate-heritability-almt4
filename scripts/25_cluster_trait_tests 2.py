#!/usr/bin/env python3
"""
PCA cluster x trait differentiation: do the two marker-PCA clusters
(Cluster 1, n = 11 PC1-high outliers; Cluster 2, n = 84 main bulk) differ in
phenotype?

Tests
-----
Per trait, compare the trait distribution between Cluster 1 and Cluster 2:
  * Mann-Whitney U (Wilcoxon rank-sum) — non-parametric, robust to skew.
  * Rank-biserial correlation (effect-size equivalent of Cohen's d,
    Kerby 2014): r_rb = 1 - 2 U / (n1 * n2). Range -1 to +1.
  * Cliff's delta (an interpretable companion effect size).
  * Benjamini-Hochberg FDR over the 13 tests.

Caveat: clusters are derived from PCA on the same SNPs, so this test
quantifies phenotype differentiation across an already-identified structure
- not independent evidence of structure. The point is to ask whether the
genetic structure has phenotypic consequences (i.e. is the small Cluster
distinct on which traits?).

Outputs (results/25_cluster_trait_tests/)
-----------------------------------------
tables/
    cluster_trait_tests.csv      one row per trait
figures/
    fig62_cluster_trait_effects.png/.pdf  effect-size forest plot
    fig63_cluster_trait_boxplots.png/.pdf  per-trait box plots side-by-side

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from _plotstyle import apply, WONG, CLUSTER_PAL
from _pheno import load_phenotypes, TRAITS_ALL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PCA_CSV = ROOT / "results" / "01_qc_pca_power" / "tables" / "pca_coords.csv"
OUT = ROOT / "results" / "25_cluster_trait_tests"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def cliffs_delta(a, b):
    """Cliff's delta = (#(a > b) - #(a < b)) / (n1 * n2)."""
    a = np.asarray(a); b = np.asarray(b)
    gt = sum((ai > bj) for ai in a for bj in b)
    lt = sum((ai < bj) for ai in a for bj in b)
    return (gt - lt) / (len(a) * len(b))


def cliffs_delta_ci(a, b, alpha=0.05):
    """Asymptotic 95 % CI on Cliff's delta (Long, Feng & Cliff 2003).

    Uses the consistent variance estimator
        var(d) = (S1 + S2 - n1 * n2 * d^2) / (n1 * n2 * (n1 - 1) * (n2 - 1))
    where S1 and S2 are derived from the dominance score for each
    element. The asymptotic normal CI is appropriate at small n where
    bootstrap of delta is unstable (and the same domain Cliff 1993
    targeted with this estimator).
    """
    a = np.asarray(a); b = np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return float("nan"), float("nan")
    d = cliffs_delta(a, b)
    # dominance score for each a_i: P(a_i > b) - P(a_i < b)
    di = np.array([np.mean(ai > b) - np.mean(ai < b) for ai in a])
    # dominance score for each b_j: P(a > b_j) - P(a < b_j) — note the
    # symmetric sign convention used by Long et al. 2003.
    dj = np.array([np.mean(a > bj) - np.mean(a < bj) for bj in b])
    s1 = np.var(di, ddof=1) if n1 > 1 else 0.0
    s2 = np.var(dj, ddof=1) if n2 > 1 else 0.0
    var_d = ((n2 - 1) * s1 + (n1 - 1) * s2) / (n1 * n2)
    se = float(np.sqrt(max(var_d, 0.0)))
    z_crit = stats.norm.ppf(1 - alpha / 2)
    # logit transformation per Long 2003 to keep CI inside [-1, 1]
    if abs(d) >= 1 or se == 0:
        return float(max(-1.0, d - z_crit * se)), float(min(1.0, d + z_crit * se))
    fl = np.arctanh(d)
    var_fl = var_d / max((1 - d * d) ** 2, 1e-12)
    se_fl = float(np.sqrt(var_fl))
    return float(np.tanh(fl - z_crit * se_fl)), float(np.tanh(fl + z_crit * se_fl))


def min_detectable_delta(n_small, n_large, power=0.80, alpha=0.05):
    """Minimum two-sided Cliff's delta detectable at the given power.

    Approximation via the Mann-Whitney U power formula expressed in
    delta. Under H0, U is normal with mean n1*n2/2 and variance
    n1*n2*(n1+n2+1)/12. Solving for the smallest |delta| that achieves
    the target power gives:
        |delta_min| = (z_alpha + z_power) * sqrt((n1+n2+1) / (3*n1*n2))
    (Lehmann 1975; conservative for tied data).
    """
    if min(n_small, n_large) < 2:
        return float("nan")
    z_a = stats.norm.ppf(1 - alpha / 2)
    z_b = stats.norm.ppf(power)
    return float((z_a + z_b) * np.sqrt(
        (n_small + n_large + 1) / (3.0 * n_small * n_large)))


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
print("[load] phenotype + cluster labels...")
pheno = load_phenotypes()
pca = pd.read_csv(PCA_CSV).set_index("sample")
# Determine the larger vs smaller cluster automatically — at 0.90 QC the
# k-means label assigned to the small subgroup is whichever has fewer lines.
cluster_sizes = pca["cluster"].value_counts().to_dict()
small_k = min(cluster_sizes, key=cluster_sizes.get)
large_k = max(cluster_sizes, key=cluster_sizes.get)
print(f"  Cluster sizes: small (k={small_k}) = {cluster_sizes[small_k]}, "
      f"large (k={large_k}) = {cluster_sizes[large_k]}")


# ---------------------------------------------------------------------------
# Per-trait Mann-Whitney + effect sizes
# ---------------------------------------------------------------------------
rows = []
for t in TRAITS_ALL:
    df = pheno[[t]].dropna().reset_index()
    df = df.merge(pca[["cluster"]], left_on="sample", right_index=True)
    a = df.loc[df["cluster"] == small_k, t].values
    b = df.loc[df["cluster"] == large_k, t].values
    if len(a) < 3 or len(b) < 3:
        continue
    U, pval = stats.mannwhitneyu(a, b, alternative="two-sided")
    rb = 1.0 - (2.0 * U) / (len(a) * len(b))   # rank-biserial
    delta = cliffs_delta(a, b)
    d_lo, d_hi = cliffs_delta_ci(a, b)
    d_min = min_detectable_delta(len(a), len(b))
    rows.append({"trait": t,
                 "n_small_cluster": int(len(a)),
                 "n_large_cluster": int(len(b)),
                 "median_small": float(np.median(a)),
                 "median_large": float(np.median(b)),
                 "median_diff_small_minus_large": float(np.median(a) - np.median(b)),
                 "mannwhitney_U": float(U),
                 "p_value": float(pval),
                 "rank_biserial_r": float(rb),
                 "cliffs_delta": float(delta),
                 "cliffs_delta_ci95_low": d_lo,
                 "cliffs_delta_ci95_high": d_hi,
                 "min_detectable_delta_80pct": d_min,
                 "detectable_below_threshold": bool(abs(delta) < d_min)})
res = pd.DataFrame(rows)
res["q_bh"] = stats.false_discovery_control(res["p_value"])
res = res.sort_values("p_value")
res.to_csv(TAB / "cluster_trait_tests.csv", index=False)

print("\n[results] (sorted by Mann-Whitney p):")
print(res[["trait", "n_small_cluster", "n_large_cluster",
           "median_small", "median_large", "p_value", "q_bh",
           "cliffs_delta", "cliffs_delta_ci95_low", "cliffs_delta_ci95_high",
           "min_detectable_delta_80pct"]].round(3).to_string(index=False))

# Headline power note: largest, smallest, and oxalate-subset detectable deltas
power_note = res[["trait", "n_small_cluster", "n_large_cluster",
                   "min_detectable_delta_80pct"]].drop_duplicates(
    subset=["n_small_cluster", "n_large_cluster"])
print("\n[power] minimum detectable Cliff's delta at 80% power, alpha = 0.05:")
print(power_note.to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 62. Forest plot of rank-biserial effect size with Cliff's delta CI
# ---------------------------------------------------------------------------
from _figstyle import publishable_axes  # noqa: E402

fig, ax = plt.subplots(figsize=(10.5, 6.5), constrained_layout=True)
df_show = res.sort_values("rank_biserial_r", ascending=True).reset_index(drop=True)
y = np.arange(len(df_show))

# Wong-palette colouring: + direction = vermillion, - direction = blue, with
# q-significant traits getting full saturation and ns traits getting reduced
# alpha so the eye is drawn to the load-bearing rows.
for i, row in df_show.iterrows():
    base_col = WONG["vermillion"] if row["rank_biserial_r"] >= 0 else WONG["blue"]
    sig = bool(row["q_bh"] < 0.10)
    alpha = 0.90 if sig else 0.55
    ax.barh(i, row["rank_biserial_r"], color=base_col, alpha=alpha,
            edgecolor="black", linewidth=0.45, zorder=3)

# Place p/q labels OUTSIDE the bars so they don't compete with the colour.
x_pad = 0.025
for i, row in df_show.iterrows():
    sig_mark = " *" if row["q_bh"] < 0.10 else ""
    label = f"p = {row['p_value']:.3f}  q = {row['q_bh']:.3f}{sig_mark}"
    if row["rank_biserial_r"] >= 0:
        ax.text(row["rank_biserial_r"] + x_pad, i, label,
                ha="left", va="center", fontsize=8, color="#222",
                zorder=5)
    else:
        ax.text(row["rank_biserial_r"] - x_pad, i, label,
                ha="right", va="center", fontsize=8, color="#222",
                zorder=5)

ax.set_yticks(y)
ax.set_yticklabels([t.replace("_", " ") for t in df_show["trait"]],
                    fontsize=10)
ax.axvline(0, color="black", linewidth=0.7, zorder=2)
ax.set_xlabel(r"Rank-biserial r (small cluster minus large cluster)",
              fontsize=11)
ax.set_xlim(-1.20, 1.20)
ax.set_title("Per-trait differentiation between marker-PCA clusters\n"
             "(positive = trait elevated in the small PC1-outlier cluster; "
             "* = BH-q < 0.10)",
             fontsize=11)
publishable_axes(ax, grid="x", grid_alpha=0.20)

# Legend explaining colour
from matplotlib.patches import Patch as _Patch
legend_handles = [
    _Patch(facecolor=WONG["vermillion"], alpha=0.90, edgecolor="black",
           linewidth=0.45, label="Elevated in small cluster (significant)"),
    _Patch(facecolor=WONG["vermillion"], alpha=0.55, edgecolor="black",
           linewidth=0.45, label="Elevated in small cluster (ns)"),
    _Patch(facecolor=WONG["blue"], alpha=0.90, edgecolor="black",
           linewidth=0.45, label="Depressed in small cluster (significant)"),
    _Patch(facecolor=WONG["blue"], alpha=0.55, edgecolor="black",
           linewidth=0.45, label="Depressed in small cluster (ns)"),
]
ax.legend(handles=legend_handles, loc="lower right",
          fontsize=8.5, framealpha=0.92, edgecolor="none",
          ncol=1)

save(fig, "fig62_cluster_trait_effects")
print("[fig] fig62_cluster_trait_effects")


# ---------------------------------------------------------------------------
# Fig 63. Per-trait box plots side-by-side, 4 x 4 grid (13 traits + spares)
# ---------------------------------------------------------------------------
n_cols = 4
n_rows = int(np.ceil(len(TRAITS_ALL) / n_cols))
fig, axes = plt.subplots(n_rows, n_cols,
                          figsize=(3.6 * n_cols, 3.0 * n_rows),
                          constrained_layout=True)
for ax, t in zip(axes.ravel(), TRAITS_ALL):
    df = pheno[[t]].dropna().reset_index()
    df = df.merge(pca[["cluster"]], left_on="sample", right_index=True)
    if df.empty:
        ax.set_visible(False); continue
    a = df.loc[df["cluster"] == small_k, t].values
    b = df.loc[df["cluster"] == large_k, t].values

    bp = ax.boxplot(
        [a, b], positions=[0, 1], widths=0.55,
        patch_artist=True,
        showmeans=True, meanline=True,
        meanprops={"color": "black", "linewidth": 1.3, "linestyle": "--"},
        medianprops={"color": "white", "linewidth": 1.5},
        whiskerprops={"color": "#333", "linewidth": 0.8},
        capprops={"color": "#333", "linewidth": 0.8},
        flierprops={"marker": "o", "markersize": 3,
                     "markerfacecolor": "#888",
                     "markeredgecolor": "none", "alpha": 0.55},
    )
    for patch, k in zip(bp["boxes"], (small_k, large_k)):
        patch.set_facecolor(CLUSTER_PAL[k])
        patch.set_alpha(0.78)
        patch.set_edgecolor("#333")
        patch.set_linewidth(0.6)

    # jittered individual points so the reader sees the raw n
    rng = np.random.default_rng(0)
    for pos, vals in [(0, a), (1, b)]:
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.10, 0.10, size=len(vals))
        ax.scatter(pos + jitter, vals, s=12,
                    color=CLUSTER_PAL[(small_k, large_k)[pos]],
                    alpha=0.55, edgecolor="white", linewidth=0.4,
                    zorder=4)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"Cluster {small_k + 1}\n(n = {len(a)})",
         f"Cluster {large_k + 1}\n(n = {len(b)})"],
        fontsize=9,
    )
    ax.set_ylabel(t.replace("_", " "), fontsize=10)
    row_match = res[res["trait"] == t]
    if not row_match.empty:
        r0 = row_match.iloc[0]
        sig_mark = " *" if r0["q_bh"] < 0.10 else ""
        ax.set_title(
            f"{t.replace('_', ' ')}\n"
            f"p = {r0['p_value']:.3f}; "
            f"r$_b$ = {r0['rank_biserial_r']:+.2f}{sig_mark}",
            fontsize=10, pad=6,
        )
    else:
        ax.set_title(t.replace("_", " "), fontsize=10, pad=6)
    publishable_axes(ax, grid="y", grid_alpha=0.20)

for ax in axes.ravel()[len(TRAITS_ALL):]:
    ax.set_visible(False)

fig.suptitle(
    "Per-trait distributions split by marker-PCA cluster "
    "(jittered points overlaid; * = BH-q < 0.10)",
    y=1.005, fontsize=13,
)
save(fig, "fig63_cluster_trait_boxplots")
print("[fig] fig63_cluster_trait_boxplots")


print(f"\nOutputs in: {OUT}")
