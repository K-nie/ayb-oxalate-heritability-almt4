#!/usr/bin/env python
"""
Regenerate Figure 13A (fig48_merit_index_top20) to show the RESTRICTED (primary)
merit index instead of the all-trait index (reviewer decision).

The restricted index averages the direction-signed z-score GEBVs over only the
traits predicted above chance (favourable-direction traits with leave-one-out CV
r > 0): Antioxidant, Mass_of_Seeds, Total/Soluble/Insoluble Oxalate. The five traits
not predicted above chance (Tannin, Phenol, Flavonoid, Seed_Coat_Tannin, Crude_Protein)
are dropped, as are the three direction-neutral seed-size dimensions.

Reads the per-trait GEBVs and LOO accuracy already on disk; recomputes only the index.
No GEBV re-fit (the primary gebv_per_accession.csv is left untouched).

Inputs : results/18_gebv/tables/gebv_per_accession.csv
         results/R1_blue_pipeline/tables/cascade_blup_vs_blue_diff.csv (LOO r)
Outputs: results/18_gebv/tables/merit_index_top20_restricted.csv
         results/18_gebv/figures/fig48_merit_index_top20.{png,pdf}  (Figure 13A, regenerated)
Author : Benjamin Narh-Madey
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
sys.path.insert(0, str(ROOT / "scripts"))
from _plotstyle import apply, WONG
apply()
TAB = ROOT / "results/18_gebv/tables"
FIG = ROOT / "results/18_gebv/figures"

DIRECTION = {"Tannin": "low", "Phenol": "high", "Flavonoid": "high",
             "Antioxidant": "high", "Mass_of_Seeds": "high",
             "Seed_Coat_Tannin": "low", "Crude_Protein": "high",
             "Total_Oxalate": "low", "Soluble_Oxalate": "low",
             "Insoluble_Oxalate": "high"}  # size dims deliberately absent

gebv = pd.read_csv(TAB / "gebv_per_accession.csv").set_index("sample")
loo = pd.read_csv(ROOT / "results/R1_blue_pipeline/tables/cascade_blup_vs_blue_diff.csv"
                  ).set_index("trait")["loo_r_BLUE_unw"]

# predictable favourable-direction traits = LOO r > 0
predictable = [t for t in DIRECTION if t in loo.index and loo[t] > 0]
expected = {"Antioxidant", "Mass_of_Seeds", "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"}
assert set(predictable) == expected, f"predictable set {set(predictable)} != {expected}"
print(f"[restricted] predictable traits (LOO r>0): {predictable}")

z = gebv[predictable].copy()
for t in predictable:
    s = z[t]
    if s.std() == 0 or s.notna().sum() < 5:
        continue
    zt = (s - s.mean()) / s.std()
    if DIRECTION[t] == "low":
        zt = -zt
    z[t] = zt
gebv["merit_index_restricted"] = z.mean(axis=1, skipna=True)
gebv["merit_index_restricted_n"] = z.notna().sum(axis=1)
r = gebv.sort_values("merit_index_restricted", ascending=False)
# Persist the restricted (primary) merit index back into the per-accession GEBV table
# so Supplementary Table S3 carries BOTH the all-trait (merit_index) and the restricted
# (merit_index_restricted) indices. Adding a column only; existing columns untouched.
gebv.to_csv(TAB / "gebv_per_accession.csv")
top20 = r.head(20).copy()
top20[["merit_index_restricted", "merit_index_restricted_n"] + predictable].to_csv(
    TAB / "merit_index_top20_restricted.csv")
print("[restricted] top-10:", list(top20.index[:10]))

# Figure 13A — regenerated with the restricted index
fig, ax = plt.subplots(figsize=(6.5, 8), constrained_layout=True)
top_show = top20.iloc[::-1]
y = np.arange(len(top_show))
colors = [WONG["blue"] if v > 0 else WONG["vermillion"] for v in top_show["merit_index_restricted"]]
ax.barh(y, top_show["merit_index_restricted"], color=colors, edgecolor="black", linewidth=0.4)
ax.set_yticks(y)
ax.set_yticklabels(top_show.index, fontsize=8)
ax.set_xlabel("Restricted merit index (avg favourable-direction z-score GEBV,\n"
              "predictable traits only: antioxidant, seed mass, 3 oxalate fractions)")
ax.axvline(0, color="black", linewidth=0.6)
ax.grid(True, axis="x")
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig48_merit_index_top20.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("[fig] regenerated fig48_merit_index_top20 (Figure 13A, restricted index)")

# ---------------------------------------------------------------------------
# Figure 13B — radar of the restricted top-10 across the five predictable traits
# ---------------------------------------------------------------------------
top10 = r.head(10)
zmat = z.loc[top10.index, predictable].fillna(0).values
zmin, zmax = z[predictable].min().values, z[predictable].max().values
rng = np.where(zmax > zmin, zmax - zmin, 1.0)
zmat_norm = np.clip((zmat - zmin) / rng, 0.0, 1.0)
theta = np.linspace(0, 2 * np.pi, len(predictable), endpoint=False)
theta_closed = np.concatenate([theta, theta[:1]])
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection="polar"),
                       constrained_layout=True)
COLORS = plt.cm.tab10.colors
NICE = {"Antioxidant": "Antioxidant", "Mass_of_Seeds": "Seed mass",
        "Total_Oxalate": "Total oxalate", "Soluble_Oxalate": "Soluble oxalate",
        "Insoluble_Oxalate": "Insoluble oxalate"}
for i in range(len(top10)):
    vals = np.concatenate([zmat_norm[i], zmat_norm[i:i + 1, 0].flatten()])
    ax.plot(theta_closed, vals, linewidth=1.5, color=COLORS[i % 10], label=top10.index[i])
    ax.fill(theta_closed, vals, alpha=0.08, color=COLORS[i % 10])
ax.set_xticks(theta)
ax.set_xticklabels([NICE[t] for t in predictable], fontsize=8)
ax.set_yticks([0.25, 0.5, 0.75, 1.0]); ax.set_yticklabels(["", "", "", ""])
ax.set_ylim(0, 1.05)
ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.05), fontsize=8, frameon=False)
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig49_merit_index_radar_top10.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
pd.DataFrame(zmat_norm, index=top10.index, columns=predictable).to_csv(
    TAB / "merit_index_top10_radar_restricted.csv")
print("[fig] regenerated fig49_merit_index_radar_top10 (Figure 13B, restricted top-10)")
