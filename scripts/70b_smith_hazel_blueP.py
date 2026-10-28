#!/usr/bin/env python
"""
Smith-Hazel selection index with P estimated from the PHENOTYPIC CORRELATIONS
AMONG STAGE-ONE BLUEs (reviewer R1.9), as a sensitivity alternative to the primary
index, whose P off-diagonals come from BLUP correlations (scripts/70_selection_index.py).

Everything else is held identical to the primary index by importing script 70's
constants and solver (INDEX_TRAITS, TRAIT_DIRECTIONS, SCENARIOS, signed_weights,
solve_smith_hazel, build_PG_matrices). Only P changes:
    P_BLUE = Cov(stage-one BLUEs over INDEX_TRAITS) + ridge
G (genetic) and the per-accession GEBV vector are unchanged, so any difference in the
elite ranking is attributable solely to the P source.

Inputs : results/18_gebv/tables/gebv_per_accession.csv (GEBVs, G + index)
         results/00_blup_pipeline/tables/phenotype_blups.csv (BLUPs, primary P + var fallback)
         results/00_blup_pipeline/tables/variance_components.csv
         results/R1_blue_pipeline/tables/phenotype_blues.csv (stage-one BLUEs, the new P)
Outputs: results/70_selection_index/tables/
           index_weights_b_blueP.csv
           per_accession_index_per_scenario_blueP.csv
           blueP_vs_blupP_comparison.csv   (Spearman + top-10 overlap per scenario)
Author : Benjamin Narh-Madey
"""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
spec = importlib.util.spec_from_file_location("sh70", ROOT / "scripts/70_selection_index.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

TAB = ROOT / "results/70_selection_index/tables"
BLUE_CSV = ROOT / "results/R1_blue_pipeline/tables/phenotype_blues.csv"
RIDGE = 1e-3

gebvs = pd.read_csv(m.GEBV_CSV).set_index("sample")
blups = pd.read_csv(m.BLUP_CSV).set_index("sample")
blues = pd.read_csv(BLUE_CSV).set_index("sample")
vc = m.load_variance_components()

common = sorted(set(gebvs.index) & set(blups.index))
gebvs, blups = gebvs.loc[common], blups.loc[common]
T = m.INDEX_TRAITS

# The BLUP-correlation P is now the SENSITIVITY estimate; the BLUE phenotypic
# covariance (script 70's primary) is the reference. G is shared.
blue_common = blues.reindex(common)
P_blup, G = m.build_PG_matrices(blups, gebvs, vc, blues=blue_common, p_source="blup")
P_blue = blue_common[T].cov().values + RIDGE * np.eye(len(T))
print(f"[blueP] P from {blue_common[T].dropna(how='all').shape[0]} BLUE rows; "
      f"P_blue diag mean = {np.mean(np.diag(P_blue)):.3g} vs P_blup {np.mean(np.diag(P_blup)):.3g}")


def index_per_scenario(P):
    wrows, idx = [], {}
    for sc, wd in m.SCENARIOS.items():
        a = m.signed_weights(wd)
        b = m.solve_smith_hazel(P, G, a)
        for t, bv in zip(T, b):
            wrows.append({"scenario": sc, "trait": t, "b": float(bv)})
        y = gebvs[T].copy()
        y = y.fillna(y.mean())
        idx[sc] = pd.Series(y.values @ b, index=common)
    return pd.DataFrame(wrows), pd.DataFrame(idx)


w_blue, I_blue = index_per_scenario(P_blue)
_, I_blup = index_per_scenario(P_blup)

w_blue.to_csv(TAB / "index_weights_b_blueP.csv", index=False)
out = I_blue.copy()
out.columns = [f"index_{s}" for s in m.SCENARIOS]
for s in m.SCENARIOS:
    out[f"rank_{s}"] = out[f"index_{s}"].rank(method="min", ascending=False).astype(int)
out.to_csv(TAB / "per_accession_index_per_scenario_blueP.csv")

# Compare BLUE-P vs BLUP-P: Spearman of the index + top-10 elite overlap, per scenario
cmp = []
for s in m.SCENARIOS:
    rho, _ = spearmanr(I_blue[s], I_blup[s])
    top_blue = set(I_blue[s].sort_values(ascending=False).head(10).index)
    top_blup = set(I_blup[s].sort_values(ascending=False).head(10).index)
    cmp.append({"scenario": s, "spearman_blueP_vs_blupP": round(rho, 3),
                "top10_overlap": len(top_blue & top_blup)})
cmp = pd.DataFrame(cmp)
cmp.to_csv(TAB / "blueP_vs_blupP_comparison.csv", index=False)
print("\nBLUE-P vs BLUP-P (primary) agreement:")
print(cmp.to_string(index=False))
