#!/usr/bin/env python
"""
24b_boundary_lrt -- calibrated boundary likelihood-ratio test of marker-based
heritability against zero, post-processed from the saved REML profile curves
(script 24).

Rationale
---------
The REML profile MLE for soluble oxalate sits at the upper grid edge (0.999)
and the one-sided lower bound (0.14) is only a descriptive edge of the 95 %
profile band. Reviewer 1 (R1.10) flagged the boundary MLE. The correct
calibrated test of H0: h^2 = 0 for a variance component that lies on the
boundary of the parameter space uses the 1/2 chi^2_0 + 1/2 chi^2_1 mixture
(Self & Liang, 1987) as the reference distribution. The LR statistic is read
directly off the profile: at the grid floor (h^2 ~= 0) the stored column
neg2_delta_logL = -2 (logL - logL_max) IS 2 (logL_max - logL_0), the LR stat.

No model is refitted here; this is a pure post-processing step on
results/24_bootstrap_h2/tables/h2_profile_curves.csv so it cannot perturb any
existing output.

Output
------
    results/24_bootstrap_h2/tables/h2_boundary_lrt.csv
        per-trait LR statistic at h^2 = 0, calibrated mixture p, naive 1-df p.

Author: Benjamin Narh-Madey
"""
from pathlib import Path
import pandas as pd
from scipy import stats

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TAB = ROOT / "results" / "24_bootstrap_h2" / "tables"
CURVES = TAB / "h2_profile_curves.csv"
OUT = TAB / "h2_boundary_lrt.csv"


def main() -> None:
    cur = pd.read_csv(CURVES)
    dcol = "neg2_delta_logL"
    rows = []
    for trait, sub in cur.groupby("trait"):
        sub = sub.sort_values("h2")
        floor = sub.iloc[0]                       # grid floor h^2 ~= 0
        lrt = float(floor[dcol])                  # = 2 (logL_max - logL_0)
        p_mix = 0.5 * stats.chi2.sf(lrt, 1)       # 1/2 chi2_0 + 1/2 chi2_1
        p_naive = float(stats.chi2.sf(lrt, 1))
        rows.append({
            "trait": trait,
            "h2_grid_floor": float(floor["h2"]),
            "lr_stat_h2_eq_0": lrt,
            "p_calibrated_mixture": p_mix,
            "p_naive_1df": p_naive,
        })
    out = pd.DataFrame(rows).sort_values("p_calibrated_mixture")
    out.to_csv(OUT, index=False)
    with pd.option_context("display.float_format", lambda v: f"{v:.4g}"):
        print(out.to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
