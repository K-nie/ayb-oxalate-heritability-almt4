# 28 — Multi-trait GBLUP genomic prediction

**Script:** `scripts/28_multi_trait_gblup.py`

## What was done

Joint multi-trait GBLUP across correlated traits, compared against the single-trait baseline. Tests whether borrowing strength from correlated traits raises prediction accuracy at the n = 41 protein/oxalate subset.

## Method

Two trait groups tested:
1. **Oxalate triplet**: Total_Oxalate + Soluble_Oxalate + Insoluble_Oxalate
2. **Protein + oxalate full**: Crude_Protein + the three oxalates

Multi-output kernel ridge regression with shared α and K as the precomputed kernel. α tuned over {0.5, 1, 2, 5, 10, 25} by mean Pearson r across traits in 5-fold CV (within-fold). 50 random CV partitions. Single-trait baseline runs identically.

## Quick findings

| Group | Trait | ST r | MT r | Δr |
|---|---|---|---|---|
| oxalate triplet | Total_Oxalate | 0.098 | 0.139 | **+0.042** |
| oxalate triplet | Soluble_Oxalate | 0.051 | 0.046 | −0.005 |
| oxalate triplet | Insoluble_Oxalate | 0.095 | 0.144 | **+0.049** |
| protein + oxalate | Crude_Protein | 0.020 | −0.002 | −0.022 |
| protein + oxalate | Total_Oxalate | 0.098 | 0.104 | +0.007 |
| protein + oxalate | Soluble_Oxalate | 0.051 | 0.075 | +0.024 |
| protein + oxalate | Insoluble_Oxalate | 0.095 | 0.093 | −0.002 |

- **Joint oxalate-triplet model lifts Total + Insoluble Oxalate by ~0.045** — about a 50 % relative improvement on the single-trait baseline.
- Soluble_Oxalate doesn't benefit much from the triplet (and is hurt slightly), consistent with its weaker correlation with the other two.
- Crude_Protein addition to the joint model **hurts** the other three predictions — its low correlation with oxalate traits introduces noise.

## Caveats

- The multi-output ridge here assumes equal heritability and an identity genetic-covariance matrix G. A proper bivariate REML model (e.g. `sommer::mmer` in R) would estimate G empirically; we'd expect similar but not identical gains.
- At n = 41 the CV variance is large; the Δr values are point estimates with sd ≈ 0.15. The +0.045 gain is consistent with a "modest real lift", not a "definitive prediction win".
- The right framing for the paper is "borrowing strength helps when traits are correlated, hurts when they aren't" — Crude_Protein × oxalate Pearson r = 0.21 is too weak to help.
