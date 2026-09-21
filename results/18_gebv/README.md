# 18 — Per-accession GEBV table + multi-trait merit index

**Script:** `scripts/18_gebv.py`

## What was done

Full-panel GBLUP refit (no cross-validation split) for each of 13 traits gives one genomic-estimated breeding value (GEBV) per accession per trait. A favourable-direction-signed z-score composite multi-trait merit index ranks accessions for selection.

## Method

For each trait, the GBLUP α (ridge regularisation) was tuned by leave-one-out within-panel Pearson r over {0.5, 1, 2, 5, 10, 25}, then refit on full data; the prediction at training points is the GBLUP-implied additive breeding value. Z-scores are computed per trait across the panel; for "low = favourable" traits the sign is flipped so positive = breeder-desirable. The merit index averages signed z-scores across the favourable-direction trait subset. Size-direction-neutral seed traits (Length, Width, Thickness) are excluded from the index.

## Quick findings

- Top-10 elites by merit index:
  1. **TSs68** — biochemistry-elite (6 traits, merit = 1.20)
  2. TSs334 (6 traits, 0.68)
  3. TSs38 (6 traits, 0.61)
  4. TSs282 (6, 0.54) — but Seed_Thickness was corrected at the loader
  5. TSs22A (6, 0.54)
  6. **TSs10** — full 10-trait elite (n = 10, merit = 0.54)
  7. TSs368 (6, 0.51)
  8. **TSs331** — full 10-trait elite (n = 10, merit = 0.51)
  9. TSs325 (10, 0.45)
  10. TSs49A (10, 0.42)

- **TSs10 and TSs331 are the cleanest balanced 10-trait elites** — both have positive GEBVs for Crude_Protein and favourable directions on Insoluble_Oxalate.

## Caveats

- LOO-tuned α gives realistic in-panel accuracy; for new-line prediction the user should validate on a held-out independent cohort (currently unavailable).
- Merit-index trait weights are equal across favourable-direction traits. A breeder-weighted scheme (e.g. doubling crude protein and oxalate) would shift the ranking but the top-5 elites are stable across reasonable reweightings.
- The 6-trait elites only score on biochem traits because they're missing from the n = 46 protein/oxalate subset. Direct comparison with the 10-trait elites is apples-to-oranges; favour the 10-trait elites for breeding decisions.
