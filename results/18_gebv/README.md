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

## Figure 13A regenerated to the restricted (primary) merit index (reviewer decision, 2026-10-27)

Fig 13A (`fig48_merit_index_top20`) now shows the RESTRICTED merit index — the mean of
direction-signed z-score GEBVs over only the traits predicted above chance (favourable-
direction traits with leave-one-out CV r > 0): Antioxidant, Mass_of_Seeds, Total_Oxalate,
Soluble_Oxalate, Insoluble_Oxalate. The five traits with LOO r < 0 (Tannin, Phenol,
Flavonoid, Seed_Coat_Tannin, Crude_Protein) and the three direction-neutral seed-size
dimensions are excluded.

- Script: `scripts/_replot_fig13a_restricted.py` (recomputes the index from the existing
  per-trait GEBVs + LOO accuracy; the primary `gebv_per_accession.csv` is NOT re-fit).
- New table: `tables/merit_index_top20_restricted.csv`.
- Restricted top-3: TSs217, TSs15, TSs131 (all-trait top-3 were TSs157A, TSs282, TSs38;
  only TSs38 is shared, at restricted rank 10 — consistent with the full-vs-restricted
  Spearman ρ = 0.33 reported in §3.13).

Fig 13B (fig49 radar) was also regenerated (2026-10-28) for the restricted top-10 over
the five predictable traits (script `_replot_fig13a_restricted.py`; table
`merit_index_top10_radar_restricted.csv`), so both Figure 13 panels are now on the
restricted basis. STILL OUTSTANDING: §3.1 and Figure 2A annotate TSs157A as "the top
merit-index accession" (the all-trait top); under the restricted-primary framing the top
accession is TSs217 -- these references/annotations should be reconciled (flagged for Benjamin).
