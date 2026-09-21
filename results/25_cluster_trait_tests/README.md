# 25 — Marker-PCA cluster × trait Mann-Whitney tests

**Script:** `scripts/25_cluster_trait_tests.py`

## What was done

Per-trait differentiation between the marker-PCA k = 2 clusters: Mann-Whitney U rank-sum test + rank-biserial effect size + Cliff's delta. Tests whether the two genetically-distinct subgroups also differ phenotypically.

## Method

For each of the 13 traits, Mann-Whitney U test with `alternative="two-sided"` between Cluster 1 (the smaller, 11-line PC1-outlier subgroup) and Cluster 2 (the larger 84-line bulk). Rank-biserial correlation r_rb = 1 − 2U/(n1·n2). Cliff's delta = (gt − lt) / (n1·n2). Benjamini-Hochberg FDR across the 13 tests.

## Quick findings

| Trait | p_MW | q_BH | r_b |
|---|---|---|---|
| **Seed_Length** | **0.013** | 0.16 | +0.47 |
| Soluble_Oxalate | 0.10 | 0.61 | +0.43 |
| Total_Oxalate | 0.17 | 0.61 | +0.36 |
| Antioxidant | 0.25 | 0.61 | +0.23 |
| Tannin | 0.31 | 0.61 | −0.20 |
| Others | > 0.32 | > 0.61 | small |

- **No trait survives FDR < 0.10** after correction.
- The strongest nominal hit is Seed_Length (small cluster has shorter seeds; median 7.93 vs 8.61 mm).
- Conclusion: the marker-PCA clusters are genetically distinct (AMOVA Φ_ST = 0.27, p = 0.001) but **not strongly phenotypically differentiated** in this trait set.

## Caveats

- Cluster sizes are very unbalanced (11 vs 84) — power is limited.
- All Mann-Whitney p-values use the larger panel where the trait is scored; the protein/oxalate subset has only 6 lines in the small cluster (limited power).
- The result reinforces the caveat that the cluster labels are an internal marker-PCA partition, not a provenance / agroecological grouping. A meaningful breeder-relevant cluster split would likely need external metadata.
