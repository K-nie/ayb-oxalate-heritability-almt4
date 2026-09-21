# 26 — Locus zoom: ALMT4 + MATE-cluster LD haplotype blocks

**Script:** `scripts/26_mate_locus_zoom.py`

## What was done

Two locus-zoom figures around the paper-bearing oxalate-transport hits: (A) the Soluble_Oxalate × ALMT4 region on Ss10, and (B) the Insoluble_Oxalate × DETOXIFICATION-cluster region on Ss04. Each shows a Manhattan p-value strip, a Funannotate gene track, and a Haploview-style triangular LD heatmap.

## Method

The two focal SNPs are hardcoded (top-hit rs IDs from the AYB-anchored scan). For each, anchored markers within ±500 kb (Ss10) or ±1,000 kb (Ss04) of the focal SNP are pulled. Pairwise r² computed on the standardised dosage. Funannotate genes in the window flagged using the pathway-keyword regex (ALMT, MATE/DETOXIFICATION, MYB, AINTEGUMENTA, PME, expansin, etc.). The LD plot is a **classic Haploview-style triangle**: each upper-triangle SNP pair (i, j) is rendered as a diamond at (cx = midpoint, cy = half-physical-distance) coloured by r², with the triangle apex pointing up toward the gene track.

## Quick findings

- **ALMT4 locus (Ss10:15.4 Mb ± 500 kb)**: 16 anchored SNPs, 88 genes, 1 pathway-flagged (ALMT4 itself). LD block is narrow and centred on the focal SNP (high r² immediately adjacent, drops sharply by ±100 kb). Suggests the hit tags a single haplotype, not a large LD region.
- **MATE_DETOX locus (Ss04:67.8 Mb ± 1,000 kb)**: 32 anchored SNPs, 231 genes, 3 pathway-flagged. Sparser LD structure consistent with the wider window. The DETOXIFICATION-family cluster sits ~30 kb from the focal SNP.

## Caveats

- LD r² in a small panel (n = 95) is noisy; off-diagonal r² values < 0.2 are within sampling error of zero.
- The 1 Mb Ss04 window is wider than the typical ±500 kb because the AYB-anchored marker density on Ss04 is sparser around the MATE-cluster region; otherwise the focal SNP would have too few neighbours for a meaningful LD plot.
- The MATE-cluster focal SNP is rank-26 (not rank-1) in the n = 41 Insoluble_Oxalate scan — kept as the paper-bearing focal because it has the mechanistic match to the published cowpea MATE-cluster orthologue.
