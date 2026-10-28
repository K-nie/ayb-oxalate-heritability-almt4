# R4 — Genome-wide LD decay and the data-driven focal window (revision Phase 4)

## What and why
Reviewer L277-279, L498-499: report the genome-wide LD half-decay distance, provide a clean
LD-decay figure, and set the candidate window around the focal SNP from the **data** (the LD
block), not an assumed ±50 kb. LD decay was already computed
(`results/91_ld_decay_genomewide_ayb`, `results/96_ld_decay_dedup_sscorr`); here we re-fit
the Hill-Weir curve (cross-check), draw the figure, and list every annotated gene inside the
LD-based window around Ss10:15,394,673.

## Inputs
- Binned r² + fitted summary: `results/96_ld_decay_dedup_sscorr/tables/`
  (dedup n=93 panel, Hill-Weir 1988 / Remington 2001 sample-size correction) — primary;
  `results/91_ld_decay_genomewide_ayb` (full n=95) — sensitivity.
- Per-chromosome summary: `.../ld_per_chrom_dedup_sscorr_summary.csv` (Ss10 radii).
- Annotation: `refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3`
  (Funannotate; Zenodo 10.5281/zenodo.13853757).

## Method / parameters
- Re-fit expected-r² decay `E[r²](d)` (Hill & Weir 1988 with finite-sample term, n = 93)
  to the binned means by `scipy.optimize.curve_fit`; **half-decay = distance where the fitted
  r² falls to half of its extrapolated d→0 value** (the project's / Remington 2001 convention).
- Focal-window radii taken from the **data**: Ss10 half-decay (99.1 kb) and Ss10 distance to
  r²<0.1 (391.3 kb). Genes within ±radius of the focal SNP extracted from the GFF.

## Reproduce
```
cd scripts
python R4_ld_decay_focal.py     # env shap-repro
```

## Outputs (`results/R4_ld_decay_focal/`)
- `tables/ld_halfdecay_summary.csv` — refit vs stored half-decay / r²<0.1; Ss10 radii.
- `tables/focal_window_genes_half_decay.csv` (±99.1 kb), `..._r2_below_0.1.csv` (±391.3 kb).
- `figures/ld_decay_genomewide.{png,pdf}` (no title on image; caption carries it).

## Key results
- **Genome-wide half-decay = 83.7 kb** (refit) vs 83.6 kb (stored, dedup n=93) — reproduced;
  full n=95 panel stores 102.3 kb. Distance to r²<0.1 = 330.3 kb (refit) = 330.3 kb (stored).
  Background (asymptotic) r² ≈ 0.04.
- **Ss10 (the focal chromosome): half-decay 99.1 kb, r²<0.1 at 391.3 kb.**
- **The LD block around the focal SNP contains ALMT4_2 plus 16 other genes.** Within
  ±99.1 kb (half-decay) there are **17 genes**; ALMT4_2 (`AYBTSS11_029726`, 3,174 bp from the
  SNP) is the 2nd-nearest. Within ±391.3 kb (r²<0.1) there are **62 genes**. ALMT4 is a
  candidate within the block, not the only gene in it — the candidate window must be read at
  LD scale (~100 kb), not at an assumed ±50 kb.

## Caveats
- Half-decay is population-average; LD varies ~2× across chromosomes (Ss07 42 kb → Ss01 113 kb,
  per `results/96_.../ld_per_chrom_dedup_sscorr_summary.csv`).
- The gene list is positional (annotation only); it does not assert function. The 17-gene
  half-decay window is the honest unit for "the candidate region," and supports
  "a candidate SNP near ALMT4," not "ALMT4 is the gene."

## Pinned versions
Python 3.9 (`shap-repro`): numpy 1.26.4, pandas 2.3.3, scipy 1.13.1, matplotlib.
