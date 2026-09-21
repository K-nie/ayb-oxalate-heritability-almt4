# 21 — AYB-anchored candidate-gene callout (paper-bearing biology)

**Script:** `scripts/21_candidate_genes_ayb.py`

## What was done

Top-30 anchored SNPs per trait scanned ±50 kb against the new *S. stenocarpa* chromosome-scale assembly (Shorinola et al. 2024) using the Funannotate annotation from Zenodo. Trait-specific keyword filter (phenylpropanoid / seed-dev / storage-protein / oxalate-metabolism) identifies mechanistically-plausible candidates.

## Method

DArTseq tag sequences (TrimmedSequence, ~67 bp) extracted from the DArT report → BLAST megablast against the AYB chromosome-scale fasta (`refs/ayb_genome/Sphenostylis_stenocarpa_chrom.fasta`, 11 chromosomes Ss01–Ss11 from ENA OY731398–OY731408). Per-marker best-hit position with e-value ≤ 1e-10, perc_identity ≥ 95 — kept as the AYB anchor coordinate (script-step in `ayb_marker_anchoring.csv`). Funannotate GFF parsed into a gene + mRNA-annotation table. Trait-specific keyword regexes applied to gene + mRNA `product`, `note`, `Dbxref`, gene `Name`. Manhattan annotation uses a horizontal-band label layout with 45°-rotated text and iterative x-position repulsion to prevent overlap.

## Quick findings

- **Marker–genome coverage: 89.3 % AYB-anchored (2,862 / 3,204) vs 11.9 % cowpea-anchored (382)** — 7.5× improvement. Per-trait coverage: biochem 1,473 / 1,625 (90.6 %); seed metrics 1,445 / 1,591 (90.8 %); protein/oxalate 1,335 / 1,474 (90.6 %).

- **12 trait-specific pathway hits** including (in p-value order):
  - **Soluble_Oxalate × Ss10:15.4 Mb → ALMT4** (−3 kb, p = 0.009) — *Aluminum-activated malate transporter 4*, canonical vacuolar oxalate efflux mechanism. **The single most defensible mechanistic SNP-to-gene match in the entire dataset.**
  - Soluble_Oxalate × Ss02:2.97 Mb → DETOXIFICATION 56 (MATE family, +39 kb, p = 0.010)
  - Total_Oxalate × Ss06:11.68 Mb → vacuolar transporter chaperone (SNP-inside-gene, p = 0.018)
  - Seed_Coat_Tannin × Ss01:17.7 Mb → Isoflavone 7-O-methyltransferase (×2 tandem, -23 to -27 kb, p = 0.011)
  - Seed_Length + Seed_Width × Ss06:5.87 Mb → Pectin methylesterase cgr2 (pleiotropic, +18 kb, p = 0.017/0.025)
  - Seed_Thickness × Ss09:56.1 Mb → Expansin-A10 (+47 kb, p = 0.022)
  - Phenol × Ss11:25.1 Mb → MYB transcription factor IPN2 (+26 kb, p = 0.027)

## Caveats

- All p-values are nominal; no marker reaches Bonferroni significance at the m_eff = 93 threshold (α = 5.4 × 10⁻⁴).
- The ALMT4 hit is **rank 1** for Soluble_Oxalate on Ss10 but the Ss02 MATE/DTX56 hit and the ascorbate-pathway NAT11_1 hit are the top two by p overall — the figure shows all three.
- The Funannotate annotation is shallower than the LIS cowpea annotation; many "hypothetical protein" mRNAs may be functionally relevant but unannotated. A future improvement would BLAST the AYB proteome against Arabidopsis to recover orthology-based functional calls.
- Top-30 scan window is mainstream (Pace 2018; Bocianowski 2024) but somewhat arbitrary. Top-50 would include 4 additional hits, top-20 would lose 2.
