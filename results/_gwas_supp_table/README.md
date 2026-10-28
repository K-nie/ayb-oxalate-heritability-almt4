# Supplementary Table S18 — GWAS per-SNP per-trait results (The Plant Genome GWAS template)

`TableS18_gwas_per_snp_per_trait.csv` — the journal-required per-SNP results
table for the GWAS, long format: **one row per SNP per trait per model**. Built
2026-10-28 by `scripts/_build_gwas_supp_table.py`.

## What the journal asks for (and where it is here)
For every SNP tested, TPG's GWAS template requires: the allele scored, a position on a
reference (or the adjacent sequence if unplaced), and the p-value. All three are columns
below; unplaced markers carry their 69-bp DArT tag in place of a position.

## Scope
- **19,525 rows**, 13 traits, 1,625 unique markers.
- Rows per trait = that trait's QC/MAF-passing marker count (biochem 1,625; seed 1,591;
  crude protein 1,585; oxalate fractions 1,161; +2 focal rows for soluble oxalate).
- `in_genome_wide_scan = TRUE` for 19,523 scan rows; `FALSE` for the 2 focal direct-test rows.

## Columns
| column | meaning |
|---|---|
| trait | phenotype tested |
| marker_id | DArTseq marker ID (encodes the SNP: `…31:T>C…` = T→C at tag position 31) |
| chromosome, position | *S. stenocarpa* reference locus (re-anchoring, `marker_placement_classification.csv`); blank for unplaced markers |
| ref_allele, alt_allele | the two alleles (from `marker_qc.csv` `alleles` = ref/alt, and the marker ID) |
| allele_scored | **the alternate allele** — the one the 0/1/2 dosage counts (confirmed in `01_qc_pca_power.py` `dosage_row`: hom-ref 0, het 1, hom-alt 2) |
| maf, n | minor-allele frequency and analysed sample size for the trait |
| beta, se, p | per-allele effect, SE and p from the reported model |
| model | `G` (kinship only) or `G + PC` (kinship + 3 PCs) |
| in_genome_wide_scan | TRUE = recovered from the scan; FALSE = focal SNP tested directly |
| tag_sequence | 69-bp DArT tag, given **only** for unplaced markers (adjacent-sequence fallback) |

## Allele-coding rule (for the legend)
DArTseq IDs encode the SNP inside the ID (`31:T>C` = T→C at tag position 31). The 0/1/2
dosages count the **alternate** allele, so `allele_scored` is the alt allele — `C` for the
focal SNP, matching the TT/TC/CC classes in Figure 5.

## Model reported per trait
- **Antioxidant, Flavonoid, Phenol, Tannin** → `G + PC` (`results/03_gwas_mlm`, M2_KQ) —
  the model the paper reports for those scans.
- **the other nine traits** (5 seed-morphometrics, crude protein, 3 oxalate fractions) →
  `G` (`results/15_gwas_new_traits`, M1_K) — the only model run for them.
- **Focal SNP** (soluble oxalate): a `G` row (β 98.65, SE 30.85, p = 0.00139) and a
  `G + PC` row (p = 0.009), both `in_genome_wide_scan = FALSE` (it fell below the 0.10
  subset MAF floor and was tested directly, not recovered from the scan). It still appears
  as an ordinary scan row for the 11 full-panel traits, where it passes the MAF filter.

## Known issues / caveats (flagged for Benjamin)
1. **Model inconsistency.** Only the 4 biochem traits have a G + PC scan table; the other
   9 have G-only. The manuscript (§2.6) states scans are reported under G + PC — so either
   that statement is specific to the biochem traits, or G + PC scans should be re-run for
   the remaining 9 traits before this table ships.
2. **Focal SNP position.** The re-anchoring places the focal SNP at **Ss10:15,394,609**;
   the manuscript text cites **Ss10:15,394,673** (a 64-bp difference). The table uses the
   anchoring value; the two should be reconciled.

## Deposit
Too large for a Word table — ships as **Supplementary Table S18** (CSV), uploaded with the
supplement and deposited in the GitHub repo and Zenodo alongside the other analysis
outputs.
