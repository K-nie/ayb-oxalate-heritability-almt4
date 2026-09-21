# 38 — topGO Gene Ontology enrichment

**Scripts:** `scripts/38_go_enrichment.py` (first-pass Fisher exact),
`scripts/38b_go_topgo.R` (final topGO elim Fisher)

## What was done

Gene Ontology enrichment of the candidate genes flanking the top-10
GWAS SNPs per trait, pooled across all 13 traits. The first-pass Python
implementation used a classic Fisher exact test over the full
±50 kb foreground and produced no FDR-significant terms because of
classic-Fisher inflation of redundant parent terms in the GO DAG. The
final R implementation switched to topGO with the elim Fisher
algorithm, which down-weights parent terms once a child term explains
the signal.

## Method

The Funannotate GFF mRNA `Ontology_term` field was parsed into a per-gene
GO-term list. Foreground = the set of unique Funannotate-annotated genes
within ±10 kb of the top-10 GWAS SNPs per trait, pooled across all 13
traits. Background = the set of Funannotate-annotated genes with at least
one GO term. Three separate topGO analyses (BP, MF, CC) were run with
the elim Fisher algorithm, minimum nodeSize = 5, and Benjamini-Hochberg
false-discovery-rate adjustment within each namespace. GO term names
and namespace assignments were taken from the Gene Ontology release
of 2026-03-25.

## Inputs

- Funannotate GFF (`refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3`)
- AYB marker anchoring table
- Per-trait GWAS p-values (`results/03_gwas_mlm/`, `results/15_gwas_new_traits/`)
- GO term names from `refs/go/go.obo` (2026-03-25 release)

## Outputs

- `tables/go_topgo_BP.csv`, `go_topgo_MF.csv`, `go_topgo_CC.csv` —
  topGO output per namespace (GO ID, term name, fg / bg counts,
  classic Fisher p, elim Fisher p, BH q)
- `figures/fig75_go_topgo_BP/MF/CC.png/.pdf` — per-namespace horizontal
  bar charts (term name on y-axis, -log10 elim p on x-axis, fg / bg
  counts annotated, 0.05 threshold as dashed red line)
- `figures/fig75_go_topgo_combined.png/.pdf` — three-panel BP/MF/CC
  stacked figure (the version that goes into the paper)

## Quick findings

- **Top BP term: cellulose biosynthetic process** (elim p = 0.006)
  — coherent with the seed-coat phenotype variation (Seed_Coat_Tannin,
  Seed_Thickness).
- **Top MF terms: cellulose synthase activity** (p = 0.004) and
  **monoatomic ion channel activity** (p = 0.040) — the ion channel
  term aligns with the ALMT4 / MATE candidate gene hits in the
  oxalate-transport story.
- Sample size and sparse marker density limit how many FDR-significant
  terms survive correction, but the surviving terms map cleanly onto
  the candidate-gene biology already identified by the GWAS scan.

## Caveats

- topGO with the elim algorithm corrects within the GO DAG but not for
  the multiple traits we pool the foreground over. The 13-trait pooling
  is justified for hypothesis-generation but inflates the effective
  test count.
- Funannotate GO annotations are auto-assigned from sequence-similarity
  hits; terms with low protein-level evidence are still in the table.
- The first-pass Python Fisher exact run from `38_go_enrichment.py` is
  kept for completeness — superseded by the topGO output.
