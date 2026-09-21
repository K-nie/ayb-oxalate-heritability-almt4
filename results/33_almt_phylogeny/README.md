# 33 — ALMT family phylogeny

**Script:** `scripts/33_almt_phylogeny.R`

## What was done

Maximum-likelihood phylogeny of the aluminum-activated malate transporter
(ALMT) family across five legumes and *Arabidopsis*, anchoring the five
AYB ALMT sequences from the Funannotate annotation in their cross-species
ortholog groups.

## Method

The five AYB ALMT protein sequences (Ss02 ALMT12, Ss02 ALMT4, Ss03 ALMT9,
Ss04 ALMT9, Ss10 ALMT4) were combined with 137 ALMT-family proteins from
UniProt across *Arabidopsis thaliana*, *Glycine max*, *Phaseolus vulgaris*,
*Vigna unguiculata*, and *Medicago truncatula* (query: protein name
"aluminum-activated malate transporter"). The 142-sequence set was
aligned with MAFFT v7.515 (Katoh & Standley 2013) under `--auto`. A
maximum-likelihood tree was inferred with IQ-TREE 3.0.1 (Minh et al.
2020) under LG+G4 with 1,000 ultrafast bootstrap replicates. The tree
was rendered in `ape::plot.phylo` as a circular fan (14 × 14 inch PDF),
with AYB tips highlighted by red-filled circles via `tiplabels()` and
internal nodes carrying bootstrap support ≥ 70 % via `nodelabels()`.

## Outputs

- `tables/almt_alignment.fasta` — 142-sequence MAFFT alignment
- `tables/almt_tree.treefile` — IQ-TREE Newick with bootstrap support
- `figures/fig72_almt_family_tree.png/.pdf` — circular fan tree

## Quick findings

- All five AYB ALMTs cluster within established legume ALMT4/9/12 clades
  rather than in *Arabidopsis*-only clades, consistent with the
  legume-lineage origin of these gene copies.
- The Ss10:15.4 Mb ALMT4 candidate from the Soluble_Oxalate GWAS hit
  sits in the ALMT4 clade with common bean and soybean orthologs,
  supporting the functional annotation.

## Caveats

- 137 UniProt ALMT sequences include duplicates and isoforms; the tree
  reflects sequence-level clustering, not species-level orthology.
- LG+G4 was selected over more elaborate ModelFinder candidates for
  reproducibility; the topology around the AYB nodes is unchanged under
  WAG+G4 and JTT+G4 (spot-checked).
