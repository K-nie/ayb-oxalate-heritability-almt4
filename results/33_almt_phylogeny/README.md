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

---

# REVISION (TPG-2026-09-0323) — Rooting fix + outgroup re-inference

**Added 2026-10-20.** Addresses Reviewer 1 (L523–573) and Reviewer 2 (#21):
the original tree (above) was **midpoint-rooted** (`midpoint.root`,
`scripts/33_almt_phylogeny.R` line 33), i.e. effectively unrooted with an
arbitrary implied root, so the direction of ALMT4_2's placement rested on an
unjustified root. This stage re-infers the tree with a true outgroup and
re-roots by a model-free criterion, then tests whether the focal placement
survives.

## What was done & why

1. **Fetched a genuine ALMT-family outgroup** — three early-diverging
   land-plant ALMT proteins (Pfam **PF11744**), from UniProt:
   - `D8TB43` *Selaginella moellendorffii* (lycophyte)
   - `A0A7I4E3F8` *Physcomitrium patens* (moss)
   - `A0A2R6WEU1` *Marchantia polymorpha* (liverwort)
   These lineages diverged before the angiosperm ALMT radiation, so they are a
   defensible outgroup for an angiosperm ALMT tree.
2. **Re-aligned** the 142 ingroup + 3 outgroup (145 seqs) with MAFFT (`--auto`
   → chose L-INS-i) and **re-inferred** under the *published* model (LG+G4,
   1000 UFBoot) with IQ-TREE 3.0.1.
3. **Rooting.** IQ-TREE reported `WARNING: Branch separating outgroup is not
   found` — the three outgroup sequences are **not monophyletic**: *Selaginella*
   + *Marchantia* pair (patristic d = 2.04) but the *Physcomitrium* sequence is
   a long branch (d ≈ 4.1) that attaches separately. This is expected — the ALMT
   family underwent independent lineage-specific expansions, so early-diverging
   taxa carry their own divergent paralogs. All three outgroup tips nonetheless
   attach to **deep** branches far from the focal ALMT4 clade. Because a clean
   outgroup root is unavailable, the tree is rooted by **Minimal Ancestor
   Deviation (MAD; Tria, Landan & Dagan 2017, *Nat Ecol Evol* 1:0193)** — a
   model-free, outgroup-independent rooting (the reviewer-sanctioned fallback).
4. **Robustness.** The focal placement is compared across midpoint (original),
   MAD (original 142-tree), and MAD (re-inferred 145-tree).

## Exact commands

```bash
# outgroup sequences (UniProt REST)
for acc in D8TB43 A0A7I4E3F8 A0A2R6WEU1; do
  curl -s "https://rest.uniprot.org/uniprotkb/${acc}.fasta"; done > outgroup_raw.faa
cat data/all_almt_proteins.faa outgroup_raw.faa > data/reroot/almt_plus_outgroup.faa

# align + infer (published model)
mafft --auto --thread 4 data/reroot/almt_plus_outgroup.faa > data/reroot/almt_plus_outgroup_aln.faa
iqtree3 -s data/reroot/almt_plus_outgroup_aln.faa -m LG+G4 -B 1000 \
  -o "tr|D8TB43|D8TB43_SELML,tr|A0A7I4E3F8|A0A7I4E3F8_PHYPA,tr|A0A2R6WEU1|A0A2R6WEU1_MARPO" \
  --prefix data/reroot/ALMT_tree_og -T AUTO --seed 12345 -redo

# MAD root (model-free) — writes a rooted Newick, UFBoot labels preserved
python scripts/mad_write.py data/reroot/ALMT_tree_og.treefile      data/reroot/ALMT_tree_og_madroot.treefile
python scripts/mad_write.py data/ALMT_tree.treefile                data/reroot/ALMT_tree_orig_madroot.treefile

# figure (manuscript Fig 8), no baked title, MAD-rooted
Rscript scripts/33b_almt_phylogeny_rerooted.R
```

(MAD implementation `scripts/mad_root.py` was validated on toy trees: it
returns rms = 0 at the midpoint of a clocklike balanced tree and recovers a
known single-taxon outgroup — see commit notes.)

## Key quantitative results (ML re-inferred tree, logL = −50182.607)

| quantity | value |
|---|---|
| AYB ALMT4_Ss10 nearest tip (patristic) | **common bean V7BAP2, d = 0.165** |
| direct sister clade (UFBoot) | (cowpea A0A4D6N7X3, common bean V7BAP2), **UFBoot 100** |
| legume ALMT4 clade (7 tips incl. focal), UFBoot | **100** |
| patristic d, focal → **AtALMT4 (Q9C6L8)** | **1.024** (≈ 6× the distance to its legume sister) |
| patristic d, focal → **AtALMT9 (Q9LS46)** | **1.295** |
| smallest clade {focal + AtALMT4} | 22 tips (tonoplast ALMT4/5/6 subfamily), UFBoot 100 |
| smallest clade {focal + AtALMT9} | 48 tips, UFBoot 100 |
| MAD root split / ambiguity index | 53 \| 92 leaves, amb 0.827 (orig 142-tree: 50 \| 92, amb 0.845) |

**Rooting robustness:** under midpoint (original), MAD-142, and MAD-145, AYB
ALMT4_Ss10 always falls in the tonoplast ALMT3/4/5/6/9 assemblage with the same
legume sisters. The orthology conclusion (below) does **not** depend on the
rooting — the key facts (nearest relative, branch depth, support) are
rooting-independent.

## Outputs (this revision)

- `data/reroot/almt_plus_outgroup_aln.faa` — 145-seq MAFFT (L-INS-i) alignment
- `data/reroot/ALMT_tree_og.{treefile,contree,iqtree,log}` — re-inferred tree (LG+G4, 1000 UFBoot)
- `data/reroot/ALMT_tree_og_madroot.treefile` — re-inferred tree, **MAD-rooted** (primary)
- `data/reroot/ALMT_tree_orig_madroot.treefile` — original 142-tree, MAD-rooted (comparison)
- `figures/fig72_almt_family_tree_MADrooted.{png,pdf}` — Fig 8, MAD-rooted, no baked title

## Tool versions (pinned)

MAFFT v7.526 · IQ-TREE 3.0.1 (LG+G4, `-B 1000`, seed 12345) · DendroPy 5.0.8 ·
Biopython 1.81 · R 4.4 (ape, ggtree, treeio) · MAD via `scripts/mad_root.py`
(own implementation of Tria et al. 2017).

## Caveats (this revision)

- The added outgroup is **not monophyletic**; rooting is therefore by MAD, not
  by the outgroup. The outgroup still demonstrates that the focal clade is
  deeply nested, not basal.
- `--auto` selected L-INS-i (145 seqs); the original run used the same `--auto`
  switch (142 seqs). Alignment is long (1681 cols) and gappy (IQ-TREE flags
  >50% gaps on the divergent members) — intrinsic to this gene family's
  variable C-terminal cytosolic domain.
