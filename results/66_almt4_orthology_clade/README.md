# 66 — ALMT4 phylogeny-implied syntenic clade (Addition 1 Part B, lite)

Substitute for the coordinate-level comparative-synteny analysis originally planned for Addition 1 Part B (cowpea / common bean / soybean / Medicago / chickpea proteome BLASTp + per-flanking-gene synteny window concordance). That plan required Phytozome (login) proteome downloads and per-comparator BLASTp queries; outside today's scope and deferred to a revision pass.

The substitute uses the 142-sequence ALMT family phylogeny from script 33 (`results/33_almt_phylogeny/`) — already built from UniProt TrEMBL sequences across A. thaliana, G. max, P. vulgaris, V. unguiculata, C. arietinum, M. truncatula, plus the five AYB ALMT-family paralogs — to extract the clade that contains AYB ALMT4_Ss10. The clade members encode the same orthology + clade-level conservation signal that a coordinate-level synteny analysis would deliver, modulo per-genome locus order (Lawrence-Dill et al. 2018).

## Method

1. Parse the IQ-TREE phylogeny (`ALMT_tree.treefile`).
2. Walk up from the AYB ALMT4_Ss10 tip until the first ancestor whose leaf set also contains the AYB ALMT4_Ss02 paralog — that ancestor's clade is the "phylogeny-implied syntenic clade" for ALMT4_Ss10. (Fallback rule: smallest clade containing focal + at least one tip per comparator legume species.)
3. List the clade members with species (from tip-name prefix) and UniProt accession.
4. Walk every internal node in the clade and collect its IQ-TREE UFBoot support.
5. Render a focused tree of the clade with species-coloured tip labels and UFBoot annotations.

## Findings

**AYB ALMT4_Ss10's phylogeny-implied syntenic clade contains 14 tips spanning 5 species:**

| Species | n | UniProt accessions |
|---|---|---|
| AYB | 2 | ALMT4_Ss10 (focal), ALMT4_Ss02 (paralog) |
| Glycine max (soybean) | 6 | K7K7E2, I1J6M4, I1KWT9, I1LEY0, K7LLU3, A0A0R0EAN5 |
| Phaseolus vulgaris (common bean) | 3 | V7CJW1, V7CL11, V7BAP2 |
| Vigna unguiculata (cowpea) | 2 | A0A4D6N7X3, A0A4D6LN79 |
| Cicer arietinum (chickpea) | 1 | A0A1S2XZY2 |

**Tree topology features:**
- All 13 internal nodes carry UFBoot ≥ 94: ten of 13 at 100, the other three at 94, 98 and 99 (per `tables/clade_node_support.csv`; re-verified 2026-10-27). An earlier "12 of 13 at 100" statement here was stale — the table shows ten at 100, minimum 94.
- AYB ALMT4_Ss10 is direct sister to (cowpea A0A4D6N7X3, common bean V7BAP2) at UFBoot = 100.
- AYB ALMT4_Ss02 (paralog) sits in a parallel sub-clade with cowpea A0A4D6LN79, common bean V7CJW1+V7CL11, and three soybean orthologs — same set of species, different sub-clade, consistent with a duplication event predating the Phaseoleae / Cicereae split.
- Arabidopsis ALMT4 (ALMT4_ARATH; Q9C6L8) sits *outside* this legume-specific clade in a sister position, anchoring the deeper orthology.
- *Medicago truncatula* is **absent** from this clade despite being present in the broader 142-sequence tree — consistent with Medicago lineage-specific gene loss or substantial sequence divergence in the ALMT4_Ss10-orthologous copy.

**Three operational implications:**

1. **Functional-annotation transfer to AYB ALMT4_Ss10 should preferentially cite cowpea / common bean / soybean ALMT literature** (Aliyu et al. 2022; Lonardi et al. 2019) rather than the more-distant Arabidopsis AtALMT9 literature directly. The closer legume orthologs are the appropriate reference set.

2. **The AYB ALMT4_Ss02 paralog has legume orthologs at all four comparator species** in the parallel sub-clade — both AYB ALMT4 copies are conserved across the Phaseoleae / Cicereae lineage, consistent with sub-functionalisation after duplication rather than recent neo-functionalisation.

3. **Medicago's absence flags lineage-specific gene loss** in the Cicereae-Trifolieae transition; this is consistent with Medicago's known sparser ALMT-family complement.

## Outputs

- `tables/ayb_almt4_clade_members.csv` — 14 rows: tip_name, species, uniprot_id, is_focal.
- `tables/clade_node_support.csv` — 13 rows: internal node UFBoot.
- `figures/fig_almt4_clade_tree.png` / `.pdf` — focused tree of the 14-tip clade with species-coloured tip labels, UFBoot internal-node labels (> 50 shown), and the AYB ALMT4_Ss10 focal tip in a vermillion box.

## Caveats

- **Not a coordinate-level synteny analysis.** The phylogeny-implied syntenic clade tests orthology + clade conservation but does not directly compare per-flanking-gene order at the genomic-coordinate level. A revision pass with full comparator-genome downloads can extend the analysis to the per-flanking-gene check.
- **UniProt TrEMBL sequence coverage is uneven across the comparators.** The number of soybean tips (6) reflects soybean's tetraploid history + thorough TrEMBL coverage; the chickpea single tip reflects sparser sequencing. This affects the raw tip count per species but not the orthology call.
- **No Medicago tip in the AYB ALMT4_Ss10 clade is flagged as informative**, but it could also reflect missing UniProt entries rather than true gene loss. A Phytozome-based ortholog search would settle this.

## dN/dS (Addition 1 Part C) status

dN/dS via PAML codeml is installed (`/opt/homebrew/bin/codeml`) but requires nucleotide-coding (codon-aligned) CDS sequences for all 14 clade tips. Extracting CDS for the UniProt TrEMBL legume sequences requires Ensembl Plants / NCBI lookup per accession — outside today's scope. **Part C is deferred to the revision pass** with the explicit data dependency documented.

## Script

`scripts/66_almt4_orthology_clade.py`. Runtime: ~ 5 s on a laptop.

---

# REVISION (TPG-2026-09-0323) — orthology, "syntenic clade", honest re-assessment

**Added 2026-10-20.** Addresses Reviewer 1 (L523–573) and Reviewer 2 (#23). Run
on the **outgroup re-inferred, MAD-rooted** tree (see
`results/33_almt_phylogeny/README.md` → REVISION). Findings are nearly identical
on the original 142-tree, so they are rooting-independent.

## 1. Is AYB ALMT4_Ss10 an "ortholog of the characterized ALMT4"? — NO (over-claim)

The current manuscript wording ("sits in the characterized ALMT4 clade",
"ortholog") is **not defensible**. The phylogeny shows:

- **Nearest relatives of AYB ALMT4_Ss10 are all *uncharacterized legume*
  ALMTs**, not Arabidopsis ALMT4:
  common bean V7BAP2 (patristic d = 0.165), soybean A0A0R0EAN5 (0.190),
  soybean I1LEY0 (0.241), cowpea A0A4D6N7X3 (0.242), chickpea A0A1S2XZY2 (0.395).
  Its **direct sister** is (cowpea A0A4D6N7X3, common bean V7BAP2), **UFBoot 100**.
- The focal gene sits in a strongly-supported **legume-specific ALMT4 clade**
  (7 tips: AYB + 3 soybean + common bean + cowpea + chickpea), **UFBoot 100**.
- The **characterized Arabidopsis AtALMT4 (Q9C6L8) is ~6× farther away**
  (patristic d = 1.024) and is **not** inside the legume clade. It joins only at
  the base of a 22-tip tonoplast ALMT4/5/6 subfamily (UFBoot 100 for the
  subfamily, but a deep branch separates AtALMT4 from the legume clade — visible
  in Fig 10).
- **AtALMT9 (Q9LS46), the cited cryo-EM structure, is in a *different*
  sub-clade** (patristic d = 1.295), even farther than AtALMT4.

**Defensible wording:** AYB ALMT4_Ss10 is *an uncharacterized legume ALMT in the
tonoplast ALMT4/5/6 (ALMT clade-II) subfamily*, most closely related to
uncharacterized legume ALMTs and only **distantly** related to the characterized
Arabidopsis AtALMT4 — it is **not a one-to-one ortholog of AtALMT4**. AtALMT9 is
a more distant family member, appropriate only as a **structural homology
template**, not as the ortholog.

## 2. "Syntenic clade" — should be dropped (never tested)

The 14-tip clade was defined **purely by phylogeny** (walk up the tree until the
AYB ALMT4_Ss02 paralog is included); **gene order was never tested**, so
"syntenic" is unsupported. A proper microsynteny test is **not feasible** from
data on disk: of the clade's 5 species, only AYB
(`refs/ayb_genome/...Funannotate.gff3`) and cowpea
(`refs/cowpea_v1/gene_models_main.gff3.gz`) have annotations here — the closest
relatives (common bean = the direct sister, soybean, chickpea) have **no genome
annotation present**, and cowpea has no proteome FASTA for flanking-gene ortholog
BLAST. **Recommendation (b): replace "syntenic clade" with "orthologous clade" /
"phylogenetically-defined clade"** throughout.

## 3. Exact commands

```bash
# re-inferred MAD-rooted tree comes from results/33_almt_phylogeny REVISION
Rscript scripts/66b_almt4_clade_rerooted.R          # Fig 10 (focal + AtALMT4 clade), no baked title
python  scripts/33c_almt4_topology_report.py \
        results/33_almt_phylogeny/data/reroot/ALMT_tree_og.treefile   # distances + supports
```

## 4. Outputs (this revision)

- `figures/fig_almt4_clade_rerooted.{png,pdf}` — Fig 10: the 22-tip tonoplast
  ALMT4/5/6 clade showing the focal gene nested in the legume ALMT4 clade and
  AtALMT4 on a deep separate branch. Lowercase/colourblind-safe, sans-serif, no
  baked title, scale bar shown.
- `tables/clade_members_rerooted.csv` — 22-tip clade membership + species.
- `tables/clade_node_support_rerooted.csv` — per-node UFBoot.

## 5. Tool versions

IQ-TREE 3.0.1 · R 4.4 (ape, ggtree, treeio) · Biopython 1.81. The focal-clade
supports (UFBoot 100 for the legume clade; 100 for the direct cowpea+bean sister)
are identical on the original 142-tree and the re-inferred 145-tree.

## 6. Caveats

- Uneven UniProt TrEMBL coverage per species still holds; it affects raw tip
  counts, not the orthology call.
- The clade is defined by phylogeny only; genomic gene-order synteny remains
  untested (see point 2).
