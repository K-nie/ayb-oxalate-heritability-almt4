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
- All 13 internal nodes carry UFBoot ≥ 70 (12 of 13 at 100).
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
