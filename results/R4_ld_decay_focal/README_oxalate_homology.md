# Oxalate-enzyme exclusion in the ALMT4_2 candidate window

**Question:** does any of the 62 genes in the ±391 kb LD-extinction window around
the focal SNP encode an oxalate-metabolism enzyme (an alternative candidate to ALMT4_2)?

**Why beyond product names:** the Funannotate reference annotation leaves 50/62 (81%)
of the window genes as "hypothetical protein", and several named genes carry
fungal-ortholog Swiss-Prot names (top-hit artefacts). The annotation DOES resolve
plant-specific genes elsewhere (e.g. chloroplastic proteins), so it is not fungal-
lineage mis-annotation, but the exclusion claim is verified at three levels so it
does not rest on product names alone.

**Methods / commands** (`scripts/_window_oxalate_homology.py`, 2026-10-27):
1. Product name — from the Funannotate GFF (`product=`).
2. Pfam/InterPro domain — from the GFF `Dbxref=` (Funannotate ran InterProScan =
   HMMER profile search); scanned for PF01072 (germin/oxalate oxidase), PF07883/
   PF00190 (cupin / oxalate decarboxylase), PF00501 (AMP-binding / AAE3), PF01554
   (MATE), PF01565 (glycolate oxidase), PF11744 (ALMT).
3. Direct BLASTP (BLAST+, -evalue 100 -max_target_seqs 50) of all 61 protein-
   coding window proteins (the 62nd window gene, AYBTSS11_029708, is a non-coding
   tRNA-Cys with no protein) (`refs/ayb_genome/
   AYB_proteins.faa`) against three reference enzymes fetched from UniProt —
   Arabidopsis AAE3 oxalyl-CoA synthetase (Q9SMT7), wheat oxalate oxidase / germin
   (P15290), Bacillus oxalate decarboxylase (O34714) — with AYB ALMT4_2 as a
   positive control.

**Result** (`tables/window_oxalate_homology.csv`):
- Pfam/InterPro: the only oxalate-relevant domain in the window is ALMT (PF11744 +
  FUSC PF13515) in ALMT4_2; **0** window genes carry a non-ALMT oxalate-enzyme domain.
- BLASTP (max_target_seqs 50, evalue 100): best hit per enzyme — oxalate
  decarboxylase e = 0.003 (46-aa cupin scrap), AAE3 e = 0.068, oxalate oxidase
  e = 0.180; smallest e anywhere 0.003, **none significant** (threshold 1e-5). Positive control
  POSCTRL_AYB_ALMT4_2 recovers ALMT4_2 at 100% id, e = 0.
- Conclusion: no oxalate-metabolism enzyme in the window by name, domain, or
  homology; ALMT4_2 is the only oxalate-relevant (organic-acid transporter) gene.

**Annotation provenance caveat:** the reference annotation is the Funannotate v1.8
models distributed with the Waweru et al. (2024) genome (Palmer & Stajich 2020).
The exact Funannotate run settings (evidence, BUSCO lineage) are a property of that
genome release; the exclusion above is robust to them because it uses domain and
sequence homology, not product names. Confirm annotation provenance with the
genome authors before submission.

## HMMER confirmation (added 2026-10-27, `scripts/_window_oxalate_hmmer.py`)
Independent of InterProScan: `hmmsearch --cut_ga` of the oxalate-enzyme Pfam HMM
profiles against the 61 window proteins, using each family's curated GATHERING
threshold (a hit above --cut_ga IS a family member). Output
`tables/window_oxalate_hmmer.csv`.
- PF00190 Cupin_1 (germin/oxalate oxidase): none above threshold.
- PF07883 Cupin_2 (oxalate decarboxylase): none.
- PF00501 / PF13193 AMP-binding (AAE3 oxalyl-CoA synthetase): none.
- PF01554 MatE (MATE): none.
- PF11744 ALMT (positive control): AYBTSS11_029726 = ALMT4_2, e = 1.6e-175 (passes).
- phmmer of AAE3 / oxalate oxidase / oxalate decarboxylase vs the window proteins: no hit.
Conclusion stands across four methods (name, domain, BLASTP, HMMER): no oxalate-
metabolism enzyme in the window; ALMT4_2 is the only oxalate-relevant gene.

## Whole-proteome positive control (reviewer request, added 2026-10-27)

To show that the "no oxalate enzyme in the candidate window" result is a true local
absence and not a search that cannot detect these enzymes, the same BLASTP and HMMER
searches were run against the WHOLE AYB proteome (30,840 proteins).

- Script: `scripts/_proteome_oxalate_positive_control.py`
- Output: `tables/proteome_oxalate_positive_control.csv`
- Inputs: `refs/ayb_genome/AYB_proteins.faa`; reference enzymes Q9SMT7 / P15290 / O34714 from UniProt; Pfam profiles extracted from the local `~/.myconote/dbs/pfam/Pfam-A.hmm` (e.g. PF11744.15, PF00501.35).
- Settings identical to the window run: BLASTP significance e < 1e-5; hmmsearch at the Pfam gathering threshold (--cut_ga).

Result — the methods detect oxalate-metabolism genes/families across the proteome:

| search | query | proteome hits | note |
|---|---|---|---|
| BLASTP | AAE3 oxalyl-CoA synthetase (Q9SMT7) | 47 (e<1e-5) | best AYBTSS11_015480, e=0.0 (76% id) |
| BLASTP | oxalate oxidase / germin (P15290) | 27 (e<1e-5) | best AYBTSS11_002440, e=2.5e-62 |
| BLASTP | oxalate decarboxylase (O34714, *Bacillus*) | 0 (e<1e-5) | best e=5.8e-5; bacterial enzyme, no true plant ortholog |
| hmmsearch | PF00190 Cupin_1 (germin/oxalate oxidase) | 37 | |
| hmmsearch | PF07883 Cupin_2 (oxalate decarboxylase) | 32 | |
| hmmsearch | PF00501 AMP-binding (AAE3) | 53 | |
| hmmsearch | PF13193 AMP-binding_C (AAE3) | 32 | |
| hmmsearch | PF01554 MatE (MATE) | 61 | |
| hmmsearch | PF11744 ALMT (internal control) | 15 | **ALMT4_2 recovered** |

The AYB proteome carries dozens of members of every oxalate-enzyme family searched (47 AAE3 and 27 germin BLASTP orthologs; counts uncapped, max_target_seqs 100000), and
the ALMT control recovers ALMT4_2, so the search is sensitive. The 62-gene candidate
window containing none of them is therefore a genuine local absence. (The bacterial
oxalate decarboxylase has no plant BLASTP ortholog, as expected — plants degrade oxalate
by other routes — but its Cupin families PF00190/PF07883 are both present proteome-wide.)
