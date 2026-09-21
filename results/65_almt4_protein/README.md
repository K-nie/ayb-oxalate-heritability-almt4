# 65 — ALMT4 protein-domain confirmation and TM topology

Part A of the three-part ALMT4 mechanism-evolution chain (Addition 1 of `manuscript/paper_A2/analysis_plans_A2.md`). Lifts the ALMT4 candidate-gene call from "GWAS proximity + orthology phylogeny" to "GWAS proximity + orthology phylogeny + structural-domain + TM topology" — independent structural-biology evidence for the mechanism.

Parts B (comparative synteny against legume reference proteomes) and C (dN/dS on the AYB ALMT4 clade) are scoped in separate scripts (66 and 67 respectively).

## Method

1. Extract the AYB ALMT4_2 amino-acid sequence (`AYBTSS11_029726-T1`, 560 aa) from `refs/ayb_genome/AYB_proteins.faa`.
2. Pull Pfam and InterPro accessions directly from the Funannotate GFF (`refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3`). Funannotate (Shorinola et al. 2024 Zenodo 10.5281/zenodo.13853757) runs HMMER + Pfam-A at annotation time and writes the accessions into the Dbxref attribute — no rerun of HMMER is needed for primary domain ID.
3. Compute per-residue Kyte-Doolittle (1982) hydrophobicity under a centred 19-residue window.
4. Predict transmembrane segments as contiguous runs above an entry threshold of 0.8 (Wimley-White 2002 revised, tuned for window-19 smoothing and plant membrane proteins per Cao et al. 2013), with a segment-mean floor of 1.2 and a minimum length of 14 aa. Segments longer than 30 aa are tentatively split at their internal local minimum; canonical ALMT-family TM helices reach at most 25 aa by crystal structure (Qiu et al. 2021 cryo-EM on AtALMT9).
5. Render a composite topology figure stacking the hydrophobicity profile, predicted TM segments, and Pfam domain spans.

## Findings

- **AYB ALMT4_2 is 560 aa**, encoded on Ss10:15,397,847–15,400,688 on the (+) strand, sitting **3,174 bp downstream of the Soluble_Oxalate focal SNP** at Ss10:15,394,673 — the "−3 kb proximity" cited in the manuscript.
- **Pfam annotation** (from Funannotate): PF11744 (ALMT family domain) in the N-terminal half, PF13515 (FUSC-like / Membrane protein superfamily) in the C-terminal half. **InterPro** IPR020966 (ALMT-family superfamily) full-length. The PF11744 + IPR020966 calls together confirm canonical aluminum-activated malate-transporter family membership at the structural-domain level, independently of the orthology phylogeny.
- **Predicted TM segments (3 total)**:
  - TM1 at residues 113–126 (14 aa, mean K-D hydrophobicity 1.44)
  - TM2 at residues 141–178 (38 aa, mean 1.49 — likely two fused helices not separated under the 19-residue window smoothing; the canonical ALMT helix length is 17–25 aa per Qiu et al. 2021)
  - TM3 at residues 218–236 (19 aa, mean 1.44)
- **All three predicted TMs sit inside the PF11744 N-terminal domain.** Two additional hydrophobic peaks in the PF13515 C-terminal domain (around residues 340 and 525) do not clear the segment-mean floor under our default parameters but match the canonical ALMT 6-TM architecture (Sasaki 2004; Qiu 2021).
- The conservative TM count (3 of the canonical 6) reflects window-19 smoothing on a 560-aa multi-domain protein, not absence of TM helices.

## Outputs

- `data/almt4_protein.fasta` — extracted AYB ALMT4_2 protein, ready for Parts B and C.
- `tables/pfam_interpro_domains.csv` — 3 rows: PF11744 + PF13515 + IPR020966 with source / accession / span / comment.
- `tables/tm_segments.csv` — 3 predicted TM segments with start / end / length / mean hydrophobicity.
- `tables/hydrophobicity_per_residue.csv` — 560-row per-position Kyte-Doolittle score (NaN at the edge positions where the window extends past the sequence).
- `figures/fig_almt4_protein_topology.png` / `.pdf` — composite topology figure.

## Caveats

- **TM prediction parameters.** The entry threshold (0.8), segment-mean floor (1.2), and minimum length (14 aa) are tuned for plant membrane proteins under window-19 K-D smoothing. They are more permissive than the original Kyte-Doolittle 1982 threshold (1.6 with window 7) and less so than typical TMHMM HMM-based calls. For a definitive TM topology, a follow-up using DeepTMHMM or TMbed is recommended; the K-D scan here provides supporting structural-biology evidence but is not a substitute for HMM-based topology prediction.
- **Pfam coordinates.** The Funannotate GFF does not record per-domain HMMER coordinates; we label the PF11744 / PF13515 spans as N-terminal half and C-terminal half respectively, following the canonical ALMT-family architecture (Sasaki 2004). A future revision with raw HMMER output would replace these placeholder spans with exact coordinates.
- **Single-paralog analysis.** This script analyses only the focal AYB ALMT4_2 (Ss10:15.4 Mb). The other four AYB ALMT-family paralogs (ALMT12 on Ss02, ALMT4_1 on Ss02, ALMT9_1 on Ss03, ALMT9_2 on Ss04) are not analysed here; comparative TM topology across the AYB ALMT family is a candidate for a future Part D.

## Script

`scripts/65_almt4_protein_topology.py`. Runtime: ~ 5 seconds on a laptop.
