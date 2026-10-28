#!/usr/bin/env python
"""
Combine the four oxalate-enzyme-exclusion methods for the 62 ALMT4_2-window genes
into one per-gene table (Supplementary Table S13).

Methods (all in results/R4_ld_decay_focal/tables/):
  1. product name            - Funannotate product (window_oxalate_homology.csv)
  2. domain                  - Pfam/InterPro domains + oxalate-domain flag (same)
  3. BLASTP                  - best hit vs 3 reference enzymes (same)
  4. HMMER                   - oxalate-family Pfam HMM hit at --cut_ga + phmmer
                               (window_oxalate_hmmer.csv), mapped back per gene

Output: results/R4_ld_decay_focal/tables/focal_window_genes_4method.csv
Author: Benjamin Narh-Madey
"""
from pathlib import Path
import pandas as pd

T = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
         "/results/R4_ld_decay_focal/tables")
hom = pd.read_csv(T / "window_oxalate_homology.csv")          # methods 1-3, per gene
hmm = pd.read_csv(T / "window_oxalate_hmmer.csv")             # method 4, per profile
win = pd.read_csv(T / "focal_window_genes_annotated.csv")[
        ["gene_id", "chr", "start", "end", "strand"]]

# method 4 per gene: which window gene was hit by an oxalate-family HMM (--cut_ga)?
hmm_pass = hmm[(hmm["passes_threshold"] == True) &
               (hmm["window_gene"].str.startswith("AYBTSS", na=False))]
hmm_map = {r["window_gene"]: f"{r['profile']} (e={r['evalue']})"
           for _, r in hmm_pass.iterrows()}

df = hom.merge(win, on="gene_id", how="left")
df["m4_hmmer_oxalate_family_hit"] = df["gene_id"].map(hmm_map).fillna("none")
# overall call: oxalate-relevant only if a genuine oxalate signal by ANY method
df["oxalate_relevant"] = df["gene_id"].eq("AYBTSS11_029726")

out = df.rename(columns={
    "product": "m1_product_name",
    "pfam_domains": "m2_pfam_domains",
    "oxalate_enzyme_domain": "m2_oxalate_enzyme_domain",
    "blast_best_oxalate_enzyme": "m3_blast_best_enzyme",
    "blast_best_evalue": "m3_blast_best_evalue",
    "blast_verdict": "m3_blast_verdict",
})[[
    "gene_id", "name", "feature_type", "chr", "start", "end", "strand",
    "distance_from_focal_kb", "in_half_decay_99kb", "is_candidate",
    "m1_product_name", "m2_pfam_domains", "m2_oxalate_enzyme_domain",
    "m3_blast_best_enzyme", "m3_blast_best_evalue", "m3_blast_verdict",
    "m4_hmmer_oxalate_family_hit", "oxalate_relevant",
]].sort_values("distance_from_focal_kb")

OUT = T / "focal_window_genes_4method.csv"
out.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(out)} genes; oxalate-relevant by any method: "
      f"{int(out['oxalate_relevant'].sum())} = {out[out.oxalate_relevant]['name'].tolist()})")
print("\nmethod summary (non-candidate genes with ANY oxalate signal):")
hit = out[(~out.is_candidate) & ((out.m2_oxalate_enzyme_domain.fillna('').str.strip() != '')
          | (out.m3_blast_verdict == 'significant oxalate hit')
          | (out.m4_hmmer_oxalate_family_hit != 'none'))]
print("  rows:", len(hit), "(0 = clean across all four methods)")