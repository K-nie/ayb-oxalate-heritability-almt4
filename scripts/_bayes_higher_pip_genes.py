#!/usr/bin/env python
"""
Reviewer item (TPG-2026-09-0323, §3.4): which SNPs carry a higher BayesB/BayesCpi
posterior inclusion probability (PIP) for soluble oxalate than the ALMT4-proximal
focal SNP, on which chromosomes, and what gene sits at each?

The point of the result is the honest statement that the Bayesian analysis does NOT
single out the ALMT4-proximal SNP: several SNPs on other chromosomes carry higher
PIP, and the genes at those SNPs are hypothetical / housekeeping, none oxalate-related.

Inputs (all read, nothing recomputed from raw):
  - results/63_bayes_regression/tables/bayesB_bayesC_pip_per_snp_per_trait.csv
      BGLR per-SNP PIPs (BayesB pi=0.05; BayesCpi pi~Beta, mean 0.05).
  - results/89_blast_uniqueness_reanchor/tables/marker_placement_classification.csv
      DArT marker -> (top_chrom, top_pos) megablast anchoring.
  - refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3
      Funannotate v1.8 gene models + products (gene-in / nearest-gene lookup).

Focal SNP: rs 100033542 (Ss10:15,394,609 anchored; ~3.2 kb 5' of ALMT4_2
AYBTSS11_029726). Focal PIP: BayesB 0.0664, BayesCpi 0.0802.

Output: results/63_bayes_regression/tables/soluble_oxalate_higher_pip_than_focal.csv
Author: Benjamin Narh-Madey
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PIP = ROOT / "results/63_bayes_regression/tables/bayesB_bayesC_pip_per_snp_per_trait.csv"
ANC = ROOT / "results/89_blast_uniqueness_reanchor/tables/marker_placement_classification.csv"
GFF = ROOT / "refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3"
OUT = ROOT / "results/63_bayes_regression/tables/soluble_oxalate_higher_pip_than_focal.csv"

FOCAL_RS = "100033542|F|0-31:T>C-31:T>C"
FOCAL_B, FOCAL_C = 0.0664, 0.0802  # focal-SNP PIPs (read from this same table)

pip = pd.read_csv(PIP)
pip = pip[pip.trait == "Soluble_Oxalate"].copy()
anc = pd.read_csv(ANC)[["rs", "top_chrom", "top_pos"]]
m = pip.merge(anc, on="rs", how="left").dropna(subset=["top_pos"])

# gene models + products
genes = []
with open(GFF) as fh:
    for ln in fh:
        if ln.startswith("#") or "\t" not in ln:
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        gid = re.search(r"ID=([^;]+)", f[8])
        genes.append((f[0], int(f[3]), int(f[4]), gid.group(1) if gid else "?"))
g = pd.DataFrame(genes, columns=["chrom", "start", "end", "gene"])

# gene products from the mRNA lines
prod = {}
with open(GFF) as fh:
    for ln in fh:
        if ln.startswith("#") or "\t" not in ln:
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        gid = re.search(r"Parent=([^;]+)", f[8])
        pm = re.search(r"product=([^;]+)", f[8])
        if gid and gid.group(1) not in prod:
            prod[gid.group(1)] = pm.group(1) if pm else "hypothetical protein"


def nearest(chrom, pos):
    gg = g[g.chrom == chrom]
    if gg.empty:
        return ("none", np.nan)
    inside = gg[(gg.start <= pos) & (gg.end >= pos)]
    if len(inside):
        return (inside.iloc[0].gene, 0)
    d = np.minimum((gg.start - pos).abs(), (gg.end - pos).abs())
    i = d.idxmin()
    return (g.loc[i, "gene"], int(d.loc[i]))


hi = m[m.pip_bayesB > FOCAL_B].copy().sort_values("pip_bayesB", ascending=False)
hi["nearest_gene"], hi["dist_bp"] = zip(*[nearest(c, int(p)) for c, p in zip(hi.top_chrom, hi.top_pos)])
hi["gene_product"] = hi.nearest_gene.map(lambda x: prod.get(x, ""))
hi["higher_under_both_priors"] = hi.pip_bayesC > FOCAL_C
hi["top_pos"] = hi.top_pos.astype(int)
hi = hi[["rs", "top_chrom", "top_pos", "pip_bayesB", "pip_bayesC",
         "higher_under_both_priors", "nearest_gene", "dist_bp", "gene_product"]]
hi.to_csv(OUT, index=False)

print(f"focal PIP: BayesB {FOCAL_B}, BayesCpi {FOCAL_C}")
print(f"SNPs with BayesB PIP > focal: {len(hi)}")
print("by chromosome:\n", hi.top_chrom.value_counts().to_string())
print(f"higher under BOTH priors: {int(hi.higher_under_both_priors.sum())}")
print("\ntop rows:")
print(hi.head(8).to_string(index=False))
print(f"\nwrote {OUT}")
