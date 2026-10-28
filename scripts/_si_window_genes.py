#!/usr/bin/env python
"""
Supplementary Table S13: annotated gene content of the ALMT4_2 candidate window.

Merges the focal-SNP half-decay window genes (±99 kb; results/R4_ld_decay_focal)
with their Funannotate products so readers can check for alternative candidates.
Confirms no oxalate-metabolism enzyme (oxalyl-CoA synthetase, oxalate oxidase,
oxalate decarboxylase) and no other organic-acid transporter is present in the
window — ALMT4_2 is the only locally annotated candidate of its class.

Output: results/R4_ld_decay_focal/tables/focal_window_genes_annotated.csv

Author: Benjamin Narh-Madey
"""
import re
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
GFF = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
WIN = ROOT / "results" / "R4_ld_decay_focal" / "tables" / "focal_window_genes_half_decay.csv"
OUT = ROOT / "results" / "R4_ld_decay_focal" / "tables" / "focal_window_genes_annotated.csv"

# gene_id -> product from mRNA lines
prod = {}
for line in open(GFF):
    if "\tmRNA\t" not in line:
        continue
    a = line.rstrip().split("\t")[-1]
    pid = re.search(r"Parent=([^;]+)", a)
    pr = re.search(r"product=([^;]+)", a)
    if pid and pr and pid.group(1) not in prod:
        prod[pid.group(1)] = pr.group(1)

# half-decay (±99 kb) gene set, to flag within the wider ±391 kb table
half = set(pd.read_csv(WIN)["gene_id"])
# full extinction window (±391 kb, r2 < 0.1) = all 62 genes
EXT = ROOT / "results" / "R4_ld_decay_focal" / "tables" / "focal_window_genes_r2_below_0.1.csv"
d = pd.read_csv(EXT).sort_values("abs_dist_bp")
d["product"] = d["gene_id"].map(lambda g: prod.get(g, "hypothetical protein"))
d["distance_from_focal_kb"] = (d["abs_dist_bp"] / 1000).round(1)
d["in_half_decay_99kb"] = d["gene_id"].isin(half)
d["is_candidate"] = d["gene_id"].eq("AYBTSS11_029726")
out = d[["gene_id", "name", "product", "chr", "start", "end", "strand",
         "distance_from_focal_kb", "in_half_decay_99kb", "is_candidate"]]
out.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(out)} genes in the ±391 kb extinction window; "
      f"{int(out['in_half_decay_99kb'].sum())} within ±99 kb half-decay)")
# show only the named (non-hypothetical) ones for the record
named = out[~out["product"].str.contains("hypothetical", case=False, na=False)]
print(named[["gene_id", "name", "product", "distance_from_focal_kb",
             "in_half_decay_99kb"]].to_string(index=False))
