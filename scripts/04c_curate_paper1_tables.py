#!/usr/bin/env python3
"""
Assemble the canonical Paper 1 table set into
`results/04_publication_plots/tables/`.

Reads from the per-analysis directories (the authoritative copies) and
writes paper-numbered duplicates so the manuscript can reference a single
Table 1..5 (main) + Table S1..S4 (supplementary) series.

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
RES  = ROOT / "results"
DST  = RES / "04_publication_plots" / "tables"
DST.mkdir(parents=True, exist_ok=True)

# (paper-table label, source CSV basename, source per-analysis dir,
#  short slug used for the destination filename).
TBLS = [
    # main paper tables
    ("Table1",  "amova_table.csv",
                "02_amova",                "amova_partition"),
    ("Table2",  "cophenetic_correlation_matrix.csv",
                "19_phylogenetic_trees",   "tree_method_cophenetic"),
    ("Table3",  "candidate_genes_pathway_hits.csv",
                "21_candidate_genes_ayb",  "candidate_genes_top30"),
    ("Table4",  "cross_pairs_least_related.csv",
                "13_grm_crosspairs",       "least_related_cross_pairs"),
    ("Table5",  "breeder_bivariate_pairs.csv",
                "17_trait_extremes_13",    "phenotype_extremes"),
    # supplementary tables
    ("TableS1", "trait_model_summary.csv",
                "03_gwas_mlm",             "trait_model_summary"),
    ("TableS2", "gblup_prediction_accuracy.csv",
                "03_gwas_mlm",             "gblup_prediction_accuracy"),
    ("TableS3", "gebv_per_accession.csv",
                "18_gebv",                 "gebv_per_accession"),
    ("TableS4", "core_collection_top20_stratified.csv",
                "06_breeding",             "core_collection_top20_stratified"),
]


print(f"[curate] target = {DST}")
manifest_rows = []
missing = []
for label, src_csv, src_dir, slug in TBLS:
    src = RES / src_dir / "tables" / src_csv
    if not src.exists():
        missing.append(f"  {label}: {src}")
        continue
    dst = DST / f"{label}_{slug}.csv"
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    manifest_rows.append((label, slug, src_csv, src_dir,
                          str(dst.relative_to(ROOT))))

manifest_path = DST.parent / "table_manifest.csv"
with open(manifest_path, "w") as fh:
    fh.write("table_label,slug,source_csv,source_dir,destination_path\n")
    for row in manifest_rows:
        fh.write(",".join(str(v) for v in row) + "\n")

print(f"[curate] copied {len(manifest_rows)} tables")
print(f"[curate] manifest -> {manifest_path}")
if missing:
    print("[!] missing sources:")
    for m in missing:
        print(m)
