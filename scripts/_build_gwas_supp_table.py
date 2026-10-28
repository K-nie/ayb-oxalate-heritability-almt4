#!/usr/bin/env python
"""
Build the journal-required per-SNP GWAS results table (The Plant Genome GWAS template),
long format: one row per SNP per trait per model.

For every tested SNP it records the allele scored (the alternate allele the 0/1/2 dosage
counts), its position on the S. stenocarpa reference (or the DArT tag sequence when
unanchored), MAF, n, effect (beta, se) and p-value, the model, and whether the row comes
from the genome-wide scan. The ALMT4-proximal focal SNP is added as its own rows
(G-only and G + PC direct tests), flagged in_genome_wide_scan = FALSE because it fell
below the subset MAF floor and was tested directly, not recovered from the scan.

Model reported per trait (matching the paper):
  - the four bulk-biochemistry traits (Antioxidant, Flavonoid, Phenol, Tannin): G + PC
    (results/03_gwas_mlm, M2_KQ) -- the model the paper reports for those scans;
  - the nine remaining traits (five seed-morphometrics, crude protein, three oxalate
    fractions): G (results/15_gwas_new_traits, M1_K) -- the only model run for them.
The focal SNP carries both a G row and a G + PC row.

Inputs : results/03_gwas_mlm/tables/gwas_<trait>_M2_KQ.csv
         results/15_gwas_new_traits/tables/gwas_<trait>_M1_K.csv
         results/01_qc_pca_power/tables/marker_qc.csv              (alleles ref/alt)
         results/89_blast_uniqueness_reanchor/.../marker_placement_classification.csv (position)
         refs/ayb_genome/dartseq_tags.fasta                        (tag sequence)
         results/60_almt4_pleiotropy_scan/tables/almt4_per_trait_effects.csv (per-trait n)
         results/R2_focal_snp_honest/tables/exact_tests_sommer.csv (focal G / G+PC)
Output : results/_gwas_supp_table/SupplementaryData_S1_gwas_per_snp_per_trait.csv
Author : Benjamin Narh-Madey
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUTDIR = ROOT / "results/_gwas_supp_table"
OUTDIR.mkdir(exist_ok=True)
OUT = OUTDIR / "TableS18_gwas_per_snp_per_trait.csv"

GPC_TRAITS = ["Antioxidant", "Flavonoid", "Phenol", "Tannin"]          # G + PC (dir 03)
G_TRAITS = ["Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
            "Seed_Coat_Tannin", "Crude_Protein", "Total_Oxalate",
            "Soluble_Oxalate", "Insoluble_Oxalate"]                     # G (dir 15)
FOCAL = "100033542|F|0-31:T>C-31:T>C"

# --- allele map (ref/alt) + per-trait n + position + tag sequence ----------
qc = pd.read_csv(ROOT / "results/01_qc_pca_power/tables/marker_qc.csv")
qc[["ref_allele", "alt_allele"]] = qc["alleles"].str.split("/", expand=True)
allele = qc.set_index("rs")[["ref_allele", "alt_allele"]]

anc = pd.read_csv(ROOT / "results/89_blast_uniqueness_reanchor/tables/"
                  "marker_placement_classification.csv").set_index("rs")
pos = anc[["top_chrom", "top_pos"]].rename(
    columns={"top_chrom": "chromosome", "top_pos": "position"})

nmap = pd.read_csv(ROOT / "results/60_almt4_pleiotropy_scan/tables/"
                   "almt4_per_trait_effects.csv").set_index("trait")["n"].to_dict()

# Tag sequences keyed by CLONE id (the field before "|"): one DArT clone = one 69-bp
# tag, so clone-level matching is robust to the SNP-position label differing between
# the GWAS marker id and the FASTA header. The longest sequence per clone is kept.
tag = {}
with open(ROOT / "refs/ayb_genome/dartseq_tags.fasta") as fh:
    h, s = None, []

    def _flush(header, seq):
        if header is None:
            return
        clone = header.split("|")[0]
        if clone not in tag or len(seq) > len(tag[clone]):
            tag[clone] = seq
    for ln in fh:
        if ln.startswith(">"):
            _flush(h, "".join(s))
            h, s = ln[1:].strip(), []
        else:
            s.append(ln.strip())
    _flush(h, "".join(s))


def rows_for(trait, path, model):
    g = pd.read_csv(path)
    mafcol = "maf" if "maf" in g.columns else "maf_subset"
    r = pd.DataFrame({
        "trait": trait, "marker_id": g["rs"],
        "maf": g[mafcol], "n": nmap.get(trait, np.nan),
        "beta": g["beta"], "se": g["se"], "p": g["p"],
        "model": model, "in_genome_wide_scan": True,
    })
    return r


frames = []
for t in GPC_TRAITS:
    frames.append(rows_for(t, ROOT / f"results/03_gwas_mlm/tables/gwas_{t}_M2_KQ.csv", "G + PC"))
for t in G_TRAITS:
    # The dir-15 "M1_K" output is MLM with K + intercept + PC1-PC3 (script 15,
    # emmax scan design matrix), i.e. G + PC -- the filename is legacy-misleading.
    frames.append(rows_for(t, ROOT / f"results/15_gwas_new_traits/tables/gwas_{t}_M1_K.csv", "G + PC"))

# --- focal SNP direct-test rows (soluble oxalate), not in the scan ---------
ex = pd.read_csv(ROOT / "results/R2_focal_snp_honest/tables/exact_tests_sommer.csv")
w = ex.iloc[0]   # 'exact Wald: Soluble ~ dose'
focal_maf = float(qc.set_index("rs").loc[FOCAL, "maf"])
focal = pd.DataFrame([
    {"trait": "Soluble_Oxalate", "marker_id": FOCAL, "maf": focal_maf, "n": 41,
     "beta": w["beta"], "se": w["se"], "p": float(w["EMMAX_Konly_p"]),
     "model": "G", "in_genome_wide_scan": False},
    {"trait": "Soluble_Oxalate", "marker_id": FOCAL, "maf": focal_maf, "n": 41,
     "beta": np.nan, "se": np.nan, "p": float(w["EMMAX_KPC_p"]),
     "model": "G + PC", "in_genome_wide_scan": False},
])
frames.append(focal)

df = pd.concat(frames, ignore_index=True)

# --- join allele, position, tag sequence -----------------------------------
df = df.join(allele, on="marker_id").join(pos, on="marker_id")
df["ref_allele"] = df["ref_allele"].fillna(df["marker_id"].str.extract(r":([ACGT])>")[0])
df["alt_allele"] = df["alt_allele"].fillna(df["marker_id"].str.extract(r">([ACGT])")[0])
df["allele_scored"] = df["alt_allele"]      # dosage counts the alternate allele
df["tag_sequence"] = df["marker_id"].str.split("|").str[0].map(tag)   # match by clone id
# tag sequence is the fallback locator; keep it only where there is no position
df.loc[df["position"].notna(), "tag_sequence"] = ""

df["position"] = df["position"].astype("Int64")
df = df[["trait", "marker_id", "chromosome", "position", "ref_allele", "alt_allele",
         "allele_scored", "maf", "n", "beta", "se", "p", "model",
         "in_genome_wide_scan", "tag_sequence"]]
df = df.sort_values(["trait", "in_genome_wide_scan", "p"], ascending=[True, True, True])
df.to_csv(OUT, index=False)

# --- report ---------------------------------------------------------------
print(f"wrote {OUT}")
print(f"rows: {len(df)}  | traits: {df.trait.nunique()}  | unique markers: {df.marker_id.nunique()}")
print("rows per trait:")
print(df.groupby("trait").size().to_string())
unplaced = df[df.position.isna()]
print(f"\nrows with no position (tag sequence given instead): {len(unplaced)} "
      f"({unplaced.marker_id.nunique()} markers)")
print("tag sequence present for every unplaced row:",
      (unplaced.tag_sequence.str.len() > 0).all())
print("\nfocal SNP rows:")
print(df[df.marker_id == FOCAL][["trait", "chromosome", "position", "allele_scored",
                                 "p", "model", "in_genome_wide_scan"]].to_string(index=False))
