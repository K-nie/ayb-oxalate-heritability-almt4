#!/usr/bin/env python
"""
Reviewer R2.16: besides the ALMT4-proximal focal SNP, is any SNP flagged by several
methods for soluble oxalate? For each SNP, collect the four per-SNP signals used in the
paper -- single-marker MLM p, FarmCPU p, BayesB PIP, BayesCpi PIP -- take each method's
top-10, and count how many methods flag each SNP.

Inputs : results/15_gwas_new_traits/tables/gwas_Soluble_Oxalate_M1_K.csv (MLM)
         results/61_farmcpu_mrmlm/tables/farmcpu_hits_per_trait.csv (FarmCPU)
         results/63_bayes_regression/tables/bayesB_bayesC_pip_per_snp_per_trait.csv (PIP)
         results/89_blast_uniqueness_reanchor/.../marker_placement_classification.csv (positions)
Output : results/61_farmcpu_mrmlm/tables/soluble_multimethod_concordance.csv
Author : Benjamin Narh-Madey
"""
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
FOCAL = "100033542|F|0-31:T>C-31:T>C"
OUT = ROOT / "results/61_farmcpu_mrmlm/tables/soluble_multimethod_concordance.csv"

mlm = pd.read_csv(ROOT / "results/15_gwas_new_traits/tables/gwas_Soluble_Oxalate_M1_K.csv")
pcol = [c for c in mlm.columns if c.lower() in ("p", "p.value", "pvalue", "p_value")][0]
rscol = [c for c in mlm.columns if c.lower() in ("snp", "rs", "marker")][0]
mlm = mlm[[rscol, pcol]].rename(columns={rscol: "rs", pcol: "mlm_p"})
fc = pd.read_csv(ROOT / "results/61_farmcpu_mrmlm/tables/farmcpu_hits_per_trait.csv")
fc = fc[fc.trait == "Soluble_Oxalate"][["SNP", "P.value"]].rename(columns={"SNP": "rs", "P.value": "farmcpu_p"})
pip = pd.read_csv(ROOT / "results/63_bayes_regression/tables/bayesB_bayesC_pip_per_snp_per_trait.csv")
pip = pip[pip.trait == "Soluble_Oxalate"][["rs", "pip_bayesB", "pip_bayesC"]]
anc = pd.read_csv(ROOT / "results/89_blast_uniqueness_reanchor/tables/marker_placement_classification.csv")[["rs", "top_chrom", "top_pos"]]

m = mlm.merge(fc, on="rs", how="outer").merge(pip, on="rs", how="outer").merge(anc, on="rs", how="left")
sets = {"MLM": set(m.dropna(subset=["mlm_p"]).nsmallest(10, "mlm_p").rs),
        "FarmCPU": set(m.dropna(subset=["farmcpu_p"]).nsmallest(10, "farmcpu_p").rs),
        "BayesB": set(m.dropna(subset=["pip_bayesB"]).nlargest(10, "pip_bayesB").rs),
        "BayesC": set(m.dropna(subset=["pip_bayesC"]).nlargest(10, "pip_bayesC").rs)}
m["n_methods_top10"] = m.rs.map(lambda r: sum(r in s for s in sets.values()))
m["is_focal"] = m.rs == FOCAL
multi = m[m.n_methods_top10 >= 2].sort_values("n_methods_top10", ascending=False)
multi.to_csv(OUT, index=False)

print("focal SNP in each method's top-10?:", {k: FOCAL in s for k, s in sets.items()})
print(f"non-focal SNPs flagged by >=2 methods: {len(multi)} ; by all 4: "
      f"{(multi.n_methods_top10 == 4).sum()}")
print(multi[["rs", "top_chrom", "top_pos", "n_methods_top10",
             "mlm_p", "farmcpu_p", "pip_bayesB", "pip_bayesC"]].head(6).to_string(index=False))
print(f"wrote {OUT}")
