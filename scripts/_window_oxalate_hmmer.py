#!/usr/bin/env python
"""
HMMER confirmation: do any of the ALMT4_2-window proteins belong to an
oxalate-metabolism enzyme family? Independent of Funannotate's InterProScan.

Two HMMER searches on the 61 protein-coding window genes (the 62nd is a tRNA):
  1. hmmsearch of curated Pfam HMM profiles with the family GATHERING THRESHOLD
     (--cut_ga): a hit above --cut_ga IS a member of that family. Profiles:
       PF01072 Germin / oxalate oxidase
       PF07883 Cupin_2 (oxalate decarboxylase)
       PF00190 Cupin_1
       PF00501 AMP-binding (AAE3 oxalyl-CoA synthetase family)
       PF01554 MatE (MATE transporter)
       PF11744 ALMT  (positive control -> must hit ALMT4_2)
  2. phmmer of the three reference enzymes (AAE3 Q9SMT7, oxalate oxidase P15290,
     oxalate decarboxylase O34714) against the window proteins.

Output: results/R4_ld_decay_focal/tables/window_oxalate_hmmer.csv
Tools : hmmsearch, phmmer (HMMER 3.4); network for Pfam HMM + UniProt fetch.
Author: Benjamin Narh-Madey
"""
import re, subprocess, tempfile, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PROT = ROOT / "refs/ayb_genome/AYB_proteins.faa"
WIN = ROOT / "results/R4_ld_decay_focal/tables/focal_window_genes_annotated.csv"
OUT = ROOT / "results/R4_ld_decay_focal/tables/window_oxalate_hmmer.csv"

HMMS = {"PF00190": "Cupin_1 (germin / oxalate oxidase)",
        "PF07883": "Cupin_2 (oxalate decarboxylase)",
        "PF00501": "AMP-binding (AAE3 oxalyl-CoA synthetase)",
        "PF13193": "AMP-binding_C (AAE3)",
        "PF01554": "MatE (MATE)", "PF11744": "ALMT (positive control)"}
REFS = {"AAE3_oxalyl_CoA_synthetase": "Q9SMT7", "oxalate_oxidase_germin": "P15290",
        "oxalate_decarboxylase": "O34714"}

tmp = Path(tempfile.mkdtemp())
ids = list(pd.read_csv(WIN)["gene_id"])
tx = {f"{g}-T1" for g in ids}
w = open(tmp / "win.faa", "w"); keep = False; n = 0
for line in open(PROT):
    if line.startswith(">"):
        keep = line[1:].split()[0] in tx; n += keep
    if keep:
        w.write(line)
w.close()
print(f"[hmmer] window proteins: {n}")

import gzip
def fetch_hmm(url, dest):
    try:
        data = urllib.request.urlopen(url, timeout=40).read()
    except Exception as e:
        print(f"  [warn] fetch error: {e}")
        return False
    if data[:2] == b"\x1f\x8b":          # gzip magic
        data = gzip.decompress(data)
    if data[:7] != b"HMMER3/":
        return False
    open(dest, "wb").write(data)
    return True

rows = []
# --- 1. Pfam HMM profiles with gathering thresholds -------------------------
for pf, label in HMMS.items():
    hmm = tmp / f"{pf}.hmm"
    if not fetch_hmm(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}?annotation=hmm", str(hmm)):
        print(f"  [warn] could not fetch HMM for {pf}; skipping")
        rows.append({"method": "hmmsearch --cut_ga", "profile": f"{pf} {label}",
                     "window_gene": "(HMM fetch failed)", "evalue": "", "score": "",
                     "passes_threshold": False})
        continue
    dom = tmp / f"{pf}.domtbl"
    # --cut_ga uses the curated family membership threshold
    subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", str(dom), str(hmm), str(tmp / "win.faa")],
                   capture_output=True)
    hits = [l.split() for l in open(dom) if not l.startswith("#")]
    if hits:
        for h in hits:
            rows.append({"method": "hmmsearch --cut_ga", "profile": f"{pf} {label}",
                         "window_gene": h[0].replace("-T1", ""), "evalue": h[6], "score": h[7],
                         "passes_threshold": True})
    else:
        rows.append({"method": "hmmsearch --cut_ga", "profile": f"{pf} {label}",
                     "window_gene": "(none above gathering threshold)", "evalue": "", "score": "",
                     "passes_threshold": False})

# --- 2. phmmer reference enzymes vs window proteins -------------------------
for name, acc in REFS.items():
    fa = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30).read().decode()
    q = tmp / f"{name}.faa"; open(q, "w").write(fa)
    tbl = tmp / f"{name}.tbl"
    subprocess.run(["phmmer", "--tblout", str(tbl), "-E", "10", str(q), str(tmp / "win.faa")],
                   capture_output=True)
    hits = [l.split() for l in open(tbl) if not l.startswith("#")]
    if hits:
        best = min(hits, key=lambda h: float(h[4]))
        rows.append({"method": "phmmer", "profile": name, "window_gene": best[0].replace("-T1", ""),
                     "evalue": best[4], "score": best[5], "passes_threshold": float(best[4]) < 1e-5})
    else:
        rows.append({"method": "phmmer", "profile": name, "window_gene": "(no hit)",
                     "evalue": "", "score": "", "passes_threshold": False})

res = pd.DataFrame(rows)
res.to_csv(OUT, index=False)
print("\n=== HMMER results ===")
print(res.to_string(index=False))
oxpass = res[(res["passes_threshold"]) & (~res["profile"].str.contains("ALMT"))]
print(f"\noxalate-family HMMER hits (excluding the ALMT positive control): {len(oxpass)}")
print(f"ALMT positive control hit: "
      f"{res[(res.profile.str.contains('ALMT'))&(res.passes_threshold)]['window_gene'].tolist()}")
print(f"wrote {OUT}")
