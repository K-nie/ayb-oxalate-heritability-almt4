#!/usr/bin/env python
"""
Homology check: does any gene in the ALMT4_2 candidate window (62 genes,
±391 kb extinction window) encode an oxalate-metabolism enzyme?

Motivation: the Funannotate reference annotation leaves 81% of the window genes
as "hypothetical protein", and several named genes carry fungal-ortholog names
(the top Swiss-Prot hit). The "no oxalate enzyme in the window" statement should
therefore not rest on product names alone. Here we check three ways:
  1. Product name (Funannotate).
  2. Pfam/InterPro domain assignment (Funannotate ran InterProScan = HMMER
     profile search) — catches domains even in unnamed proteins.
  3. Direct BLASTP of the window proteins against the three alternatives named
     in the Introduction (Arabidopsis AAE3 oxalyl-CoA synthetase Q9SMT7; wheat
     oxalate oxidase / germin P15290; Bacillus oxalate decarboxylase O34714),
     with AYB ALMT4_2 as a positive control.

Inputs : refs/ayb_genome/AYB_proteins.faa, the Funannotate GFF3, the window gene
         list (results/R4_ld_decay_focal). Reference enzymes fetched from UniProt.
Output : results/R4_ld_decay_focal/tables/window_oxalate_homology.csv
Tools  : blastp, makeblastdb (BLAST+); network for the UniProt fetch.

Author: Benjamin Narh-Madey
"""
import re, subprocess, tempfile, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
GFF = ROOT / "refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3"
PROT = ROOT / "refs/ayb_genome/AYB_proteins.faa"
WIN = ROOT / "results/R4_ld_decay_focal/tables/focal_window_genes_annotated.csv"
OUT = ROOT / "results/R4_ld_decay_focal/tables/window_oxalate_homology.csv"

OX_PFAM = {"PF01072": "Germin/oxalate oxidase", "PF07883": "Cupin_2 (oxalate decarboxylase)",
           "PF00190": "Cupin_1", "PF00501": "AMP-binding (AAE3 oxalyl-CoA synthetase)",
           "PF01554": "MatE (MATE)", "PF01565": "FAD-oxidase (glycolate oxidase)",
           "PF11744": "ALMT"}
REFS = {"AAE3_oxalyl_CoA_synthetase": "Q9SMT7",
        "oxalate_oxidase_germin": "P15290",
        "oxalate_decarboxylase": "O34714"}

win = pd.read_csv(WIN)
ids = list(win["gene_id"])
tx = {f"{g}-T1" for g in ids}

# --- 1/2: feature type, product, and Pfam domains for ALL 62 genes ----------
# Parse any transcript line (mRNA, tRNA, rRNA, ...) whose Parent is a window gene.
feat, prod, dom = {}, {}, {}
for line in open(GFF):
    c = line.rstrip().split("\t")
    if len(c) < 9 or c[2] in ("gene", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
        continue
    a = c[-1]
    pid = re.search(r"Parent=([^;]+)", a)
    if not pid:
        continue
    g = pid.group(1)
    if g not in ids or g in feat:
        continue
    feat[g] = c[2]
    pr = re.search(r"product=([^;]+)", a)
    dx = re.search(r"Dbxref=([^;]+)", a)
    prod[g] = pr.group(1) if pr else ""
    dom[g] = dx.group(1) if dx else ""

def pfams(d):
    return ",".join(re.findall(r"PFAM:(PF\d+)", d))

def ox_domains(d):
    return "; ".join(f"{p}:{OX_PFAM[p]}" for p in OX_PFAM if p in d) or ""

# --- 3: extract the protein-coding window proteins, BLASTP vs 3 oxalate enzymes
tmp = Path(tempfile.mkdtemp())
w = open(tmp / "win.faa", "w"); keep = False; have_prot = set()
for line in open(PROT):
    if line.startswith(">"):
        h = line[1:].split()[0]
        keep = h in tx
        if keep:
            have_prot.add(h.replace("-T1", ""))
    if keep:
        w.write(line)
w.close()
# reference enzyme DB
r = open(tmp / "refs.faa", "w")
for name, acc in REFS.items():
    fa = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30).read().decode()
    r.write(f">{name}\n" + "".join(fa.splitlines(keepends=True)[1:]))
r.close()
subprocess.run(["makeblastdb", "-in", str(tmp / "refs.faa"), "-dbtype", "prot",
                "-out", str(tmp / "refdb")], check=True, capture_output=True)
# query = all window proteins; subject = the 3 oxalate enzymes
raw = subprocess.run(["blastp", "-query", str(tmp / "win.faa"), "-db", str(tmp / "refdb"),
                      "-evalue", "100", "-max_target_seqs", "50",
                      "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                                 "qstart qend sstart send evalue bitscore"],
                     capture_output=True, text=True).stdout
RAW = OUT.parent / "window_oxalate_blast_raw.csv"
cols = ["window_protein", "oxalate_enzyme", "pident", "aln_len", "mismatch", "gapopen",
        "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
rawrows = [ln.split("\t") for ln in raw.splitlines()]
pd.DataFrame(rawrows, columns=cols).to_csv(RAW, index=False)
# per-gene best oxalate hit
best = {}
for ln in rawrows:
    g = ln[0].replace("-T1", ""); e = float(ln[10])
    if g not in best or e < best[g][1]:
        best[g] = (ln[1], e, float(ln[2]), int(ln[3]))

# --- assemble the complete per-gene table for all 62 ------------------------
SIG = 1e-5
out_rows = []
for g in ids:
    w_ = win.set_index("gene_id").loc[g]
    coding = g in have_prot
    if not coding:
        bstat, benz, bev = "non-coding (no protein)", "", ""
    elif g in best:
        enz, e, pid_, l_ = best[g]
        benz, bev = enz, e
        bstat = "significant oxalate hit" if e < SIG else "no significant oxalate hit"
    else:
        benz, bev, bstat = "", "", "no oxalate hit"
    out_rows.append({
        "gene_id": g, "name": w_["name"], "feature_type": feat.get(g, "?"),
        "product": prod.get(g, ""), "pfam_domains": pfams(dom.get(g, "")),
        "oxalate_enzyme_domain": ox_domains(dom.get(g, "")),
        "blast_best_oxalate_enzyme": benz, "blast_best_evalue": bev,
        "blast_verdict": bstat, "distance_from_focal_kb": w_["distance_from_focal_kb"],
        "in_half_decay_99kb": w_["in_half_decay_99kb"], "is_candidate": g == "AYBTSS11_029726"})
res = pd.DataFrame(out_rows)
res.to_csv(OUT, index=False)

n_coding = int((res["feature_type"] != "tRNA").sum())
n_oxdom = int(((res["oxalate_enzyme_domain"] != "") & (~res["oxalate_enzyme_domain"].str.contains("ALMT"))).sum())
n_sig = int((res["blast_verdict"] == "significant oxalate hit").sum())
print(f"[homology] 62 window genes: {int((res.feature_type=='mRNA').sum())} protein-coding blasted, "
      f"{int((res.feature_type=='tRNA').sum())} tRNA (non-coding).")
print(f"  non-ALMT oxalate-enzyme domains: {n_oxdom}")
print(f"  significant (e<1e-5) BLASTP oxalate hits: {n_sig}")
print(f"  ALMT4_2 oxalate_enzyme_domain: {res.set_index('gene_id').loc['AYBTSS11_029726','oxalate_enzyme_domain']}")
print(f"wrote {OUT}")
print(f"wrote {RAW}")
