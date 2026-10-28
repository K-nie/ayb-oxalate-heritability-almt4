#!/usr/bin/env python
"""
Positive control for the Section 3.6 "no oxalate enzyme in the candidate window"
claim (reviewer request): run the SAME oxalate-enzyme BLASTP and HMMER searches used
on the window genes, but on the WHOLE AYB proteome (30,840 proteins). If the method
recovers oxalate-metabolism genes elsewhere in the genome, then the window negative
is a true absence, not a search that cannot detect these enzymes.

Two searches, identical settings to the window run:
  1. BLASTP of three reference oxalate enzymes vs the whole proteome
       AAE3 oxalyl-CoA synthetase  Q9SMT7
       oxalate oxidase / germin    P15290
       Bacillus oxalate decarboxylase O34714
     significance threshold e < 1e-5 (as in scripts/_window_oxalate_homology.py).
  2. hmmsearch of the oxalate-enzyme Pfam families at their curated GATHERING
     thresholds (--cut_ga), plus PF11744 (ALMT) as an internal positive control:
       PF00190 Cupin_1 (germin / oxalate oxidase)
       PF07883 Cupin_2 (oxalate decarboxylase)
       PF00501 AMP-binding (AAE3)        PF13193 AMP-binding_C (AAE3)
       PF01554 MatE (MATE)               PF11744 ALMT (control -> must hit ALMT4_2)
     Profiles are hmmfetch'd from the local Pfam-A.hmm (no network).

Inputs : refs/ayb_genome/AYB_proteins.faa ; ~/.myconote/dbs/pfam/Pfam-A.hmm ; UniProt (3 refs)
Output : results/R4_ld_decay_focal/tables/proteome_oxalate_positive_control.csv
Author : Benjamin Narh-Madey
"""
import os, re, subprocess, tempfile, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PROT = ROOT / "refs/ayb_genome/AYB_proteins.faa"
PFAM = Path.home() / ".myconote/dbs/pfam/Pfam-A.hmm"
OUT = ROOT / "results/R4_ld_decay_focal/tables/proteome_oxalate_positive_control.csv"
WORK = Path("/private/tmp/claude-501/-Users-black-einstein-Desktop-Other-Projects-"
            "Prof-Adewale-s-Data-manuscript-paper-B/9625fe92-1972-43ee-ab20-"
            "ffd2006b7b78/scratchpad/poscontrol")
WORK.mkdir(parents=True, exist_ok=True)

REFS = {"AAE3_oxalyl_CoA_synthetase": "Q9SMT7",
        "oxalate_oxidase_germin": "P15290",
        "oxalate_decarboxylase": "O34714"}
HMMS = {"PF00190": "Cupin_1 (germin/oxalate oxidase)",
        "PF07883": "Cupin_2 (oxalate decarboxylase)",
        "PF00501": "AMP-binding (AAE3)",
        "PF13193": "AMP-binding_C (AAE3)",
        "PF01554": "MatE (MATE)",
        "PF11744": "ALMT (positive control)"}
EVAL = 1e-5
rows = []


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def extract_profile(pfam_hmm, acc, dest):
    """Pull one HMM record (text) from Pfam-A.hmm by accession prefix."""
    keep, buf, hit = False, [], False
    with open(pfam_hmm) as fh:
        for ln in fh:
            if ln.startswith("HMMER3"):
                buf, hit = [ln], False
            buf.append(ln) if not ln.startswith("HMMER3") else None
            if ln.startswith("ACC") and ln.split()[1].startswith(acc):
                hit = True
            if ln.strip() == "//":
                if hit:
                    Path(dest).write_text("".join(buf))
                    return True
                buf = []
    return False


# --- proteome BLAST db -------------------------------------------------------
db = WORK / "ayb_prot"
if not (WORK / "ayb_prot.pin").exists():
    r = sh(f"makeblastdb -in '{PROT}' -dbtype prot -out '{db}'")
    assert (WORK / "ayb_prot.pin").exists(), f"makeblastdb failed: {r.stderr}"

# --- 1. BLASTP the three reference enzymes vs the whole proteome ------------
for name, acc in REFS.items():
    q = WORK / f"{acc}.fa"
    if not q.exists():
        data = urllib.request.urlopen(
            f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=40).read().decode()
        q.write_text(data)
    r = sh(f"blastp -query '{q}' -db '{db}' -evalue 10 "
           f"-outfmt '6 sseqid pident length evalue bitscore' -max_target_seqs 100000")
    hits = [ln.split("\t") for ln in r.stdout.strip().splitlines() if ln]
    sig = [h for h in hits if float(h[3]) < EVAL]
    best = min(hits, key=lambda h: float(h[3])) if hits else None
    rows.append({"method": "blastp", "query": f"{name} ({acc})",
                 "n_hits_evalue<1e-5": len(sig),
                 "best_hit": best[0] if best else "none",
                 "best_evalue": best[3] if best else "",
                 "best_pident": best[1] if best else ""})
    print(f"[blastp] {name}: {len(sig)} proteome hits e<1e-5; "
          f"best {best[0] if best else '-'} e={best[3] if best else '-'}")

# --- 2. hmmsearch each oxalate Pfam profile (--cut_ga) vs the proteome -------
for pf, label in HMMS.items():
    prof = WORK / f"{pf}.hmm"
    if not prof.exists() or prof.stat().st_size == 0:
        ok = extract_profile(PFAM, pf, prof)
        assert ok, f"profile {pf} not found in {PFAM}"
    tbl = WORK / f"{pf}.tblout"
    sh(f"hmmsearch --cut_ga --tblout '{tbl}' '{prof}' '{PROT}'")
    hit_ids = []
    if tbl.exists():
        for ln in tbl.read_text().splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            hit_ids.append(ln.split()[0])
    almt4 = any("029726" in h for h in hit_ids)
    rows.append({"method": "hmmsearch --cut_ga", "query": f"{pf} {label}",
                 "n_hits_evalue<1e-5": len(hit_ids), "best_hit": hit_ids[0] if hit_ids else "none",
                 "best_evalue": "", "best_pident": ("ALMT4_2 recovered" if almt4 else "")})
    print(f"[hmmsearch] {pf} {label}: {len(hit_ids)} proteome members at --cut_ga"
          + ("  (ALMT4_2 recovered)" if almt4 else ""))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f"\nwrote {OUT}")
print(f"Pfam-A profiles from {PFAM} (local release)")
