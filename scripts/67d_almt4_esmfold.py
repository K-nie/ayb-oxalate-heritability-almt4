#!/usr/bin/env python3
"""De-novo ESMFold model of AYB ALMT4_2 (single-sequence, no MSA).

TPG-2026-09-0323 revision, Fig 9A. STRUCTURAL CONTEXT ONLY — confirms the
canonical six-TM ALMT fold by homology/structure; makes NO functional or
oxalate-transport claim.

Uses HuggingFace `facebook/esmfold_v1` (ESMFold; Lin et al. 2023, Science
379:1123) locally on CPU. No MSA, no 400-aa cap (that cap is only the public
web API). Writes a PDB with per-residue pLDDT in the B-factor column.

Input : results/67_almt4_structural_context/data/AYB_ALMT4_2_AYBTSS11_029726.faa
Output: results/67_almt4_structural_context/data/AYB_ALMT4_2_esmfold.pdb
        results/67_almt4_structural_context/tables/esmfold_plddt.json

Author: Benjamin Narh-Madey
"""
from pathlib import Path
import json
import numpy as np
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUT = ROOT / "results/67_almt4_structural_context"
FAA = OUT / "data/AYB_ALMT4_2_AYBTSS11_029726.faa"
PDB_OUT = OUT / "data/AYB_ALMT4_2_esmfold.pdb"

seq = "".join(l.strip() for l in open(FAA) if not l.startswith(">"))
print(f"[seq] AYB ALMT4_2 length {len(seq)}")

torch.set_num_threads(max(1, torch.get_num_threads()))
print("[load] facebook/esmfold_v1 (this downloads ~2.8 GB on first run)")
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
# standard load (fast) -- 34 GB RAM comfortably holds the 3B-param model + load buffer
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
model = model.eval()
# ESMFold casts the ESM-2 trunk to fp16 by default; CPU has no Half LayerNorm
# kernel, so force the whole model to float32 for CPU inference.
model.esm = model.esm.float()
model = model.float()
# memory: chunk the folding trunk so 560 residues fit comfortably in RAM
model.trunk.set_chunk_size(64)
import transformers
print(f"[versions] transformers {transformers.__version__} torch {torch.__version__}")

ids = tok([seq], return_tensors="pt", add_special_tokens=False)["input_ids"]
print("[fold] running ESMFold forward on CPU ...")
with torch.no_grad():
    out = model(ids)

pdb = model.output_to_pdb(out)[0]
PDB_OUT.write_text(pdb)
print(f"[write] {PDB_OUT}")

# per-residue pLDDT: output_to_pdb already wrote pLDDT into B-factors.
# read them back from CA atoms for a consistent summary.
from Bio.PDB import PDBParser
import warnings; warnings.simplefilter("ignore")
s = PDBParser(QUIET=True).get_structure("ayb", str(PDB_OUT))
ca = [a for a in s.get_atoms() if a.name == "CA"]
b = np.array([a.bfactor for a in ca])
tmd = b[:340]            # N-terminal TM domain region
summary = {
    "esmfold": "facebook/esmfold_v1 (HuggingFace transformers)",
    "transformers_version": transformers.__version__,
    "torch_version": torch.__version__,
    "n_residues": int(len(b)),
    "mean_pLDDT": round(float(b.mean()), 1),
    "TMdomain_1_340_pLDDT": round(float(tmd.mean()), 1),
    "frac_pLDDT_gt70": round(float((b > 70).mean()), 2),
    "TMdomain_frac_gt70": round(float((tmd > 70).mean()), 2),
}
json.dump(summary, open(OUT / "tables/esmfold_plddt.json", "w"), indent=2)
print("[plddt]", summary)
