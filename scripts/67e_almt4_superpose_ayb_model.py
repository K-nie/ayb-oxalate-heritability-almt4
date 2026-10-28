"""Superpose the de-novo AYB ALMT4_2 ESMFold model on GmALMT12/QUAC1 (7W6K)
and regenerate Fig 9A. Structural context only; no functional claim.

Run with: pymol -cq scripts/67e_almt4_superpose_ayb_model.py
Reuses the logic of scripts/67b but uses the AYB-specific ESMFold model as the
primary structure (soybean proxy kept only as a cross-check file).
"""
from pymol import cmd, stored
import json, os, urllib.request
import numpy as np

ROOT = "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
OUT = ROOT + "/results/67_almt4_structural_context"
AYB_MODEL = OUT + "/data/AYB_ALMT4_2_ESMFold.pdb"
TMPL_A = OUT + "/data/7W6K_chainA.pdb"

# ensure template chain A exists (download 7W6K + extract chain A if needed)
if not os.path.exists(TMPL_A):
    raw = OUT + "/data/7W6K.pdb"
    if not os.path.exists(raw):
        req = urllib.request.Request("https://files.rcsb.org/download/7W6K.pdb",
                                     headers={"User-Agent": "Mozilla/5.0"})
        open(raw, "wb").write(urllib.request.urlopen(req).read())
    cmd.load(raw, "full7w6k")
    cmd.save(TMPL_A, "full7w6k and chain A and polymer")
    cmd.delete("full7w6k")

cmd.load(TMPL_A, "tmpl")
cmd.load(AYB_MODEL, "ayb")

r = cmd.cealign("tmpl", "ayb")
print("CEALIGN_FULL_RMSD %.3f NALN %d" % (r["RMSD"], r["alignment_length"]))
cmd.create("tmpl_tm", "tmpl and resi 36-330")
cmd.create("ayb_tm", "ayb and resi 1-345")
rt = cmd.cealign("tmpl_tm", "ayb_tm")
print("CEALIGN_TMcore_RMSD %.3f NALN %d" % (rt["RMSD"], rt["alignment_length"]))
cmd.delete("tmpl_tm"); cmd.delete("ayb_tm")

# six experimentally-defined TM helices of 7W6K -> map to AYB residues (nearest CA)
TMPL_TM = [(55, 78), (80, 99), (104, 131), (135, 160), (168, 188), (189, 215)]
def ca(sel):
    stored.xyz = {}
    cmd.iterate_state(1, sel + " and name CA", "stored.xyz[int(resi)]=(x,y,z)")
    return dict(stored.xyz)
tca, aca = ca("tmpl"), ca("ayb")
ak = list(aca); aarr = np.array([aca[i] for i in ak])
ayb_tm = []
for (a, b) in TMPL_TM:
    mapped = []
    for ri in range(a, b + 1):
        if ri not in tca:
            continue
        d = np.linalg.norm(aarr - np.array(tca[ri]), axis=1)
        j = int(d.argmin())
        if d[j] < 4.0:
            mapped.append(ak[j])
    if mapped:
        ayb_tm.append((min(mapped), max(mapped)))
print("AYB_TM_RANGES", ayb_tm, "N=%d" % len(ayb_tm))

json.dump({"model": "de-novo ESMFold (facebook/esmfold_v1) of AYB ALMT4_2 (560 aa)",
           "template": "GmALMT12/QUAC1 PDB 7W6K chain A (Qin et al. 2022 Sci Adv 8:eabm3238)",
           "cealign_full_rmsd_A": round(r["RMSD"], 3), "cealign_full_naln": r["alignment_length"],
           "tmcore_cealign_rmsd_A": round(rt["RMSD"], 3), "tmcore_cealign_naln": rt["alignment_length"],
           "n_TM_helices_mapped_to_template": len(ayb_tm),
           "ayb_TM_helix_ranges": ayb_tm,
           "note": "TM count established by 1:1 superposition onto the six experimentally-defined TM helices of the template."},
          open(OUT + "/tables/superpose_metrics.json", "w"), indent=2)

# ---- Fig 9A: AYB TM helices Okabe-Ito on grey 7W6K ----
cmd.hide("everything"); cmd.bg_color("white")
cmd.show("cartoon", "tmpl and resi 36-225")
cmd.show("cartoon", "ayb and resi 40-235")
cmd.set("cartoon_transparency", 0.6, "tmpl")
cmd.color("grey70", "tmpl")
cmd.color("grey88", "ayb")
okabe = ["0x0072B2", "0xD55E00", "0x009E73", "0xCC79A7", "0xE69F00", "0x56B4E9"]
for k, (a, b) in enumerate(ayb_tm[:6]):
    cmd.color(okabe[k % len(okabe)], "ayb and resi %d-%d" % (a, b))
cmd.set("cartoon_transparency", 0.0, "ayb")
cmd.orient("ayb and resi 40-235")
cmd.turn("z", 90)
cmd.zoom("ayb and resi 40-235", 3)
cmd.set("ray_opaque_background", 0); cmd.set("antialias", 2); cmd.set("ray_shadows", 0)
cmd.ray(1700, 2100)
cmd.png(OUT + "/figures/fig_almt4_structure_superposition.png", dpi=300)
# also a PDB-session for later vector/figure edits
cmd.save(OUT + "/data/fig9A_superposition_session.pse")
print("WROTE Fig 9A (AYB ESMFold model) N_TM=%d" % len(ayb_tm))
