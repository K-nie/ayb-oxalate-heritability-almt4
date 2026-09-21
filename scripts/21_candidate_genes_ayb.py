#!/usr/bin/env python3
"""
Candidate-gene callout for all 13 AYB traits, re-anchored to the new
Sphenostylis stenocarpa chromosome-scale assembly (Shorinola, Domelevo
Entfellner, Uauy et al. 2024, Sci Data 11:1389; assembly accessions ENA
OY731398-OY731408; Funannotate annotation at Zenodo 10.5281/zenodo.13853757).

This replaces the cowpea (Vigna unguiculata IT97K-499-35 v1.0) proxy
anchoring used in scripts 07 and 16. Coverage rises from 382 cowpea-
anchored markers (11.9 %) to 2,862 AYB-anchored markers (89.3 %).

Inputs
------
* AYB SNP anchoring table:
    refs/ayb_genome/ayb_marker_anchoring.csv
  built from megablast of all DArTseq tags (Report_DAf18-2580_SNP_2.csv,
  TrimmedSequence column) against the AYB chromosome-scale fasta
  (Sphenostylis_stenocarpa_chrom.fasta, 11 chr, 650 Mb).
* Funannotate GFF3:
    refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3
  Ss01-Ss11 chromosome IDs; ~46k gene records.
* GWAS tables from scripts 03 (4 biochem traits) and 15 (9 new traits).

Method
------
For each trait:
  1. Restrict GWAS results to AYB-anchored markers.
  2. Take the top-N (default 30) by p-value.
  3. Scan +/- WIN_KB (default 50) of each top SNP against the GFF gene
     records.
  4. Match gene Name / Note / cross-reference fields against trait-specific
     keyword sets:
       biochem traits  -> phenylpropanoid / flavonoid pathway
       seed-metric     -> seed-development regulators
       Crude_Protein   -> seed-storage proteins
       Oxalate         -> oxalate metabolism + MATE transporters

Outputs (results/21_candidate_genes_ayb/)
-----------------------------------------
tables/
    snp_top30_anchored_<TRAIT>.csv          AYB-anchored top SNPs
    candidate_genes_<TRAIT>.csv             all genes within +/- WIN_KB
    candidate_genes_pathway_hits.csv        trait-keyword filtered hits
    ayb_anchoring_summary.csv               coverage stats per trait
figures/
    fig46_manhattan_ayb_<TRAIT>.png/.pdf    annotated Manhattan, AYB coords

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG, TRAIT_PAL
from _pheno import TRAITS_BIOCHEM, TRAITS_SEED, TRAITS_PROT_OX
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
ANCHOR_CSV = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
GFF_PATH = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
GWAS_DIR_BIOCHEM = ROOT / "results" / "03_gwas_mlm" / "tables"
GWAS_DIR_NEW = ROOT / "results" / "15_gwas_new_traits" / "tables"
OUT = ROOT / "results" / "21_candidate_genes_ayb"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

ALL_TRAITS = TRAITS_BIOCHEM + TRAITS_SEED + TRAITS_PROT_OX
TOP_N = 30
WIN_KB = 50
WIN_BP = WIN_KB * 1_000
CHR_ORDER = [f"Ss{i:02d}" for i in range(1, 12)]


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Trait-specific keyword sets (same templates as scripts 07 + 16; the AYB
# annotation strings in the Funannotate GFF are different but the keyword
# regexes still match: products / notes mention "MATE", "MYB", "chalcone",
# "auxin response factor", "AINTEGUMENTA", "cinnamate-4-hydroxylase", etc.).
# ---------------------------------------------------------------------------
PHENYL_KEYS = [
    # Arabidopsis-style names (preserved for backwards-compat if AYB GFF
    # ever picks them up via InterPro descriptions)
    (r"\bphenylalanine ammonia[- ]?lyase\b", "PAL"),
    (r"\bcinnamate[- ]?4[- ]?hydroxylase\b", "C4H"),
    (r"\bCYP73A\b", "C4H"),
    (r"\b4[- ]?coumarate.*CoA ligase\b", "4CL"),
    # Funannotate / InterPro / Pfam strings as actually emitted in the
    # Sphenostylis_stenocarpa annotation
    (r"\bchalcone synthase\b", "CHS"),
    (r"\bchalcone[- ]?(flavanone)? isomerase\b", "CHI"),
    (r"\bnaringenin\b", "CHI/F3H"),
    (r"\bflavanone 3[- ]?hydroxylase\b", "F3H"),
    (r"\bdihydroflavonol[- ]?4[- ]?reductase\b", "DFR"),
    (r"\banthocyanidin synthase\b", "ANS"),
    (r"\banthocyanidin reductase\b", "ANR"),
    (r"\bleuco?anthocyanidin reductase\b", "LAR"),
    (r"\btransparent testa\b", "TT-regulator"),
    (r"\bMyb\b", "MYB"),
    (r"\bmyb[- ]?like\b", "MYB"),
    (r"\btranscription factor myb\b", "MYB"),
    (r"\bflavonol synthase\b", "FLS"),
    (r"\bglucosyltransferase\b", "UGT"),
    (r"\bhydroxycinnamoyl\b", "HCT/HQT"),
    (r"\bproanthocyanidin\b", "PA-biosynth"),
    (r"\bO-methyltransferase\b", "OMT"),
    (r"\bcaffeoyl\b", "caffeoyl-OMT"),
    (r"\bshikimate\b", "shikimate-path"),
    (r"\bchorismate\b", "chorismate-path"),
]

SEED_DEV_KEYS = [
    # Arabidopsis-style names retained
    (r"\bAINTEGUMENTA\b", "ANT"),
    (r"\bAP2[- ]?like\b", "AP2"),
    (r"\bSHATTERPROOF\b", "SHP"),
    (r"\bSEEDSTICK\b", "STK"),
    (r"\bHAIKU\b", "IKU"),
    (r"\bMINI3\b", "MINI3"),
    (r"\bBIG BROTHER\b", "BB"),
    (r"\bDA1\b", "DA1"),
    # Funannotate-style strings actually emitted by AYB GFF
    (r"\bauxin[- ]?response factor\b", "ARF"),
    (r"\bARF[ \-]?[0-9]+\b", "ARF"),
    (r"\bADP[- ]?ribosylation factor\b", "ARF-related"),
    (r"\bAUX[/_]?IAA\b", "AUX/IAA"),
    (r"\bbrassinosteroid\b", "BR-pathway"),
    (r"\bexpansin\b", "expansin"),
    (r"\bpectin methylesterase\b", "PME"),
    (r"\bcell wall\b", "cell-wall-mod"),
    (r"\bintegument\b", "integument-dev"),
    (r"\bSHORT INTEGUMENTS\b", "SIN"),
    (r"\bovule\b", "ovule-dev"),
    (r"\bseed coat\b", "seed-coat-dev"),
    (r"\bCYP78A\b", "KLU"),
    (r"\bcytochrome P450 78\b", "KLU"),
    (r"\bAGAMOUS\b", "AG"),
    (r"\bMADS[- ]?box\b", "MADS-box"),
    (r"\btranscription factor [A-Z]{2,}\b", "TF-other"),
]

STORAGE_PROTEIN_KEYS = [
    (r"\bcupin\b", "cupin"),
    (r"\blegumin\b", "legumin"),
    (r"\bvicilin\b", "vicilin"),
    (r"\bconglutin\b", "conglutin"),
    (r"\bconglycinin\b", "conglycinin"),
    (r"\b11S globulin\b", "11S-globulin"),
    (r"\b7S globulin\b", "7S-globulin"),
    (r"\bseed storage protein\b", "SSP"),
    (r"\bglutelin\b", "glutelin"),
    (r"\balbumin\b", "albumin"),
    (r"\bNF[- ]?Y[ABC][0-9]?\b", "NF-Y"),
    (r"\bnuclear factor Y\b", "NF-Y"),
    (r"\bLEAFY COTYLEDON\b", "LEC"),
    (r"\bFUSCA3\b", "FUS3"),
    (r"\bABI3\b", "ABI3"),
    (r"\bABSCISIC ACID INSENSITIVE\b", "ABI"),
    (r"\bglutamine synthetase\b", "GS"),
    (r"\basparagine synthetase\b", "AS"),
    (r"\bglutamate dehydrogenase\b", "GDH"),
    (r"\bglutamate synthase\b", "GOGAT"),
    (r"\bnitrogen\b", "N-metabolism"),
    (r"\bproteinase inhibitor\b", "PI"),
    (r"\bKunitz\b", "Kunitz-PI"),
    (r"\bBowman[- ]?Birk\b", "BBI"),
]

OXALATE_KEYS = [
    # canonical names
    (r"\boxalate oxidase\b", "OXO"),
    (r"\bgermin\b", "germin-OXO"),
    (r"\boxalate decarboxylase\b", "OXDC"),
    (r"\boxalyl[- ]?CoA\b", "oxalyl-CoA"),
    (r"\bAAE3\b", "AAE3"),
    (r"\bisocitrate lyase\b", "ICL"),
    (r"\bglycolate oxidase\b", "GLO"),
    (r"\bglyoxylate\b", "glyoxylate"),
    # Funannotate naming: MATE -> "DETOXIFICATION" family
    (r"\bDETOXIFICATION\b", "MATE/DTX"),
    (r"\bDTX[0-9]+\b", "MATE/DTX"),
    (r"\bMulti[- ]?Antimicrobial Extrusion\b", "MATE"),
    (r"\bMATE\b", "MATE"),
    # ALMT / aluminum-activated malate transporter (organic acid efflux)
    (r"\baluminum[- ]?activated malate transporter\b", "ALMT"),
    (r"\bALMT[0-9]+\b", "ALMT"),
    # calcium transporters
    (r"\bcalcium transporter\b", "Ca-transporter"),
    (r"\bCAX[0-9]?\b", "CAX"),
    (r"\bCa2\+[- ]?ATPase\b", "Ca-ATPase"),
    (r"\bP[- ]?type.*ATPase\b", "P-ATPase"),
    # oxalate-pathway adjacent
    (r"\bascorbate\b", "ascorbate-pathway"),
    (r"\b2-oxoglutarate.*oxygenase\b", "2OG-oxygenase"),
    (r"\bvacuolar\b", "vacuolar"),
    (r"\befflux\b", "efflux"),
    (r"\btonoplast\b", "tonoplast"),
]

TRAIT_KEYWORDS = {t: PHENYL_KEYS for t in TRAITS_BIOCHEM + ["Seed_Coat_Tannin"]}
TRAIT_KEYWORDS.update({t: SEED_DEV_KEYS for t in TRAITS_SEED if t != "Seed_Coat_Tannin"})
TRAIT_KEYWORDS["Crude_Protein"] = STORAGE_PROTEIN_KEYS
for t in ["Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"]:
    TRAIT_KEYWORDS[t] = OXALATE_KEYS
COMPILED = {t: [(re.compile(p, re.IGNORECASE), lab) for p, lab in v]
            for t, v in TRAIT_KEYWORDS.items()}


def classify(text: str, trait: str):
    if not isinstance(text, str) or not text:
        return None
    hits = []
    for pat, lab in COMPILED[trait]:
        if pat.search(text):
            hits.append(lab)
    if not hits:
        return None
    seen, out = set(), []
    for h in hits:
        if h not in seen:
            seen.add(h); out.append(h)
    return ",".join(out)


# ---------------------------------------------------------------------------
# Load AYB anchoring + GFF gene table
# ---------------------------------------------------------------------------
print("[load] AYB anchoring table...")
anc = pd.read_csv(ANCHOR_CSV)
anc_ok = anc.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc_ok["snp_pos_ayb"] = anc_ok["snp_pos_ayb"].astype(int)
print(f"  {len(anc_ok):,} of {len(anc):,} markers AYB-anchored "
      f"({100 * len(anc_ok) / len(anc):.1f}%)")

print("[gff] parsing Funannotate GFF...")
genes = []
attr_re = re.compile(r"(\w+)=([^;]+)")
with open(GFF_PATH) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        seqid = f[0]
        if seqid not in CHR_ORDER:
            continue
        start, end = int(f[3]), int(f[4])
        strand = f[6]
        attrs = dict(attr_re.findall(f[8]))
        gene_id = attrs.get("ID", "")
        name = attrs.get("Name", "")
        note = attrs.get("Note", "") or attrs.get("note", "")
        product = attrs.get("product", "")
        genes.append({"chr": seqid, "start": start, "end": end, "strand": strand,
                      "gene_id": gene_id, "name": name, "note": note,
                      "product": product})
genes = pd.DataFrame(genes)
print(f"  {len(genes):,} gene records on chromosomes Ss01-Ss11")

# Try to recover annotation from mRNA-level rows where the gene row lacks
# product/note (Funannotate sometimes puts annotation only on mRNA features).
print("[gff] backfilling product/note from mRNA records...")
mrna_ann = []
with open(GFF_PATH) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        attrs = dict(attr_re.findall(f[8]))
        parent = attrs.get("Parent", "")
        product = attrs.get("product", "")
        note = attrs.get("note", "") or attrs.get("Note", "")
        dbxref = attrs.get("Dbxref", "")
        if parent and (product or note or dbxref):
            mrna_ann.append({"gene_id": parent, "mrna_product": product,
                             "mrna_note": note, "Dbxref": dbxref})
mrna_df = (pd.DataFrame(mrna_ann)
           .drop_duplicates(subset=["gene_id"], keep="first"))
genes = genes.merge(mrna_df, on="gene_id", how="left")
genes[["product", "note", "mrna_product", "mrna_note", "Dbxref"]] = (
    genes[["product", "note", "mrna_product", "mrna_note", "Dbxref"]].fillna("")
)
genes["_search"] = (genes["name"] + " | " + genes["product"]
                    + " | " + genes["note"]
                    + " | " + genes["mrna_product"]
                    + " | " + genes["mrna_note"]
                    + " | " + genes["Dbxref"])
print(f"  Genes with non-empty annotation: "
      f"{(genes['_search'].str.strip() != '|||||').sum():,}")


# ---------------------------------------------------------------------------
# Load all per-trait GWAS tables; merge AYB anchoring
# ---------------------------------------------------------------------------
def gwas_path(trait):
    if trait in TRAITS_BIOCHEM:
        return GWAS_DIR_BIOCHEM / f"gwas_{trait}_M1_K.csv"
    return GWAS_DIR_NEW / f"gwas_{trait}_M1_K.csv"


gwas = {}
for t in ALL_TRAITS:
    p = gwas_path(t)
    if not p.exists():
        print(f"  {t}: no GWAS table"); continue
    g = pd.read_csv(p)
    g = g.drop(columns=[c for c in ("chrom", "pos") if c in g.columns])
    g = g.merge(anc_ok[["rs", "chr_ayb", "snp_pos_ayb"]],
                left_on="rs", right_on="rs", how="left")
    gwas[t] = g

# anchoring summary
print("[anch] per-trait anchored marker counts:")
summary = []
for t in ALL_TRAITS:
    g = gwas.get(t)
    if g is None: continue
    n_total = len(g)
    n_anc = g["chr_ayb"].notna().sum()
    print(f"  {t:>18}: {n_anc:>5d} / {n_total:>5d} ({100*n_anc/n_total:.1f}%) AYB-anchored")
    summary.append({"trait": t, "n_total_markers": n_total,
                    "n_ayb_anchored": int(n_anc),
                    "pct_anchored": float(100*n_anc/n_total)})
pd.DataFrame(summary).to_csv(TAB / "ayb_anchoring_summary.csv", index=False)


# top-N anchored per trait
top = {}
for t, g in gwas.items():
    anchored = g.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
    anchored = anchored.sort_values("p").head(TOP_N).reset_index(drop=True)
    anchored.to_csv(TAB / f"snp_top{TOP_N}_anchored_{t}.csv", index=False)
    top[t] = anchored


# ---------------------------------------------------------------------------
# +/- WIN_KB scan with trait-specific keywords
# ---------------------------------------------------------------------------
print(f"[hit] +/- {WIN_KB} kb scan per trait...")
all_hits = []
for t, snps in top.items():
    rows = []
    for _, snp in snps.iterrows():
        c = snp["chr_ayb"]
        p = int(snp["snp_pos_ayb"])
        sub = genes[(genes["chr"] == c) &
                    (genes["end"] >= p - WIN_BP) &
                    (genes["start"] <= p + WIN_BP)].copy()
        if sub.empty:
            continue
        sub["snp_rs"] = snp["rs"]
        sub["snp_chr"] = c
        sub["snp_pos"] = p
        sub["snp_p"] = snp["p"]
        sub["snp_beta"] = snp["beta"]
        sub["snp_neglog10p"] = -np.log10(snp["p"])
        def signed_dist(row):
            if row["start"] <= p <= row["end"]:
                return 0
            return (p - row["start"]) if p < row["start"] else (p - row["end"])
        sub["dist_bp"] = sub.apply(signed_dist, axis=1)
        sub["abs_dist_bp"] = sub["dist_bp"].abs()
        sub["pathway_label"] = sub["_search"].apply(lambda s: classify(s, t))
        sub["trait"] = t
        rows.append(sub)
    if rows:
        cand = pd.concat(rows, ignore_index=True)
    else:
        cand = pd.DataFrame()
    if cand.empty:
        print(f"  {t}: no hits"); continue
    cand = cand[["trait", "snp_rs", "snp_chr", "snp_pos", "snp_p",
                 "snp_neglog10p", "snp_beta", "gene_id", "name",
                 "start", "end", "strand", "dist_bp", "abs_dist_bp",
                 "pathway_label", "product", "mrna_product",
                 "note", "mrna_note", "Dbxref"]].sort_values(
        ["snp_p", "abs_dist_bp"])
    cand.to_csv(TAB / f"candidate_genes_{t}.csv", index=False)
    n_path = cand["pathway_label"].notna().sum()
    print(f"  {t:>18}: {len(cand):>4d} gene hits, {n_path:>2d} pathway-flagged")
    all_hits.append(cand[cand["pathway_label"].notna()])

if all_hits:
    pathway_hits = pd.concat(all_hits, ignore_index=True)
else:
    pathway_hits = pd.DataFrame()
pathway_hits = pathway_hits.sort_values(["trait", "snp_p", "abs_dist_bp"])
pathway_hits.to_csv(TAB / "candidate_genes_pathway_hits.csv", index=False)
print(f"[hit] union pathway-flagged hits: {len(pathway_hits)} rows")


# ---------------------------------------------------------------------------
# Annotated Manhattans (AYB chromosome layout)
# ---------------------------------------------------------------------------
print("[fig] annotated Manhattans (AYB)...")
N_MARKERS_TESTED = int(gwas[ALL_TRAITS[0]].shape[0])  # all markers


def manhattan_ayb(trait, ax):
    g = gwas[trait].dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
    if g.empty:
        ax.set_title(f"{trait} — no AYB-anchored SNPs")
        return
    g["pos"] = g["snp_pos_ayb"].astype(int)
    g["chr"] = g["chr_ayb"]
    g = g.sort_values(["chr", "pos"]).reset_index(drop=True)
    chrs = [c for c in CHR_ORDER if c in set(g["chr"])]
    offsets, mids, x_cursor = {}, [], 0
    xs = []
    for c in chrs:
        sub = g[g["chr"] == c]
        offsets[c] = x_cursor
        xs.extend((sub["pos"].values + x_cursor).tolist())
        mids.append(x_cursor + (sub["pos"].max() - sub["pos"].min()) / 2)
        x_cursor += sub["pos"].max() + 5_000_000
    g["x"] = xs
    g["nlp"] = -np.log10(g["p"])
    base_col = TRAIT_PAL.get(trait, WONG["blue"])
    palette = [base_col, "#7c7c7c"]
    for i, c in enumerate(chrs):
        sub = g[g["chr"] == c]
        ax.scatter(sub["x"], sub["nlp"], s=12,
                   color=palette[i % 2], alpha=0.85, edgecolor="none")
    alpha_bonf = 0.05 / len(g)
    bonf = -np.log10(alpha_bonf)
    ax.axhline(bonf, ls="--", color="black", linewidth=0.8,
               label=f"Bonferroni $\\alpha$ = {alpha_bonf:.2e}")

    # Collect one annotation per focal SNP. Cap to single best label per
    # SNP (the trait-keyword "label" + the most-informative gene name) so
    # we don't stack multi-line strings on top of each other.
    sub = pathway_hits[pathway_hits["trait"] == trait]
    grouped = []
    for snp_rs, gsub in sub.groupby("snp_rs", sort=False):
        # pick the single best representative gene = the closest one to
        # the SNP among the pathway-flagged hits for this SNP
        gbest = gsub.sort_values("abs_dist_bp").iloc[0]
        sym = gbest["name"] if gbest["name"] else ""
        full = f"{gbest['pathway_label']} ({sym})" if sym else gbest["pathway_label"]
        grouped.append({
            "snp_chr": gsub["snp_chr"].iloc[0],
            "snp_pos": int(gsub["snp_pos"].iloc[0]),
            "snp_nlp": gsub["snp_neglog10p"].iloc[0],
            "label": full,
        })

    # Rotate labels 45 deg so they take less horizontal space, then place in
    # a band above the data and repel only as much as needed. Labels are
    # constrained to stay within the data x-range so leader lines never run
    # to the figure edge.
    if grouped:
        bonf_local = -np.log10(0.05 / len(g))
        band_y = max(bonf_local + 0.6, g["nlp"].max() + 1.5)
        grouped = sorted(grouped, key=lambda d: d["snp_pos"]
                         + offsets.get(d["snp_chr"], 0))
        x_init = np.array([d["snp_pos"] + offsets.get(d["snp_chr"], 0)
                           for d in grouped], dtype=float)
        x_min = min(offsets.values())
        x_max = x_cursor
        # at 45 deg the horizontal projection of a label of n chars is
        # roughly n * 0.7 * char_w; estimate char_w as ~0.4 % of axis span
        # (matplotlib at fontsize 7, dpi 300, with 11 in wide figure).
        ax_span = x_max - x_min
        char_w = ax_span * 0.0045
        widths = np.array([len(d["label"]) * char_w * 0.7
                           for d in grouped])
        half_widths = widths / 2.0
        x_lab = x_init.copy()
        # iterative repel; bound to data range
        for _ in range(200):
            moved = False
            for i in range(len(x_lab) - 1):
                gap = (x_lab[i + 1] - x_lab[i]) - (half_widths[i]
                                                   + half_widths[i + 1])
                if gap < 0:
                    push = (-gap) / 2 + char_w * 0.3
                    x_lab[i] -= push
                    x_lab[i + 1] += push
                    moved = True
            x_lab = np.clip(x_lab, x_min, x_max)
            if not moved:
                break
        for i, h in enumerate(grouped):
            x_snp = h["snp_pos"] + offsets.get(h["snp_chr"], 0)
            ax.annotate(h["label"],
                        xy=(x_snp, h["snp_nlp"]),
                        xytext=(x_lab[i], band_y),
                        textcoords="data",
                        rotation=45, rotation_mode="anchor",
                        fontsize=7, ha="left", va="bottom",
                        color="black",
                        arrowprops=dict(arrowstyle="-",
                                        connectionstyle="arc3,rad=0",
                                        color="black", lw=0.5, alpha=0.6))
    ax.set_xticks(mids)
    ax.set_xticklabels(chrs, rotation=0, fontsize=8)
    ax.set_ylabel(r"$-\log_{10}\,p$")
    ax.set_title(f"{trait}   AYB-anchored ({len(g):,} SNPs); "
                 f"top {TOP_N} scanned +/- {WIN_KB} kb", fontsize=10)
    # Headroom: 45-degree rotated labels rise above band_y; longest label
    # is ~25 chars which adds approximately the same vertical extent as
    # horizontal so we add ~2 to band_y.
    headroom = (band_y + 2.0) if grouped else (bonf + 0.5)
    ax.set_ylim(0, max(headroom, g["nlp"].max() + 2.0))
    ax.grid(True, axis="y")
    ax.legend(loc="upper left", fontsize=8)


for t in ALL_TRAITS:
    if t not in gwas: continue
    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    manhattan_ayb(t, ax)
    save(fig, f"fig46_manhattan_ayb_{t}")
    print(f"  fig46_manhattan_ayb_{t}")


# ---------------------------------------------------------------------------
# Composite multi-panel Manhattans (replacements for legacy fig41/fig42 which
# were cowpea-anchored). Each composite stacks the per-trait AYB Manhattans
# in a single figure so reviewers can read the seed-trait or protein/oxalate
# story in a single view rather than as separate supplementary figures.
# ---------------------------------------------------------------------------
SEED_PANEL_ORDER = ["Mass_of_Seeds", "Seed_Length", "Seed_Width",
                    "Seed_Thickness", "Seed_Coat_Tannin"]
PROTOX_PANEL_ORDER = ["Crude_Protein", "Soluble_Oxalate",
                      "Insoluble_Oxalate", "Total_Oxalate"]


def composite_manhattan(traits, fname, panel_letter_offset=0):
    traits_present = [t for t in traits if t in gwas]
    n = len(traits_present)
    if n == 0:
        print(f"  [skip] {fname}: no traits with data")
        return
    # Taller per-panel allocation + extra hspace so 45-degree rotated
    # candidate-gene callout labels stay within their own panel and don't
    # intrude into the panel above.
    fig, axes = plt.subplots(n, 1, figsize=(11, 3.8 * n),
                              sharex=False)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.97, bottom=0.05,
                         hspace=0.95)
    if n == 1:
        axes = [axes]
    for i, (ax, t) in enumerate(zip(axes, traits_present)):
        manhattan_ayb(t, ax)
        letter = chr(ord("A") + panel_letter_offset + i)
        # Replace the per-trait title with one that carries the panel letter.
        current = ax.get_title()
        ax.set_title(f"({letter}) {current}", fontsize=10, loc="center")
    save(fig, fname)
    print(f"  {fname}: {n}-panel composite")


print("[fig] composite AYB-anchored Manhattans...")
composite_manhattan(SEED_PANEL_ORDER,   "fig41_manhattan_panels_seed_ayb")
composite_manhattan(PROTOX_PANEL_ORDER, "fig42_manhattan_panels_protox_ayb")


print(f"\nOutputs in: {OUT}")
print(f"  union pathway-flagged hits: {len(pathway_hits)} rows")
