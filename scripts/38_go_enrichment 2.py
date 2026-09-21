#!/usr/bin/env python3
"""
Functional enrichment (GO terms + Pfam domains) for genes flanking the
top-suggestive SNPs in the AYB GWAS scans.

Two refinements over the first-pass script:
  1. Tighter foreground: top-10 SNPs per trait (vs top-30) within +/- 10 kb
     (vs +/- 50 kb). Reduces foreground from ~2.4 k genes to ~300, raising
     signal-to-noise for the Fisher's-exact tests.
  2. GO term names from the current OBO file (downloaded once to
     refs/go/go.obo) so the figure labels are human-readable.
  3. Three independent enrichment runs split by GO namespace
     (biological_process, molecular_function, cellular_component) plus a
     Pfam-domain run, each with its own BH-FDR.

Outputs (results/38_go_enrichment/)
-----------------------------------
tables/
    go_enrichment_full.csv     all GO terms tested, with names + namespace
    go_enrichment_top.csv      BH q < 0.10 set
    pfam_enrichment_full.csv   Pfam domain enrichment
figures/
    fig75_go_enrichment_bubble.png/.pdf      top GO terms with names
    fig76_pfam_enrichment_bubble.png/.pdf    top Pfam domains

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control

from _plotstyle import apply, WONG
from _pheno import TRAITS_BIOCHEM, TRAITS_SEED, TRAITS_PROT_OX, TRAITS_ALL
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
GFF = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
OBO = ROOT / "refs" / "go" / "go.obo"
GWAS_BIOCHEM = ROOT / "results" / "03_gwas_mlm" / "tables"
GWAS_NEW = ROOT / "results" / "15_gwas_new_traits" / "tables"
OUT = ROOT / "results" / "38_go_enrichment"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

TOP_N_PER_TRAIT = 10
WIN_KB = 10
WIN_BP = WIN_KB * 1_000

NAMESPACE_LABEL = {"biological_process": "BP",
                   "molecular_function": "MF",
                   "cellular_component": "CC"}
NAMESPACE_COLOR = {"biological_process": WONG["blue"],
                   "molecular_function": WONG["vermillion"],
                   "cellular_component": WONG["green"]}


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Parse OBO for GO term names + namespace
# ---------------------------------------------------------------------------
go_meta = {}
if OBO.exists():
    print("[obo] parsing GO ontology...")
    cur = {}
    with open(OBO) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                if cur.get("id") and cur.get("name"):
                    go_meta[cur["id"]] = {"name": cur["name"],
                                          "namespace": cur.get("namespace", "")}
                cur = {}
            elif line.startswith("id: "):
                cur["id"] = line[4:]
            elif line.startswith("name: "):
                cur["name"] = line[6:]
            elif line.startswith("namespace: "):
                cur["namespace"] = line[11:]
        if cur.get("id") and cur.get("name"):
            go_meta[cur["id"]] = {"name": cur["name"],
                                  "namespace": cur.get("namespace", "")}
    print(f"[obo] loaded {len(go_meta):,} term records")
else:
    print("[obo] OBO file missing; GO IDs will be reported without names")


# ---------------------------------------------------------------------------
# Parse Funannotate GFF for gene coordinates + per-gene GO and Pfam
# ---------------------------------------------------------------------------
print("[gff] parsing Funannotate gene table...")
attr_re = re.compile(r"(\w+)=([^;]+)")
genes = {}
mrna_go = {}
mrna_pfam = {}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] == "gene":
            attrs = dict(attr_re.findall(f[8]))
            gid = attrs.get("ID", "")
            genes[gid] = {"chr": f[0], "start": int(f[3]), "end": int(f[4]),
                          "name": attrs.get("Name", "")}
        elif f[2] == "mRNA":
            attrs = dict(attr_re.findall(f[8]))
            pid = attrs.get("Parent", "")
            ont = attrs.get("Ontology_term", "")
            dbx = attrs.get("Dbxref", "")
            if pid:
                if ont:
                    gos = [g.strip() for g in ont.split(",") if g.startswith("GO:")]
                    if gos:
                        mrna_go.setdefault(pid, set()).update(gos)
                # Pfam IDs come through Dbxref like "PFAM:PF18029,InterPro:IPR029068"
                if dbx:
                    pfs = [d.split(":")[1] for d in dbx.split(",")
                           if d.startswith("PFAM:")]
                    if pfs:
                        mrna_pfam.setdefault(pid, set()).update(pfs)

gene_go   = {gid: sorted(mrna_go.get(gid, set())) for gid in genes}
gene_pfam = {gid: sorted(mrna_pfam.get(gid, set())) for gid in genes}
n_with_go = sum(1 for v in gene_go.values() if v)
n_with_pf = sum(1 for v in gene_pfam.values() if v)
print(f"  {len(genes):,} genes; {n_with_go:,} with >=1 GO; "
      f"{n_with_pf:,} with >=1 Pfam")


# ---------------------------------------------------------------------------
# Build foreground gene set: top-N SNPs per trait, +/- WIN_KB
# ---------------------------------------------------------------------------
print(f"[fg] foreground = +/- {WIN_KB} kb of top-{TOP_N_PER_TRAIT} SNPs per trait...")
anc = pd.read_csv(ANCHOR)
anc = anc.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc["snp_pos_ayb"] = anc["snp_pos_ayb"].astype(int)
anc_lookup = anc.set_index("rs")[["chr_ayb", "snp_pos_ayb"]].to_dict("index")

fg_gene_ids = set()
for trait in TRAITS_ALL:
    if trait in TRAITS_BIOCHEM:
        gp = GWAS_BIOCHEM / f"gwas_{trait}_M1_K.csv"
    else:
        gp = GWAS_NEW / f"gwas_{trait}_M1_K.csv"
    if not gp.exists(): continue
    g = pd.read_csv(gp)
    g["chr_ayb"] = g["rs"].map(lambda r: anc_lookup.get(r, {}).get("chr_ayb"))
    g["snp_pos_ayb"] = g["rs"].map(lambda r: anc_lookup.get(r, {}).get("snp_pos_ayb"))
    top = (g.dropna(subset=["chr_ayb"])
           .sort_values("p").head(TOP_N_PER_TRAIT))
    for _, snp in top.iterrows():
        c = snp["chr_ayb"]; p = int(snp["snp_pos_ayb"])
        for gid, gd in genes.items():
            if gd["chr"] == c and gd["end"] >= p - WIN_BP and gd["start"] <= p + WIN_BP:
                fg_gene_ids.add(gid)
print(f"[fg] foreground gene set: {len(fg_gene_ids):,} unique genes")


# ---------------------------------------------------------------------------
# Per-term Fisher's exact + BH-FDR — one run per annotation table
# ---------------------------------------------------------------------------
def enrichment(annot_map, fg_ids, label):
    """One-sided Fisher's exact test for over-representation in foreground.

    annot_map : dict gene_id -> list of term IDs
    """
    fg_terms = Counter()
    bg_terms = Counter()
    for gid, ts in annot_map.items():
        for t in ts:
            bg_terms[t] += 1
            if gid in fg_ids:
                fg_terms[t] += 1
    n_fg = sum(1 for gid in fg_ids if annot_map.get(gid))
    N_bg = sum(1 for v in annot_map.values() if v)
    rows = []
    for term, K in bg_terms.items():
        k = fg_terms.get(term, 0)
        if k == 0: continue
        a = k
        b = n_fg - k
        c = K - k
        d = N_bg - n_fg - K + k
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        fold = ((k / n_fg) / (K / N_bg)) if N_bg > 0 else np.nan
        rows.append({"term": term, "k_in_fg": k, "K_in_bg": K,
                     "n_fg": n_fg, "N_bg": N_bg,
                     "fold_enrichment": float(fold), "p_value": float(p)})
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_bh"] = false_discovery_control(df["p_value"])
        df = df.sort_values("p_value").reset_index(drop=True)
    print(f"  {label}: n_fg = {n_fg}, N_bg = {N_bg}, "
          f"{len(df):,} terms tested, "
          f"{(df['q_bh'] < 0.10).sum() if 'q_bh' in df else 0} BH q < 0.10")
    return df


print("\n[enrich] running GO enrichment (split by namespace)...")
go_df = enrichment(gene_go, fg_gene_ids, "GO all")
# annotate with names + namespace
go_df["name"] = go_df["term"].map(lambda t: go_meta.get(t, {}).get("name", ""))
go_df["namespace"] = go_df["term"].map(
    lambda t: go_meta.get(t, {}).get("namespace", ""))
# per-namespace FDR (more honest than pooled FDR)
go_df["q_bh_within_namespace"] = np.nan
for ns in go_df["namespace"].unique():
    if not ns: continue
    mask = go_df["namespace"] == ns
    go_df.loc[mask, "q_bh_within_namespace"] = (
        false_discovery_control(go_df.loc[mask, "p_value"].values))
go_df = go_df[["term", "name", "namespace", "k_in_fg", "K_in_bg",
               "n_fg", "N_bg", "fold_enrichment", "p_value",
               "q_bh", "q_bh_within_namespace"]]
go_df.to_csv(TAB / "go_enrichment_full.csv", index=False)

go_top = go_df[(go_df["q_bh_within_namespace"] < 0.10)].copy()
go_top.to_csv(TAB / "go_enrichment_top.csv", index=False)
print(f"[enrich] GO terms passing q_bh_within_namespace < 0.10: "
      f"{len(go_top)} (split by BP/MF/CC)")
for ns in ("biological_process", "molecular_function", "cellular_component"):
    n_sig = ((go_df["namespace"] == ns) &
             (go_df["q_bh_within_namespace"] < 0.10)).sum()
    print(f"    {NAMESPACE_LABEL[ns]}: {n_sig} significant terms")


print("\n[enrich] running Pfam-domain enrichment...")
pf_df = enrichment(gene_pfam, fg_gene_ids, "Pfam")
pf_df = pf_df.rename(columns={"term": "Pfam_id"})
pf_df.to_csv(TAB / "pfam_enrichment_full.csv", index=False)
pf_top = pf_df[pf_df["q_bh"] < 0.10] if "q_bh" in pf_df else pd.DataFrame()
print(f"[enrich] Pfam domains passing BH q < 0.10: {len(pf_top)}")

print("\n[top GO]:")
print(go_df.head(15)[["term", "name", "namespace",
                       "k_in_fg", "K_in_bg",
                       "fold_enrichment", "p_value",
                       "q_bh_within_namespace"]].to_string(index=False))
print("\n[top Pfam]:")
print(pf_df.head(15)[["Pfam_id", "k_in_fg", "K_in_bg",
                       "fold_enrichment", "p_value", "q_bh"]].to_string(index=False))


# ---------------------------------------------------------------------------
# Fig 75. GO enrichment bubble plot, split by namespace, with names
# ---------------------------------------------------------------------------
import matplotlib.patches as mpatches

# Take top-8 per namespace by p
plot_blocks = []
for ns in ("biological_process", "molecular_function", "cellular_component"):
    sub = go_df[go_df["namespace"] == ns].head(8)
    plot_blocks.append(sub)
plot_df = pd.concat(plot_blocks, ignore_index=True)
plot_df = plot_df.dropna(subset=["fold_enrichment", "p_value"])

fig, ax = plt.subplots(figsize=(11, 7.5), constrained_layout=True)
for ns, sub in plot_df.groupby("namespace"):
    color = NAMESPACE_COLOR.get(ns, WONG["blue"])
    ax.scatter(sub["fold_enrichment"], -np.log10(sub["p_value"]),
               s=sub["k_in_fg"] * 50,
               c=color, alpha=0.7, edgecolor="black", linewidth=0.4,
               label=f"{NAMESPACE_LABEL[ns]} (n = {len(sub)})")
# labels
for _, row in plot_df.iterrows():
    label = row["name"] if row["name"] else row["term"]
    if len(label) > 40: label = label[:37] + "..."
    ax.text(row["fold_enrichment"] * 1.04,
            -np.log10(row["p_value"]),
            label, fontsize=7,
            va="center", ha="left", color="black")
ax.set_xlabel("Fold enrichment (foreground vs background)")
ax.set_ylabel(r"$-\log_{10}\,p$ (Fisher's exact, one-sided)")
ax.set_title(f"GO enrichment top-8 per namespace "
             f"(foreground: top-{TOP_N_PER_TRAIT} SNPs × 13 traits, "
             f"+/- {WIN_KB} kb)")
# add per-namespace q < 0.10 threshold lines if any term passes
for ns in ("biological_process", "molecular_function", "cellular_component"):
    sub_full = go_df[(go_df["namespace"] == ns) &
                     (go_df["q_bh_within_namespace"] < 0.10)]
    if not sub_full.empty:
        thresh_p = sub_full["p_value"].max()
        ax.axhline(-np.log10(thresh_p),
                   color=NAMESPACE_COLOR[ns], linewidth=0.7,
                   linestyle="--", alpha=0.5,
                   label=f"{NAMESPACE_LABEL[ns]} q < 0.10")
ax.legend(loc="lower right", fontsize=8)
ax.grid(True, alpha=0.4)
ax.set_xscale("log")
save(fig, "fig75_go_enrichment_bubble")
print("[fig] fig75_go_enrichment_bubble")


# ---------------------------------------------------------------------------
# Fig 76. Pfam-domain enrichment bubble
# ---------------------------------------------------------------------------
pf_show = pf_df.head(20).copy()
if not pf_show.empty:
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    ax.scatter(pf_show["fold_enrichment"], -np.log10(pf_show["p_value"]),
               s=pf_show["k_in_fg"] * 50,
               c=WONG["blue"], alpha=0.7, edgecolor="black", linewidth=0.4)
    for _, row in pf_show.iterrows():
        ax.text(row["fold_enrichment"] * 1.04,
                -np.log10(row["p_value"]),
                row["Pfam_id"], fontsize=7, va="center", color="black")
    ax.set_xlabel("Fold enrichment (foreground vs background)")
    ax.set_ylabel(r"$-\log_{10}\,p$ (Fisher's exact, one-sided)")
    ax.set_title(f"Pfam domain enrichment, top 20 "
                 f"(foreground: top-{TOP_N_PER_TRAIT} SNPs × 13 traits, "
                 f"+/- {WIN_KB} kb)")
    if (pf_df["q_bh"] < 0.10).any():
        thresh_p = pf_df.loc[pf_df["q_bh"] < 0.10, "p_value"].max()
        ax.axhline(-np.log10(thresh_p), color=WONG["vermillion"],
                   linestyle="--", linewidth=0.7,
                   label=f"BH q < 0.10")
        ax.legend()
    ax.grid(True, alpha=0.4)
    ax.set_xscale("log")
    save(fig, "fig76_pfam_enrichment_bubble")
    print("[fig] fig76_pfam_enrichment_bubble")

print(f"\nOutputs in: {OUT}")
