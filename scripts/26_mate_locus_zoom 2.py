#!/usr/bin/env python3
"""
Locus-zoom: LD haplotype block + gene neighbourhood around the strongest
oxalate-transport hits in the AYB panel.

Configurable via the TARGETS list: each entry picks the lowest-p AYB-anchored
SNP on a given chromosome for a given trait, builds an LD-r^2 heatmap around
it, and overlays the Funannotate gene track + a -log10 p Manhattan strip.

Default targets:
  Soluble_Oxalate  Ss10  -> ALMT4 locus (best AYB-anchored hit, p = 0.009)
  Insoluble_Oxalate Ss04 -> MATE-cluster locus (cowpea orthologue map; DTX
                            paralogs in the same region)

Outputs (results/26_mate_locus_zoom/)
-------------------------------------
tables/
    zoom_<label>_snps.csv
    zoom_<label>_ld_matrix.csv
    zoom_<label>_genes.csv
figures/
    fig63_zoom_<label>.png/.pdf

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
HAPMAP = ROOT / "AYB_SNP_Result_Report-DAf18-2580" / "Report_DAf18-2580_SNP_HapMap.csv"
ANCHOR = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
GFF = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
GWAS_BIOCHEM_DIR = ROOT / "results" / "03_gwas_mlm" / "tables"
GWAS_NEW_DIR = ROOT / "results" / "15_gwas_new_traits" / "tables"
OUT = ROOT / "results" / "26_mate_locus_zoom"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

SAMP_CR = 0.90
MARK_CR = 0.90
MIN_MAF = 0.05
WIN_KB = 500
WIN_BP = WIN_KB * 1_000

TARGETS = [
    # Soluble_Oxalate ALMT4 locus -- top AYB-anchored hit, p = 0.009, n = 41.
    # Window expanded to +/- 500 kb to span the ALMT4 gene + neighbours.
    {"trait": "Soluble_Oxalate", "chrom": "Ss10",
     "label": "ALMT4_locus",
     "focal_rs": "100033542|F|0-31:T>C-31:T>C", "win_bp": 500_000},
    # Ss04 DETOXIFICATION / MATE-family cluster — picked at the rank-26
    # Insoluble_Oxalate SNP that the AYB-anchored top-30 scan flagged for
    # four tandem DETOX-family genes (cowpea ortholog of Vigun04g028900-
    # 029300 MATE cluster). Window widened to +/- 1 Mb because Ss04 has
    # sparse anchored markers in this region.
    {"trait": "Insoluble_Oxalate", "chrom": "Ss04",
     "label": "MATE_DETOX_locus",
     "focal_rs": "100017924|F|0-5:G>C-5:G>C", "win_bp": 1_000_000},
]

HIGHLIGHT_PATTERNS = [
    (re.compile(r"\bDETOXIFICATION\b|\bMATE\b", re.IGNORECASE), "MATE"),
    (re.compile(r"\baluminum[- ]?activated malate transporter\b|\bALMT[0-9]*\b",
                re.IGNORECASE), "ALMT"),
    (re.compile(r"\bformate dehydrogenase\b", re.IGNORECASE), "FDH"),
    (re.compile(r"\boxalate\b|\bglyoxylate\b", re.IGNORECASE), "OXO"),
    (re.compile(r"\bAINTEGUMENTA\b|\bintegument\b", re.IGNORECASE), "ANT"),
    (re.compile(r"\bchalcone\b|\bisoflavone\b", re.IGNORECASE), "phenyl"),
    (re.compile(r"\bMyb\w*\b", re.IGNORECASE), "MYB"),
    (re.compile(r"\bpectin methylesterase\b", re.IGNORECASE), "PME"),
    (re.compile(r"\bexpansin\b", re.IGNORECASE), "expansin"),
    (re.compile(r"\bauxin\b|\bArf\b", re.IGNORECASE), "ARF"),
]
HL_COLORS = {"MATE": WONG["vermillion"], "ALMT": WONG["blue"],
             "FDH": WONG["orange"], "OXO": "#d95f02",
             "ANT": WONG["green"], "phenyl": WONG["purple"],
             "MYB": "#666666", "PME": "#b15928",
             "expansin": "#33a02c", "ARF": "#e7298a"}


def label_gene(row):
    blob = f"{row['name']} | {row['product']} | {row['note']}"
    for pat, lab in HIGHLIGHT_PATTERNS:
        if pat.search(blob):
            return lab
    return None


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Load anchoring + filtered dosage
# ---------------------------------------------------------------------------
print("[load] HapMap + filtered dosage...")
anc = pd.read_csv(ANCHOR)
anc_ok = anc.dropna(subset=["chr_ayb", "snp_pos_ayb"]).copy()
anc_ok["snp_pos_ayb"] = anc_ok["snp_pos_ayb"].astype(int)

hm = pd.read_csv(HAPMAP, low_memory=False)
META = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode"]
sample_cols = [c for c in hm.columns if c not in META]
calls = hm[sample_cols].astype(str)
ref_alt = hm["alleles"].str.split("/", expand=True)
ref_alt.columns = ["ref", "alt"]
dosage = np.full((len(hm), len(sample_cols)), np.nan)
for i in range(len(hm)):
    ref, alt = ref_alt.iloc[i]
    arr = calls.iloc[i].values
    dosage[i, arr == ref + ref] = 0.0
    dosage[i, (arr == ref + alt) | (arr == alt + ref)] = 1.0
    dosage[i, arr == alt + alt] = 2.0
dosage = pd.DataFrame(dosage, index=hm["rs#"].values, columns=sample_cols)
call_rate_m = dosage.notna().sum(axis=1) / dosage.shape[1]
maf = np.minimum(dosage.mean(axis=1, skipna=True) / 2.0,
                 1.0 - dosage.mean(axis=1, skipna=True) / 2.0)
call_rate_s = dosage.notna().sum(axis=0) / dosage.shape[0]
keep_m = ((call_rate_m >= MARK_CR) & (maf >= MIN_MAF)).values
keep_s = (call_rate_s >= SAMP_CR).values
dos_filt = dosage.loc[keep_m, keep_s].copy()
# per-marker mean impute
row_mean = dos_filt.mean(axis=1, skipna=True)
for r in dos_filt.index:
    row = dos_filt.loc[r]
    if row.isna().any():
        dos_filt.loc[r] = row.fillna(row_mean.loc[r])
samples = dos_filt.columns.tolist()
print(f"[load] dosage = {dos_filt.shape}")


# ---------------------------------------------------------------------------
# Parse GFF once (gene + mRNA tables)
# ---------------------------------------------------------------------------
print("[gff] parsing Funannotate GFF (Ss01-Ss11)...")
attr_re = re.compile(r"(\w+)=([^;]+)")
gene_rows = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        if not re.match(r"^Ss\d+$", f[0]): continue
        attrs = dict(attr_re.findall(f[8]))
        gene_rows.append({"chr": f[0], "start": int(f[3]), "end": int(f[4]),
                          "strand": f[6], "gene_id": attrs.get("ID", ""),
                          "name": attrs.get("Name", "")})
mrna_rows = {}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "mRNA": continue
        if not re.match(r"^Ss\d+$", f[0]): continue
        attrs = dict(attr_re.findall(f[8]))
        pid = attrs.get("Parent", "")
        if pid and pid not in mrna_rows:
            mrna_rows[pid] = {"product": attrs.get("product", ""),
                              "note": attrs.get("note", "")}
all_genes = pd.DataFrame(gene_rows)
all_genes["product"] = all_genes["gene_id"].map(
    lambda x: mrna_rows.get(x, {}).get("product", ""))
all_genes["note"] = all_genes["gene_id"].map(
    lambda x: mrna_rows.get(x, {}).get("note", ""))
all_genes["highlight"] = all_genes.apply(label_gene, axis=1)
print(f"[gff] {len(all_genes)} genes total, "
      f"{all_genes['highlight'].notna().sum()} pathway-flagged genome-wide")


# ---------------------------------------------------------------------------
# Per-target zoom
# ---------------------------------------------------------------------------
def pick_focal(trait, chrom, focal_rs=None):
    if trait in ("Tannin", "Phenol", "Flavonoid", "Antioxidant"):
        gp = GWAS_BIOCHEM_DIR / f"gwas_{trait}_M1_K.csv"
    else:
        gp = GWAS_NEW_DIR / f"gwas_{trait}_M1_K.csv"
    g_ = pd.read_csv(gp)
    g_ = g_.drop(columns=[c for c in ("chrom", "pos") if c in g_.columns])
    g_ = g_.merge(anc_ok[["rs", "chr_ayb", "snp_pos_ayb"]],
                  on="rs", how="left")
    if focal_rs is not None:
        sub = g_[g_["rs"] == focal_rs]
        if sub.empty:
            # Requested SNP not in current GWAS (e.g. dropped by QC
            # re-run). Fall back to the lowest-p anchored SNP on chrom
            # within +/- 200 kb of the originally-published focal.
            anc_focal = anc_ok[anc_ok["rs"] == focal_rs]
            if anc_focal.empty:
                return None
            anchor_pos = int(anc_focal.iloc[0]["snp_pos_ayb"])
            sub = (g_.dropna(subset=["chr_ayb"])
                   .query("chr_ayb == @chrom")
                   .query("snp_pos_ayb >= @anchor_pos - 200_000")
                   .query("snp_pos_ayb <= @anchor_pos + 200_000")
                   .sort_values("p").head(1))
            if sub.empty:
                return None
            print(f"  [fallback] requested focal {focal_rs} not in current "
                  f"GWAS -- using nearest anchored hit on {chrom}")
    else:
        sub = (g_.dropna(subset=["chr_ayb"])
               .query("chr_ayb == @chrom")
               .sort_values("p").head(1))
        if sub.empty: return None
    r = sub.iloc[0]
    return g_, str(r["rs"]), int(r["snp_pos_ayb"]), float(r["p"])


for target in TARGETS:
    trait = target["trait"]; chrom = target["chrom"]; label = target["label"]
    focal_rs_req = target.get("focal_rs")
    win_bp = int(target.get("win_bp", WIN_BP))
    win_kb = win_bp // 1000
    print(f"\n========== {label}: {trait} on {chrom} (+/- {win_kb} kb) ==========")
    picked = pick_focal(trait, chrom, focal_rs=focal_rs_req)
    if picked is None:
        print(f"  no anchored SNPs; skip"); continue
    g, focal_rs, focal_pos, focal_p = picked
    print(f"  focal: {focal_rs}  {chrom}:{focal_pos:,}  p = {focal_p:.4f}")

    rs_in_gwas = set(g["rs"].astype(str))
    win_snps = anc_ok[(anc_ok["chr_ayb"] == chrom) &
                      (anc_ok["snp_pos_ayb"] >= focal_pos - win_bp) &
                      (anc_ok["snp_pos_ayb"] <= focal_pos + win_bp) &
                      (anc_ok["rs"].isin(rs_in_gwas))].copy()
    win_snps = win_snps.merge(dos_filt.reset_index().rename(
        columns={"index": "rs"}), on="rs", how="inner")
    win_snps = win_snps.sort_values("snp_pos_ayb").reset_index(drop=True)
    print(f"  {len(win_snps)} markers in +/- {win_kb} kb window")
    if len(win_snps) < 2:
        print("  too few markers for an LD heatmap"); continue
    g_sub = g.set_index("rs").loc[win_snps["rs"], ["p", "beta", "se", "fdr_bh"]].reset_index()
    out_snps = win_snps[["rs", "chr_ayb", "snp_pos_ayb",
                          "pident", "evalue", "bitscore"]].merge(
        g_sub, on="rs", how="left")
    out_snps.to_csv(TAB / f"zoom_{label}_snps.csv", index=False)

    # LD r^2
    X = win_snps[samples].values.astype(float)
    Xz = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True,
                                                     ddof=1)
    Xz = np.nan_to_num(Xz, nan=0.0)
    n_samp = len(samples)
    LD = (Xz @ Xz.T) / (n_samp - 1)
    LD = LD ** 2
    np.fill_diagonal(LD, 1.0)
    LD_df = pd.DataFrame(LD, index=win_snps["rs"], columns=win_snps["rs"])
    LD_df.to_csv(TAB / f"zoom_{label}_ld_matrix.csv")

    # Genes in window
    win_genes = all_genes[(all_genes["chr"] == chrom) &
                          (all_genes["end"] >= focal_pos - win_bp) &
                          (all_genes["start"] <= focal_pos + win_bp)].copy()
    win_genes = win_genes.sort_values("start").reset_index(drop=True)
    win_genes.to_csv(TAB / f"zoom_{label}_genes.csv", index=False)
    n_path = win_genes["highlight"].notna().sum()
    print(f"  {len(win_genes)} genes in window, {n_path} pathway-flagged")

    # Composite figure
    fig = plt.figure(figsize=(12, 8), constrained_layout=False)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.7, 3.4], hspace=0.05)
    ax_top = fig.add_subplot(gs[0, 0])
    positions = win_snps["snp_pos_ayb"].values
    p_vals = g.set_index("rs").loc[win_snps["rs"], "p"].values
    nlp = -np.log10(p_vals)
    focal_idx = int(np.where(win_snps["rs"].values == focal_rs)[0][0])
    colors = ["#7c7c7c"] * len(win_snps)
    colors[focal_idx] = WONG["vermillion"]
    ax_top.scatter(positions / 1e6, nlp, s=22, color=colors,
                   edgecolor="black", linewidth=0.3, zorder=3)
    # Pin the focal-SNP label to the top-left corner of the panel so it
    # never collides with the legend or the title text. Arrow points to the
    # actual focal SNP location.
    ax_top.annotate(f"focal SNP\n{focal_rs}\np = {focal_p:.4f}",
                    xy=(focal_pos / 1e6, nlp[focal_idx]),
                    xytext=(0.02, 0.92), textcoords="axes fraction",
                    ha="left", va="top", fontsize=7.5,
                    color=WONG["vermillion"],
                    arrowprops=dict(arrowstyle="-",
                                    color=WONG["vermillion"], lw=0.6,
                                    alpha=0.85,
                                    connectionstyle="arc3,rad=-0.05"))
    ax_top.set_ylabel(r"$-\log_{10}\,p$" "\n" f"({trait})")
    ax_top.set_xlim((focal_pos - win_bp) / 1e6,
                    (focal_pos + win_bp) / 1e6)
    ax_top.set_xticklabels([])
    ax_top.grid(True, axis="y", alpha=0.4)
    ax_top.set_title(f"{label}: {trait} on {chrom} — "
                     f"+/- {win_kb} kb around focal "
                     f"({chrom}:{focal_pos:,}, p = {focal_p:.4f})")

    ax_mid = fig.add_subplot(gs[1, 0])
    for _, gd in win_genes.iterrows():
        color = HL_COLORS.get(gd["highlight"], "#cccccc")
        ax_mid.add_patch(plt.Rectangle(
            (gd["start"] / 1e6, 0.3),
            (gd["end"] - gd["start"]) / 1e6, 0.4,
            facecolor=color, edgecolor="black", linewidth=0.3,
            alpha=0.85 if gd["highlight"] else 0.45))
        if gd["highlight"]:
            ax_mid.annotate(gd["highlight"],
                            xy=((gd["start"] + gd["end"]) / 2 / 1e6, 0.95),
                            ha="center", fontsize=6, color=color)
    ax_mid.set_xlim((focal_pos - win_bp) / 1e6,
                    (focal_pos + win_bp) / 1e6)
    ax_mid.set_ylim(0, 1.1)
    ax_mid.set_yticks([])
    ax_mid.set_xticklabels([])
    ax_mid.set_ylabel("Funannotate\ngenes")
    ax_mid.axvline(focal_pos / 1e6, color=WONG["vermillion"],
                   linewidth=0.8, alpha=0.6)

    # Triangular (Haploview-style) LD heatmap: each upper-triangle pair (i, j)
    # is plotted as a diamond at center x = midpoint(x_i, x_j), y = half the
    # physical distance between them. The result is a 45-degree-rotated
    # triangle pointing up toward the gene track. Colour = r^2.
    ax_ld = fig.add_subplot(gs[2, 0])
    import matplotlib.collections as mcoll
    n_w = len(positions)
    pos_mb = positions / 1e6
    diamonds = []
    colors = []
    # mean SNP spacing for diamond size (so they tile rather than overlap)
    if n_w > 1:
        spacing = np.diff(pos_mb).mean()
    else:
        spacing = 1.0
    half = spacing / 2.0
    for i in range(n_w):
        for j in range(i + 1, n_w):
            cx = (pos_mb[i] + pos_mb[j]) / 2.0
            cy = (pos_mb[j] - pos_mb[i]) / 2.0
            # diamond verts: top, right, bottom, left
            verts = [(cx, cy + half), (cx + half, cy),
                     (cx, cy - half), (cx - half, cy)]
            diamonds.append(verts)
            colors.append(LD[i, j])
    coll = mcoll.PolyCollection(diamonds, array=np.asarray(colors),
                                 cmap="Reds",
                                 edgecolors="white", linewidths=0.15)
    coll.set_clim(0, 1)
    ax_ld.add_collection(coll)
    ax_ld.set_xlim((focal_pos - win_bp) / 1e6,
                   (focal_pos + win_bp) / 1e6)
    # triangle height = max half-distance between any pair in the window
    pair_max_d = (pos_mb.max() - pos_mb.min()) / 2.0 if n_w > 1 else 0.5
    ax_ld.set_ylim(0, pair_max_d + half)
    ax_ld.set_xlabel(f"{chrom} position (Mb)")
    ax_ld.set_ylabel("pair separation / 2 (Mb)")
    ax_ld.axvline(focal_pos / 1e6, color=WONG["vermillion"],
                  linewidth=0.6, alpha=0.5)
    # SNP-position ticks along the bottom
    ax_ld.scatter(pos_mb, np.zeros_like(pos_mb), s=12, marker="v",
                  color="black", zorder=4)
    cb_ax = fig.add_axes([0.92, 0.08, 0.015, 0.35])
    fig.colorbar(coll, cax=cb_ax, label="r$^2$")

    # Gene-category legend anchored beneath the LD heatmap's x-axis so it
    # cannot collide with either the focal-SNP label up top or the heatmap
    # itself. Two-row layout keeps the figure footprint compact.
    handles = [mpatches.Patch(color=c, label=k)
               for k, c in HL_COLORS.items()]
    ax_ld.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=(0.5, -0.18),
                  fontsize=7.5, frameon=False, ncol=min(len(handles), 5))

    save(fig, f"fig63_zoom_{label}")
    print(f"  [fig] fig63_zoom_{label}")

print(f"\nOutputs in: {OUT}")
