#!/usr/bin/env python3
"""
Permutation null on the candidate-gene-callout enrichment.

The top-30 candidate-gene callout (script 21; Table 2 / Table 3 of the
A2 draft) is heuristic: top-30 SNPs swept +/- 50 kb against trait-specific
keyword sets, with no formal multiple-testing correction on the gene-
overlap statistic. Reviewers will ask whether the recovered pathway-gene
hits are real or whether the same number would surface from 30 random
SNPs sampled from the same AYB-anchored marker pool.

Method
------
For each trait:
  1. Observed enrichment metric = the number of top-30 SNPs (out of 30)
     whose +/- 50 kb candidate-gene set contains at least one keyword-
     matched gene under the trait's keyword set.
  2. Null distribution: 1,000 permutation replicates. Each replicate
     samples 30 SNPs from the full AYB-anchored marker pool (2,862 SNPs),
     stratified per chromosome so the chromosome-density bias of the AYB
     annotation is preserved.
  3. Per-replicate null metric: same +/- 50 kb sweep, same keyword set,
     same "at least one matched gene" rule.
  4. Empirical p-value = (1 + #{null >= observed}) / (1 + B), with B = 1000.
  5. Benjamini-Hochberg FDR across the 13 traits.

The chromosome-stratified null is essential. An unstratified null would
inflate p-values because chromosome-density-rich chromosomes (Ss10 with
the ALMT4 cluster; Ss03 with the AYB seed-protein cluster) are over-
represented among the top-30 SNPs by gene-density, not by trait
association.

Outputs (results/69_callout_permutation/)
-----------------------------------------
tables/
    per_trait_enrichment_p.csv          observed + null mean / SD + p + BH q
figures/
    fig_callout_enrichment_forest.png    forest of observed vs null per trait

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
ANCHOR_CSV = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
GFF_PATH = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
TOP30_DIR = ROOT / "results" / "21_candidate_genes_ayb" / "tables"
OUT = ROOT / "results" / "69_callout_permutation"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

TRAITS = [
    "Tannin", "Phenol", "Flavonoid", "Antioxidant",
    "Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
    "Seed_Coat_Tannin", "Crude_Protein",
    "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate",
]
WIN_BP = 50_000   # +/- 50 kb candidate-gene window
N_PERM = 1_000
SEED = 1234

# Trait-specific keyword sets, identical to script 21 conventions.
# Patterns are case-insensitive when matched.
PHENYL_KEYS = [
    (r"\bphenylalanine ammonia[- ]?lyase\b", "PAL"),
    (r"\bcinnamate[- ]?4[- ]?hydroxylase\b", "C4H"),
    (r"\bCYP73A\b", "C4H"),
    (r"\b4[- ]?coumarate.*CoA ligase\b", "4CL"),
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
    (r"\bhydroxycinnamoyl\b", "HCT/HQT"),
    (r"\bproanthocyanidin\b", "PA-biosynth"),
    (r"\bisoflavone\b", "isoflavone"),
    (r"\bO[- ]?methyltransferase\b", "OMT"),
]

SEED_DEV_KEYS = [
    (r"\bAINTEGUMENTA\b", "ANT"),
    (r"\bAP2[- ]?like\b", "AP2"),
    (r"\bauxin[- ]?response factor\b", "ARF"),
    (r"\bARF[ \-]?[0-9]+\b", "ARF"),
    (r"\bAUX[/_]?IAA\b", "AUX/IAA"),
    (r"\bbrassinosteroid\b", "BR-pathway"),
    (r"\bexpansin\b", "expansin"),
    (r"\bpectin methylesterase\b", "PME"),
    (r"\bcell wall\b", "cell-wall-mod"),
    (r"\bintegument\b", "integument-dev"),
    (r"\bovule\b", "ovule-dev"),
    (r"\bseed coat\b", "seed-coat-dev"),
    (r"\bCYP78A\b", "KLU"),
    (r"\bcytochrome P450 78\b", "KLU"),
    (r"\bMADS\b", "MADS"),
    (r"\bSEEDSTICK\b", "STK"),
    (r"\bSHATTERPROOF\b", "SHP"),
    (r"\bAGAMOUS\b", "AG/STK"),
    (r"\bGT[- ]?1\b", "TF-other"),
    (r"\btrihelix\b", "TF-other"),
    (r"\bsite[- ]?1 protease\b", "TF-other"),
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
    (r"\bnuclear factor Y\b", "NF-Y"),
    (r"\bLEAFY COTYLEDON\b", "LEC"),
    (r"\bFUSCA3\b", "FUS3"),
    (r"\bABI3\b", "ABI3"),
    (r"\bglutamine synthetase\b", "GS"),
    (r"\basparagine synthetase\b", "AS"),
    (r"\bglutamate dehydrogenase\b", "GDH"),
    (r"\bglutamate synthase\b", "GOGAT"),
    (r"\bproteinase inhibitor\b", "PI"),
    (r"\bKunitz\b", "Kunitz-PI"),
]

OXALATE_KEYS = [
    (r"\boxalate oxidase\b", "OXO"),
    (r"\bgermin\b", "germin-OXO"),
    (r"\boxalate decarboxylase\b", "OXDC"),
    (r"\boxalyl[- ]?CoA\b", "oxalyl-CoA"),
    (r"\bAAE3\b", "AAE3"),
    (r"\bisocitrate lyase\b", "ICL"),
    (r"\bglycolate oxidase\b", "GLO"),
    (r"\bglyoxylate\b", "glyoxylate"),
    (r"\bDETOXIFICATION\b", "MATE/DTX"),
    (r"\bDTX[0-9]+\b", "MATE/DTX"),
    (r"\bMATE\b", "MATE"),
    (r"\baluminum[- ]?activated malate transporter\b", "ALMT"),
    (r"\bALMT[0-9]+\b", "ALMT"),
    (r"\bcalcium transporter\b", "Ca-transporter"),
    (r"\bvacuolar transporter\b", "vacuolar"),
    (r"\bvacuolar\b", "vacuolar"),
    (r"\btonoplast\b", "tonoplast"),
    (r"\bnucleobase[- ]?ascorbate transporter\b", "ascorbate-pathway"),
    (r"\bascorbate\b", "ascorbate-pathway"),
]

TRAIT_KEYWORDS: dict[str, list[tuple[str, str]]] = {}
for t in ("Tannin", "Phenol", "Flavonoid", "Antioxidant", "Seed_Coat_Tannin"):
    TRAIT_KEYWORDS[t] = PHENYL_KEYS
for t in ("Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds"):
    TRAIT_KEYWORDS[t] = SEED_DEV_KEYS
TRAIT_KEYWORDS["Crude_Protein"] = STORAGE_PROTEIN_KEYS
for t in ("Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate"):
    TRAIT_KEYWORDS[t] = OXALATE_KEYS

COMPILED = {t: [(re.compile(p, re.IGNORECASE), label) for p, label in v]
            for t, v in TRAIT_KEYWORDS.items()}

# A single compound regex per trait that fires on ANY of the keywords; the
# per-pattern compiled regexes above are kept for labelling but never used
# in the hot path.
COMPOUND_RE = {
    t: re.compile("|".join(f"(?:{p})" for p, _label in v), re.IGNORECASE)
    for t, v in TRAIT_KEYWORDS.items()
}


# ---------------------------------------------------------------------------
# GFF parsing -- build a gene table with chrom / start / end / search string
# ---------------------------------------------------------------------------
def parse_gff() -> dict[str, dict]:
    """Return a per-chromosome dict of arrays for fast windowed lookup.

    Each chromosome's entry has:
      start  : np.array of mRNA start positions (sorted)
      end    : np.array of mRNA end positions (parallel)
      text   : np.array of search_text strings (parallel, dtype=object)
      kw_hit : per-trait np.array of bool, True if the gene text matches
               the trait's compound keyword regex (pre-computed once for
               every (chrom, trait) pair).

    Pre-computing the per-gene keyword hit at load time lets the per-
    permutation hot path skip regex entirely -- it only does a window
    selection + boolean reduce."""
    print(f"[gff] reading {GFF_PATH.name}")
    rows = []
    with open(GFF_PATH) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "mRNA":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]
            search_text = attrs.replace(";", " ").replace("=", " ")
            rows.append((chrom, start, end, search_text))

    print(f"  parsed {len(rows)} mRNA records; pre-computing per-trait keyword hits")
    by_chrom: dict[str, dict] = {}
    for chrom, start, end, text in rows:
        d = by_chrom.setdefault(chrom, {"start": [], "end": [], "text": []})
        d["start"].append(start)
        d["end"].append(end)
        d["text"].append(text)

    for chrom, d in by_chrom.items():
        order = np.argsort(d["start"])
        d["start"] = np.asarray(d["start"], dtype=np.int64)[order]
        d["end"] = np.asarray(d["end"], dtype=np.int64)[order]
        d["text"] = np.asarray(d["text"], dtype=object)[order]
        # Per-trait pre-scan: bool array, True where the gene's text
        # matches the trait's compound regex.
        d["kw_hit"] = {}
        for trait, compound in COMPOUND_RE.items():
            d["kw_hit"][trait] = np.asarray(
                [bool(compound.search(t)) for t in d["text"]],
                dtype=bool,
            )
    print(f"  chromosomes indexed: {sorted(by_chrom.keys())}")
    return by_chrom


# ---------------------------------------------------------------------------
# Per-SNP keyword hit -- vectorised window lookup
# ---------------------------------------------------------------------------
def snp_has_keyword_match(chrom: str, pos: int,
                          by_chrom: dict[str, dict],
                          trait: str,
                          window_bp: int = WIN_BP) -> bool:
    """True if any gene within +/- window_bp of (chrom, pos) was tagged as
    a keyword hit for this trait during the pre-scan in parse_gff."""
    d = by_chrom.get(chrom)
    if d is None:
        return False
    starts = d["start"]
    ends = d["end"]
    # A gene overlaps the window [pos - W, pos + W] iff start <= pos + W
    # AND end >= pos - W. Use boolean masks on the sorted start array.
    in_window = (starts <= pos + window_bp) & (ends >= pos - window_bp)
    if not np.any(in_window):
        return False
    return bool(np.any(d["kw_hit"][trait][in_window]))


def count_hits_for_snp_set(snp_chrom: np.ndarray,
                            snp_pos: np.ndarray,
                            by_chrom: dict[str, dict],
                            trait: str) -> int:
    return int(sum(
        snp_has_keyword_match(str(c), int(p), by_chrom, trait)
        for c, p in zip(snp_chrom, snp_pos)
    ))


# ---------------------------------------------------------------------------
# Per-trait permutation
# ---------------------------------------------------------------------------
def run_per_trait(trait: str,
                  by_chrom: dict[str, dict],
                  anchored_by_chrom: dict[str, np.ndarray],
                  rng: np.random.Generator,
                  n_perm: int = N_PERM) -> dict:
    # Observed: top-30 anchored SNPs for this trait, from script 21.
    top30_path = TOP30_DIR / f"snp_top30_anchored_{trait}.csv"
    top30 = pd.read_csv(top30_path)
    observed = count_hits_for_snp_set(
        top30["chr_ayb"].astype(str).values,
        top30["snp_pos_ayb"].astype(float).astype(int).values,
        by_chrom, trait,
    )

    # Chromosome-stratified null: preserves the top-30's per-chromosome
    # count so chromosome-density bias does not inflate p-values.
    chrom_counts = top30["chr_ayb"].value_counts().to_dict()
    null = np.zeros(n_perm, dtype=int)
    for b in range(n_perm):
        sampled_chrom = []
        sampled_pos = []
        for chrom, k in chrom_counts.items():
            pool_pos = anchored_by_chrom.get(chrom)
            if pool_pos is None or len(pool_pos) == 0:
                continue
            if len(pool_pos) < k:
                idx = rng.integers(len(pool_pos), size=k)
            else:
                idx = rng.choice(len(pool_pos), size=k, replace=False)
            sampled_chrom.extend([chrom] * k)
            sampled_pos.extend(pool_pos[idx].tolist())
        null[b] = count_hits_for_snp_set(
            np.asarray(sampled_chrom), np.asarray(sampled_pos),
            by_chrom, trait,
        )

    emp_p = float((1 + np.sum(null >= observed)) / (1 + n_perm))
    return {
        "trait": trait,
        "observed_hits": observed,
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
        "null_max": int(null.max()),
        "n_top30_chroms": len(chrom_counts),
        "empirical_p": emp_p,
    }


# ---------------------------------------------------------------------------
# Benjamini-Hochberg
# ---------------------------------------------------------------------------
def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    m = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    adj_sorted = np.zeros(m)
    running = 1.0
    for k in range(m - 1, -1, -1):
        candidate = sorted_p[k] * m / (k + 1)
        running = min(running, candidate)
        adj_sorted[k] = running
    adj = np.zeros(m)
    adj[order] = adj_sorted
    return adj


# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------
def plot_enrichment_forest(results: pd.DataFrame, out_path: Path) -> None:
    results = results.copy().sort_values("empirical_p").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    y = np.arange(len(results))

    # Null mean +/- 1.96 SD as a horizontal interval per trait.
    for i, (_, row) in enumerate(results.iterrows()):
        lo = row["null_mean"] - 1.96 * row["null_sd"]
        hi = row["null_mean"] + 1.96 * row["null_sd"]
        ax.plot([lo, hi], [i, i], color="grey", linewidth=2.2, alpha=0.55,
                solid_capstyle="round")
        ax.scatter([row["null_mean"]], [i], color="grey", s=35, zorder=2)
        # Observed point coloured by significance.
        sig_q = row["bh_q"] < 0.10
        sig_p = row["empirical_p"] < 0.05
        if sig_q:
            colour = WONG["vermillion"]
        elif sig_p:
            colour = WONG["orange"]
        else:
            colour = WONG["blue"]
        ax.scatter([row["observed_hits"]], [i], color=colour, s=110,
                   edgecolor="black", linewidth=0.7, zorder=4,
                   label=None)
        ax.text(max(row["observed_hits"], row["null_mean"]) + 0.4, i,
                f"p={row['empirical_p']:.3f}, q={row['bh_q']:.2f}",
                fontsize=8, va="center")

    ax.set_yticks(y)
    ax.set_yticklabels(results["trait"])
    ax.set_xlabel("Number of top-30 SNPs hitting >=1 keyword-matched gene within +/- 50 kb")
    ax.set_title(
        "Permutation null on candidate-gene callout enrichment\n"
        f"observed (colour) vs null mean +/- 1.96 SD (grey) "
        f"across {N_PERM} chromosome-stratified replicates"
    )
    ax.invert_yaxis()
    ax.set_xlim(left=-0.5)
    ax.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=WONG["vermillion"],
                    markersize=10, label="BH q < 0.10"),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=WONG["orange"],
                    markersize=10, label="empirical p < 0.05 only"),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=WONG["blue"],
                    markersize=10, label="n.s."),
        plt.Line2D([0], [0], color='grey', linewidth=2.2, label="null mean +/- 1.96 SD"),
    ], loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    rng = np.random.default_rng(SEED)

    by_chrom = parse_gff()

    print(f"[anchor] reading {ANCHOR_CSV.name}")
    anchored = pd.read_csv(ANCHOR_CSV)
    anchored = anchored.loc[anchored["ayb_anchored"] == True].copy()
    anchored = anchored.dropna(subset=["chr_ayb", "snp_pos_ayb"])
    print(f"  AYB-anchored pool: {len(anchored):,} markers across {anchored['chr_ayb'].nunique()} chromosomes")

    # Pre-group anchored marker positions by chromosome for fast sampling.
    anchored_by_chrom: dict[str, np.ndarray] = {}
    for chrom, g in anchored.groupby("chr_ayb"):
        anchored_by_chrom[str(chrom)] = g["snp_pos_ayb"].astype(float).astype(int).values

    results = []
    for i, trait in enumerate(TRAITS, start=1):
        print(f"[{i:2d}/{len(TRAITS)}] {trait} ...", flush=True)
        rec = run_per_trait(trait, by_chrom, anchored_by_chrom, rng, n_perm=N_PERM)
        results.append(rec)
        print(f"      observed = {rec['observed_hits']}, "
              f"null mean = {rec['null_mean']:.2f} +/- {rec['null_sd']:.2f}, "
              f"emp p = {rec['empirical_p']:.4f}", flush=True)

    results_df = pd.DataFrame(results)
    results_df["bh_q"] = bh_fdr(results_df["empirical_p"].values)
    results_df = results_df.sort_values("empirical_p").reset_index(drop=True)

    out_csv = TAB / "per_trait_enrichment_p.csv"
    results_df.to_csv(out_csv, index=False)
    print(f"\n[result] wrote {out_csv}")
    print(results_df.to_string(index=False))

    print("\n[plot] enrichment forest")
    plot_enrichment_forest(results_df, FIG / "fig_callout_enrichment_forest")

    print("\n=== manuscript-text summary ===")
    sig_q = results_df.loc[results_df["bh_q"] < 0.10, "trait"].tolist()
    sig_p_only = results_df.loc[(results_df["empirical_p"] < 0.05)
                                & (results_df["bh_q"] >= 0.10), "trait"].tolist()
    print(f"  traits with BH q < 0.10                : {sig_q}")
    print(f"  traits with raw p < 0.05 (not q < 0.10): {sig_p_only}")
    sol = results_df.loc[results_df["trait"] == "Soluble_Oxalate"]
    if not sol.empty:
        r = sol.iloc[0]
        print(f"  Soluble_Oxalate observed enrichment    : "
              f"{int(r['observed_hits'])} of 30 SNPs hit, null mean "
              f"{r['null_mean']:.2f}, emp p = {r['empirical_p']:.4f}, "
              f"BH q = {r['bh_q']:.4f}")


if __name__ == "__main__":
    main()
