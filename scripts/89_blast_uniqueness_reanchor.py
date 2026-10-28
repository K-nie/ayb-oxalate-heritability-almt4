#!/usr/bin/env python3
"""
89 - BLAST placement uniqueness and paralogy -> heterozygosity test.

Reviewer 1 major point 3: the original anchoring used megablast with
`-max_target_seqs 1`, which cannot establish that a marker maps to a
*single* genomic position - it simply reports one hit regardless of how
many exist (Shah et al. 2019; NCBI's own guidance). The "2,862 of 3,204
(89.3%) mapped to a single position" claim therefore has to be re-derived
from a run that reports *all* qualifying hits per tag, so uniqueness can
actually be counted.

This script:
  1. reads a megablast run that reported up to 50 hits per tag under the
     same thresholds used in the paper (-perc_identity 95 -evalue 1e-10),
  2. collapses hits to distinct genomic loci per marker (merging HSPs and
     the two tag representations of each marker within a 1 kb window),
  3. classifies each of the 3,204 study markers as uniquely placed
     (exactly one locus), multi-mapping (>=2 loci), or unplaced (0 hits),
  4. tests whether multi-mapping markers - the expected signature of
     collapsed paralogues / repeats - carry inflated observed
     heterozygosity relative to uniquely placed markers (Mann-Whitney U),
     which speaks directly to the panel's H_obs ~ 0.30 anomaly.

Inputs
  - refs/ayb_genome/dartseq_tags.fasta       (6,200 tags; 2 per marker)
  - /tmp/ayb_blast multihit BLAST tsv        (outfmt 6, up to 50 hits/tag)
  - refs/ayb_genome/ayb_marker_anchoring.csv (3,204 study markers)
  - results/01_qc_pca_power/tables/marker_qc.csv (per-marker het)

Outputs land in results/89_blast_uniqueness_reanchor/.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
MULTIHIT = Path("/tmp/ayb_blast/dartseq_to_ayb.multihit.tsv")
ANCHOR_CSV = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
MARKER_QC = ROOT / "results" / "01_qc_pca_power" / "tables" / "marker_qc.csv"

OUTDIR = ROOT / "results" / "89_blast_uniqueness_reanchor"
TAB = OUTDIR / "tables"
TAB.mkdir(parents=True, exist_ok=True)

# BLAST significance thresholds (identical to the paper's anchoring run).
PIDENT_MIN = 95.0
EVALUE_MAX = 1e-10
# Two hits on the same chromosome within this distance are treated as the
# same genomic locus (one 67 bp tag cannot resolve finer than this, and
# the reference vs SNP-allele tag of a marker land at the same spot).
LOCUS_MERGE_BP = 1000


def marker_of(tag: str) -> str:
    """Base DArT clone id = substring before the first '|'."""
    return tag.split("|", 1)[0]


def count_loci(hits: pd.DataFrame) -> int:
    """Number of distinct genomic loci among a marker's qualifying hits.
    Hits on the same chromosome within LOCUS_MERGE_BP bp collapse to one."""
    n = 0
    for _chrom, sub in hits.groupby("sseqid"):
        starts = np.sort(sub["hit_start"].to_numpy())
        clusters = 1
        last = starts[0]
        for s in starts[1:]:
            if s - last > LOCUS_MERGE_BP:
                clusters += 1
            last = s
        n += clusters
    return n


def best_second_gap(hits: pd.DataFrame) -> float:
    """Bitscore gap between the best and the best *other-locus* hit.
    NaN when the marker has a single locus. Large gap => the top locus is
    unambiguous even if secondary hits exist."""
    if hits.empty:
        return np.nan
    ordered = hits.sort_values("bitscore", ascending=False)
    top = ordered.iloc[0]
    for _, h in ordered.iloc[1:].iterrows():
        same_locus = (h["sseqid"] == top["sseqid"]
                      and abs(h["hit_start"] - top["hit_start"]) <= LOCUS_MERGE_BP)
        if not same_locus:
            return float(top["bitscore"] - h["bitscore"])
    return np.nan


def main() -> None:
    print("[load] multihit BLAST")
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    bl = pd.read_csv(MULTIHIT, sep="\t", names=cols)
    print(f"  raw HSPs: {len(bl)}, tags with >=1 HSP: {bl['qseqid'].nunique()}")

    # Enforce the paper thresholds explicitly (belt and braces).
    bl = bl.loc[(bl["pident"] >= PIDENT_MIN) & (bl["evalue"] <= EVALUE_MAX)].copy()
    bl["hit_start"] = bl[["sstart", "send"]].min(axis=1)
    bl["marker"] = bl["qseqid"].map(marker_of)
    print(f"  HSPs passing pident>={PIDENT_MIN} & evalue<={EVALUE_MAX}: {len(bl)}")

    print("[markers] study set")
    anchor = pd.read_csv(ANCHOR_CSV)
    study_rs = anchor["rs"].astype(str).unique()
    n_rs = len(study_rs)
    n_clone = pd.Index([marker_of(r) for r in study_rs]).nunique()
    print(f"  study SNP-markers (rs): {n_rs}; unique DArT clones: {n_clone}")

    # Placement is a property of the clone/tag, so classify each clone once
    # and propagate to every SNP-marker (rs) carried on that clone. This
    # keeps the manuscript's 3,204 SNP-marker denominator while counting
    # genomic loci at the correct (clone) resolution.
    hits_by_marker = {m: sub for m, sub in bl.groupby("marker")}
    clone_class: dict[str, dict] = {}

    def classify_clone(m: str) -> dict:
        if m in clone_class:
            return clone_class[m]
        h = hits_by_marker.get(m)
        if h is None or h.empty:
            rec = {"n_loci": 0, "placement": "unplaced", "top_chrom": "",
                   "top_pos": np.nan, "top_bitscore": np.nan,
                   "best_second_gap": np.nan}
        else:
            nloci = count_loci(h)
            top = h.sort_values("bitscore", ascending=False).iloc[0]
            rec = {"n_loci": nloci,
                   "placement": "unique" if nloci == 1 else "multi",
                   "top_chrom": top["sseqid"], "top_pos": top["hit_start"],
                   "top_bitscore": float(top["bitscore"]),
                   "best_second_gap": best_second_gap(h)}
        clone_class[m] = rec
        return rec

    rows = []
    for rs in study_rs:
        rec = classify_clone(marker_of(rs))
        rows.append({"rs": rs, "clone": marker_of(rs), **rec})
    clip = pd.DataFrame(rows)

    counts = clip["placement"].value_counts()
    n_total = len(clip)
    n_unique = int(counts.get("unique", 0))
    n_multi = int(counts.get("multi", 0))
    n_unplaced = int(counts.get("unplaced", 0))
    print("\n=== placement classification (n = %d SNP-markers) ===" % n_total)
    print(f"  uniquely placed (1 locus)     : {n_unique} ({100*n_unique/n_total:.1f}%)")
    print(f"  multi-mapping   (>=2 loci)    : {n_multi} ({100*n_multi/n_total:.1f}%)")
    print(f"  unplaced        (0 hits)      : {n_unplaced} ({100*n_unplaced/n_total:.1f}%)")
    placed = n_unique + n_multi
    if placed:
        print(f"  among PLACED markers, unique  : {n_unique}/{placed} "
              f"({100*n_unique/placed:.1f}%)")

    # Clone-level view (deduplicated on clone id).
    clone_df = clip.drop_duplicates("clone")
    cc = clone_df["placement"].value_counts()
    print(f"  [clone-level] unique={int(cc.get('unique',0))}, "
          f"multi={int(cc.get('multi',0))}, unplaced={int(cc.get('unplaced',0))} "
          f"of {len(clone_df)} clones")

    # Ambiguity of the multi-mappers: how many have a clear best locus?
    multi = clip.loc[clip["placement"] == "multi"]
    n_ambiguous = 0
    n_confident = n_unique
    if len(multi):
        n_ambiguous = int((multi["best_second_gap"] < 2.0).sum())
        n_clear_best = len(multi) - n_ambiguous  # secondary hit, but clear top locus
        n_confident = n_unique + n_clear_best
        print(f"  multi-mappers with best-vs-second bitscore gap < 2: {n_ambiguous} "
              f"({100*n_ambiguous/len(multi):.1f}% of multi-mapping SNP-markers)")
        print(f"  multi-mappers with an unambiguous best locus (gap >= 2): {n_clear_best}")
    print(f"  --> confidently placed to a best locus "
          f"(single locus OR clear best): {n_confident} "
          f"({100*n_confident/n_total:.1f}%)")
    print(f"  --> genuinely ambiguous (>=2 near-equal loci): {n_ambiguous} "
          f"({100*n_ambiguous/n_total:.1f}%)")

    # -------- paralogy -> heterozygosity test --------
    print("\n[het] paralogy -> heterozygosity enrichment test")
    qc = pd.read_csv(MARKER_QC)
    qc["rs"] = qc["rs"].astype(str)
    clip = clip.merge(qc[["rs", "het"]].rename(columns={"het": "het_obs"}),
                      on="rs", how="left")
    clip.to_csv(TAB / "marker_placement_classification.csv", index=False)

    uni_het = clip.loc[(clip["placement"] == "unique"), "het_obs"].dropna()
    mul_het = clip.loc[(clip["placement"] == "multi"), "het_obs"].dropna()
    print(f"  unique-marker het: n={len(uni_het)}, "
          f"median={uni_het.median():.3f}, mean={uni_het.mean():.3f}")
    print(f"  multi-marker  het: n={len(mul_het)}, "
          f"median={mul_het.median():.3f}, mean={mul_het.mean():.3f}")
    mwu_p = np.nan
    if len(uni_het) and len(mul_het):
        U, p = mannwhitneyu(mul_het, uni_het, alternative="greater")
        mwu_p = float(p)
        n1, n2 = len(mul_het), len(uni_het)
        rbc = 1.0 - (2.0 * U) / (n1 * n2)  # rank-biserial (negative => mul>uni)
        print(f"  Mann-Whitney U (multi > unique het): U={U:.0f}, "
              f"p={p:.3e}, rank-biserial={-rbc:.3f}")
        # How much of the panel's H_obs excess sits in multi-mappers?
        print(f"  panel mean het (all markers): {clip['het_obs'].mean():.3f}")

    # Summary table for the manuscript / Author Response.
    summ = pd.DataFrame([
        {"metric": "SNP-markers (rs)", "value": n_total},
        {"metric": "unique DArT clones", "value": n_clone},
        {"metric": "uniquely placed (1 locus)", "value": n_unique},
        {"metric": "uniquely placed pct of all", "value": round(100*n_unique/n_total, 1)},
        {"metric": "confidently placed (1 locus OR clear best)", "value": n_confident},
        {"metric": "confidently placed pct of all", "value": round(100*n_confident/n_total, 1)},
        {"metric": "multi-mapping (>=2 loci)", "value": n_multi},
        {"metric": "multi-mapping pct of all", "value": round(100*n_multi/n_total, 1)},
        {"metric": "genuinely ambiguous (>=2 near-equal loci)", "value": n_ambiguous},
        {"metric": "genuinely ambiguous pct of all", "value": round(100*n_ambiguous/n_total, 1)},
        {"metric": "unplaced (0 hits)", "value": n_unplaced},
        {"metric": "unique pct of placed", "value": round(100*n_unique/placed, 1) if placed else np.nan},
        {"metric": "unique-marker median het", "value": round(float(uni_het.median()), 4)},
        {"metric": "multi-marker median het", "value": round(float(mul_het.median()), 4)},
        {"metric": "MWU multi>unique het p", "value": mwu_p},
    ])
    summ.to_csv(TAB / "placement_summary.csv", index=False)
    print("\n[write] tables/marker_placement_classification.csv, placement_summary.csv")


if __name__ == "__main__":
    main()
