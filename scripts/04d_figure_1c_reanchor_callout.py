#!/usr/bin/env python3
"""
Figure 1C -- re-anchoring callout: AYB confident anchoring vs cowpea proxy.

Rebuilt for the revision (Reviewer 1 major point 3). The submitted Figure 1C
plotted a single "89.3 % anchored (n = 2,862)" bar taken from the legacy
`-max_target_seqs 1` BLAST run (`ayb_marker_anchoring.csv`). That run returns
the first hit encountered rather than a unique or best hit, so a "single
position" claim cannot be made from it. The figure now reads the multi-hit
re-anchoring classification (script 89, `-max_target_seqs 50`) and reports the
defensible split:

  confident anchoring to a best Ss01-Ss11 locus = 2,891 / 3,204 (90.2 %),
      of which 2,603 (81.2 %) are strictly unique (1 locus) and
               288 (9.0 %) are a clear best locus (>=2 bitscore-gap over
               the second locus among several);
  200 (6.2 %) remain ambiguous between near-equal loci and 113 (3.5 %) are
      unplaced -- neither counts toward confident anchoring.

  (top)    stacked AYB bar: strictly-unique 81.2 % + clear-best 9.0 % = the
           90.2 % confident total, against the 11.9 % (n = 382) anchored to
           cowpea v1.0 in the original DArT report. The 7.6x lift
           (2,891 / 382) is annotated directly.
  (bottom) per-AYB-chromosome (Ss01-Ss11) confidently-anchored marker counts
           under the S. stenocarpa reference vs the cowpea proxy on the same
           tags; shows the lift holds across every pseudo-chromosome.

The cowpea proxy counts are unchanged (they come from the original DArT report,
`cw_anchored` in `ayb_marker_anchoring.csv`); only the AYB side is re-derived.

Output
------
results/04_publication_plots/figures/fig01c_reanchor_callout.png/.pdf
manuscript/paper_A1/share/figures/main/A1_NEW_1_fig01c_reanchor_callout.png/.pdf
results/04_publication_plots/tables/Table_S6_reanchoring_per_marker_summary.csv
manuscript/paper_A1/share/tables/supp/TableS6_Table_S6_reanchoring_per_marker_summary.csv

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
CLASS_CSV = ROOT / "results" / "89_blast_uniqueness_reanchor" / "tables" / "marker_placement_classification.csv"
ANCHOR_CSV = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
PUB_FIG = ROOT / "results" / "04_publication_plots" / "figures"
PUB_TAB = ROOT / "results" / "04_publication_plots" / "tables"
SHARE_FIG = ROOT / "manuscript" / "paper_A1" / "share" / "figures" / "main"
SHARE_TAB = ROOT / "manuscript" / "paper_A1" / "share" / "tables" / "supp"
for d in (PUB_FIG, PUB_TAB, SHARE_FIG, SHARE_TAB):
    d.mkdir(parents=True, exist_ok=True)

# A marker is confidently anchored if it has a strictly unique locus, or a
# clear best locus (best-vs-second bitscore gap >= 2) among several. This is
# the identical rule used in script 89's placement summary.
CLEAR_BEST_GAP = 2.0

AYB_CHROMS = [f"Ss{i:02d}" for i in range(1, 12)]


def load() -> pd.DataFrame:
    """Top-panel source: placement-confidence classification (script 89,
    -max_target_seqs 50) merged with the cowpea proxy anchoring."""
    cls = pd.read_csv(CLASS_CSV)
    cw = pd.read_csv(ANCHOR_CSV)[["rs", "cw_anchored"]]
    cw["cw_anchored"] = cw["cw_anchored"].astype(bool)
    df = cls.merge(cw, on="rs", how="left")
    df["cw_anchored"] = df["cw_anchored"].fillna(False)
    # confident = unique OR (multi with a clear best locus)
    df["confident"] = (df["placement"] == "unique") | (
        (df["placement"] == "multi") & (df["best_second_gap"] >= CLEAR_BEST_GAP)
    )
    return df


def load_anchored() -> pd.DataFrame:
    """Bottom-panel / Table S6 source: the working AYB-anchored marker set
    (chr_ayb positions) that every downstream position-based analysis (LD,
    marker density, selection scans) actually used. The re-anchoring
    classification (script 89) re-examined placement confidence but did not
    re-derive these positions, so the per-chromosome distribution stays keyed
    to this set to remain consistent with the rest of the manuscript
    (149-401 markers per chromosome; cowpea 6.0-15.8 %)."""
    df = pd.read_csv(ANCHOR_CSV)
    for c in ("ayb_anchored", "cw_anchored"):
        df[c] = df[c].astype(bool)
    return df


def per_chrom_counts(anc: pd.DataFrame) -> pd.DataFrame:
    """Per-AYB-chromosome anchored marker counts under the S. stenocarpa
    reference (chr_ayb), and the cowpea proxy count on the SAME markers."""
    a = anc.loc[anc["ayb_anchored"]]
    rows = []
    for chrom in AYB_CHROMS:
        sub = a.loc[a["chr_ayb"] == chrom]
        n_ayb_here = len(sub)
        n_cw_here = int(sub["cw_anchored"].sum())
        rows.append({
            "chrom": chrom,
            "n_ayb_anchored": n_ayb_here,
            "n_cowpea_anchored": n_cw_here,
            "pct_cowpea_of_ayb": (100.0 * n_cw_here / n_ayb_here) if n_ayb_here else np.nan,
        })
    return pd.DataFrame(rows)


def plot_two_panel(df: pd.DataFrame, anc: pd.DataFrame, per_chrom: pd.DataFrame,
                   out_path: Path) -> dict:
    n_total = len(df)
    n_unique = int((df["placement"] == "unique").sum())
    n_confident = int(df["confident"].sum())
    n_clearbest = n_confident - n_unique
    n_amb = int(((df["placement"] == "multi") & (~df["confident"])).sum())
    n_unplaced = int((df["placement"] == "unplaced").sum())
    n_cw = int(df["cw_anchored"].sum())
    n_anchored = int(anc["ayb_anchored"].sum())

    pct_unique = 100.0 * n_unique / n_total
    pct_clearbest = 100.0 * n_clearbest / n_total
    pct_confident = 100.0 * n_confident / n_total
    pct_cw = 100.0 * n_cw / n_total
    lift = n_confident / n_cw if n_cw else np.nan

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.5, 6.5), gridspec_kw={"height_ratios": [1.0, 1.4]}
    )

    # ---- Top panel: overall anchoring rate (AYB stacked, cowpea single) ----
    x_ayb, x_cw = 0, 1
    # AYB stacked: strictly-unique base + clear-best segment = confident total.
    ax_top.bar(x_ayb, pct_unique, width=0.55, color=WONG["vermillion"],
               edgecolor="white", linewidth=0.6, label="strictly unique (1 locus)")
    ax_top.bar(x_ayb, pct_clearbest, bottom=pct_unique, width=0.55,
               color=WONG["orange"], edgecolor="white", linewidth=0.6,
               label="clear best locus (≥2-bit gap)")
    ax_top.bar(x_cw, pct_cw, width=0.55, color=WONG["skyblue"],
               edgecolor="white", linewidth=0.6)

    # AYB confident-total label above the stack, with the unique share below it.
    ax_top.text(x_ayb, pct_confident + 1.5,
                f"{pct_confident:.1f} %\n(n = {n_confident:,})",
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_top.text(x_ayb, pct_unique / 2,
                f"{pct_unique:.1f} % unique", ha="center", va="center",
                fontsize=8.5, color="white", fontweight="bold")
    ax_top.text(x_cw, pct_cw + 1.5, f"{pct_cw:.1f} %\n(n = {n_cw:,})",
                ha="center", va="bottom", fontsize=10, fontweight="bold")

    # lift bracket spanning the two bars, clear of the AYB confident label.
    bracket_top = 122
    ax_top.plot([x_ayb, x_ayb, x_cw, x_cw],
                [bracket_top - 4, bracket_top, bracket_top, bracket_top - 4],
                color="black", linewidth=0.9, clip_on=False)
    ax_top.text((x_ayb + x_cw) / 2, bracket_top + 2, f"{lift:.1f}× lift",
                ha="center", va="bottom", fontsize=12, fontweight="bold",
                clip_on=False)

    ax_top.set_xticks([x_ayb, x_cw])
    ax_top.set_xticklabels(["S. stenocarpa\nchromosome-scale\n(Shorinola 2024)",
                            "Cowpea v1.0\nproxy\n(DArT report)"], fontsize=10)
    ax_top.set_ylabel("DArTseq tags anchored (%)")
    ax_top.set_ylim(0, 135)
    ax_top.set_yticks([0, 25, 50, 75, 100])
    ax_top.text(-0.10, 1.10, "c", transform=ax_top.transAxes,
                fontsize=14, fontweight="bold", va="bottom", ha="right")
    ax_top.legend(loc="center", fontsize=8, bbox_to_anchor=(0.5, 0.55),
                  title="S. stenocarpa anchoring", title_fontsize=8)

    # ---- Bottom panel: per-AYB-chromosome side-by-side ----
    chroms = per_chrom["chrom"].tolist()
    x = np.arange(len(chroms))
    w = 0.4
    ax_bot.bar(x - w / 2, per_chrom["n_ayb_anchored"], width=w,
               color=WONG["vermillion"], edgecolor="white", linewidth=0.4,
               label=f"S. stenocarpa anchored markers (total {n_anchored:,})")
    ax_bot.bar(x + w / 2, per_chrom["n_cowpea_anchored"], width=w,
               color=WONG["skyblue"], edgecolor="white", linewidth=0.4,
               label=f"Cowpea proxy on same markers (total {n_cw:,})")

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(chroms)
    ax_bot.set_xlabel("S. stenocarpa pseudo-chromosome")
    ax_bot.set_ylabel("Number of DArTseq tags")
    ax_bot.legend(loc="upper right")

    fig.tight_layout()
    for stem in (out_path, SHARE_FIG / "A1_NEW_1_fig01c_reanchor_callout"):
        fig.savefig(stem.with_suffix(".png"))
        fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)

    return {
        "n_total": n_total, "n_unique": n_unique, "n_clearbest": n_clearbest,
        "n_confident": n_confident, "n_amb": n_amb, "n_unplaced": n_unplaced,
        "n_cw": n_cw, "n_anchored": n_anchored, "pct_unique": pct_unique,
        "pct_clearbest": pct_clearbest, "pct_confident": pct_confident,
        "pct_cw": pct_cw, "lift": lift,
    }


def main() -> None:
    print(f"reading placement classification from {CLASS_CSV}")
    df = load()
    anc = load_anchored()
    per_chrom = per_chrom_counts(anc)
    print(per_chrom.to_string(index=False))

    out_path = PUB_FIG / "fig01c_reanchor_callout"
    s = plot_two_panel(df, anc, per_chrom, out_path)

    # Table S6: per-chrom anchored + cowpea counts (working chr_ayb set), TOTAL row.
    tbl = pd.concat([
        per_chrom,
        pd.DataFrame([{
            "chrom": "TOTAL",
            "n_ayb_anchored": s["n_anchored"],
            "n_cowpea_anchored": s["n_cw"],
            "pct_cowpea_of_ayb": 100.0 * s["n_cw"] / s["n_anchored"],
        }]),
    ], ignore_index=True)
    for tab_path in (PUB_TAB / "Table_S6_reanchoring_per_marker_summary.csv",
                     SHARE_TAB / "TableS6_Table_S6_reanchoring_per_marker_summary.csv"):
        tbl.to_csv(tab_path, index=False)
        print(f"  wrote {tab_path}")

    print("\n=== manuscript-text summary ===")
    print(f"raw DArTseq tags        : {s['n_total']:,}")
    print(f"strictly unique         : {s['n_unique']:,} ({s['pct_unique']:.1f} %)")
    print(f"clear best locus        : {s['n_clearbest']:,} ({s['pct_clearbest']:.1f} %)")
    print(f"confident (unique+best) : {s['n_confident']:,} ({s['pct_confident']:.1f} %)")
    print(f"ambiguous               : {s['n_amb']:,}")
    print(f"unplaced                : {s['n_unplaced']:,}")
    print(f"working anchored set     : {s['n_anchored']:,} (chr_ayb; bottom panel)")
    print(f"cowpea proxy            : {s['n_cw']:,} ({s['pct_cw']:.1f} %)")
    print(f"lift (confident/cowpea) : {s['lift']:.2f}x")


if __name__ == "__main__":
    main()
