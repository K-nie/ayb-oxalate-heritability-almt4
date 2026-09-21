#!/usr/bin/env python3
"""
Figure 1C — re-anchoring callout: AYB vs cowpea coverage.

The paper-A1 Figure 1 currently has marker QC and sample QC panels but
no anchoring callout. Figure 1C closes that gap with a two-panel summary:

  (top)    overall anchoring: 89.3 % of the 3,204 DArTseq tags re-mapped
           to the S. stenocarpa chromosome-scale assembly (Shorinola
           2024) vs the 11.9 % anchored to cowpea v1.0 in the original
           DArT report. The 7.5x lift is annotated directly.
  (bottom) per-AYB-chromosome (Ss01-Ss11) marker counts under each
           reference, side-by-side; shows that the lift holds across
           every pseudo-chromosome rather than concentrating on one.

Output
------
results/04_publication_plots/figures/
    fig01c_reanchor_callout.png/.pdf
results/04_publication_plots/tables/
    Table_S6_reanchoring_per_marker_summary.csv

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
ANCHOR_CSV = ROOT / "refs" / "ayb_genome" / "ayb_marker_anchoring.csv"
PUB_FIG = ROOT / "results" / "04_publication_plots" / "figures"
PUB_TAB = ROOT / "results" / "04_publication_plots" / "tables"
for d in (PUB_FIG, PUB_TAB):
    d.mkdir(parents=True, exist_ok=True)

# The 11 AYB pseudo-chromosomes in canonical order. Anything outside this
# set (e.g. unscaffolded contigs that briefly leaked into the BLAST output)
# is grouped as "other".
AYB_CHROMS = [f"Ss{i:02d}" for i in range(1, 12)]


def load_anchoring() -> pd.DataFrame:
    df = pd.read_csv(ANCHOR_CSV)
    # ayb_anchored / cw_anchored come in as bool strings on read; the
    # explicit cast keeps the downstream sum() honest.
    for col in ("ayb_anchored", "cw_anchored"):
        df[col] = df[col].astype(bool)
    return df


def per_chrom_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Per-AYB-chromosome counts: total markers anchored to that chrom under
    AYB, and the matched count under cowpea on the SAME markers (i.e. how
    many of the markers that anchored to AYB chrom X were ALSO anchored to
    cowpea)."""
    rows = []
    for chrom in AYB_CHROMS:
        sub = df.loc[df["chr_ayb"] == chrom]
        n_ayb_here = int(sub["ayb_anchored"].sum())
        # Of the markers anchored to this AYB chrom, how many had a cowpea
        # anchor in the original DArT report. The fraction is interpretable
        # as the cowpea anchoring rate stratified by AYB chromosome.
        n_cw_here = int(sub["cw_anchored"].sum())
        rows.append({
            "chrom": chrom,
            "n_ayb_anchored": n_ayb_here,
            "n_cowpea_anchored": n_cw_here,
            "pct_cowpea_of_ayb": (100.0 * n_cw_here / n_ayb_here) if n_ayb_here else np.nan,
        })
    return pd.DataFrame(rows)


def plot_two_panel(df: pd.DataFrame, per_chrom: pd.DataFrame, out_path: Path) -> None:
    n_total = len(df)
    n_ayb = int(df["ayb_anchored"].sum())
    n_cw = int(df["cw_anchored"].sum())
    pct_ayb = 100.0 * n_ayb / n_total
    pct_cw = 100.0 * n_cw / n_total
    lift = n_ayb / n_cw if n_cw else np.nan

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(8.5, 6.5), gridspec_kw={"height_ratios": [1.0, 1.4]}
    )

    # ---- Top panel: overall anchoring rate ----
    bar_x = np.array([0, 1])
    bar_h = np.array([pct_ayb, pct_cw])
    bar_n = np.array([n_ayb, n_cw])
    bar_col = [WONG["vermillion"], WONG["skyblue"]]
    ax_top.bar(bar_x, bar_h, color=bar_col, width=0.55,
               edgecolor="white", linewidth=0.6)

    # Per-bar labels positioned above each bar; offset chosen so the AYB
    # bar's two-line label (89.3 % + n) fits below the bracket without
    # overlapping the title.
    for x, h, n in zip(bar_x, bar_h, bar_n):
        ax_top.text(x, h + 1.5, f"{h:.1f} %\n(n = {n:,})",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    # 7.5x lift annotation as an arc-bracket spanning the two bars. The
    # bracket sits well clear of the AYB-bar label (which ends near 100)
    # and the title is suppressed in favour of an axis label so the
    # bracket-plus-text region has room.
    bracket_top = 122
    ax_top.plot([bar_x[0], bar_x[0], bar_x[1], bar_x[1]],
                [bracket_top - 4, bracket_top, bracket_top, bracket_top - 4],
                color="black", linewidth=0.9, clip_on=False)
    ax_top.text((bar_x[0] + bar_x[1]) / 2, bracket_top + 2,
                f"{lift:.1f}× lift",
                ha="center", va="bottom", fontsize=12, fontweight="bold",
                clip_on=False)

    ax_top.set_xticks(bar_x)
    ax_top.set_xticklabels(["S. stenocarpa\nchromosome-scale\n(Shorinola 2024)",
                            "Cowpea v1.0\nproxy\n(DArT report)"],
                           fontsize=10)
    ax_top.set_ylabel("DArTseq tags anchored (%)")
    ax_top.set_ylim(0, 135)
    ax_top.set_yticks([0, 25, 50, 75, 100])
    ax_top.set_title(f"Overall re-anchoring of {n_total:,} DArTseq tags",
                     pad=22)

    # ---- Bottom panel: per-AYB-chromosome side-by-side ----
    chroms = per_chrom["chrom"].tolist()
    x = np.arange(len(chroms))
    w = 0.4
    ax_bot.bar(x - w/2, per_chrom["n_ayb_anchored"], width=w,
               color=WONG["vermillion"], edgecolor="white", linewidth=0.4,
               label=f"S. stenocarpa anchor (total {n_ayb:,})")
    ax_bot.bar(x + w/2, per_chrom["n_cowpea_anchored"], width=w,
               color=WONG["skyblue"], edgecolor="white", linewidth=0.4,
               label=f"Cowpea proxy on same tags (total {n_cw:,})")

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(chroms)
    ax_bot.set_xlabel("S. stenocarpa pseudo-chromosome")
    ax_bot.set_ylabel("Number of DArTseq tags")
    ax_bot.set_title("Per-chromosome anchoring under each reference (same tags, both estimators)")
    ax_bot.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)

    return {
        "n_total": n_total,
        "n_ayb": n_ayb,
        "n_cw": n_cw,
        "pct_ayb": pct_ayb,
        "pct_cw": pct_cw,
        "lift": lift,
    }


def main() -> None:
    print(f"reading anchoring table from {ANCHOR_CSV}")
    df = load_anchoring()
    print(f"  {len(df)} markers in the raw DArT report")

    print("computing per-chromosome counts")
    per_chrom = per_chrom_counts(df)
    print(per_chrom.to_string(index=False))

    print("rendering Figure 1C")
    out_path = PUB_FIG / "fig01c_reanchor_callout"
    summary = plot_two_panel(df, per_chrom, out_path)

    summary_csv = PUB_TAB / "Table_S6_reanchoring_per_marker_summary.csv"
    pd.concat([
        per_chrom,
        pd.DataFrame([{
            "chrom": "TOTAL",
            "n_ayb_anchored": summary["n_ayb"],
            "n_cowpea_anchored": summary["n_cw"],
            "pct_cowpea_of_ayb": 100.0 * summary["n_cw"] / summary["n_ayb"],
        }]),
    ], ignore_index=True).to_csv(summary_csv, index=False)
    print(f"  wrote {summary_csv}")

    print("\n=== manuscript-text summary ===")
    print(f"raw DArTseq tags                : {summary['n_total']:,}")
    print(f"AYB-anchored                    : {summary['n_ayb']:,} ({summary['pct_ayb']:.1f} %)")
    print(f"cowpea-anchored (DArT report)   : {summary['n_cw']:,} ({summary['pct_cw']:.1f} %)")
    print(f"lift                            : {summary['lift']:.2f}x")


if __name__ == "__main__":
    main()
