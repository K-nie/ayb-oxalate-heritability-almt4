#!/usr/bin/env python3
"""
REVISION Phase 4 — genome-wide LD decay and the data-driven candidate
window around the focal ALMT4 SNP (reviewer L277-279, L498-499).

LD decay was already computed (results/91_ld_decay_genomewide_ayb,
results/96_ld_decay_dedup_sscorr). Here we:
  1. Re-fit the Hill-Weir (1988) / Remington (2001) expected-r² decay to
     the binned r² (primary = the dedup, sample-size-corrected panel,
     results/96) and report the genome-wide half-decay distance and the
     distance at which r² falls below 0.1, cross-checking the stored
     values.
  2. Draw a clean, publication-grade LD-decay figure (Nature style:
     lowercase panel letters handled at assembly; no title on the image).
  3. List EVERY annotated gene within the LD-based window around the focal
     SNP Ss10:15,394,673 — using the data-driven half-decay distance and
     the r²<0.1 distance as the two window radii, NOT an assumed +/-50 kb.

Inputs : results/96_.../ld_genomewide_dedup_sscorr_bins.csv + _summary.csv
         results/96_.../ld_per_chrom_dedup_sscorr_summary.csv (Ss10 radius)
         refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3
Output : results/R4_ld_decay_focal/
Author : Benjamin Narh-Madey
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from _plotstyle import apply, WONG
apply()
warnings.filterwarnings("ignore")

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
BINS = ROOT / "results" / "96_ld_decay_dedup_sscorr" / "tables" / "ld_genomewide_dedup_sscorr_bins.csv"
GSUM = ROOT / "results" / "96_ld_decay_dedup_sscorr" / "tables" / "ld_genomewide_dedup_sscorr_summary.csv"
CSUM = ROOT / "results" / "96_ld_decay_dedup_sscorr" / "tables" / "ld_per_chrom_dedup_sscorr_summary.csv"
GFF = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"
OUT = ROOT / "results" / "R4_ld_decay_focal"
TAB = OUT / "tables"; FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

FOCAL_CHR = "Ss10"
FOCAL_POS = 15_394_673
ALMT4_START, ALMT4_END = 15_397_847, 15_400_688
N_PANEL = 93   # dedup panel size for the Hill-Weir n term

# ---------------------------------------------------------------------------
bins = pd.read_csv(BINS)
gsum = pd.read_csv(GSUM).iloc[0]
csum = pd.read_csv(CSUM)
d = bins["mid_bp"].values.astype(float)
r2 = bins["mean"].values.astype(float)


def hill_weir(d_bp, rho):
    """Expected r² under drift-recombination equilibrium (Hill & Weir 1988;
    Remington et al. 2001), with the finite-sample correction, n = N_PANEL.
    C = rho * d (d in Mb); returns E[r²]."""
    C = rho * (d_bp / 1e6)
    n = N_PANEL
    term1 = (10 + C) / ((2 + C) * (11 + C))
    term2 = 1 + ((3 + C) * (12 + 12 * C + C ** 2)) / (n * (2 + C) * (11 + C))
    return term1 * term2


popt, _ = curve_fit(hill_weir, d, r2, p0=[1.0], maxfev=20000)
rho_fit = float(popt[0])
grid = np.linspace(1.0, d.max(), 40000)
fit = hill_weir(grid, rho_fit)
baseline = fit[-1]                       # asymptotic r² (background)
# Half-decay = distance at which the fitted r² falls to half of its
# extrapolated value at distance 0 (the standard convention used by the
# project's LD-decay scripts and by Remington et al. 2001).
r2_at_zero = float(hill_weir(np.array([0.0]), rho_fit)[0])
half_level = r2_at_zero / 2.0
below_half = np.where(fit < half_level)[0]
half_decay_bp = float(grid[below_half[0]]) if len(below_half) else float("nan")
# distance at which fitted r² < 0.1
below01 = np.where(fit < 0.1)[0]
dist_r2_01_bp = float(grid[below01[0]]) if len(below01) else float("nan")

ss10 = csum.loc[csum["chr"] == FOCAL_CHR].iloc[0]
summary = {
    "panel": str(gsum.get("panel", "dedup n=93")),
    "n_markers": int(gsum["n_markers"]),
    "n_pairs": int(gsum["n_pairs"]),
    "rho_hat_refit": rho_fit,
    "rho_hat_stored": float(gsum["rho_hat"]),
    "half_decay_kb_refit": half_decay_bp / 1e3,
    "half_decay_kb_stored": float(gsum["half_decay_kb"]),
    "dist_r2_below_0.1_kb_refit": dist_r2_01_bp / 1e3,
    "dist_r2_below_0.1_kb_stored": float(gsum["dist_r2_below_0.1_kb"]),
    "Ss10_half_decay_kb": float(ss10["half_decay_kb"]),
    "Ss10_dist_r2_below_0.1_kb": float(ss10["dist_r2_below_0.1_kb"]),
    "baseline_r2": float(baseline),
}
pd.DataFrame([summary]).to_csv(TAB / "ld_halfdecay_summary.csv", index=False)
print("[ld] genome-wide half-decay (refit) = %.1f kb (stored %.1f kb)"
      % (half_decay_bp / 1e3, float(gsum["half_decay_kb"])))
print("[ld] r2<0.1 distance (refit) = %.1f kb (stored %.1f kb)"
      % (dist_r2_01_bp / 1e3, float(gsum["dist_r2_below_0.1_kb"])))
print("[ld] Ss10 half-decay = %.1f kb ; Ss10 r2<0.1 = %.1f kb"
      % (float(ss10["half_decay_kb"]), float(ss10["dist_r2_below_0.1_kb"])))

# ---------------------------------------------------------------------------
# LD-decay figure (no title; single panel)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.6, 4.2), constrained_layout=True)
ax.scatter(d / 1e3, r2, s=26, color=WONG["skyblue"], edgecolor="black",
           linewidth=0.3, alpha=0.85, zorder=2, label="Binned mean $r^2$")
ax.plot(grid / 1e3, fit, color=WONG["vermillion"], lw=1.8, zorder=3,
        label="Hill-Weir fit")
ax.axhline(0.1, color="grey", ls=":", lw=0.8)
ax.axvline(half_decay_bp / 1e3, color="black", ls="--", lw=0.9)
ax.text(half_decay_bp / 1e3 + 6, ax.get_ylim()[1] * 0.82,
        f"half-decay\n{half_decay_bp/1e3:.0f} kb", fontsize=8)
ax.set_xlabel("Distance between SNPs (kb)")
ax.set_ylabel("Linkage disequilibrium ($r^2$)")
ax.set_xlim(0, min(1500, d.max() / 1e3))
ax.legend(frameon=False, fontsize=9)
fig.savefig(FIG / "ld_decay_genomewide.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "ld_decay_genomewide.pdf", bbox_inches="tight")
plt.close(fig)
print(f"[fig] ld_decay_genomewide -> {FIG}")

# ---------------------------------------------------------------------------
# Genes within the LD-based window around the focal SNP
# ---------------------------------------------------------------------------
def parse_attr(s, key):
    for part in s.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return ""


genes = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or f[0] != FOCAL_CHR:
            continue
        start, end = int(f[3]), int(f[4])
        attrs = f[8]
        gid = parse_attr(attrs, "ID")
        # Name field may appear twice (locus name + functional name); take last
        names = [parse_attr(p + ";", "Name") for p in attrs.split(";") if p.startswith("Name=")]
        name = names[-1] if names else ""
        genes.append({"chr": f[0], "start": start, "end": end, "strand": f[6],
                      "gene_id": gid, "name": name})
gdf = pd.DataFrame(genes)
# signed distance from focal SNP to nearest edge of each gene (0 if overlapping)
def sdist(row):
    if row["start"] <= FOCAL_POS <= row["end"]:
        return 0
    return int(min(abs(row["start"] - FOCAL_POS), abs(row["end"] - FOCAL_POS)))
gdf["abs_dist_bp"] = gdf.apply(sdist, axis=1)

for radius_kb, tag in [(float(ss10["half_decay_kb"]), "half_decay"),
                       (float(ss10["dist_r2_below_0.1_kb"]), "r2_below_0.1")]:
    rad = radius_kb * 1e3
    win = gdf[gdf["abs_dist_bp"] <= rad].sort_values("abs_dist_bp").copy()
    win["window_radius_kb"] = round(radius_kb, 1)
    win["window_def"] = tag
    win["is_ALMT4_2"] = win["gene_id"] == "AYBTSS11_029726"
    win.to_csv(TAB / f"focal_window_genes_{tag}.csv", index=False)
    print(f"\n[window:{tag}] radius +/-{radius_kb:.0f} kb -> {len(win)} genes "
          f"on {FOCAL_CHR} ({FOCAL_POS-rad:,.0f}-{FOCAL_POS+rad:,.0f})")
    print(f"  ALMT4_2 present: {bool(win['is_ALMT4_2'].any())}; "
          f"nearest 6 genes:")
    print(win.head(6)[["gene_id", "name", "start", "end", "abs_dist_bp"]].to_string(index=False))

print(f"\n[done] Phase-4 outputs in {OUT}")
