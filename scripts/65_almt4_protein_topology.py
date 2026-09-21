#!/usr/bin/env python3
"""
ALMT4 protein-domain confirmation and transmembrane topology.

Part A of the three-part ALMT4 mechanism-evolution chain detailed in
manuscript/paper_A2/analysis_plans_A2.md Addition 1. Parts B (comparative
synteny) and C (dN/dS on the ALMT4 clade) are scoped in separate scripts.

The case for ALMT4 as the Soluble_Oxalate candidate gene currently rests
on (i) GWAS proximity at -3 kb of the focal SNP, (ii) Bayesian + REML h^2
triangulation, (iii) GBLUP holdout r = 0.54, and (iv) the 142-sequence
ALMT-family phylogeny. The phylogeny call is orthology only. This script
adds the structural-biology layer: a Pfam-domain annotation (already
present in the Funannotate GFF and extracted here as the primary domain
ID) plus a Kyte & Doolittle 1982 hydrophobicity-based transmembrane scan
that predicts the canonical 4-6 TM segments characteristic of ALMT-family
vacuolar transporters (Sasaki et al. 2004; Hoekenga et al. 2006).

The output topology figure stacks (i) per-residue Kyte-Doolittle
hydrophobicity profile, (ii) predicted TM segments as shaded bars, and
(iii) the PF11744 ALMT family domain span and PF13515 FUSC-like /
Membrane-protein-superfamily span pulled from the Funannotate GFF.

Outputs (results/65_almt4_protein/)
-----------------------------------
data/
    almt4_protein.fasta             extracted AYB ALMT4_2 protein
tables/
    pfam_interpro_domains.csv        domain ID, source, span (from GFF)
    tm_segments.csv                   predicted TM start/end + mean score
    hydrophobicity_per_residue.csv   per-position Kyte-Doolittle score
figures/
    fig_almt4_protein_topology.png/.pdf  composite topology figure

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
PROT_FAA = ROOT / "refs" / "ayb_genome" / "AYB_proteins.faa"
GFF = ROOT / "refs" / "ayb_genome" / "Sphenostylis_stenocarpa_Funannotate.gff3"

OUT = ROOT / "results" / "65_almt4_protein"
DATA = OUT / "data"
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (DATA, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

# Focal ALMT4_2 gene -- the AYB gene that the Soluble_Oxalate top SNP
# at Ss10:15,394,673 lies -3,174 bp upstream of.
GENE_ID = "AYBTSS11_029726"
TRANSCRIPT_ID = "AYBTSS11_029726-T1"
TARGET_CHROM = "Ss10"

# ---------------------------------------------------------------------------
# Kyte-Doolittle hydrophobicity values for the 20 amino acids.
# Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132. Higher = more
# hydrophobic. TM segments are predicted as contiguous windowed scores
# above the threshold.
# ---------------------------------------------------------------------------
KD = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
    "X":  0.0, "U":  2.5, "B": -3.5, "Z": -3.5,
}
TM_WINDOW = 19
# Threshold tuned for window = 19. The Kyte-Doolittle (1982) original
# threshold of 1.6 was for a 7-residue window; smoothing over 19 residues
# dampens the per-position peak. We use a permissive 0.8 entry rule
# (catches plant-protein TM helices, which often have more polar residues
# than animal counterparts; Cao et al. 2013 review on plant membrane
# proteins) and require the segment mean to clear 1.2 so the call is not
# driven by isolated high-hydrophobicity peaks. TM helices in the ALMT
# family are typically 17-25 residues (Sasaki 2004; Hoekenga 2006) but
# can be as short as 14 in flexible regions, so we use a 14 aa floor.
TM_ENTRY_THRESHOLD = 0.8
TM_SEGMENT_MEAN_MIN = 1.2
TM_MIN_LENGTH = 14


def load_protein() -> str:
    """Return the AYB ALMT4_2 amino-acid sequence as a string."""
    for rec in SeqIO.parse(PROT_FAA, "fasta"):
        if TRANSCRIPT_ID in rec.id or TRANSCRIPT_ID in rec.description:
            return str(rec.seq)
    raise RuntimeError(f"Protein {TRANSCRIPT_ID} not found in {PROT_FAA}")


def parse_gff_attributes(attr_string: str) -> dict[str, str]:
    """Parse a GFF3 attribute column into a dict. Handles repeated keys
    (e.g. Name=...,Name=...) by joining with ';'."""
    out: dict[str, list[str]] = {}
    for kv in attr_string.strip().rstrip(";").split(";"):
        kv = kv.strip()
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        out.setdefault(k, []).append(v)
    return {k: ";".join(vs) for k, vs in out.items()}


def find_almt4_gff_record() -> dict:
    """Locate the ALMT4_2 mRNA / CDS in the Funannotate GFF and return its
    parsed attributes plus chrom / position."""
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "mRNA":
                continue
            adict = parse_gff_attributes(attrs)
            mrna_id = adict.get("ID", "")
            parent = adict.get("Parent", "")
            if (TRANSCRIPT_ID == mrna_id
                    or GENE_ID == parent
                    or TRANSCRIPT_ID in mrna_id
                    or GENE_ID in parent):
                return {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attributes": adict,
                }
    raise RuntimeError(f"ALMT4_2 (gene {GENE_ID}) mRNA not found in GFF")


def extract_pfam_interpro(adict: dict[str, str], seq_length: int) -> pd.DataFrame:
    """Pull Pfam and InterPro accessions from the Funannotate GFF attribute
    block. Funannotate runs HMMER+Pfam-A at annotation time and writes the
    hit accessions into the Dbxref field, so we do not re-run HMMER here.

    Notes on span. Funannotate emits accessions without per-domain coords;
    the canonical Pfam architecture for ALMT-family proteins is the PF11744
    N-terminal ALMT domain spanning roughly the first half of the protein
    (Sasaki et al. 2004), with PF13515 FUSC-like domain in the C-terminal
    half. We label the spans as "full-length annotated" and flag the
    canonical N / C-terminal halves on the figure rather than mis-claiming
    HMMER coordinates we do not have.
    """
    dbxref = adict.get("Dbxref", "")
    ontology = adict.get("Ontology_term", "")

    rows = []
    # GFF3 spec: multiple values for one attribute key (here Dbxref) are
    # comma-separated within the value string; the semicolon separates
    # different attribute keys, not values inside one key.
    for token in dbxref.split(","):
        token = token.strip()
        if not token:
            continue
        # token like "PFAM:PF11744" or "InterPro:IPR020966"
        m = re.match(r"^(PFAM|Pfam|InterPro):(\w+)$", token)
        if not m:
            continue
        source, accession = m.group(1), m.group(2)
        # Canonical ALMT-family domain architecture: PF11744 N-terminal,
        # PF13515 C-terminal. Use protein length as a span proxy on the
        # figure; flag this explicitly in the description.
        if accession == "PF11744":
            span_start, span_end = 1, seq_length // 2
            label = "PF11744 (ALMT family domain; N-terminal half)"
        elif accession == "PF13515":
            span_start, span_end = seq_length // 2, seq_length
            label = "PF13515 (FUSC-like / Membrane protein superfamily; C-terminal half)"
        elif accession.startswith("IPR"):
            span_start, span_end = 1, seq_length
            label = f"{accession} (InterPro: ALMT-family superfamily)"
        else:
            span_start, span_end = 1, seq_length
            label = f"{accession} (full-length annotated)"
        rows.append({
            "source": source,
            "accession": accession,
            "label": label,
            "span_start_aa": span_start,
            "span_end_aa": span_end,
            "comment": "Funannotate HMMER+Pfam-A; canonical ALMT architecture used for figure spans"
                       if source.upper() == "PFAM" else "",
        })
    return pd.DataFrame(rows)


def kyte_doolittle_profile(seq: str, window: int = TM_WINDOW) -> np.ndarray:
    """Centred-window mean Kyte-Doolittle hydrophobicity per residue. Edge
    positions where the window would extend past the sequence are filled
    with NaN so the figure visibly tapers at the ends."""
    half = window // 2
    n = len(seq)
    raw = np.array([KD.get(a.upper(), 0.0) for a in seq])
    profile = np.full(n, np.nan)
    for i in range(half, n - half):
        profile[i] = np.mean(raw[i - half:i + half + 1])
    return profile


def predict_tm_segments(profile: np.ndarray,
                        entry_threshold: float = TM_ENTRY_THRESHOLD,
                        segment_mean_min: float = TM_SEGMENT_MEAN_MIN,
                        min_length: int = TM_MIN_LENGTH,
                        split_above_length: int = 30) -> pd.DataFrame:
    """Predict transmembrane segments as contiguous runs of positions
    where the centred-window hydrophobicity exceeds entry_threshold, of
    minimum length min_length, AND whose mean hydrophobicity clears
    segment_mean_min.

    Segments longer than split_above_length residues are split at their
    internal local minimum. ALMT-family TM helices reach at most ~ 25 aa
    by crystal structure (Qiu et al. 2021 on AtALMT9); a 30+ aa
    hydrophobic stretch with a soft interior dip almost always reflects
    two fused helices whose connector did not drop below the entry
    threshold under a 19-residue window."""
    above = np.where(profile > entry_threshold, 1, 0)
    above = np.nan_to_num(above, nan=0).astype(int)

    raw_segments = []
    i = 0
    n = len(above)
    while i < n:
        if above[i] == 1:
            j = i
            while j < n and above[j] == 1:
                j += 1
            length = j - i
            if length >= min_length:
                segment_mean = float(np.nanmean(profile[i:j]))
                if segment_mean >= segment_mean_min:
                    raw_segments.append((i, j))
            i = j
        else:
            i += 1

    # Split any segment longer than the canonical ALMT TM length at its
    # interior local minimum.
    refined = []
    for start, end in raw_segments:
        if end - start > split_above_length:
            interior_slice = profile[start + 5:end - 5]
            if len(interior_slice) > 0:
                rel_min = int(np.nanargmin(interior_slice))
                split_at = start + 5 + rel_min
                # Only accept the split if both halves still meet the
                # min-length floor; otherwise leave the segment intact.
                left_len = split_at - start
                right_len = end - split_at
                if left_len >= min_length and right_len >= min_length:
                    refined.append((start, split_at))
                    refined.append((split_at, end))
                    continue
        refined.append((start, end))

    rows = []
    for k, (start, end) in enumerate(refined, start=1):
        rows.append({
            "tm_segment_number": k,
            "start_aa": start + 1,
            "end_aa": end,
            "length_aa": end - start,
            "mean_hydrophobicity": float(np.nanmean(profile[start:end])),
        })
    return pd.DataFrame(rows)


def plot_topology(seq_length: int,
                  profile: np.ndarray,
                  tm_segments: pd.DataFrame,
                  domains: pd.DataFrame,
                  out_path: Path) -> None:
    fig, (ax_hp, ax_tm, ax_dom) = plt.subplots(
        3, 1, figsize=(11.5, 6.5),
        gridspec_kw={"height_ratios": [3.0, 0.8, 0.8]},
        sharex=True,
    )

    # Top: hydrophobicity profile.
    x = np.arange(1, seq_length + 1)
    ax_hp.plot(x, profile, color=WONG["blue"], linewidth=1.0)
    ax_hp.fill_between(x, 0, profile,
                       where=(profile > TM_ENTRY_THRESHOLD),
                       color=WONG["vermillion"], alpha=0.30,
                       interpolate=True,
                       label=f"hydrophobic stretch (> {TM_ENTRY_THRESHOLD})")
    ax_hp.axhline(0, color="grey", linewidth=0.5)
    ax_hp.axhline(TM_ENTRY_THRESHOLD, color=WONG["vermillion"],
                  linewidth=0.8, linestyle="--",
                  label=f"entry threshold = {TM_ENTRY_THRESHOLD}")
    ax_hp.axhline(TM_SEGMENT_MEAN_MIN, color=WONG["vermillion"],
                  linewidth=0.6, linestyle=":",
                  label=f"segment-mean floor = {TM_SEGMENT_MEAN_MIN}")
    ax_hp.set_ylabel("Kyte-Doolittle\nhydrophobicity\n(window = 19)")
    ax_hp.set_title(
        f"AYB ALMT4_2 ({GENE_ID}, transcript {TRANSCRIPT_ID}, {seq_length} aa) -- "
        f"predicted topology"
    )
    ax_hp.legend(loc="upper right", fontsize=8)

    # Middle: TM segments as filled bars.
    ax_tm.set_yticks([])
    ax_tm.set_ylim(0, 1)
    ax_tm.set_ylabel("Predicted\nTM segments", rotation=0,
                      labelpad=55, ha="right", va="center")
    for _, row in tm_segments.iterrows():
        ax_tm.add_patch(plt.Rectangle((row["start_aa"], 0.2),
                                      row["end_aa"] - row["start_aa"], 0.6,
                                      facecolor=WONG["vermillion"],
                                      edgecolor="black", linewidth=0.5))
        ax_tm.text((row["start_aa"] + row["end_aa"]) / 2, 0.5,
                   f"TM{int(row['tm_segment_number'])}",
                   ha="center", va="center", fontsize=8, color="white",
                   fontweight="bold")
    ax_tm.spines["left"].set_visible(False)

    # Bottom: Pfam / InterPro domain bar.
    pfam_rows = domains.loc[domains["source"].str.upper() == "PFAM"]
    ax_dom.set_yticks([])
    ax_dom.set_ylim(0, len(pfam_rows) + 0.5)
    ax_dom.set_ylabel("Pfam\ndomains", rotation=0,
                      labelpad=55, ha="right", va="center")
    pfam_colours = [WONG["green"], WONG["purple"], WONG["orange"]]
    for idx, (_, row) in enumerate(pfam_rows.iterrows()):
        y = len(pfam_rows) - idx
        ax_dom.add_patch(plt.Rectangle(
            (row["span_start_aa"], y - 0.35),
            row["span_end_aa"] - row["span_start_aa"], 0.7,
            facecolor=pfam_colours[idx % len(pfam_colours)],
            edgecolor="black", linewidth=0.5, alpha=0.85,
        ))
        ax_dom.text((row["span_start_aa"] + row["span_end_aa"]) / 2, y,
                    row["label"], ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
    ax_dom.set_xlabel("Residue position (amino acid)")
    ax_dom.spines["left"].set_visible(False)
    ax_dom.set_xlim(1, seq_length)

    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    print(f"[load] reading AYB ALMT4_2 protein ({TRANSCRIPT_ID})")
    seq = load_protein()
    seq_length = len(seq)
    print(f"  length: {seq_length} aa")

    # Save the extracted sequence for the synteny + dN/dS follow-ups.
    fasta_out = DATA / "almt4_protein.fasta"
    with open(fasta_out, "w") as fh:
        fh.write(f">{TRANSCRIPT_ID} AYB ALMT4_2 Spste.TSs11.10G297260.1 "
                 f"Aluminum-activated malate transporter 4 (chr Ss10)\n")
        for i in range(0, seq_length, 60):
            fh.write(seq[i:i + 60] + "\n")
    print(f"  wrote {fasta_out}")

    print(f"[gff] locating ALMT4_2 in {GFF.name}")
    rec = find_almt4_gff_record()
    print(f"  Ss10:{rec['start']:,}-{rec['end']:,} ({rec['strand']} strand)")
    print(f"  product: {rec['attributes'].get('product', 'n/a')}")
    print(f"  Dbxref:  {rec['attributes'].get('Dbxref', 'n/a')}")

    print("[domain] extracting Pfam / InterPro from GFF attributes")
    domains = extract_pfam_interpro(rec["attributes"], seq_length)
    domains.to_csv(TAB / "pfam_interpro_domains.csv", index=False)
    print(f"  {len(domains)} domain hits")
    print(domains[["source", "accession", "label"]].to_string(index=False))

    print(f"[tm] computing Kyte-Doolittle hydrophobicity (window={TM_WINDOW})")
    profile = kyte_doolittle_profile(seq, window=TM_WINDOW)
    pd.DataFrame({
        "position": np.arange(1, seq_length + 1),
        "residue": list(seq),
        "kyte_doolittle_window19": profile,
    }).to_csv(TAB / "hydrophobicity_per_residue.csv", index=False)

    print(f"[tm] predicting TM segments (entry_threshold={TM_ENTRY_THRESHOLD}, "
          f"segment_mean_min={TM_SEGMENT_MEAN_MIN}, min_length={TM_MIN_LENGTH})")
    tm = predict_tm_segments(profile,
                              entry_threshold=TM_ENTRY_THRESHOLD,
                              segment_mean_min=TM_SEGMENT_MEAN_MIN,
                              min_length=TM_MIN_LENGTH)
    tm.to_csv(TAB / "tm_segments.csv", index=False)
    print(f"  predicted TM segments: {len(tm)}")
    print(tm.to_string(index=False))

    print("[plot] rendering composite topology figure")
    plot_topology(seq_length, profile, tm, domains,
                   FIG / "fig_almt4_protein_topology")

    print("\n=== manuscript-text summary ===")
    print(f"  AYB ALMT4_2 protein length             : {seq_length} aa")
    print(f"  Pfam domains (Funannotate)             : "
          f"{', '.join(domains.loc[domains['source'].str.upper() == 'PFAM', 'accession'])}")
    print(f"  predicted TM segments (Kyte-Doolittle) : {len(tm)}")
    if len(tm) > 0:
        print(f"  TM segment length range                : "
              f"{tm['length_aa'].min()}-{tm['length_aa'].max()} aa")
        print(f"  mean hydrophobicity in TM segments     : "
              f"{tm['mean_hydrophobicity'].mean():.2f}")


if __name__ == "__main__":
    main()
