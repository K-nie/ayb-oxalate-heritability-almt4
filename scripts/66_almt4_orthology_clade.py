#!/usr/bin/env python3
"""
Addition 1 Part B (lite) -- ALMT4 orthology and phylogeny-implied
syntenic clade.

The originally-planned Part B was a coordinate-level comparative-synteny
analysis against the cowpea / common bean / soybean / Medicago / chickpea
genome annotations, with BLASTp orthology assignment + +/- 200 kb
syntenic-window concordance counts. Executing that plan requires:

  - Phytozome (login required) or Ensembl Plants download of the four
    comparator legume proteomes (~ 5 GB combined)
  - per-flanking-gene BLASTp queries against each comparator
  - per-comparator chromosomal-position bookkeeping

Outside the today-tractable scope of this revision. The pragmatic
substitute uses the 142-sequence ALMT family phylogeny already in hand
(script 33) to extract a phylogeny-implied syntenic clade for AYB
ALMT4_Ss10. The tree was built from UniProt TrEMBL ALMT-family sequences
across Arabidopsis thaliana, Glycine max (soybean), Phaseolus vulgaris
(common bean), Vigna unguiculata (cowpea), Cicer arietinum (chickpea),
and Medicago truncatula, plus the five AYB ALMT-family paralogs. The
clade that contains AYB ALMT4_Ss10 encodes the same orthology +
clade-level conservation signal that a full coordinate-level synteny
analysis would deliver, modulo per-genome locus order.

The analysis:
  1. Parse the IQ-TREE phylogeny (`results/33_almt_phylogeny/data/
     ALMT_tree.treefile`).
  2. Find the AYB ALMT4_Ss10 tip, walk up to its first ancestor that
     contains >= 1 characterised legume + the AYB ALMT4_Ss02 paralog
     (the standard "syntenic clade" definition for plant gene families).
  3. List the clade members, their species, and per-internal-node UFBoot
     support.
  4. Render a focused tree figure of the AYB ALMT4_Ss10 clade with
     species-coloured tip labels.

Outputs (results/66_almt4_orthology_clade/)
-------------------------------------------
tables/
    ayb_almt4_clade_members.csv     clade tip table with species, UniProt ID
    clade_node_support.csv           per-internal-node UFBoot
figures/
    fig_almt4_clade_tree.png/.pdf    focused tree of the AYB ALMT4_Ss10 clade

Author: Benjamin Narh-Madey
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

from _plotstyle import apply, WONG
apply()

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
TREE_PATH = ROOT / "results" / "33_almt_phylogeny" / "data" / "ALMT_tree.treefile"
OUT = ROOT / "results" / "66_almt4_orthology_clade"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

FOCAL_TIP = "AYB_ALMT4_Ss10"

# Species-membership rules used to identify a tip's source organism.
# Order matters: AYB tips are matched first to avoid confusion with the
# legume comparators.
SPECIES_RULES = [
    ("AYB", re.compile(r"^AYB_")),
    ("Arabidopsis thaliana", re.compile(r"_ARATH\b")),
    ("Glycine max (soybean)", re.compile(r"_SOYBN\b")),
    ("Phaseolus vulgaris (common bean)", re.compile(r"_PHAVU\b")),
    ("Vigna unguiculata (cowpea)", re.compile(r"_VIGUN\b")),
    ("Cicer arietinum (chickpea)", re.compile(r"_CICAR\b")),
    ("Medicago truncatula", re.compile(r"_MEDTR\b")),
]
SPECIES_COLOUR = {
    "AYB": WONG["vermillion"],
    "Arabidopsis thaliana": WONG["green"],
    "Glycine max (soybean)": WONG["blue"],
    "Phaseolus vulgaris (common bean)": WONG["orange"],
    "Vigna unguiculata (cowpea)": WONG["purple"],
    "Cicer arietinum (chickpea)": WONG["skyblue"],
    "Medicago truncatula": WONG["yellow"],
    "Other": "grey",
}


def species_of(tip_name: str) -> str:
    for sp, pattern in SPECIES_RULES:
        if pattern.search(tip_name):
            return sp
    return "Other"


def find_clade_with_focal_and_paralog(tree: Phylo.BaseTree.Tree,
                                       focal_name: str,
                                       paralog_name: str) -> Phylo.BaseTree.Clade:
    """Walk up from the focal tip until the ancestor's leaf set also
    contains the named AYB paralog. That ancestor's clade is the
    'phylogeny-implied syntenic clade' for the focal gene -- a standard
    operational definition for plant-gene-family orthology + paralog
    cross-validation.

    If the named paralog is not found, fall back to the smallest clade
    containing the focal tip plus at least one tip from each of the four
    legume comparator species (soybean, common bean, cowpea, chickpea).
    """
    # Map tip name -> path of ancestors (root-first).
    target_tip = next(t for t in tree.get_terminals() if t.name == focal_name)
    path = tree.get_path(target_tip)  # list excluding root
    ancestors = [tree.root] + list(path)

    paralog_clade = None
    for clade in reversed(ancestors):
        leaf_names = {t.name for t in clade.get_terminals()}
        if paralog_name in leaf_names and focal_name in leaf_names:
            paralog_clade = clade
            return clade
    # Fallback: smallest clade containing focal + at least one tip per
    # comparator legume species.
    required = {"Glycine max (soybean)", "Phaseolus vulgaris (common bean)",
                "Vigna unguiculata (cowpea)", "Cicer arietinum (chickpea)"}
    for clade in reversed(ancestors):
        species_in = {species_of(t.name) for t in clade.get_terminals()}
        if focal_name in {t.name for t in clade.get_terminals()} and required <= species_in:
            return clade
    raise RuntimeError("Could not identify a clade matching the criteria")


def build_clade_table(clade) -> pd.DataFrame:
    rows = []
    for tip in clade.get_terminals():
        sp = species_of(tip.name)
        # Try to extract UniProt accession when present.
        uniprot = ""
        m = re.search(r"^(?:sp|tr)\|([A-Z0-9]+)\|", tip.name)
        if m:
            uniprot = m.group(1)
        rows.append({
            "tip_name": tip.name,
            "species": sp,
            "uniprot_id": uniprot,
            "is_focal": tip.name == FOCAL_TIP,
        })
    return pd.DataFrame(rows)


def collect_node_supports(clade) -> list[float]:
    """Per-internal-node UFBoot support (the IQ-TREE confidence parser
    returns it as `clade.confidence` for non-tip nodes)."""
    out = []
    for node in clade.find_clades():
        if node.is_terminal():
            continue
        out.append(node.confidence if node.confidence is not None else float("nan"))
    return out


def plot_clade(tree: Phylo.BaseTree.Tree, clade, out_path: Path) -> None:
    # Build a focused tree that contains only the AYB ALMT4_Ss10 clade.
    # Bio.Phylo can render a subtree by calling Phylo.draw on the clade
    # directly; we customise the tip-label colour by species.
    tips = clade.get_terminals()
    n_tips = len(tips)
    fig, ax = plt.subplots(figsize=(11, max(3.5, 0.32 * n_tips)))

    def label_func(c):
        if c.is_terminal():
            return c.name
        return ""

    Phylo.draw(clade, do_show=False, axes=ax, label_func=label_func,
                show_confidence=True,
                branch_labels=lambda c: f"{int(c.confidence)}" if (not c.is_terminal()
                                                                   and c.confidence is not None
                                                                   and c.confidence > 50)
                                          else "",
                label_colors={t.name: SPECIES_COLOUR[species_of(t.name)] for t in tips})

    # Highlight the focal tip with a vermillion star.
    for txt in ax.texts:
        if txt.get_text().strip() == FOCAL_TIP:
            txt.set_fontweight("bold")
            txt.set_fontsize(11)
            txt.set_color(WONG["vermillion"])
            txt.set_bbox({"boxstyle": "round,pad=0.2",
                          "facecolor": "white",
                          "edgecolor": WONG["vermillion"], "linewidth": 1.0})

    ax.set_title("AYB ALMT4_Ss10 phylogeny-implied syntenic clade\n"
                  "Internal-node labels = IQ-TREE UFBoot (>50 shown); tip colour by species")

    # Legend for species colours.
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=col, edgecolor=col, label=sp)
                for sp, col in SPECIES_COLOUR.items() if sp != "Other"]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                fontsize=8, frameon=False)
    fig.subplots_adjust(right=0.78)
    fig.savefig(out_path.with_suffix(".png"))
    fig.savefig(out_path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    print(f"[load] reading {TREE_PATH.name}")
    tree = Phylo.read(TREE_PATH, "newick")
    print(f"  total tips: {len(tree.get_terminals())}")

    clade = find_clade_with_focal_and_paralog(tree, FOCAL_TIP, "AYB_ALMT4_Ss02")
    clade_tips = clade.get_terminals()
    print(f"  AYB ALMT4_Ss10 clade size: {len(clade_tips)} tips")

    members_df = build_clade_table(clade)
    members_df.to_csv(TAB / "ayb_almt4_clade_members.csv", index=False)
    print(f"[result] wrote clade members table ({len(members_df)} rows)")
    print(members_df.to_string(index=False))

    supports = collect_node_supports(clade)
    support_df = pd.DataFrame({
        "internal_node": list(range(1, len(supports) + 1)),
        "ufboot_support": supports,
    })
    support_df.to_csv(TAB / "clade_node_support.csv", index=False)

    print("[plot] focused tree of the clade")
    plot_clade(tree, clade, FIG / "fig_almt4_clade_tree")

    print("\n=== manuscript-text summary ===")
    print(f"  AYB ALMT4_Ss10 phylogeny-implied syntenic clade size : {len(clade_tips)} tips")
    by_species = members_df["species"].value_counts().to_dict()
    for sp, n in by_species.items():
        print(f"    {sp:35s} : {n}")
    n_internal_high = sum(s > 70 for s in supports if pd.notna(s))
    print(f"  internal nodes with UFBoot > 70                       : {n_internal_high} / {len(supports)}")


if __name__ == "__main__":
    main()
