#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Render the OUTGROUP-ROOTED ALMT-family phylogeny (manuscript Fig 8).
#
# Revision fix (TPG-2026-09-0323, Reviewer 1 L523-573, Reviewer 2 #21):
# the original tree (script 33) was MIDPOINT-rooted, i.e. effectively
# unrooted with an arbitrary implied root. Here the 142 ingroup ALMT
# proteins were re-aligned together with three early-diverging land-plant
# ALMT-family outgroup proteins (Selaginella moellendorffii D8TB43,
# Physcomitrium patens A0A7I4E3F8, Marchantia polymorpha A0A2R6WEU1; all
# Pfam PF11744) and re-inferred under the published model (LG+G4, 1000
# UFBoot). The three outgroup sequences did NOT form a clade (IQ-TREE:
# "Branch separating outgroup is not found") -- they are divergent ALMT
# paralogs that attach at separate deep branches (the moss sequence is
# long-branched), a known property of lineage-specific ALMT expansions.
# A monophyletic-outgroup rooting is therefore not available, so the tree
# is rooted by MINIMAL ANCESTOR DEVIATION (MAD; Tria, Landan & Dagan 2017,
# Nat Ecol Evol 1:0193) -- a model-free rooting that does not assume an
# outgroup. The MAD root is independent of the midpoint root and the focal
# placement is identical under both.
#
# Input  : results/33_almt_phylogeny/data/reroot/ALMT_tree_og_madroot.treefile
# Output : results/33_almt_phylogeny/figures/fig72_almt_family_tree_MADrooted.{png,pdf}
#
# No title is baked into the image (Benjamin figure rule); panel text is
# sans-serif; colours are Okabe-Ito colour-vision-safe.
#
# Author: Benjamin Narh-Madey
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
})

ROOT <- "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
TREE_FILE <- file.path(ROOT,
  "results/33_almt_phylogeny/data/reroot/ALMT_tree_og_madroot.treefile")
FIG <- file.path(ROOT, "results/33_almt_phylogeny/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

cat("[load] MAD-rooted tree...\n")
tree <- read.tree(TREE_FILE)   # already MAD-rooted (scripts/mad_write.py)
cat(sprintf("  %d tips, rooted=%s\n", length(tree$tip.label), is.rooted(tree)))

classify_tip <- function(label) {
  if (grepl("SELML|PHYPA|MARPO", label)) return("Outgroup (basal land plant)")
  if (grepl("^AYB_", label)) return("AYB")
  if (grepl("ARATH", label)) return("Arabidopsis")
  if (grepl("SOYBN", label, ignore.case = TRUE)) return("Soybean")
  if (grepl("PHAVU", label, ignore.case = TRUE)) return("Common bean")
  if (grepl("VIGUN", label, ignore.case = TRUE)) return("Cowpea")
  if (grepl("CICAR", label, ignore.case = TRUE)) return("Chickpea")
  if (grepl("MEDTR", label, ignore.case = TRUE)) return("Medicago")
  return("Other")
}
species <- sapply(tree$tip.label, classify_tip)

COLORS <- c(AYB = "#D55E00", Arabidopsis = "#0072B2",
            Soybean = "#009E73", `Common bean` = "#CC79A7",
            Cowpea = "#E69F00", Chickpea = "#56B4E9",
            Medicago = "#F0E442", `Outgroup (basal land plant)` = "#000000",
            Other = "#999999")

# Tidy display labels; name the characterized Arabidopsis ALMTs explicitly.
ARATH_NAMED <- c(Q9SJE9 = "AtALMT1", Q9SJE8 = "AtALMT2", Q9LPQ8 = "AtALMT3",
                 Q9C6L8 = "AtALMT4", Q93Z29 = "AtALMT5", Q9SHM1 = "AtALMT6",
                 Q9XIN1 = "AtALMT7", Q9SRM9 = "AtALMT8", Q9LS46 = "AtALMT9")
display_label <- function(lbl) {
  if (grepl("^AYB_", lbl)) return(sub("_", " ", sub("AYB_", "", lbl)))
  m <- regmatches(lbl, regexec("^[a-z]+\\|([A-Z0-9]+)\\|([A-Z0-9_]+)", lbl))[[1]]
  if (length(m) >= 3) {
    acc <- m[2]; prot <- m[3]
    if (acc %in% names(ARATH_NAMED)) return(ARATH_NAMED[[acc]])
    return(acc)
  }
  return(lbl)
}
tree$tip.label_display <- sapply(tree$tip.label, display_label)

tip_df <- data.frame(label = tree$tip.label,
                     display = tree$tip.label_display,
                     species = species,
                     is_ayb = species == "AYB",
                     stringsAsFactors = FALSE)

bs_vals <- suppressWarnings(as.numeric(tree$node.label))
n_tips <- length(tree$tip.label)
bs_df <- data.frame(
  node = (n_tips + 1):(n_tips + tree$Nnode),
  ufboot = bs_vals,
  ufboot_label = ifelse(!is.na(bs_vals) & bs_vals >= 70,
                        as.character(round(bs_vals)), NA_character_))

p <- ggtree(tree, layout = "fan", open.angle = 18,
            size = 0.55, colour = "grey30") %<+% tip_df %<+% bs_df +
  geom_tippoint(aes(colour = species, size = is_ayb), show.legend = TRUE) +
  scale_size_manual(values = c(`TRUE` = 4.6, `FALSE` = 2.0), guide = "none") +
  geom_tiplab(aes(label = display, colour = species,
                  fontface = ifelse(is_ayb, "bold", "plain")),
              size = 3.3, offset = 0.03, show.legend = FALSE, family = "sans") +
  scale_colour_manual(values = COLORS, name = "Lineage",
                      breaks = names(COLORS),
                      guide = guide_legend(override.aes = list(size = 5.5))) +
  geom_text2(aes(label = ufboot_label, subset = !is.na(ufboot_label)),
             size = 2.4, colour = "grey30", hjust = 1.15, vjust = -0.35) +
  theme(legend.position = c(0.93, 0.13),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(FIG, "fig72_almt_family_tree_MADrooted.png"), p,
       width = 18, height = 18, dpi = 300, bg = "white")
ggsave(file.path(FIG, "fig72_almt_family_tree_MADrooted.pdf"), p,
       width = 18, height = 18, bg = "white")
cat(sprintf("[fig] fig72_almt_family_tree_MADrooted written to %s\n", FIG))
