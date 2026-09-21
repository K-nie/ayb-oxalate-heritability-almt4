#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# topGO GO-term enrichment for genes flanking the top-suggestive SNPs.
#
# topGO (Alexa et al. 2006 Bioinformatics) is the academic-standard GO
# enrichment tool. It accounts for the GO DAG hierarchy via the elim and
# weight01 algorithms, which down-weight parent terms once a child term
# explains the signal. This avoids the classic-Fisher inflation of redundant
# parent terms that the script-38 Python version produced.
#
# Three separate analyses, one per GO namespace (BP, MF, CC), each using
# elim + classic. Foreground = genes within +/- 10 kb of top-10 SNPs per
# trait across all 13 traits.
#
# Input:
#   Funannotate GFF gene table -> gene_id : list of GO term IDs (mRNA-level
#   Ontology_term field, dedup to one record per gene_id)
#   AYB SNP anchoring table -> defines the foreground gene set
#
# Output (results/38_go_enrichment/):
#   tables/go_topgo_<NAMESPACE>.csv     elim + classic p per term
#   figures/fig75_go_topgo_<NAMESPACE>.png/.pdf   horizontal bar with names
#
# Author: Benjamin Narh-Madey
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(topGO)
  library(GO.db)
})

ROOT <- "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
GFF  <- file.path(ROOT, "refs/ayb_genome/Sphenostylis_stenocarpa_Funannotate.gff3")
ANCHOR <- file.path(ROOT, "refs/ayb_genome/ayb_marker_anchoring.csv")
GWAS_BIO <- file.path(ROOT, "results/03_gwas_mlm/tables")
GWAS_NEW <- file.path(ROOT, "results/15_gwas_new_traits/tables")
OUT  <- file.path(ROOT, "results/38_go_enrichment")
FIG  <- file.path(OUT, "figures")
TAB  <- file.path(OUT, "tables")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

TOP_N <- 10        # top-N SNPs per trait
WIN_BP <- 10000    # +/- 10 kb of each focal SNP

TRAITS_BIO <- c("Tannin", "Phenol", "Flavonoid", "Antioxidant")
TRAITS_NEW <- c("Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
                "Seed_Coat_Tannin", "Crude_Protein", "Total_Oxalate",
                "Soluble_Oxalate", "Insoluble_Oxalate")


# ---------------------------------------------------------------------------
# Parse GFF for gene coordinates + per-gene GO sets
# ---------------------------------------------------------------------------
cat("[gff] parsing Funannotate gene table...\n")
gff_lines <- readLines(GFF)
gff_lines <- gff_lines[!grepl("^#", gff_lines)]
gff_lines <- gff_lines[nchar(gff_lines) > 0]

split_attr <- function(attr_str) {
  pairs <- strsplit(attr_str, ";", fixed = TRUE)[[1]]
  kv <- strsplit(pairs, "=", fixed = TRUE)
  ks <- sapply(kv, `[`, 1)
  vs <- sapply(kv, function(x) if (length(x) > 1) x[2] else NA)
  setNames(vs, ks)
}

# Defensive named-list lookup that returns NA when the key is absent.
get_attr <- function(a, k) {
  if (k %in% names(a)) a[[k]] else NA_character_
}

genes_chr <- character(); genes_start <- integer(); genes_end <- integer()
genes_id  <- character()
mrna_parent <- character(); mrna_go <- character()
for (line in gff_lines) {
  f <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(f) < 9) next
  a <- split_attr(f[9])
  if (f[3] == "gene") {
    gid <- get_attr(a, "ID")
    if (!is.na(gid)) {
      genes_id    <- c(genes_id, gid)
      genes_chr   <- c(genes_chr, f[1])
      genes_start <- c(genes_start, as.integer(f[4]))
      genes_end   <- c(genes_end, as.integer(f[5]))
    }
  } else if (f[3] == "mRNA") {
    pid <- get_attr(a, "Parent")
    ont <- get_attr(a, "Ontology_term")
    if (!is.na(pid) && !is.na(ont)) {
      mrna_parent <- c(mrna_parent, pid)
      mrna_go     <- c(mrna_go, ont)
    }
  }
}

gene_tbl <- data.frame(gene_id = genes_id, chr = genes_chr,
                       start = genes_start, end = genes_end,
                       stringsAsFactors = FALSE)

# Build gene_id -> character vector of GO terms (deduplicated)
go_split <- strsplit(mrna_go, ",", fixed = TRUE)
gene2go <- list()
for (i in seq_along(mrna_parent)) {
  gid <- mrna_parent[i]
  terms <- trimws(go_split[[i]])
  terms <- terms[grepl("^GO:", terms)]
  if (length(terms) > 0) {
    prior <- gene2go[[gid]]
    gene2go[[gid]] <- unique(c(prior, terms))
  }
}
cat(sprintf("  %d genes, %d with at least one GO term\n",
            nrow(gene_tbl), length(gene2go)))


# ---------------------------------------------------------------------------
# Foreground gene set: +/- WIN_BP of top-N anchored SNPs per trait
# ---------------------------------------------------------------------------
cat(sprintf("[fg] +/- %d kb of top-%d SNPs per trait...\n",
            WIN_BP / 1000, TOP_N))
anc <- read.csv(ANCHOR)
anc <- anc[!is.na(anc$chr_ayb) & !is.na(anc$snp_pos_ayb), ]
anc$snp_pos_ayb <- as.integer(anc$snp_pos_ayb)

fg <- character()
for (trait in c(TRAITS_BIO, TRAITS_NEW)) {
  gp <- if (trait %in% TRAITS_BIO)
    file.path(GWAS_BIO, paste0("gwas_", trait, "_M1_K.csv"))
  else
    file.path(GWAS_NEW, paste0("gwas_", trait, "_M1_K.csv"))
  if (!file.exists(gp)) next
  g <- read.csv(gp)
  g <- merge(g, anc[, c("rs", "chr_ayb", "snp_pos_ayb")],
             by = "rs", all.x = TRUE)
  g <- g[!is.na(g$chr_ayb), ]
  g <- g[order(g$p), ][seq_len(min(TOP_N, nrow(g))), ]
  for (i in seq_len(nrow(g))) {
    c_ <- g$chr_ayb[i]; p_ <- g$snp_pos_ayb[i]
    hit <- gene_tbl[gene_tbl$chr == c_ &
                    gene_tbl$end >= p_ - WIN_BP &
                    gene_tbl$start <= p_ + WIN_BP, "gene_id"]
    fg <- c(fg, hit)
  }
}
fg <- unique(fg)
cat(sprintf("[fg] foreground gene set: %d unique genes\n", length(fg)))


# ---------------------------------------------------------------------------
# topGO run per namespace
# ---------------------------------------------------------------------------
# the topGO 'allGenes' vector is a factor labelling each gene as foreground
# (1) or background (0). We use only genes WITH at least one GO term as
# the testable universe (consistent with topGO best practice).
universe <- names(gene2go)
fg_in_univ <- intersect(fg, universe)
cat(sprintf("[topgo] universe = %d annotated genes; foreground in universe = %d\n",
            length(universe), length(fg_in_univ)))

all_genes <- factor(as.integer(universe %in% fg_in_univ))
names(all_genes) <- universe

run_topgo <- function(ontology, label) {
  cat(sprintf("\n=== topGO %s ===\n", label))
  go_data <- new("topGOdata",
                 description = paste0("AYB candidate-region ", label),
                 ontology    = ontology,
                 allGenes    = all_genes,
                 annot       = annFUN.gene2GO,
                 gene2GO     = gene2go,
                 nodeSize    = 5)
  test_classic <- new("classicCount", testStatistic = GOFisherTest,
                      name = "Fisher classic")
  test_elim    <- new("elimCount",    testStatistic = GOFisherTest,
                      name = "Fisher elim")
  res_classic <- getSigGroups(go_data, test_classic)
  res_elim    <- getSigGroups(go_data, test_elim)
  tab <- GenTable(go_data,
                  classic = res_classic,
                  elim    = res_elim,
                  orderBy = "elim",
                  topNodes = 100)
  tab$ontology <- label
  for (col in c("classic", "elim")) {
    tab[[col]] <- suppressWarnings(as.numeric(tab[[col]]))
  }
  # BH-FDR over the topGO-tested terms
  tab$q_classic_bh <- p.adjust(tab$classic, method = "BH")
  tab$q_elim_bh    <- p.adjust(tab$elim,    method = "BH")
  out_csv <- file.path(TAB, paste0("go_topgo_", label, ".csv"))
  write.csv(tab, out_csv, row.names = FALSE)
  cat(sprintf("  %d terms reported, %d with elim p < 0.01, %d with elim q_BH < 0.10\n",
              nrow(tab), sum(tab$elim < 0.01, na.rm = TRUE),
              sum(tab$q_elim_bh < 0.10, na.rm = TRUE)))
  return(tab)
}

bp <- run_topgo("BP", "BP")
mf <- run_topgo("MF", "MF")
cc <- run_topgo("CC", "CC")


# ---------------------------------------------------------------------------
# Figures: per-namespace horizontal bar chart of top terms by elim p
# ---------------------------------------------------------------------------
NS_PAL <- c(BP = "#0072B2", MF = "#D55E00", CC = "#009E73")

# Q1 convention: show GO term NAMES, -log10 p, and a small "k / N" count.
# No q-values on individual bars; the FDR threshold is in the figure title /
# caption only. Combined three-panel figure across BP / MF / CC.
render_combined <- function(bp_, mf_, cc_, fname) {
  topn <- 10
  top_each <- function(tab, ns_label, col_) {
    s <- head(tab[order(tab$elim), ], topn)
    s$ns <- ns_label
    s$col <- col_
    s
  }
  combined <- rbind(
    top_each(bp_, "BP", NS_PAL["BP"]),
    top_each(mf_, "MF", NS_PAL["MF"]),
    top_each(cc_, "CC", NS_PAL["CC"])
  )
  for (dev_open in list(
        function() png(file.path(FIG, paste0(fname, ".png")),
                       width = 2800, height = 2400, res = 220),
        function() pdf(file.path(FIG, paste0(fname, ".pdf")),
                       width = 13, height = 11))) {
    dev_open()
    layout(matrix(c(1, 2, 3), nrow = 3), heights = c(1, 1, 1))
    par(mar = c(4, 26, 2.5, 1))
    for (ns in c("BP", "MF", "CC")) {
      s <- combined[combined$ns == ns, ]
      s <- s[order(s$elim), ]
      nlp <- -log10(s$elim)
      # bar plot, term-name labels on the left
      bp_obj <- barplot(rev(nlp), horiz = TRUE,
                        names.arg = rev(s$Term),
                        las = 1, cex.names = 0.8,
                        col = unique(s$col), border = "black",
                        xlab = expression(-log[10] * " elim Fisher p"),
                        xlim = c(0, max(nlp) * 1.25),
                        main = sprintf("topGO %s — top %d terms",
                                       ns, topn),
                        cex.main = 1.0, font.main = 1)
      # tiny fg / bg count just past each bar
      for (i in seq_along(nlp)) {
        ii <- length(nlp) - i + 1
        text(nlp[ii], bp_obj[i],
             sprintf("(%d / %d)", s$Significant[ii], s$Annotated[ii]),
             pos = 4, cex = 0.7, col = "grey30", xpd = NA)
      }
      abline(v = -log10(0.05), col = "red", lty = "dashed", lwd = 0.6)
    }
    dev.off()
  }
  cat(sprintf("  rendered %s (PNG + PDF)\n", fname))
}

render_combined(bp, mf, cc, "fig75_go_topgo_combined")

cat("\n[done] topGO enrichment complete.\n")
