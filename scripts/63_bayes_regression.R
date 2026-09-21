#!/usr/bin/env Rscript
# BayesB / BayesC posterior inclusion probabilities per SNP per trait.
#
# Independent corroboration of the candidate-gene callout from a
# fundamentally different shrinkage prior. RR-BLUP and FarmCPU / mrMLM
# use frequentist shrinkage or variable-selection frameworks; BayesB and
# BayesC pi use a spike-and-slab prior that assigns a per-SNP probability
# of having a non-zero effect.
#
# A high posterior inclusion probability (PIP) at the ALMT4 SNP under
# both BayesB and BayesC is independent evidence for the candidate gene
# call: the data prefer the SNP to have a non-zero effect under a method
# that is willing to set most other SNPs' effects to zero.
#
# Inputs
# ------
#   - HapMap CSV (script 01 QC re-applied here)
#   - BLUP table per trait
#
# Outputs (results/63_bayes_regression/)
# --------------------------------------
#   tables/
#       bayesB_pip_per_snp_per_trait.csv    PIP per SNP per trait
#       bayesC_posterior_effects.csv         posterior mean effect per SNP per trait
#       almt4_focal_summary.csv               PIP at the ALMT4 SNP per trait
#
# Author: Benjamin Narh-Madey

suppressMessages({
  library(BGLR)
  library(data.table)
})

ROOT <- "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data"
HAPMAP <- file.path(ROOT, "AYB_SNP_Result_Report-DAf18-2580",
                    "Report_DAf18-2580_SNP_HapMap.csv")
BLUP_CSV <- file.path(ROOT, "results", "00_blup_pipeline", "tables",
                       "phenotype_blups.csv")
OUT <- file.path(ROOT, "results", "63_bayes_regression")
TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
# BGLR scratch directory keeps the chain trace files out of the working
# directory; we point it at a results subfolder so they survive.
SCRATCH <- file.path(OUT, "bglr_scratch")
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)
setwd(SCRATCH)

SAMP_CR <- 0.90
MARK_CR <- 0.90
MAF_MIN <- 0.05

NITER <- 30000
BURNIN <- 5000
THIN <- 5

TRAITS <- c("Tannin", "Phenol", "Flavonoid", "Antioxidant",
            "Seed_Length", "Seed_Width", "Seed_Thickness", "Mass_of_Seeds",
            "Seed_Coat_Tannin", "Crude_Protein",
            "Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate")

# Focal SNP at Ss10:15,394,673 (ALMT4 candidate)
FOCAL_RS <- "100033542|F|0-31:T>C-31:T>C"


# ---------------------------------------------------------------------------
# HapMap -> dosage matrix with the standard QC cascade
# ---------------------------------------------------------------------------
build_dosage <- function() {
  cat("[load] reading HapMap\n")
  hm <- fread(HAPMAP, sep = ",", header = TRUE)
  meta_cols <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                  "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
  sample_cols <- setdiff(colnames(hm), meta_cols)
  cat(sprintf("  raw: %d markers x %d samples\n", nrow(hm), length(sample_cols)))

  # Build dosage matrix: rows = markers, cols = samples; 0 = hom ref,
  # 1 = het, 2 = hom alt, NA = missing.
  alleles_split <- strsplit(hm$alleles, "/")
  ref <- vapply(alleles_split, `[`, character(1), 1)
  alt <- vapply(alleles_split, `[`, character(1), 2)
  hom_ref <- paste0(ref, ref)
  hom_alt <- paste0(alt, alt)
  het1 <- paste0(ref, alt)
  het2 <- paste0(alt, ref)

  call_mat <- as.matrix(hm[, .SD, .SDcols = sample_cols])
  dose <- matrix(NA_real_, nrow = nrow(hm), ncol = length(sample_cols),
                  dimnames = list(hm[["rs#"]], sample_cols))
  for (i in seq_len(nrow(hm))) {
    row_calls <- as.character(call_mat[i, ])
    is_ref <- row_calls == hom_ref[i]
    is_alt <- row_calls == hom_alt[i]
    is_het <- row_calls == het1[i] | row_calls == het2[i]
    dose[i, is_ref] <- 0
    dose[i, is_het] <- 1
    dose[i, is_alt] <- 2
  }

  # Sample QC.
  samp_cr <- colMeans(!is.na(dose))
  keep_samples <- names(samp_cr)[samp_cr >= SAMP_CR]
  dose <- dose[, keep_samples]
  # Marker QC.
  mark_cr <- rowMeans(!is.na(dose))
  mean_dose <- rowMeans(dose, na.rm = TRUE)
  maf <- pmin(mean_dose / 2, 1 - mean_dose / 2)
  keep_markers <- which(mark_cr >= MARK_CR & maf >= MAF_MIN)
  dose <- dose[keep_markers, ]

  cat(sprintf("  post-QC: %d markers x %d samples\n", nrow(dose), ncol(dose)))

  # Mean impute residual missingness per marker.
  mean_dose <- rowMeans(dose, na.rm = TRUE)
  for (i in seq_len(nrow(dose))) {
    na_idx <- is.na(dose[i, ])
    if (any(na_idx)) {
      dose[i, na_idx] <- mean_dose[i]
    }
  }

  list(dose = dose, samples = colnames(dose), markers = rownames(dose))
}


# ---------------------------------------------------------------------------
# Per-trait BayesB and BayesC fit
# ---------------------------------------------------------------------------
fit_bayes_one_trait <- function(trait, dose_mat, blups, kind = "BayesB") {
  # Subset to samples with both BLUP and dosage.
  y_full <- blups[[trait]]
  names(y_full) <- blups$sample
  common <- intersect(colnames(dose_mat), names(y_full))
  common <- common[!is.na(y_full[common])]
  y <- y_full[common]
  X <- t(dose_mat[, common])  # samples x markers
  X_centred <- scale(X, center = TRUE, scale = FALSE)

  cat(sprintf("  [%s] %s: n = %d, m = %d\n", kind, trait,
               length(y), ncol(X)))
  ETA <- list(list(X = X_centred, model = kind, probIn = 0.05))
  fit <- BGLR(y = y, ETA = ETA,
               nIter = NITER, burnIn = BURNIN, thin = THIN,
               saveAt = paste0(kind, "_", trait, "_"),
               verbose = FALSE)
  # Posterior inclusion probability and posterior mean effect per SNP.
  bHat <- fit$ETA[[1]]$b
  # `d` is the posterior inclusion indicator (frequency of inclusion
  # across MCMC samples) under the BayesB / BayesC parameterisation.
  dHat <- fit$ETA[[1]]$d
  list(samples = common, b = bHat, pip = dHat)
}


main <- function() {
  geno <- build_dosage()
  blups <- fread(BLUP_CSV)
  cat(sprintf("[load] BLUP table: %d accessions x %d traits\n",
               nrow(blups), ncol(blups) - 1))

  pip_rows <- list()
  effect_rows <- list()
  focal_rows <- list()

  for (trait in TRAITS) {
    cat(sprintf("[fit] %s under BayesB and BayesC\n", trait))

    # BayesB
    resB <- fit_bayes_one_trait(trait, geno$dose, blups, kind = "BayesB")
    # BayesC (probIn estimated from the data)
    resC <- fit_bayes_one_trait(trait, geno$dose, blups, kind = "BayesC")

    pip_df <- data.table(
      trait = trait,
      rs = geno$markers,
      pip_bayesB = resB$pip,
      pip_bayesC = resC$pip
    )
    pip_rows[[trait]] <- pip_df

    effect_df <- data.table(
      trait = trait,
      rs = geno$markers,
      effect_bayesB = resB$b,
      effect_bayesC = resC$b
    )
    effect_rows[[trait]] <- effect_df

    focal_row <- data.table(
      trait = trait,
      n = length(resB$samples),
      focal_rs = FOCAL_RS,
      pip_bayesB_focal = pip_df[rs == FOCAL_RS, pip_bayesB],
      pip_bayesC_focal = pip_df[rs == FOCAL_RS, pip_bayesC],
      effect_bayesB_focal = effect_df[rs == FOCAL_RS, effect_bayesB],
      effect_bayesC_focal = effect_df[rs == FOCAL_RS, effect_bayesC]
    )
    focal_rows[[trait]] <- focal_row

    cat(sprintf(
      "  [%s] focal SNP: PIP_B = %.3f, PIP_C = %.3f, beta_B = %+0.3f, beta_C = %+0.3f\n",
      trait,
      focal_row$pip_bayesB_focal, focal_row$pip_bayesC_focal,
      focal_row$effect_bayesB_focal, focal_row$effect_bayesC_focal))
  }

  fwrite(rbindlist(pip_rows),
          file.path(TAB, "bayesB_bayesC_pip_per_snp_per_trait.csv"))
  fwrite(rbindlist(effect_rows),
          file.path(TAB, "bayesB_bayesC_effects_per_snp_per_trait.csv"))
  fwrite(rbindlist(focal_rows),
          file.path(TAB, "almt4_focal_summary.csv"))

  cat("\n[result] tables written to results/63_bayes_regression/tables/\n")
}


main()
