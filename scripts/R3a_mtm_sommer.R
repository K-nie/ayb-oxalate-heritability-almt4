#!/usr/bin/env Rscript
# REVISION Phase 3a — PROPER multivariate (multi-trait) GBLUP in sommer.
#
# Reviewer (L327-361): the script-28 "multi-output ridge with a shared
# alpha" is trivially equivalent to independent single-trait fits (hence
# its delta-r == 0 is meaningless). A proper multivariate GBLUP estimates
# an unstructured genetic covariance G0 AND a residual covariance R0 across
# traits, and borrows strength through the genetic covariance. We fit it
# with sommer (Covarrubias-Pazaran 2016) and report:
#   - the REML genetic covariance / correlation and residual covariance,
#   - an HONEST multi-trait vs single-trait predictive-ability delta-r with
#     a trait-assisted (CV2) cross-validation, where MT CAN help: the
#     target trait is masked in the test fold while the OTHER traits stay
#     observed on the same test lines.
#
# Trait sets:
#   oxalate_triplet   {Total, Soluble, Insoluble}  -- EXPECTED to be
#       rank-deficient because Total = Soluble + Insoluble exactly; we
#       attempt it and record the degeneracy.
#   oxalate_pair      {Soluble, Insoluble}         -- the two independent
#       fractions (genetic r ~ -0.8).
#   protein_oxalate   {Crude_Protein, Soluble, Insoluble}
#
# Output: results/R3_breeding_tools/tables/
#   mtm_genetic_residual_cov_<set>.csv
#   mtm_vs_st_cv_<set>.csv
# Author: Benjamin Narh-Madey
suppressMessages(library(sommer))
set.seed(42)

TAB <- "/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data/results/R3_breeding_tools/tables"
N_FOLDS <- 5; N_REPS <- 20; TOL <- 1e-2

load_set <- function(tag) {
  d <- read.csv(file.path(TAB, paste0("mtm_data_", tag, ".csv")), check.names = FALSE)
  K <- as.matrix(read.csv(file.path(TAB, paste0("mtm_K_", tag, ".csv")),
                          row.names = 1, check.names = FALSE))
  colnames(K) <- rownames(K)
  d <- d[d$sample %in% rownames(K), ]
  d$sample <- factor(d$sample, levels = rownames(K))
  list(d = d, K = K)
}

# univariate GBLUP (sommer) predictive ability on masked test lines
st_predict <- function(d, K, trait, test) {
  d2 <- d; d2[[trait]][test] <- NA
  f <- as.formula(paste(trait, "~ 1"))
  m <- tryCatch(mmer(f, random = ~ vsr(sample, Gu = K), rcov = ~ units,
                     data = d2, verbose = FALSE, date.warning = FALSE, tolParInv = TOL),
                error = function(e) NULL)
  if (is.null(m)) return(rep(NA, length(test)))
  u <- m$U[[1]][[trait]]; b <- m$Beta$Estimate[1]
  as.numeric(b + u[as.character(d$sample[test])])
}

# multivariate GBLUP predictive ability on masked target (others observed)
mt_predict <- function(d, K, traits, target, test) {
  d2 <- d; d2[[target]][test] <- NA
  Tn <- length(traits)
  f <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
  m <- tryCatch(mmer(f, random = ~ vsr(sample, Gu = K, Gtc = unsm(Tn)),
                     rcov = ~ vsr(units, Gtc = diag(Tn)),
                     data = d2, verbose = FALSE, date.warning = FALSE, tolParInv = TOL),
                error = function(e) NULL)
  if (is.null(m)) return(list(pred = rep(NA, length(test)), G0 = NULL, R0 = NULL))
  u <- m$U[[1]][[target]]
  b <- m$Beta$Estimate[m$Beta$Trait == target][1]
  list(pred = as.numeric(b + u[as.character(d$sample[test])]),
       G0 = m$sigma[[1]], R0 = m$sigma[[2]])
}

run_set <- function(tag, traits) {
  cat(sprintf("\n=== %s : %s ===\n", tag, paste(traits, collapse = ", ")))
  S <- load_set(tag); d <- S$d; K <- S$K; n <- nrow(d)
  Tn <- length(traits)

  # Full-data MTM once to report G0 / R0 (unstructured genetic + diag resid)
  f <- as.formula(paste0("cbind(", paste(traits, collapse = ","), ") ~ 1"))
  mfull <- tryCatch(mmer(f, random = ~ vsr(sample, Gu = K, Gtc = unsm(Tn)),
                         rcov = ~ vsr(units, Gtc = diag(Tn)),
                         data = d, verbose = FALSE, date.warning = FALSE, tolParInv = TOL),
                    error = function(e) NULL)
  if (!is.null(mfull)) {
    G0 <- mfull$sigma[[1]]; R0 <- mfull$sigma[[2]]
    dimnames(G0) <- list(traits, traits); dimnames(R0) <- list(traits, traits)
    gcor <- cov2cor(G0)
    write.csv(cbind(stat = "genetic_cov", as.data.frame(G0)),
              file.path(TAB, paste0("mtm_genetic_cov_", tag, ".csv")), row.names = TRUE)
    write.csv(cbind(stat = "genetic_cor", as.data.frame(round(gcor, 3))),
              file.path(TAB, paste0("mtm_genetic_cor_", tag, ".csv")), row.names = TRUE)
    write.csv(cbind(stat = "residual_cov", as.data.frame(R0)),
              file.path(TAB, paste0("mtm_residual_cov_", tag, ".csv")), row.names = TRUE)
    cat("  genetic correlation matrix:\n"); print(round(gcor, 3))
    cond_G <- kappa(G0)
    cat(sprintf("  G0 condition number = %.3g (rank-deficient if huge)\n", cond_G))
  } else {
    cat("  FULL MTM failed (likely rank-deficient genetic covariance).\n")
    cond_G <- NA
  }

  # Trait-assisted CV2: per target trait, mask target in test fold, keep others.
  rows <- list()
  for (target in traits) {
    st_r <- numeric(0); mt_r <- numeric(0)
    for (rep in 1:N_REPS) {
      folds <- sample(rep(1:N_FOLDS, length.out = n))
      st_pred <- rep(NA, n); mt_pred <- rep(NA, n)
      for (k in 1:N_FOLDS) {
        test <- which(folds == k)
        st_pred[test] <- st_predict(d, K, target, test)
        mt_pred[test] <- mt_predict(d, K, traits, target, test)$pred
      }
      obs <- d[[target]]
      ok <- is.finite(st_pred) & is.finite(mt_pred) & is.finite(obs)
      if (sum(ok) > 5 && sd(st_pred[ok]) > 0 && sd(mt_pred[ok]) > 0) {
        st_r <- c(st_r, cor(st_pred[ok], obs[ok]))
        mt_r <- c(mt_r, cor(mt_pred[ok], obs[ok]))
      }
    }
    dr <- mt_r - st_r
    tt <- if (length(dr) > 2 && sd(dr) > 0) t.test(dr) else NULL
    ci_lo <- if (!is.null(tt)) tt$conf.int[1] else NA
    ci_hi <- if (!is.null(tt)) tt$conf.int[2] else NA
    rows[[target]] <- data.frame(
      set = tag, target = target, n = n, n_reps = length(dr),
      ST_mean_r = mean(st_r), MT_mean_r = mean(mt_r),
      delta_r = mean(dr), delta_r_ci_lo = ci_lo, delta_r_ci_hi = ci_hi,
      p_paired = if (!is.null(tt)) tt$p.value else NA,
      G0_condition = cond_G)
    cat(sprintf("  %-18s ST r=%+.3f  MT r=%+.3f  dr=%+.3f [%+.3f,%+.3f] p=%.3f\n",
                target, mean(st_r), mean(mt_r), mean(dr), ci_lo, ci_hi,
                if (!is.null(tt)) tt$p.value else NA))
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(TAB, paste0("mtm_vs_st_cv_", tag, ".csv")), row.names = FALSE)
  out
}

all <- list()
all[["triplet"]] <- tryCatch(run_set("oxalate_triplet",
                    c("Total_Oxalate", "Soluble_Oxalate", "Insoluble_Oxalate")),
                    error = function(e) { cat("triplet ERR:", e$message, "\n"); NULL })
all[["pair"]] <- run_set("oxalate_pair", c("Soluble_Oxalate", "Insoluble_Oxalate"))
all[["protox"]] <- run_set("protein_oxalate",
                    c("Crude_Protein", "Soluble_Oxalate", "Insoluble_Oxalate"))

cat(sprintf("\n[done] sommer %s ; MTM-vs-ST tables in %s\n",
            as.character(packageVersion("sommer")), TAB))
