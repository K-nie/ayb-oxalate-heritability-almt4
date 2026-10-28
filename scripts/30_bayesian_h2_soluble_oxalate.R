# Bayesian marker-based heritability for Soluble_Oxalate.
#
# The REML profile-likelihood estimate for Soluble_Oxalate at n = 41 lands
# at the upper boundary of the search grid (MLE h^2 = 0.999), and the
# 95 % CI [0.14, 0.999] is therefore reported as a one-sided lower bound
# (the upper edge sits at the grid limit; the likelihood is flat past
# 0.14). The audit recommendation for Q1 submission is to pair the REML
# one-sided bound with a Bayesian estimate under a weakly informative
# prior on h^2, so the posterior is regularised off the boundary and a
# two-sided 95 % HPD credible interval is reportable.
#
# Method: BGLR (de los Campos & Perez 2010) RKHS / GBLUP regression with
# a precomputed kinship matrix K = VanRaden method 1. Priors:
#   beta_intercept ~ N(0, 1e6)   (flat)
#   sigma2_g       ~ Inv-chi-sq (Beta(2, 2)-style weakly informative on h^2)
#   sigma2_e       ~ Inv-chi-sq (Beta(2, 2)-style weakly informative on h^2)
# A 20,000-iteration Gibbs sampler with 5,000 burn-in yields posterior
# samples of (sigma2_g, sigma2_e); h^2 = sigma2_g / (sigma2_g + sigma2_e).
#
# Outputs (results/30_bayesian_h2_soluble_oxalate/)
# -------------------------------------------------
# tables/
#     h2_posterior_summary.csv        n_eff, mean, median, 95 % HPD, SD
#     h2_posterior_samples.csv        the trimmed posterior chain (for plotting)
#     reml_vs_bayes_comparison.csv    REML one-sided vs Bayesian two-sided
# figures/
#     fig68_bayes_h2_posterior.png/.pdf
#
# Author: Benjamin Narh-Madey

suppressMessages({
    library(BGLR)
    library(coda)
})

setwd("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUT <- "results/30_bayesian_h2_soluble_oxalate"
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

cat("[load] phenotype BLUPs + dosage...\n")
blups <- read.csv("results/00_blup_pipeline/tables/phenotype_blups.csv",
                  row.names = 1, check.names = FALSE)
hapmap <- read.csv("AYB_SNP_Result_Report-DAf18-2580/Report_DAf18-2580_SNP_HapMap.csv",
                   check.names = FALSE, stringsAsFactors = FALSE)
meta_cols <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
sample_cols <- setdiff(colnames(hapmap), meta_cols)

# Build dosage (0/1/2) on the alternate allele
allele_split <- do.call(rbind, strsplit(hapmap$alleles, "/"))
colnames(allele_split) <- c("ref", "alt")
dosage_mat <- matrix(NA_real_, nrow = nrow(hapmap), ncol = length(sample_cols))
for (i in seq_len(nrow(hapmap))) {
    ref <- allele_split[i, "ref"]; alt <- allele_split[i, "alt"]
    calls <- as.character(unlist(hapmap[i, sample_cols]))
    dosage_mat[i, calls == paste0(ref, ref)] <- 0
    dosage_mat[i, calls %in% c(paste0(ref, alt), paste0(alt, ref))] <- 1
    dosage_mat[i, calls == paste0(alt, alt)] <- 2
}
rownames(dosage_mat) <- hapmap$`rs#`
colnames(dosage_mat) <- sample_cols

# QC filters: marker CR >= 0.90, MAF >= 0.05, sample CR >= 0.90
mark_cr <- rowMeans(!is.na(dosage_mat))
maf <- rowMeans(dosage_mat / 2, na.rm = TRUE)
maf <- pmin(maf, 1 - maf)
samp_cr <- colMeans(!is.na(dosage_mat))
keep_m <- mark_cr >= 0.90 & maf >= 0.05
keep_s <- samp_cr >= 0.90
M <- t(dosage_mat[keep_m, keep_s])
# Mean-impute remaining missing dosages
col_mean <- colMeans(M, na.rm = TRUE)
for (j in seq_len(ncol(M))) {
    idx <- is.na(M[, j])
    if (any(idx)) M[idx, j] <- col_mean[j]
}
cat(sprintf("[load] post-QC: %d samples x %d markers\n", nrow(M), ncol(M)))

# Pull Soluble_Oxalate phenotypes for n = 41 subset
y_full <- blups$Soluble_Oxalate
names(y_full) <- rownames(blups)
shared <- intersect(rownames(M), names(y_full)[!is.na(y_full)])
y <- y_full[shared]
M_sub <- M[shared, ]
# Re-MAF on the subset, then bump to MAF >= 0.10 (matches script 15 convention)
p <- colMeans(M_sub) / 2
maf_sub <- pmin(p, 1 - p)
M_sub <- M_sub[, maf_sub >= 0.10]
cat(sprintf("[Soluble_Oxalate] n = %d, m (MAF >= 0.10) = %d\n",
             length(y), ncol(M_sub)))

# Build K (VanRaden method 1)
p <- colMeans(M_sub) / 2
p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
Z <- scale(M_sub, center = 2 * p, scale = FALSE)
K <- (Z %*% t(Z)) / sum(2 * p * (1 - p))

# BGLR setup. RKHS prior gives a Bayesian GBLUP equivalent on the kinship.
# The default scale / df settings yield a weakly informative prior on h^2
# (effectively Beta(2, 2)-like, mass shifted off the boundaries).
NITER <- 20000
BURN <- 5000
THIN <- 5
SCALE_E <- var(y) * 0.5
SCALE_G <- var(y) * 0.5
DF <- 5

ETA <- list(K = list(K = K, model = "RKHS",
                       df0 = DF, S0 = SCALE_G))
set.seed(20250519)

cat(sprintf("[bglr] running RKHS Gibbs sampler: %d iter, %d burn, thin %d...\n",
             NITER, BURN, THIN))
fm <- BGLR(y = y, ETA = ETA, nIter = NITER, burnIn = BURN, thin = THIN,
            df0 = DF, S0 = SCALE_E,
            saveAt = file.path(TAB, "_bglr_"),
            verbose = FALSE)

# Read the saved variance-component chains (BGLR writes them to disk)
var_g <- scan(file.path(TAB, "_bglr_ETA_K_varU.dat"), quiet = TRUE)
var_e <- scan(file.path(TAB, "_bglr_varE.dat"), quiet = TRUE)
# Drop burn-in safety (BGLR already removed nominal burn but a margin helps)
keep <- (length(var_g) %/% 4):length(var_g)
var_g <- var_g[keep]; var_e <- var_e[keep]
h2_post <- var_g / (var_g + var_e)

# Posterior summary + HPD CI via coda
mcmc_obj <- as.mcmc(h2_post)
hpd <- HPDinterval(mcmc_obj, prob = 0.95)
n_eff <- effectiveSize(mcmc_obj)
post_summary <- data.frame(
    parameter = "h2_marker_Soluble_Oxalate",
    n_post_samples = length(h2_post),
    n_eff = as.numeric(n_eff),
    posterior_mean = mean(h2_post),
    posterior_median = median(h2_post),
    posterior_sd = sd(h2_post),
    hpd95_low = as.numeric(hpd[1, "lower"]),
    hpd95_high = as.numeric(hpd[1, "upper"]),
    prob_h2_gt_0_5 = mean(h2_post > 0.5),
    prob_h2_gt_0_1 = mean(h2_post > 0.1)
)
write.csv(post_summary, file.path(TAB, "h2_posterior_summary.csv"),
            row.names = FALSE)
write.csv(data.frame(h2 = h2_post), file.path(TAB, "h2_posterior_samples.csv"),
            row.names = FALSE)

cat(sprintf("\n[posterior] Soluble_Oxalate h^2 (Bayesian, BGLR RKHS):\n"))
cat(sprintf("  n_eff = %.1f\n", n_eff))
cat(sprintf("  mean = %.3f, median = %.3f, sd = %.3f\n",
             mean(h2_post), median(h2_post), sd(h2_post)))
cat(sprintf("  95 %% HPD = [%.3f, %.3f]\n", hpd[1, "lower"], hpd[1, "upper"]))
cat(sprintf("  P(h^2 > 0.10) = %.3f, P(h^2 > 0.50) = %.3f\n",
             mean(h2_post > 0.10), mean(h2_post > 0.50)))

# Compare REML one-sided to Bayesian two-sided
reml_table <- read.csv("results/24_bootstrap_h2/tables/h2_profile_summary.csv")
reml_row <- reml_table[reml_table$trait == "Soluble_Oxalate", ]
comparison <- data.frame(
    method = c("REML profile-likelihood (one-sided)",
                 "Bayesian RKHS (two-sided HPD)"),
    point_estimate = c(reml_row$h2_mle, median(h2_post)),
    ci_low = c(reml_row$h2_ci95_low, hpd[1, "lower"]),
    ci_high = c(reml_row$h2_ci95_high, hpd[1, "upper"]),
    ci_kind = c(reml_row$ci_kind, "two_sided"),
    notes = c("Upper edge sits at grid limit; report as h^2 >= 0.14",
                "Regularised off boundary by weakly informative prior")
)
write.csv(comparison, file.path(TAB, "reml_vs_bayes_comparison.csv"),
            row.names = FALSE)

cat("\n[compare] REML vs Bayesian:\n")
print(comparison)

# Posterior plot
WONG_BLUE <- "#0072B2"
WONG_VERMILLION <- "#D55E00"
pdf(file.path(FIG, "fig68_bayes_h2_posterior.pdf"),
     width = 7, height = 4.5)
op <- par(mar = c(4, 4, 3, 1))
hist(h2_post, breaks = 50, freq = FALSE, col = WONG_BLUE,
      border = "white", xlim = c(0, 1),
      xlab = expression(h^2),
      main = sprintf("Posterior on marker-based h^2 (Soluble_Oxalate, n = %d)",
                      length(y)))
abline(v = median(h2_post), col = WONG_VERMILLION, lwd = 2)
abline(v = hpd[1, ], col = WONG_VERMILLION, lwd = 1.4, lty = 2)
abline(v = reml_row$h2_ci95_low, col = "black", lwd = 1.4, lty = 3)
legend("topright",
        legend = c(sprintf("posterior median = %.2f", median(h2_post)),
                    sprintf("95 %% HPD = [%.2f, %.2f]",
                             hpd[1, "lower"], hpd[1, "upper"]),
                    sprintf("REML one-sided lower bound = %.2f",
                             reml_row$h2_ci95_low)),
        col = c(WONG_VERMILLION, WONG_VERMILLION, "black"),
        lty = c(1, 2, 3), lwd = c(2, 1.4, 1.4), bty = "n",
        cex = 0.85)
par(op)
dev.off()
# Same plot, PNG copy
png(file.path(FIG, "fig68_bayes_h2_posterior.png"),
     width = 7 * 300, height = 4.5 * 300, res = 300)
op <- par(mar = c(4, 4, 3, 1))
hist(h2_post, breaks = 50, freq = FALSE, col = WONG_BLUE,
      border = "white", xlim = c(0, 1),
      xlab = expression(h^2),
      main = sprintf("Posterior on marker-based h^2 (Soluble_Oxalate, n = %d)",
                      length(y)))
abline(v = median(h2_post), col = WONG_VERMILLION, lwd = 2)
abline(v = hpd[1, ], col = WONG_VERMILLION, lwd = 1.4, lty = 2)
abline(v = reml_row$h2_ci95_low, col = "black", lwd = 1.4, lty = 3)
legend("topright",
        legend = c(sprintf("posterior median = %.2f", median(h2_post)),
                    sprintf("95 %% HPD = [%.2f, %.2f]",
                             hpd[1, "lower"], hpd[1, "upper"]),
                    sprintf("REML one-sided lower bound = %.2f",
                             reml_row$h2_ci95_low)),
        col = c(WONG_VERMILLION, WONG_VERMILLION, "black"),
        lty = c(1, 2, 3), lwd = c(2, 1.4, 1.4), bty = "n",
        cex = 0.85)
par(op)
dev.off()
cat(sprintf("\n[done] outputs in %s\n", OUT))

# Clean up the BGLR scratch files so they do not clutter the share copy
bglr_scratch <- list.files(TAB, pattern = "^_bglr_", full.names = TRUE)
file.remove(bglr_scratch)
