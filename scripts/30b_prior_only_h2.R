#!/usr/bin/env Rscript
# 30b_prior_only_h2 -- prior-only (no data) draw of marker-based heritability,
# to quantify how much of the RKHS posterior for soluble oxalate is prior
# versus data (Reviewer 1, R1.10: "RKHS prior-driven").
#
# In script 30 the RKHS variance components carry inverse-chi^2 priors with
# df0 = 5 and genetic/residual scales each set to half the phenotypic variance.
# Because the two scales are equal and the df0 are equal, the common scale
# cancels in h^2 = sigma2_g / (sigma2_g + sigma2_e): the prior-implied h^2 is
# exactly Beta(df0/2, df0/2) = Beta(2.5, 2.5), symmetric about 0.5 and
# independent of the data. We confirm this by direct Monte-Carlo draw from the
# two inverse-chi^2 priors and report the prior-only median and 95 % interval
# (equal-tailed and HPD, which coincide here because the prior is symmetric).
#
# Output: results/30_bayesian_h2_soluble_oxalate/tables/h2_prior_only_summary.csv
#
# Author: Benjamin Narh-Madey

suppressMessages(library(coda))

OUT <- file.path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data",
                 "results", "30_bayesian_h2_soluble_oxalate", "tables")
DF <- 5          # df0 in script 30
N  <- 200000
set.seed(20250519)

# scaled inverse-chi^2 draw: sigma2 = S * df0 / chi2_{df0}; S cancels in h^2,
# so set S = 1 for both components.
sg <- DF / rchisq(N, DF)
se <- DF / rchisq(N, DF)
h2 <- sg / (sg + se)

q  <- quantile(h2, c(0.025, 0.5, 0.975))
hpd <- HPDinterval(as.mcmc(h2), prob = 0.95)

summ <- data.frame(
  parameter        = "h2_prior_only_Soluble_Oxalate",
  prior            = "Beta(2.5,2.5) = ratio of two inv-chi2(df0=5, equal scale)",
  n_draws          = N,
  prior_median     = as.numeric(q[2]),
  eqtail_low       = as.numeric(q[1]),
  eqtail_high      = as.numeric(q[3]),
  hpd95_low        = as.numeric(hpd[1, "lower"]),
  hpd95_high       = as.numeric(hpd[1, "upper"])
)
write.csv(summ, file.path(OUT, "h2_prior_only_summary.csv"), row.names = FALSE)
cat(sprintf("prior-only h^2: median %.3f  95%% equal-tailed [%.3f, %.3f]  HPD [%.3f, %.3f]\n",
            q[2], q[1], q[3], hpd[1, "lower"], hpd[1, "upper"]))
cat(sprintf("wrote %s\n", file.path(OUT, "h2_prior_only_summary.csv")))
