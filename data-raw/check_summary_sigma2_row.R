suppressPackageStartupMessages({
  pkgload::load_all("../glmbayesCore", quiet = TRUE)
  pkgload::load_all(".", quiet = TRUE)
})
source("tests/manual/_small5_lmerb_fixture.R")

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
ps <- Prior_Setup_lmebayes(
  fx$form, data = fx$dat, pwt = 0.01, pwt_measurement = 0.49
)
pf <- pfamily_list(ps)
m_disp <- ps$ing_prior_measurement
disp_pf <- dGamma(
  shape = m_disp$shape, rate = m_disp$rate,
  beta = matrix(0, 1, 1), Inv_Dispersion = TRUE,
  disp_lower = m_disp$disp_lower, disp_upper = m_disp$disp_upper
)
fit <- lmerb(
  fx$form, data = fx$dat, pfamily_list = pf,
  dispersion_ranef = disp_pf, n = 50L, progbar = FALSE, verbose = FALSE
)
s <- summary(fit)
cat("--- tau2_prior ---\n")
print(round(s$tau2_prior_overview, 3))
cat("\n--- tau2_overview ---\n")
print(round(s$tau2_overview, 3))
cat("\n--- tau2_pct ---\n")
print(round(s$tau2_percentiles_overview, 1))
cat("\n--- tau2_sd_pct ---\n")
print(round(s$tau2_sd_percentiles_overview, 3))
cat("\n--- print summary dispersion section ---\n")
print(s)
