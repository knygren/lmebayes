# Quick check: window_diagnostics flow ps -> dGamma_list -> lmerb -> summary_sigma2
pkg_root <- normalizePath(getwd())
setwd(pkg_root)
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

fx <- .prepare_small5_all_full_rank_manual()
ps <- Prior_Setup_lmebayes(
  fx$form,
  data = fx$dat,
  pwt = 0.01,
  pwt_measurement = 0.1,
  max_disp_perc = 0.8
)
disp <- dGamma_list(ps, warn_asymmetric = FALSE)
pf <- pfamily_list(ps)
fit <- lmerb(
  fx$form,
  data = fx$dat,
  pfamily_list = pf,
  dispersion_ranef = disp,
  n = 2L,
  progbar = FALSE
)
stopifnot(!is.null(fit$prior$window_diagnostics))
flagged <- fit$prior$window_diagnostics$group[
  fit$prior$window_diagnostics$asymmetric_window
]
stopifnot(all(c("6", "33", "41") %in% flagged))
s <- summary_sigma2(fit, type = "prior")
stopifnot(all(c("R_lo", "asymmetric_window") %in% colnames(s$prior)))
cat("test_window_diag_quick: OK; flagged:", paste(flagged, collapse = ", "), "\n")
