# Empirical check: does tightening max_disp_perc (tighter, posterior-connected
# dispersion bounds) increase the frequency of the "negligible" UB2 sign-
# violation diagnostic? If the mechanism is "candidate dispersion lands near
# the [low,upp] boundary, where two independently-implemented closed forms
# for the same theoretical RSS_j(d)/UB2_j(d) quantity disagree by floating-
# point round-off", then a WIDER window (larger max_disp_perc, e.g. 0.999)
# should produce far FEWER such warnings than a TIGHTER one (e.g. 0.90).
#
#   cd .../lmebayes
#   Rscript data-raw/check_ub2_warning_frequency.R

source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
dat  <- fx$dat
design <- fx$design
group_levels <- levels(design$group)
form_block <- score_ppvt ~ 1 + distracted_ppvt
form_lmer  <- score_ppvt ~ 1 + distracted_ppvt + (1 + distracted_ppvt || school_id)

fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
vc       <- as.data.frame(lme4::VarCorr(fit_lmer))
tau2 <- c(
  `(Intercept)`   = vc$vcov[vc$grp == "school_id"   & vc$var1 == "(Intercept)"],
  distracted_ppvt = vc$vcov[vc$grp == "school_id.1" & vc$var1 == "distracted_ppvt"]
)
Sigma_re <- diag(tau2, nrow = 2L, ncol = 2L, names = TRUE)
dimnames(Sigma_re) <- list(names(tau2), names(tau2))
mu_re <- lme4::fixef(fit_lmer)[names(tau2)]

pwt_measurement_group <- 0.25
ps_block <- stats::setNames(lapply(group_levels, function(lev) {
  dat_j <- dat[dat$school_id == lev, , drop = FALSE]
  Prior_Setup(form_block, data = dat_j, family = gaussian(), pwt = pwt_measurement_group)
}), group_levels)

run_once <- function(max_disp_perc, n) {
  pfamily_list_j <- stats::setNames(lapply(group_levels, function(lev) {
    ps <- ps_block[[lev]]
    glmbayesCore::dIndependent_Normal_Gamma(
      mu = mu_re, Sigma = Sigma_re,
      shape = ps$shape_ING, rate = ps$rate_gamma,
      max_disp_perc = max_disp_perc
    )
  }), group_levels)

  con <- textConnection("captured", "w", local = TRUE)
  sink(con)
  out <- tryCatch(
    lmbBlock(
      form_block, block = "school_id", pfamily_list = pfamily_list_j,
      data = dat, n = n, use_parallel = FALSE, verbose = FALSE
    ),
    error = function(e) e
  )
  sink()
  close(con)

  n_warn <- sum(grepl("UB2 sign violation", captured))
  n_diag <- sum(grepl("UB2 diagnostics", captured))
  list(n_warn = n_warn, n_diag = n_diag, failed = inherits(out, "error"),
       err = if (inherits(out, "error")) conditionMessage(out) else NA_character_)
}

N <- 500L
for (mdp in c(0.90, 0.95, 0.99, 0.999)) {
  res <- run_once(mdp, N)
  cat(sprintf(
    "max_disp_perc=%.3f | n=%d x 5 schools | UB2-sign-violation warnings=%d | UB2-diagnostic lines=%d | failed=%s%s\n",
    mdp, N, res$n_warn, res$n_diag, res$failed,
    if (res$failed) paste0(" (", res$err, ")") else ""
  ))
}

cat("check_ub2_warning_frequency: done\n")
