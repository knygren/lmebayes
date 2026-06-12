pkgload::load_all(export_all = FALSE, quiet = TRUE)

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

run_case <- function(label, ...) {
  ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01, ...)
  pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
  out <- capture.output(
    fit <- lmerb(form_lmer, data = dat,
                 pfamily_list = pf,
                 dispersion_ranef = ps$dispersion_ranef,
                 n = 10L, seed = 1L)
  )
  re_names <- names(ps$prior_list)
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("  pwt_disp     : %s  [%s]\n",
              paste(sprintf("%.4g", ps$pwt_dispersion), collapse = ", "),
              attr(ps$pwt_dispersion, "source")))
  cat(sprintf("  n_prior_disp : %s\n",
              paste(sprintf("%.4g", ps$n_prior_dispersion), collapse = ", ")))
  for (k in re_names) {
    pr     <- pf[[k]]$prior_list
    tau2_k <- unname(ps$prior_list[[k]]$dispersion_fixef)
    cat(sprintf(
      "  [%s] shape = %.3f, rate = %.3f, disp_lower = %.4f (tau^2 = %.4f, ratio = %.3f)\n",
      k, pr$shape, pr$rate, pr$disp_lower, tau2_k, pr$disp_lower / tau2_k
    ))
  }
  cat(sprintf("  lambda* = %.4f,  m_min = %d  (method: %s)\n",
              fit$convergence$lambda_star, fit$convergence$m_min,
              fit$convergence$method))
  invisible(fit$convergence)
}

cv1 <- run_case("Default pwt_disp (0.2)")
cv2 <- run_case("pwt_dispersion = 0.2", pwt_dispersion = 0.2)

cat(sprintf(
  "\nSummary: m_min %d -> %d  (lambda* %.4f -> %.4f)\n",
  cv1$m_min, cv2$m_min, cv1$lambda_star, cv2$lambda_star
))
