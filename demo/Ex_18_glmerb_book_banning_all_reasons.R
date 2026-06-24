## Demo: glmerb() binomial workflow on bayesrules::book_banning (all challenge reasons)
##
## Full challenge-reason model with six binary predictors and matching state
## random intercepts and slopes (lmebayes population-mean-slope requirement).
## Two formulas below (after the removal summary):
##   form_glmerb_full  -- all six challenge reasons (Ch. 18 target model)
##   form_glmerb       -- reduced model used for glmerb / Prior_Setup (edit here)
## Predictors must exist in dat as <reason>_i (built from reasons).
##
## Per-sweep Block 2 diagnostics are stored on the fit object, not printed
## during sampling. After fitting:
##   print(fit, sweep_history = TRUE, max_sweeps = 5)
##   print(fit$sweep_history$main, components = c("violent_i", "(Intercept)"))
##
##   demo("Ex_18_glmerb_book_banning_all_reasons", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

data(book_banning, package = "bayesrules", envir = environment())

## Raw challenge-reason flags in book_banning (edit to subset predictors)
reasons <- c(
  "explicit", "antifamily", "occult",
  "language", "lgbtq", "violent"
)

dat <- book_banning[, c("state", "removed", reasons)]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(
  dat$removed == TRUE | dat$removed == 1L | dat$removed == "1"
)
for (v in reasons) {
  dat[[paste0(v, "_i")]] <- as.integer(
    dat[[v]] == TRUE | dat[[v]] == 1L | dat[[v]] == "TRUE"
  )
}

cat("\n=== Challenge reasons: n and mean removal rate by level ===\n\n")
reason_summary <- do.call(
  rbind,
  lapply(reasons, function(v) {
    x <- as.integer(
      dat[[v]] == TRUE | dat[[v]] == 1L | dat[[v]] == "TRUE" | dat[[v]] == "1"
    )
    t(sapply(0:1, function(lvl) {
      idx <- which(x == lvl)
      c(
        reason       = v,
        level        = lvl,
        n            = length(idx),
        mean_removed = if (length(idx)) mean(dat$removed_i[idx]) else NA_real_
      )
    }))
  })
)
rownames(reason_summary) <- NULL
reason_summary <- as.data.frame(reason_summary, stringsAsFactors = FALSE)
reason_summary$level <- ifelse(reason_summary$level == 1L, "yes (1)", "no (0)")
reason_summary$mean_removed <- as.numeric(reason_summary$mean_removed)
print(
  reason_summary[order(reason_summary$reason, reason_summary$level), ],
  row.names = FALSE
)

print_glmer_check <- function(fit, label) {
  cat("\n--- glmer:", label, "---\n\n")
  cat("Optimizer status (0 = OK):", fit@optinfo$conv$opt, "\n")
  lme4c <- fit@optinfo$conv$lme4
  if (!is.null(lme4c$messages) && length(lme4c$messages)) {
    cat("checkConv:\n")
    for (m in lme4c$messages) cat("  ", m, sep = "")
    cat("\n")
  }
  sm <- summary(fit)
  cat("\nFixed effects (Wald z-tests):\n")
  print(sm$coefficients)
  cat("\nRandom effects (VarCorr):\n")
  print(lme4::VarCorr(fit))
}

## Side-by-side: classical glmer Wald tests vs glmerb Block 2 hyperparameters.
## glmer fixed effects map to glmerb RE components:
##   (Intercept) -> (Intercept)::(Intercept); slope_k -> k::(Intercept).
print_glmer_glmerb_fixed_compare <- function(glmerb_fit,
                                             glmer_fit = glmerb_fit$glmer,
                                             digits = 4L) {
  if (is.null(glmer_fit)) {
    stop("Need a reference glmer fit (fit$glmer).", call. = FALSE)
  }
  re_names <- glmerb_fit$model_setup$re_coef_names
  smb      <- summary(glmerb_fit)
  gtab     <- summary(glmer_fit)$coefficients

  rows <- do.call(rbind, lapply(re_names, function(k) {
    fe_nm <- if (identical(k, "(Intercept)")) "(Intercept)" else k
    if (!fe_nm %in% rownames(gtab)) {
      stop("glmer fixed effect '", fe_nm, "' not found.", call. = FALSE)
    }
    g  <- gtab[fe_nm, , drop = TRUE]
    pt <- smb$fixef[[k]]$coefficients["(Intercept)", , drop = TRUE]
    pm <- smb$fixef[[k]]$coefficients1["(Intercept)", "Prior Mean"]
    data.frame(
      parameter    = fe_nm,
      glmer_est    = unname(g["Estimate"]),
      glmer_se     = unname(g["Std. Error"]),
      glmer_z      = unname(g["z value"]),
      glmer_p      = unname(g["Pr(>|z|)"]),
      glmer_p_1s   = unname(g["Pr(>|z|)"]) / 2,
      prior_mean   = unname(pm),
      post_mean    = unname(pt["Post.Mean"]),
      post_sd      = unname(pt["Post.Sd"]),
      prior_tail   = unname(pt["Pr(Prior_tail)"]),
      is_slope     = !identical(k, "(Intercept)"),
      stringsAsFactors = FALSE
    )
  }))

  cat("\n=== Fixed effects: glmer (Wald) vs glmerb (Block 2 hyperparameters) ===\n\n")
  cat(
    "  glmer: Wald Pr(>|z|) is two-sided; glmer_p/2 is the one-sided counterpart.\n",
    "  glmerb: Pr(Prior_tail) is one-sided (posterior mass on one side of the\n",
    "    prior mean).  For slopes (prior mean 0), compare glmer_p/2 to Pr(Prior_tail).\n",
    "  (Intercept) uses a null random-intercept prior mean, not 0 -- not comparable.\n\n",
    sep = ""
  )

  w_par <- max(nchar(rows$parameter), nchar("parameter"))
  fmt_n <- function(x) formatC(x, digits = digits, format = "f")
  fmt_p <- function(x) {
    if (is.na(x)) "           -"
    else if (x < 0.001) formatC(x, digits = 2, format = "e")
    else formatC(x, digits = digits, format = "f")
  }

  hdr <- sprintf(
    paste0(
      "  %-*s  %12s  %12s  %10s  %10s  |  %12s  %12s  %12s  %10s"
    ),
    w_par, "parameter",
    "glmer_est", "glmer_se", "glmer_p/2", "Pr(prior)",
    "prior_mean", "post_mean", "post_sd", "glmer_p"
  )
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(hdr, "\n", sep, "\n", sep = "")
  for (i in seq_len(nrow(rows))) {
    p1s <- if (rows$is_slope[i]) rows$glmer_p_1s[i] else NA_real_
    pt  <- rows$prior_tail[i]
    cat(sprintf(
      paste0(
        "  %-*s  %12s  %12s  %10s  %10s  |  %12s  %12s  %12s  %10s"
      ),
      w_par, rows$parameter[i],
      fmt_n(rows$glmer_est[i]), fmt_n(rows$glmer_se[i]),
      fmt_p(p1s), fmt_p(pt),
      fmt_n(rows$prior_mean[i]),
      fmt_n(rows$post_mean[i]), fmt_n(rows$post_sd[i]),
      fmt_p(rows$glmer_p[i])
    ), "\n", sep = "")
  }
  cat("\n  glmer_p column (right): two-sided Wald p-value for reference.\n\n")
  invisible(rows)
}

## --- Full glmerb formula (six reasons) ----------------------------------------
form_glmerb_full <- removed_i ~
  explicit_i + antifamily_i + occult_i + language_i + lgbtq_i + violent_i +
  (1 + explicit_i + antifamily_i + occult_i + language_i +
     lgbtq_i + violent_i || state)

print_glmer_check(
  lme4::glmer(form_glmerb_full, data = dat, family = binomial()),
  "FULL formula (six reasons)"
)

## --- Reduced glmerb formula (edit manually; used for glmerb below) ------------
form_glmerb <- removed_i ~
  explicit_i + language_i + violent_i +
  (1 + explicit_i + language_i + violent_i || state)

print_glmer_check(
  fit_glmer_reduced <- lme4::glmer(form_glmerb, data = dat, family = binomial()),
  "REDUCED formula (used for model_setup / Prior_Setup / glmerb)"
)

## Reference glmer (Ch. 18 book: fixed effects only, random intercept)
form_book <- removed ~
  explicit + antifamily + occult + language + lgbtq + violent +
  (1 | state)

design <- model_setup(form_glmerb, data = dat, family = binomial())
cat("\n=== model_setup (X_hyper per RE component) ===\n\n")
for (k in design$re_coef_names) {
  cat(k, ":\n")
  print(colnames(design$X_hyper[[k]]))
}

ps <- tryCatch(
  Prior_Setup_lmebayes(
    form_glmerb,
    data   = dat,
    family = binomial(),
    pwt    = 0.01
  ),
  error = function(e) {
    cat("\nPrior_Setup_lmebayes() stopped:\n  ",
        conditionMessage(e), "\n\n", sep = "")
    NULL
  }
)

if (!is.null(ps)) {
  print(ps)

  set.seed(42L)
  fit <- glmerb(
    form_glmerb,
    data         = dat,
    family       = binomial(),
    pfamily_list = pfamily_list(ps),
    n            = 1000L,
    mode_gap_max = 1.0,
    progbar      = FALSE
  )

  cat("\nm_convergence used:", fit$convergence$m_convergence, "\n")
  cat(sprintf(
    "Pilot vs mode (chi-squared): p = %.4g (n_pilot = %d)\n",
    fit$pilot_chisq$p_value,
    fit$pilot_chisq$n_pilot
  ))

  lmebayes:::print_coef_means(fit)
  print(fit)
  summary(fit)

  print_glmer_glmerb_fixed_compare(fit, glmer_fit = fit_glmer_reduced)
}

cat("\n--- Reference glmer (random intercept; all six reasons) ---\n")
cat(
  "Note: lme4 may warn that max|grad| exceeds 0.002 (~0.048 here).\n",
  "  Overlapping/rare reason flags make the fixed-effect block ill-conditioned;\n",
  "  this fit is for Ch. 18 comparison only (not used for glmerb priors).\n\n"
)
fit_book <- lme4::glmer(form_book, data = dat, family = binomial())
cat(
  "Reference glmer optimizer status (0 = OK):",
  fit_book@optinfo$conv$opt, "\n"
)
print(lme4::fixef(fit_book))
