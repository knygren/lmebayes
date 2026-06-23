## Demo: glmerb() binomial workflow on bayesrules::book_banning (all challenge reasons)
##
## Full challenge-reason model with six binary predictors and matching state
## random intercepts and slopes (lmebayes population-mean-slope requirement).
## Two formulas below (after the removal summary):
##   form_glmerb_full  -- all six challenge reasons (Ch. 18 target model)
##   form_glmerb       -- reduced model used for glmerb / Prior_Setup (edit here)
## Predictors must exist in dat as <reason>_i (built from reasons).##
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
  }
  cat("\nFixed effects:\n")
  print(lme4::fixef(fit))
  cat("\nRandom effects (VarCorr):\n")
  print(lme4::VarCorr(fit))
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
  lme4::glmer(form_glmerb, data = dat, family = binomial()),
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
