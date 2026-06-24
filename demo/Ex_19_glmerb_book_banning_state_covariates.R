## Demo: glmerb() binomial workflow on bayesrules::book_banning (state covariates)
##
## Extends demo/Ex_18_glmerb_book_banning_all_reasons.R: state-level covariates
## enter as intercept-only predictors of the state random intercept (X_hyper for
## "(Intercept)" only), not as moderators of random slopes (cf. Ex_17).
##
## Covariates (Bayes Rules book_banning; hs_grad_rate / college_grad_rate are
## already expressed relative to the national state average in the data):
##   political_value_index  -- negative = Republican, 0 = neutral, positive = Democrat
##   median_income          -- state median income (scaled at state level below)
##   hs_grad_rate           -- high-school graduation rate (relative in raw data)
##   college_grad_rate      -- college graduation rate (relative in raw data)
##
## All four covariates together fail lme4 checkConv on the reference glmer fit
## (max|grad| ~ 0.029).  This demo uses the largest passing subset:
##   political_value_index, median_income, hs_grad_rate
## (college_grad_rate omitted so Prior_Setup_lmebayes() can calibrate from glmer).
##
##   demo("Ex_19_glmerb_book_banning_state_covariates", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

data(book_banning, package = "bayesrules", envir = environment())

reasons <- c("explicit", "language", "violent")
state_covars <- c(
  "political_value_index", "median_income",
  "hs_grad_rate", "college_grad_rate"
)
## Subset used for glmerb (reference glmer passes checkConv)
state_covars_fit <- c(
  "political_value_index", "median_income", "hs_grad_rate"
)

dat <- book_banning[, c("state", "removed", reasons, state_covars)]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(
  dat$removed == TRUE | dat$removed == 1L | dat$removed == "1"
)
for (v in reasons) {
  dat[[paste0(v, "_i")]] <- as.integer(
    dat[[v]] == TRUE | dat[[v]] == 1L | dat[[v]] == "TRUE"
  )
}

## One row per state; z-score covariates for numerical stability in glmer
st_cov <- unique(dat[, c("state", state_covars)])
for (v in state_covars) {
  st_cov[[paste0(v, "_c")]] <- as.numeric(scale(st_cov[[v]]))
}
dat <- merge(
  dat[, setdiff(names(dat), state_covars)],
  st_cov[, c("state", paste0(state_covars, "_c"))],
  by = "state",
  all.x = TRUE
)

covar_rhs <- function(covars) {
  paste0(covars, "_c")
}

re_rhs <- "(1 + explicit_i + language_i + violent_i || state)"

build_form <- function(covars = state_covars_fit) {
  fe <- c("explicit_i", "language_i", "violent_i", covar_rhs(covars))
  stats::as.formula(
    paste("removed_i ~", paste(fe, collapse = " + "), "+", re_rhs)
  )
}

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
## Population-mean slopes map to k::(Intercept); intercept-block columns of
## X_hyper[["(Intercept)"]] map to the same-named glmer fixed effects.
print_glmer_glmerb_fixed_compare <- function(glmerb_fit,
                                             glmer_fit = glmerb_fit$glmer,
                                             digits = 4L) {
  if (is.null(glmer_fit)) {
    stop("Need a reference glmer fit (fit$glmer).", call. = FALSE)
  }
  re_names <- glmerb_fit$model_setup$re_coef_names
  smb      <- summary(glmerb_fit)
  gtab     <- summary(glmer_fit)$coefficients

  row_for <- function(k, col) {
    fe_nm <- if (identical(k, "(Intercept)")) col else k
    if (!fe_nm %in% rownames(gtab)) {
      stop("glmer fixed effect '", fe_nm, "' not found.", call. = FALSE)
    }
    g  <- gtab[fe_nm, , drop = TRUE]
    pt <- smb$fixef[[k]]$coefficients[col, , drop = TRUE]
    pm <- smb$fixef[[k]]$coefficients1[col, "Prior Mean"]
    is_global_int <- identical(k, "(Intercept)") && identical(col, "(Intercept)")
    is_slope <- !identical(k, "(Intercept)")
    data.frame(
      parameter  = fe_nm,
      glmer_est  = unname(g["Estimate"]),
      glmer_se   = unname(g["Std. Error"]),
      glmer_z    = unname(g["z value"]),
      glmer_p    = unname(g["Pr(>|z|)"]),
      glmer_p_1s = unname(g["Pr(>|z|)"]) / 2,
      prior_mean = unname(pm),
      post_mean  = unname(pt["Post.Mean"]),
      post_sd    = unname(pt["Post.Sd"]),
      prior_tail = unname(pt["Pr(Prior_tail)"]),
      is_slope   = is_slope,
      is_global_int = is_global_int,
      stringsAsFactors = FALSE
    )
  }

  rows <- do.call(rbind, lapply(re_names, function(k) {
    cols <- rownames(smb$fixef[[k]]$coefficients)
    do.call(rbind, lapply(cols, row_for, k = k))
  }))

  cat("\n=== Fixed effects: glmer (Wald) vs glmerb (Block 2 hyperparameters) ===\n\n")
  cat(
    "  glmer: Wald Pr(>|z|) is two-sided; glmer_p/2 is the one-sided counterpart.\n",
    "  glmerb: Pr(Prior_tail) is one-sided (posterior mass on one side of the\n",
    "    prior mean).  For slopes and state covariates (prior mean 0), compare\n",
    "    glmer_p/2 to Pr(Prior_tail).  Global (Intercept) uses a null-model\n",
    "    prior mean, not 0 -- not comparable on that scale.\n\n",
    sep = ""
  )

  w_par <- max(nchar(rows$parameter), nchar("parameter"))
  fmt_n <- function(x) formatC(x, digits = digits, format = "f")
  fmt_p <- function(x) {
    if (is.na(x)) "           —"
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
    p1s <- if (!rows$is_global_int[i]) rows$glmer_p_1s[i] else NA_real_
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

## --- Reference glmer: all four covariates (fails checkConv) -------------------
form_all_covars <- build_form(state_covars)
print_glmer_check(
  lme4::glmer(form_all_covars, data = dat, family = binomial()),
  "all four state covariates (checkConv failure expected)"
)

## --- glmerb formula (three covariates; used below) ----------------------------
form_glmerb <- build_form(state_covars_fit)

print_glmer_check(
  fit_glmer <- lme4::glmer(form_glmerb, data = dat, family = binomial()),
  "three state covariates (Prior_Setup / glmerb reference)"
)

design <- model_setup(form_glmerb, data = dat, family = binomial())
cat("\n=== model_setup (X_hyper per RE component) ===\n\n")
for (k in design$re_coef_names) {
  cat(k, ":\n")
  print(colnames(design$X_hyper[[k]]))
}
cat(
  "\nState covariates appear only under (Intercept); slope components have\n",
  "intercept-only X_hyper (population-mean slopes for challenge reasons).\n\n",
  sep = ""
)

ps <- Prior_Setup_lmebayes(
  form_glmerb,
  data   = dat,
  family = binomial(),
  pwt    = 0.01
)

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

print_glmer_glmerb_fixed_compare(fit, glmer_fit = fit_glmer)
