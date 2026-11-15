## Check that rglmerb stores the Corollary 5.1 ridge sweep count and uses
## it as the main-stage length. The pilot stays on the Theorem 3 schedule.
##
## Eight-state book banning, binomial, random intercept and violent slope.
## At rwt = 0.01 the main run is a few hundred sweeps. At rwt = 1 it is
## a handful. The likelihood draw is still the original binomial likelihood.
##
## Run from the lmebayes package root, after installing lmebayesCore:
##   Rscript data-raw/check_ridge_sweep_count.R

suppressPackageStartupMessages({
  library(glmbayesCore)
  library(lmebayesCore)
})

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

tau <- 0.01

## Same inversion as .rGLMM_k_cor51(). q is the dimension of gamma.
k_cor51 <- function(rho, q, tau) {
  if (!is.finite(rho) || !is.finite(q) || q < 1 || rho >= 1) {
    return(NA_integer_)
  }
  if (rho <= 0) {
    return(1L)
  }
  z <- 2 * stats::qnorm((1 + tau) / 2)
  D <- sqrt(q / (1 - rho))
  if (!is.finite(D) || !is.finite(z) || D <= z) {
    return(1L)
  }
  as.integer(ceiling(log(D / z) / log(1 / rho)))
}

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

## Independent of the stored ridge list: mode weights, then rwt times those
## weights, then Corollary 5.1. Matches .rGLMM_ridge_sweep_bound().
independent_ridge <- function(fit, rwt, tau) {
  b_mode <- fit$groupef.mode
  design <- fit$design
  mode_w <- lmebayesCore:::two_block_mode_weights(
    x            = design$D,
    group        = design$group,
    b_mode       = b_mode,
    family       = fit$family,
    group_levels = rownames(b_mode)
  )
  rate <- lmebayesCore::two_block_rate_from_pfamily_list(
    x                 = design$D,
    group             = design$group,
    x_hyper           = design$W,
    prior_list_block1 = fit$dispprior_list,
    pfamily_list      = fit$pfamily_list,
    weights           = rwt * mode_w$weights,
    family            = fit$family,
    group_levels      = rownames(b_mode),
    warn_slow         = FALSE
  )
  rho <- rate$lambda_star
  q <- rate$dims$q
  k <- k_cor51(rho, q, tau)
  k_installed <- lmebayesCore:::.rGLMM_k_cor51(rho, q, tau)
  list(
    lambda_star = rho,
    q           = q,
    k           = k,
    k_installed = k_installed,
    m_ridge     = if (is.na(k)) NA_integer_ else as.integer(k + 1L)
  )
}

book_eight <- local({
  data(book_banning, package = "bayesrules", envir = environment())
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  dat$state <- factor(dat$state)
  keep <- names(sort(table(dat$state), decreasing = TRUE))[seq_len(8L)]
  droplevels(subset(dat, state %in% keep))
})

form <- removed_i ~ violent_i + (1 + violent_i || state)
design <- glmbayesCore::model_setup(form, data = book_eight, family = binomial())
ps <- glmbayesCore::Prior_Setup_GLMM(
  form, data = book_eight, family = binomial(), pop.pwt = 0.01
)
pf <- glmbayesCore::pfamily_list(ps)

check_one <- function(rwt) {
  set.seed(1)
  logged <- utils::capture.output(
    fit <- lmebayesCore::rglmerb(
      n            = 1L,
      design       = design,
      pfamily_list = pf,
      family       = binomial(),
      rwt          = rwt,
      verbose      = TRUE,
      progbar      = FALSE
    ),
    type = "output"
  )
  ridge <- fit$convergence_info$ridge
  indep <- independent_ridge(fit, rwt, tau)
  ran <- fit$m_convergence
  n_sweeps <- fit$sweep_history$n_sweeps

  cat(sprintf(
    "\nrwt = %g\n  ran m_convergence = %s (%s main sweeps)\n  stored ridge    k = %s, m = %s, lambda* = %.6f, q = %s\n  independent     k = %s, m = %s, lambda* = %.6f, q = %s\n",
    rwt,
    ran, n_sweeps,
    ridge$k, ridge$m_convergence, ridge$lambda_star, ridge$q,
    indep$k, indep$m_ridge, indep$lambda_star, indep$q
  ))

  if (!isTRUE(all.equal(fit$rwt, rwt))) {
    fail("fit$rwt is ", fit$rwt, ", expected ", rwt)
  }
  if (is.null(ridge)) {
    fail("convergence_info$ridge is missing")
  }
  if (!isTRUE(all.equal(ridge$lambda_star, indep$lambda_star))) {
    fail("stored lambda*_ridge does not match the independent rate")
  }
  if (!identical(as.integer(ridge$q), as.integer(indep$q))) {
    fail("stored q does not match the independent rate")
  }
  if (!identical(as.integer(ridge$k), as.integer(indep$k))) {
    fail("stored k does not match the local Corollary 5.1 inversion")
  }
  if (!identical(as.integer(indep$k), as.integer(indep$k_installed))) {
    fail("installed .rGLMM_k_cor51() does not match the local inversion")
  }
  if (!identical(as.integer(ridge$m_convergence), as.integer(indep$m_ridge))) {
    fail("stored m_ridge is not k + 1")
  }
  if (!identical(as.integer(n_sweeps), as.integer(ran))) {
    fail("main sweep_history has ", n_sweeps, " sweeps; m_convergence is ", ran)
  }
  if (!identical(as.integer(ran), as.integer(ridge$m_convergence))) {
    fail(
      "the main run used ", ran,
      " sweeps; the ridge count is ", ridge$m_convergence
    )
  }
  text <- paste(logged, collapse = "\n")
  needle <- sprintf("m_ridge = %d", ridge$m_convergence)
  if (!grepl("Corollary 5.1", text, fixed = TRUE) ||
      !grepl("main stage", text, fixed = TRUE) ||
      !grepl(needle, text, fixed = TRUE)) {
    fail("verbose output did not report '", needle, "' as the main stage")
  }
  hit <- logged[grepl("Corollary 5.1", logged, fixed = TRUE) |
                  grepl("m_ridge =", logged, fixed = TRUE)]
  cat("  print: ", paste(hit, collapse = "\n         "), "\n", sep = "")
  invisible(fit)
}

cat("Eight-state book banning. The main run uses the ridge count. The pilot stays on Theorem 3.\n")
check_one(0.01)
check_one(1)
cat("\nPASS: the main run used the stored ridge count.\n")
