## Sweep counts for omega-ridge group likelihoods on the book-banning and
## Airbnb GLMMs.
##
## The rate is rho^G(omega * C^Q): the Gaussian two-block rate at omega times
## the GLM curvature weights evaluated at the joint ICM mode. Corollary 5.1
## of the log-concave note (Core 30) turns that rate into a sweep count for
## the modified posterior. It is not a bound on the distance to the true
## posterior; that still needs the posterior total-variation term.
##
## The reference point is groupef.mode from rglmerb(simulate = FALSE), the
## joint ICM mode. The note's reference point is the conditional mode of
## each group effect at the marginal mode of gamma.
##
## Same formulas as:
##   demo/Ex_16_glmerb_book_banning.R
##   lmebayesCore inst/examples/Ex_rglmerb.R
##   demo/Ex_14_glmerb_airbnb_small.R
##   demo/Ex_13_glmerb_Airbnb.R   (only with argument "full")
##
## Run from the lmebayes package root:
##   Rscript data-raw/ridge_sweep_bounds_book_airbnb.R
##   Rscript data-raw/ridge_sweep_bounds_book_airbnb.R full

suppressPackageStartupMessages(library(lmebayes))

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

tau <- 0.01
omegas <- c(0.01, 0.05, 0.10, 1)
run_full_airbnb <- "full" %in% commandArgs(trailingOnly = TRUE)

## Corollary 5.1, start in the admissible set (the mode is in that set).
## Returns the sweep index k for gamma. The stored group-effect draw lags
## one block update, so the engine count is k + 1.
k_cor51 <- function(rho, q, tau) {
  if (!is.finite(rho) || rho >= 1) {
    return(NA_integer_)
  }
  if (rho <= 0) {
    return(1L)
  }
  z <- 2 * stats::qnorm((1 + tau) / 2)
  D <- sqrt(q / (1 - rho))
  if (D <= z) {
    return(1L)
  }
  as.integer(ceiling(log(D / z) / log(1 / rho)))
}

## Call glmbayesCore directly. The installed lmebayes::model_setup /
## Prior_Setup_GLMM copies still look up helpers inside lmebayesCore, and
## .lmebayes_normalize_family now lives in glmbayesCore.
fit_mode <- function(form, data, family) {
  design <- glmbayesCore::model_setup(form, data = data, family = family)
  ps <- glmbayesCore::Prior_Setup_GLMM(
    form, data = data, family = family, pop.pwt = 0.01
  )
  lmebayesCore::rglmerb(
    n            = 1L,
    design       = design,
    pfamily_list = glmbayesCore::pfamily_list(ps),
    family       = family,
    simulate     = FALSE,
    verbose      = FALSE,
    progbar      = FALSE
  )
}

ridge_table <- function(label, fit) {
  design <- fit$design
  b_mode <- fit$groupef.mode
  re_names <- colnames(b_mode)
  group_levels <- rownames(b_mode)
  pf <- fit$prior$pfamily_list
  P <- glmbayesCore:::.rLMM_P_from_pfamily_list(pf, re_names)
  prior_b1 <- list(P = P)

  mw <- glmbayesCore:::two_block_mode_weights(
    x            = design$D,
    group        = design$group,
    b_mode       = b_mode,
    family       = fit$family,
    wt           = if (is.null(design$weights)) 1 else design$weights,
    offset       = if (is.null(design$offset)) 0 else design$offset,
    group_levels = group_levels
  )

  rate_at <- function(weights) {
    glmbayesCore::two_block_rate_from_pfamily_list(
      x                 = design$D,
      group             = design$group,
      x_hyper           = design$W,
      prior_list_block1 = prior_b1,
      pfamily_list      = pf,
      weights           = weights,
      family            = fit$family,
      group_levels      = group_levels,
      warn_slow         = FALSE
    )
  }

  ## What rglmerb uses today: unscaled mode weights, Nygren Theorem 3.
  ## A local-Gaussian heuristic for the true logit or Poisson posterior.
  rate_mode <- rate_at(mw$weights)
  m_heuristic <- tryCatch(
    glmbayesCore:::two_block_l_for_tv(
      rate_mode, tau, method = "theorem3", warn = FALSE
    ) + 1L,
    error = function(e) NA_integer_
  )

  cat("\n=== ", label, " ===\n", sep = "")
  cat(
    "J = ", length(group_levels),
    ", q = ", rate_mode$dims$q,
    ", p_re = ", rate_mode$dims$p_re,
    ", n_obs = ", rate_mode$dims$l2,
    "\n",
    sep = ""
  )
  cat(sprintf(
    "curvature weights at ICM mode: min = %.4g, median = %.4g, max = %.4g\n",
    min(mw$weights), stats::median(mw$weights), max(mw$weights)
  ))
  cat(sprintf(
    "current heuristic (unscaled mode weights, Theorem 3): rho = %.6f, m = %s\n",
    rate_mode$lambda_star,
    if (is.na(m_heuristic)) "NA" else as.character(m_heuristic)
  ))
  cat(
    "Ridge rows use Corollary 5.1 at weights = omega * curvature. ",
    "m = k + 1. omega = 1 is the tangent-quadratic (Gaussian) model; ",
    "its Theorem 3 count is printed beside Corollary 5.1.\n",
    sep = ""
  )

  rows <- lapply(omegas, function(omega) {
    rate <- rate_at(omega * mw$weights)
    rho <- rate$lambda_star
    k <- k_cor51(rho, rate$dims$q, tau)
    ## omega = 1 is the unscaled mode-weight model, so its Theorem 3
    ## count is the heuristic already computed above.
    m_t3 <- if (abs(as.numeric(omega) - 1) < 1e-12) {
      m_heuristic
    } else {
      NA_integer_
    }
    data.frame(
      omega = omega,
      rho = rho,
      k = k,
      m_ridge = if (is.na(k)) NA_integer_ else k + 1L,
      m_theorem3 = m_t3,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE, digits = 6)
  invisible(tab)
}

## --- book banning, full data, intercept and violent slope -----------------
book_full <- local({
  data(book_banning, package = "bayesrules", envir = environment())
  dat <- book_banning[, c("state", "removed", "violent", "antifamily", "language")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  form <- removed_i ~ violent_i + (1 + violent_i || state)
  fit_mode(form, dat, binomial())
})

## --- book banning, eight largest states, same random-slope formula --------
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
  dat <- droplevels(subset(dat, state %in% keep))
  form <- removed_i ~ violent_i + (1 + violent_i || state)
  fit_mode(form, dat, binomial())
})

## --- Airbnb small, Poisson, neighborhood intercept and rating slope -------
airbnb_small_fit <- local({
  data(airbnb_small, package = "bayesrules", envir = environment())
  dat <- airbnb_small
  dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
  dat$walk_c <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
  dat <- dat[complete.cases(dat[, c(
    "reviews", "rating_c", "walk_c", "neighborhood"
  )]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))
  form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)
  fit_mode(form, dat, poisson())
})

ridge_table(
  "book_banning, binomial, violent slope, all states",
  book_full
)
ridge_table(
  "book_banning, binomial, violent slope, 8 largest states",
  book_eight
)
ridge_table(
  "airbnb_small, Poisson, rating slope",
  airbnb_small_fit
)

if (isTRUE(run_full_airbnb)) {
  airbnb_full <- local({
    data(airbnb, package = "bayesrules", envir = environment())
    dat <- airbnb
    dat$rating_c <- dat$rating - mean(dat$rating)
    dat$log_price_c <- scale(log(dat$price + 1))[, 1]
    dat$walk_c <- dat$walk_score - mean(dat$walk_score)
    dat$transit_c <- dat$transit_score - mean(dat$transit_score)
    dat <- dat[complete.cases(dat[, c(
      "reviews", "rating", "rating_c", "price", "log_price_c",
      "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
    )]), ]
    dat$neighborhood <- droplevels(factor(dat$neighborhood))
    form <- reviews ~
      walk_c + transit_c +
      rating_c + log_price_c +
      walk_c:rating_c + transit_c:log_price_c +
      (1 + rating_c + log_price_c || neighborhood)
    fit_mode(form, dat, poisson())
  })
  ridge_table(
    "airbnb, Poisson, rating and log-price slopes",
    airbnb_full
  )
} else {
  cat(
    "\nFull airbnb model skipped. Re-run with argument 'full' to include it.\n"
  )
}
