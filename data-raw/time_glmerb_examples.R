## Wall-clock time for several sampling demos at the n those demos use.
##
##   Rscript data-raw/time_glmerb_examples.R
##
## inst/examples/Ex_glmerb.R uses simulate = FALSE, so
## data-raw/time_examples.R does not measure this cost.
## The all-state book-banning fits are the long ones: n = 3000 and,
## at the default ridge weight, a main stage of about 1380 sweeps.

suppressPackageStartupMessages({
  library(glmbayesCore)
  library(lmebayes)
})

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

rows <- list()

time_one <- function(label, n, expr) {
  cat(sprintf("\n--- %s (n = %d) ---\n", label, n))
  set.seed(1)
  t0 <- proc.time()
  fit <- force(expr)
  seconds <- (proc.time() - t0)[["elapsed"]]
  m <- fit$convergence$m_convergence
  if (is.null(m)) m <- fit$m_convergence
  ridge_m <- fit$convergence$ridge$m_convergence
  n_sweeps <- fit$sweep_history$main$n_sweeps
  if (is.null(n_sweeps)) n_sweeps <- fit$sweep_history$n_sweeps
  cat(sprintf(
    "  %.1f s (%.2f h)   main sweeps = %s   m_convergence = %s   m_ridge = %s\n",
    seconds, seconds / 3600,
    if (is.null(n_sweeps)) NA else n_sweeps,
    if (is.null(m)) NA else m,
    if (is.null(ridge_m)) NA else ridge_m
  ))
  rows[[length(rows) + 1L]] <<- data.frame(
    example = label,
    n = as.integer(n),
    seconds = round(seconds, 1),
    hours = round(seconds / 3600, 2),
    m_convergence = if (is.null(m)) NA_integer_ else as.integer(m),
    stringsAsFactors = FALSE
  )
  invisible(fit)
}

## Ex_14_glmerb_airbnb_small. Poisson, n = 1000.
air <- local({
  data(airbnb_small, package = "bayesrules")
  dat <- airbnb_small
  dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
  dat$walk_c <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
  dat <- dat[complete.cases(dat[, c(
    "reviews", "rating_c", "walk_c", "neighborhood"
  )]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))
  dat
})
form_air <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)
ps_air <- glmbayesCore::Prior_Setup_GLMM(
  form_air, data = air, family = poisson(), pop.pwt = 0.01
)
time_one("airbnb_small Poisson", 1000L, glmerb(
  form_air,
  data = air,
  family = poisson(),
  pfamily_list = glmbayesCore::pfamily_list(ps_air),
  n = 1000L,
  progbar = TRUE
))

book_dat <- local({
  data(book_banning, package = "bayesrules", envir = environment())
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
  dat$violent_i <- as.integer(
    dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
  )
  dat$state <- factor(dat$state)
  dat
})
form_book <- removed_i ~ violent_i + (1 + violent_i || state)

## Ex_16. All states, n = 3000. Main stage is about 1380 sweeps.
ps_all <- glmbayesCore::Prior_Setup_GLMM(
  form_book, data = book_dat, family = binomial(), pop.pwt = 0.01
)
time_one("book banning, all states", 3000L, glmerb(
  form_book,
  data = book_dat,
  family = binomial(),
  pfamily_list = glmbayesCore::pfamily_list(ps_all),
  n = 3000L,
  mode_gap_max = 1.0,
  progbar = TRUE
))

## Ex_22. Same model with sampled tau^2, n = 3000.
ps_ing <- glmbayesCore::Prior_Setup_GLMM(
  form_book,
  data = book_dat,
  family = binomial(),
  pop.pwt = 0.01,
  pop.dispersion.pwt = 0.2
)
time_one("book banning ING, all states", 3000L, glmerb(
  form_book,
  data = book_dat,
  family = binomial(),
  pfamily_list = glmbayesCore::pfamily_list(
    ps_ing, ptypes = "dIndependent_Normal_Gamma"
  ),
  n = 3000L,
  mode_gap_max = 1.0,
  progbar = TRUE
))

cat("\n===== elapsed =====\n")
print(do.call(rbind, rows), row.names = FALSE)
