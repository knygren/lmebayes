# Regression test: pfamily_list() method for lmebayes_prior_setup objects.
#
# Builds named lists of glmbayesCore pfamily objects (one per random-effect
# coefficient) from Prior_Setup_lmebayes() output.  Checks:
#   - dNormal path: mu/Sigma/dispersion match prior_list fields exactly.
#   - dIndependent_Normal_Gamma path: mu/Sigma match; shape/rate follow the
#     glmbayesCore default calibration (shape_ING with b_0 = tau2*(shape-1),
#     n0 = J * pwt/(1-pwt)) so the inverse-Gamma prior mean of the
#     dispersion is exactly tau^2_k.
#   - ptypes recycling (single string), positional vectors, named vectors
#     (any order), and list input.
#   - validation errors: bad type strings, wrong length, wrong names.
#
#   Rscript data-raw/test_pfamily_list.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

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

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
re_names <- names(ps$prior_list)
stopifnot(length(re_names) == 3L)

J       <- nlevels(ps$design$groups)
## Default dispersion-prior weight is a flat 0.2 (decoupled from pwt), so
## n_prior_dispersion defaults to J * 0.2 / 0.8 = J/4 per component.
n_prior <- J * 0.2 / (1 - 0.2)

stopifnot(
  is.numeric(ps$n_prior_dispersion),
  identical(names(ps$n_prior_dispersion), re_names),
  isTRUE(all.equal(as.vector(ps$n_prior_dispersion), rep(n_prior, 3L)))
)

## --- 1. default: all dNormal ------------------------------------------------
pf <- pfamily_list(ps)
stopifnot(
  is.list(pf),
  identical(names(pf), re_names),
  all(vapply(pf, inherits, logical(1L), what = "pfamily")),
  all(vapply(pf, function(p) p$pfamily, character(1L)) == "dNormal")
)
for (k in re_names) {
  pl <- ps$prior_list[[k]]
  pr <- pf[[k]]$prior_list
  stopifnot(
    isTRUE(all.equal(as.numeric(pr$mu), unname(pl$mu_fixef))),
    isTRUE(all.equal(unname(as.matrix(pr$Sigma)), unname(pl$Sigma_fixef))),
    isTRUE(all.equal(pr$dispersion, unname(pl$dispersion_fixef))),
    identical(pr$ddef, FALSE)
  )
}
cat("dNormal default path: OK\n")

## --- 2. all dIndependent_Normal_Gamma ---------------------------------------
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
stopifnot(
  identical(names(pf_ing), re_names),
  all(vapply(pf_ing, function(p) p$pfamily, character(1L)) ==
        "dIndependent_Normal_Gamma")
)
for (k in re_names) {
  pl  <- ps$prior_list[[k]]
  pr  <- pf_ing[[k]]$prior_list
  p_k <- length(pl$mu_fixef)
  shape_exp <- (n_prior + 1) / 2 + p_k / 2
  ## glmbayesCore default rate b_0 = tau2 * (n_prior + p_k - 1)/2
  ##                              = tau2 * (shape_ING - 1): mean-matched.
  rate_exp  <- unname(pl$dispersion_fixef) * (n_prior + p_k - 1) / 2
  ## Default truncation window: central 98% prior-mass interval for the
  ## inverse-Gamma dispersion (0.01 and 0.99 quantiles); both bounds are
  ## required for sampling so the tau^2 window is fixed across Gibbs sweeps.
  disp_lower_exp <- 1 / stats::qgamma(0.99, shape = shape_exp, rate = rate_exp)
  disp_upper_exp <- 1 / stats::qgamma(0.01, shape = shape_exp, rate = rate_exp)
  stopifnot(
    isTRUE(all.equal(as.numeric(pr$mu), unname(pl$mu_fixef))),
    isTRUE(all.equal(unname(as.matrix(pr$Sigma)), unname(pl$Sigma_fixef))),
    isTRUE(all.equal(pr$shape, shape_exp)),
    isTRUE(all.equal(pr$rate, rate_exp)),
    isTRUE(all.equal(pr$disp_lower, disp_lower_exp)),
    isTRUE(all.equal(pr$disp_upper, disp_upper_exp)),
    pr$disp_lower > 0,
    pr$disp_upper > pr$disp_lower
  )
  ## Mean-matching identity: rate = tau2 * (shape - 1), so the implied
  ## inverse-Gamma prior on the dispersion has mean exactly tau^2_k for
  ## every n0 and p_k (and the 98% window brackets tau^2_k).
  tau2_k <- unname(pl$dispersion_fixef)
  stopifnot(
    isTRUE(all.equal(pr$rate / (pr$shape - 1), tau2_k)),
    pr$disp_lower < tau2_k,
    pr$disp_upper > tau2_k
  )
  ## Single source of truth: the pfamily must carry exactly the calibration
  ## stored on the setup object by Prior_Setup_lmebayes() (ing_prior field).
  ing_k <- ps$ing_prior[[k]]
  stopifnot(
    !is.null(ing_k),
    isTRUE(all.equal(pr$shape,      ing_k$shape)),
    isTRUE(all.equal(pr$rate,       ing_k$rate)),
    isTRUE(all.equal(pr$disp_lower, ing_k$disp_lower)),
    isTRUE(all.equal(pr$disp_upper, ing_k$disp_upper))
  )
}
cat("dIndependent_Normal_Gamma path: OK\n")

## --- 3. mixed ptypes: positional vector -------------------------------------
mix <- c("dNormal", "dIndependent_Normal_Gamma", "dNormal")
pf_mix <- pfamily_list(ps, ptypes = mix)
stopifnot(identical(
  unname(vapply(pf_mix, function(p) p$pfamily, character(1L))),
  mix
))
cat("Positional mixed ptypes: OK\n")

## --- 4. named vector in scrambled order -------------------------------------
named <- stats::setNames(
  c("dIndependent_Normal_Gamma", "dNormal", "dIndependent_Normal_Gamma"),
  re_names[c(3L, 1L, 2L)]
)
pf_named <- pfamily_list(ps, ptypes = named)
stopifnot(
  identical(names(pf_named), re_names),
  all(vapply(re_names, function(k) {
    pf_named[[k]]$pfamily == named[[k]]
  }, logical(1L)))
)
cat("Named (scrambled) ptypes: OK\n")

## --- 5. list input -----------------------------------------------------------
pf_list <- pfamily_list(
  ps,
  ptypes = list("dNormal", "dNormal", "dIndependent_Normal_Gamma")
)
stopifnot(identical(
  unname(vapply(pf_list, function(p) p$pfamily, character(1L))),
  c("dNormal", "dNormal", "dIndependent_Normal_Gamma")
))
cat("List ptypes: OK\n")

## --- 6. validation errors ----------------------------------------------------
expect_error <- function(expr, pattern) {
  err <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  if (is.null(err)) stop("Expected an error matching: ", pattern)
  if (!grepl(pattern, err)) {
    stop("Error message ", sQuote(err), " does not match ", sQuote(pattern))
  }
  invisible(TRUE)
}

expect_error(pfamily_list(ps, ptypes = "dBanana"),
             "Invalid 'ptypes'")
expect_error(pfamily_list(ps, ptypes = c("dNormal", "dNormal")),
             "length 2 but the prior setup has 3")
expect_error(
  pfamily_list(ps, ptypes = stats::setNames(
    rep("dNormal", 3L), c("a", "b", "c"))),
  "Names of 'ptypes'"
)
expect_error(pfamily_list(ps, ptypes = 1L),
             "character vector or list")
expect_error(pfamily_list(ps, ptypes = list("dNormal", 2, "dNormal")),
             "single string")

## Prior-vs-data guard: pwt_dispersion > 0.5 implies n_prior_dispersion > J,
## which the ING dispersion envelope cannot support (sampler caps the
## log-tilt at J/2).  Building such a pfamily must fail early.
ps_heavy <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01,
                                 pwt_dispersion = 0.9)
expect_error(
  pfamily_list(ps_heavy, ptypes = "dIndependent_Normal_Gamma"),
  "n_prior_dispersion <= J"
)
## dNormal is unaffected by the dispersion-prior guard.
pf_heavy_norm <- pfamily_list(ps_heavy)
stopifnot(length(pf_heavy_norm) == 3L)
cat("Validation errors: OK\n")

cat("\nAll pfamily_list() tests passed.\n")
