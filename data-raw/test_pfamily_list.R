# Regression test: pfamily_list() method for lmebayes_prior_setup objects.
#
# Builds named lists of glmbayesCore pfamily objects (one per random-effect
# coefficient) from Prior_Setup_lmebayes() output.  Checks:
#   - dNormal path: mu/Sigma/dispersion match prior_list fields exactly.
#   - dIndependent_Normal_Gamma path: mu/Sigma match; shape/rate follow the
#     shape_ING convention with n0 = J * pwt/(1-pwt) so the Gamma prior mean
#     of the dispersion is ~ tau^2_k.
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
n_prior <- (ps$pwt / (1 - ps$pwt)) * J

## The setup object now carries per-component n_prior_dispersion; with a
## scalar pwt it must equal the classic derivation used below.
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
  rate_exp  <- unname(pl$dispersion_fixef) * (n_prior / 2)
  ## Default disp_lower: 0.01 quantile of the inverse-Gamma dispersion prior
  ## (reciprocal of the 99th percentile of the Gamma precision prior).
  disp_lower_exp <- 1 / stats::qgamma(0.99, shape = shape_exp, rate = rate_exp)
  stopifnot(
    isTRUE(all.equal(as.numeric(pr$mu), unname(pl$mu_fixef))),
    isTRUE(all.equal(unname(as.matrix(pr$Sigma)), unname(pl$Sigma_fixef))),
    isTRUE(all.equal(pr$shape, shape_exp)),
    isTRUE(all.equal(pr$rate, rate_exp)),
    isTRUE(all.equal(pr$disp_lower, disp_lower_exp)),
    pr$disp_lower > 0
  )
  ## Prior mean of the precision: shape/rate = (n0 + 1 + p_k) / (n0 * tau^2_k).
  ## For weak priors (small n0) this is inflated relative to 1/tau^2_k by the
  ## factor (n0 + 1 + p_k)/n0; verify the identity rather than closeness.
  prec_mean_exp <- (n_prior + 1 + p_k) / (n_prior * unname(pl$dispersion_fixef))
  stopifnot(isTRUE(all.equal(pr$shape / pr$rate, prec_mean_exp)))
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
cat("Validation errors: OK\n")

cat("\nAll pfamily_list() tests passed.\n")
