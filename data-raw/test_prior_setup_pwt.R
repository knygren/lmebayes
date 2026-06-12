# Regression test: flexible pwt and decoupled dispersion-prior arguments in
# Prior_Setup_lmebayes().
#
# Checks:
#   - scalar pwt: unchanged behavior; $pwt stays scalar; n_prior_dispersion
#     and pwt_dispersion derived per component (n_k = J*w/(1-w)).
#   - list pwt (per-component scalars, named/scrambled and positional):
#     Sigma_fixef scaled per component by (1-w_k)/w_k.
#   - vector pwt within a component (named, scrambled): elementwise
#     sqrt(s_i)*sqrt(s_j) scaling.
#   - pwt_dispersion / n_prior_dispersion arguments: consistency
#     w_k = n_k/(n_k + J); mutual exclusivity; source attribute.
#   - pfamily_list() ING calibration uses object$n_prior_dispersion.
#   - validation errors for malformed inputs.
#
#   Rscript data-raw/test_prior_setup_pwt.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

expect_error <- function(expr, pattern) {
  err <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  if (is.null(err)) stop("Expected an error matching: ", pattern)
  if (!grepl(pattern, err)) {
    stop("Error message ", sQuote(err), " does not match ", sQuote(pattern))
  }
  invisible(TRUE)
}

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

## --- 1. scalar pwt baseline --------------------------------------------------
w0 <- 0.01
ps0 <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = w0)
re_names <- names(ps0$prior_list)
J <- nlevels(ps0$design$groups)
stopifnot(length(re_names) == 3L)

stopifnot(
  is.numeric(ps0$pwt), length(ps0$pwt) == 1L, ps0$pwt == w0,
  is.numeric(ps0$n_prior_dispersion),
  identical(names(ps0$n_prior_dispersion), re_names),
  is.numeric(ps0$pwt_dispersion),
  identical(names(ps0$pwt_dispersion), re_names),
  identical(attr(ps0$n_prior_dispersion, "source"), "derived from pwt")
)
n_exp <- J * w0 / (1 - w0)
stopifnot(
  isTRUE(all.equal(as.vector(ps0$n_prior_dispersion), rep(n_exp, 3L))),
  isTRUE(all.equal(as.vector(ps0$pwt_dispersion), rep(w0, 3L))),
  ## consistency identity w = n / (n + J)
  isTRUE(all.equal(
    as.vector(ps0$pwt_dispersion),
    as.vector(ps0$n_prior_dispersion / (ps0$n_prior_dispersion + J))
  ))
)
cat("Scalar pwt baseline: OK\n")

## --- 2. list pwt of all-equal scalars reproduces scalar result ---------------
ps_eq <- Prior_Setup_lmebayes(
  form_lmer, data = dat,
  pwt = stats::setNames(as.list(rep(w0, 3L)), re_names)
)
stopifnot(
  is.list(ps_eq$pwt),
  identical(names(ps_eq$pwt), re_names)
)
for (k in re_names) {
  stopifnot(
    isTRUE(all.equal(
      ps_eq$prior_list[[k]]$Sigma_fixef,
      ps0$prior_list[[k]]$Sigma_fixef
    )),
    isTRUE(all.equal(
      ps_eq$prior_list[[k]]$mu_fixef,
      ps0$prior_list[[k]]$mu_fixef
    ))
  )
}
stopifnot(isTRUE(all.equal(
  as.vector(ps_eq$n_prior_dispersion), as.vector(ps0$n_prior_dispersion)
)))
cat("List pwt (all equal) == scalar pwt: OK\n")

## --- 3. per-component scalar weights (named, scrambled order) ----------------
w_by_comp <- stats::setNames(list(0.05, 0.01, 0.2), re_names[c(2L, 1L, 3L)])
ps_c <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = w_by_comp)
stopifnot(identical(names(ps_c$pwt), re_names))
s0 <- (1 - w0) / w0
for (k in re_names) {
  w_k <- w_by_comp[[k]]
  stopifnot(all(ps_c$pwt[[k]] == w_k))
  s_k <- (1 - w_k) / w_k
  stopifnot(isTRUE(all.equal(
    ps_c$prior_list[[k]]$Sigma_fixef,
    ps0$prior_list[[k]]$Sigma_fixef * (s_k / s0)
  )))
}
## dispersion prior derived from per-component mean weight
stopifnot(isTRUE(all.equal(
  as.vector(ps_c$n_prior_dispersion),
  vapply(re_names, function(k) {
    w_k <- w_by_comp[[k]]
    J * w_k / (1 - w_k)
  }, numeric(1L), USE.NAMES = FALSE)
)))
cat("Per-component scalar pwt: OK\n")

## --- 4. per-predictor vector weights within a component ----------------------
k1     <- re_names[1L]
cols_1 <- colnames(ps0$design$X_hyper[[k1]])
p_1    <- length(cols_1)
w_vec  <- stats::setNames(seq(0.02, 0.3, length.out = p_1), cols_1)
w_vec_scrambled <- w_vec[rev(seq_len(p_1))]

pwt_mixed <- stats::setNames(vector("list", 3L), re_names)
pwt_mixed[[k1]] <- w_vec_scrambled
for (k in re_names[-1L]) pwt_mixed[[k]] <- w0

ps_v <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = pwt_mixed)
stopifnot(
  identical(names(ps_v$pwt[[k1]]), cols_1),
  isTRUE(all.equal(ps_v$pwt[[k1]], w_vec))
)
s_vec <- sqrt((1 - w_vec) / w_vec)
Sigma_exp <- (ps0$prior_list[[k1]]$Sigma_fixef / s0) * outer(s_vec, s_vec)
dimnames(Sigma_exp) <- dimnames(ps0$prior_list[[k1]]$Sigma_fixef)
stopifnot(isTRUE(all.equal(ps_v$prior_list[[k1]]$Sigma_fixef, Sigma_exp)))
## other components untouched
for (k in re_names[-1L]) {
  stopifnot(isTRUE(all.equal(
    ps_v$prior_list[[k]]$Sigma_fixef, ps0$prior_list[[k]]$Sigma_fixef
  )))
}
## n_prior_dispersion for the vector component uses the mean weight
stopifnot(isTRUE(all.equal(
  unname(ps_v$n_prior_dispersion[[k1]]),
  J * mean(w_vec) / (1 - mean(w_vec))
)))
cat("Per-predictor vector pwt: OK\n")

## --- 5. pwt_dispersion argument ----------------------------------------------
wd <- 0.5
ps_wd <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = w0,
                              pwt_dispersion = wd)
stopifnot(
  isTRUE(all.equal(as.vector(ps_wd$pwt_dispersion), rep(wd, 3L))),
  isTRUE(all.equal(as.vector(ps_wd$n_prior_dispersion),
                   rep(J * wd / (1 - wd), 3L))),
  identical(attr(ps_wd$n_prior_dispersion, "source"),
            "user-supplied (pwt_dispersion)"),
  ## coefficient priors unaffected
  isTRUE(all.equal(ps_wd$prior_list, ps0$prior_list))
)

## per-component list, named and scrambled
wd_list <- stats::setNames(list(0.3, 0.5, 0.7), re_names[c(3L, 1L, 2L)])
ps_wd2 <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = w0,
                               pwt_dispersion = wd_list)
for (k in re_names) {
  w_k <- wd_list[[k]]
  stopifnot(
    isTRUE(all.equal(unname(ps_wd2$pwt_dispersion[[k]]), w_k)),
    isTRUE(all.equal(unname(ps_wd2$n_prior_dispersion[[k]]),
                     J * w_k / (1 - w_k)))
  )
}
cat("pwt_dispersion argument: OK\n")

## --- 6. n_prior_dispersion argument ------------------------------------------
nd <- 10
ps_nd <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = w0,
                              n_prior_dispersion = nd)
stopifnot(
  isTRUE(all.equal(as.vector(ps_nd$n_prior_dispersion), rep(nd, 3L))),
  isTRUE(all.equal(as.vector(ps_nd$pwt_dispersion), rep(nd / (nd + J), 3L))),
  identical(attr(ps_nd$n_prior_dispersion, "source"),
            "user-supplied (n_prior_dispersion)")
)

## per-component numeric vector, named and scrambled
nd_vec <- stats::setNames(c(5, 10, 20), re_names[c(2L, 3L, 1L)])
ps_nd2 <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = w0,
                               n_prior_dispersion = nd_vec)
for (k in re_names) {
  n_k <- unname(nd_vec[[k]])
  stopifnot(
    isTRUE(all.equal(unname(ps_nd2$n_prior_dispersion[[k]]), n_k)),
    isTRUE(all.equal(unname(ps_nd2$pwt_dispersion[[k]]), n_k / (n_k + J)))
  )
}
cat("n_prior_dispersion argument: OK\n")

## --- 7. pfamily_list ING uses n_prior_dispersion -----------------------------
pf_ing <- pfamily_list(ps_nd2, ptypes = "dIndependent_Normal_Gamma")
for (k in re_names) {
  pl  <- ps_nd2$prior_list[[k]]
  pr  <- pf_ing[[k]]$prior_list
  n_k <- unname(ps_nd2$n_prior_dispersion[[k]])
  p_k <- length(pl$mu_fixef)
  stopifnot(
    isTRUE(all.equal(pr$shape, (n_k + 1) / 2 + p_k / 2)),
    isTRUE(all.equal(pr$rate, unname(pl$dispersion_fixef) * (n_k / 2)))
  )
}
cat("pfamily_list ING uses n_prior_dispersion: OK\n")

## --- 8. print methods run without error --------------------------------------
invisible(capture.output(print(ps0)))
invisible(capture.output(print(ps_v)))
invisible(capture.output(print(ps_nd2)))
cat("print methods: OK\n")

## --- 9. validation errors -----------------------------------------------------
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 1.5),
  "'pwt' must be a scalar in \\(0, 1\\)"
)
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat, pwt = list(0.1, 0.2)),
  "length 2 but there are 3"
)
expect_error(
  Prior_Setup_lmebayes(
    form_lmer, data = dat,
    pwt = stats::setNames(list(0.1, 0.2, 0.3), c("a", "b", "c"))
  ),
  "Names of 'pwt'"
)
expect_error(
  Prior_Setup_lmebayes(
    form_lmer, data = dat,
    pwt = stats::setNames(as.list(c(0.1, 1.2, 0.3)), re_names)
  ),
  "must be numeric with all values in \\(0, 1\\)"
)
bad_len <- stats::setNames(as.list(rep(w0, 3L)), re_names)
bad_len[[k1]] <- rep(w0, p_1 + 1L)
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat, pwt = bad_len),
  sprintf("length 1 or %d", p_1)
)
bad_nms <- stats::setNames(as.list(rep(w0, 3L)), re_names)
bad_nms[[k1]] <- stats::setNames(rep(w0, p_1), paste0("z", seq_len(p_1)))
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat, pwt = bad_nms),
  "must match the Block 2 predictors"
)
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat,
                       pwt_dispersion = 0.5, n_prior_dispersion = 10),
  "at most one of"
)
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat, pwt_dispersion = 1.1),
  "must be in \\(0, 1\\)"
)
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat, n_prior_dispersion = -2),
  "positive and finite"
)
expect_error(
  Prior_Setup_lmebayes(form_lmer, data = dat,
                       n_prior_dispersion = c(1, 2)),
  "length 1 or 3"
)
expect_error(
  Prior_Setup_lmebayes(
    form_lmer, data = dat,
    n_prior_dispersion = stats::setNames(c(1, 2, 3), c("a", "b", "c"))
  ),
  "must match the random-effect coefficient names"
)
cat("Validation errors: OK\n")

cat("\nAll Prior_Setup_lmebayes pwt/dispersion tests passed.\n")
