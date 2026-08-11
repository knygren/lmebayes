# Manual check: fixed per-group dispersion vector for dispersion_ranef
# (dispersion_mode = "fixed_vector"), on lmerb() and glmerb(family = gaussian()).
#
# Covers, using lme4::sleepstudy (no bayesrules dependency):
#   1. lmerb() with a named numeric dispersion_ranef vector + dispformula = ~Subject:
#      dispersion_mode dispatch, per-group sigma2 attachment (matches the input,
#      reordered to group_levels), no glmmTMB reference fit, summary_sigma2() output.
#   2. Validation errors: wrong length, mismatched group names, dispformula = ~1
#      rejecting a per-group vector.
#   3. glmerb(family = gaussian()) parity with lmerb() (shared code path).
#   4. Light sanity check that Block~2 fixef posterior means stay close to the
#      lmer() reference fit, to catch any accidental scalar-dispersion assumption
#      left over in the rate-calibration / ICM code paths for this mode.
#
#   cd .../lmebayes
#   Rscript tests/manual/test_lmerb_fixed_vector_dispersion.R

source("tests/manual/_load.R")
.manual_test_load(require_bayesrules = FALSE)

dat <- lme4::sleepstudy
dat$Subject <- factor(dat$Subject)
form <- Reaction ~ Days + (Days || Subject)

grp_levels <- levels(dat$Subject)
J <- length(grp_levels)

ps <- Prior_Setup_GLMM(form, data = dat, pwt = 0.01)
pooled_sigma2 <- as.numeric(ps$dispersion_ranef)

set.seed(42L)
disp_vec <- stats::setNames(
  pooled_sigma2 * stats::runif(J, 0.5, 1.5),
  sample(grp_levels)
)

cat("\n=== 1. lmerb() with a fixed per-group dispersion_ranef vector ===\n\n")
cat("Input dispersion_ranef (unordered names):\n")
print(disp_vec)

fit <- lmerb(
  form,
  data             = dat,
  pfamily_list     = pfamily_list(ps),
  dispersion_ranef = disp_vec,
  dispformula      = ~Subject,
  n                = 500L,
  progbar          = FALSE
)

stopifnot(identical(fit$prior$dispersion_mode, "fixed_vector"))
stopifnot(is.null(fit$dispersion_fit))
stopifnot(identical(names(fit$group.dispersion), grp_levels))
stopifnot(isTRUE(all.equal(as.numeric(fit$group.dispersion), as.numeric(disp_vec[grp_levels]))))
stopifnot(identical(fit$group.dispersion, fit$group.dispersion.mean))

cat("\nfit$group.dispersion (reordered to group_levels, matches input exactly): OK\n")
cat("fit$prior$dispersion_mode:", fit$prior$dispersion_mode, "\n")
cat("fit$dispersion_fit (should be NULL, no glmmTMB fit needed):",
    if (is.null(fit$dispersion_fit)) "NULL" else "NOT NULL -- FAIL", "\n")

cat("\n--- summary_sigma2(fit) ---\n\n")
summ <- summary_sigma2(fit)
print(summ)
stopifnot(identical(summ$mode, "fixed_vector"))
stopifnot(is.null(summ$overview))
stopifnot(is.null(summ$percentiles))
stopifnot(isTRUE(all.equal(
  summ$prior[grp_levels, "Fixed sigma2"],
  as.numeric(disp_vec[grp_levels])
)))

cat("\n--- summary(fit) (pooled Residual row should be absent) ---\n\n")
s <- summary(fit)
stopifnot(!"Residual" %in% rownames(s$tau2_overview))
print(s)

cat("\n=== 2. Validation errors ===\n\n")

err <- tryCatch({
  lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = stats::setNames(rep(1, 3), grp_levels[1:3]),
    dispformula      = ~Subject,
    n                = 2L,
    progbar          = FALSE
  )
  NULL
}, error = function(e) conditionMessage(e))
cat("Wrong length ->", err, "\n")
stopifnot(!is.null(err), grepl("length", err))

bad_vec <- stats::setNames(rep(1, J), grp_levels)
names(bad_vec)[1] <- "not_a_subject"
err <- tryCatch({
  lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = bad_vec,
    dispformula      = ~Subject,
    n                = 2L,
    progbar          = FALSE
  )
  NULL
}, error = function(e) conditionMessage(e))
cat("Mismatched names ->", err, "\n")
stopifnot(!is.null(err), grepl("group levels", err))

good_vec <- stats::setNames(rep(1, J), grp_levels)
err <- tryCatch({
  lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = good_vec,
    dispformula      = ~1,
    n                = 2L,
    progbar          = FALSE
  )
  NULL
}, error = function(e) conditionMessage(e))
cat("dispformula = ~1 with per-group vector ->", err, "\n")
stopifnot(!is.null(err), grepl("dispformula", err))

cat("\n=== 3. glmerb(family = gaussian()) parity with lmerb() ===\n\n")

fit_glm <- glmerb(
  form,
  data             = dat,
  family           = gaussian(),
  pfamily_list     = pfamily_list(ps),
  dispersion_ranef = disp_vec,
  dispformula      = ~Subject,
  n                = 500L,
  progbar          = FALSE
)

stopifnot(identical(fit_glm$prior$dispersion_mode, "fixed_vector"))
stopifnot(is.null(fit_glm$dispersion_fit))
stopifnot(identical(names(fit_glm$group.dispersion), grp_levels))
stopifnot(isTRUE(all.equal(as.numeric(fit_glm$group.dispersion), as.numeric(disp_vec[grp_levels]))))

cat("glmerb(family = gaussian()) matches lmerb(): OK\n")

cat("\n=== 4. Sanity check: Block~2 fixef posterior mean vs lmer() reference ===\n\n")

lmer_fit <- lme4::lmer(form, data = dat, REML = TRUE)
lmer_fe  <- lme4::fixef(lmer_fit)

post_mean_intercept <- mean(fit$fixef.means[["(Intercept)"]])
post_mean_days       <- mean(fit$fixef.means[["Days"]])

cat(sprintf(
  paste0(
    "  (Intercept): lmer = %.4f, lmerb posterior mean = %.4f\n",
    "  Days       : lmer = %.4f, lmerb posterior mean = %.4f\n"
  ),
  lmer_fe[["(Intercept)"]], post_mean_intercept,
  lmer_fe[["Days"]], post_mean_days
))

## Loose sanity bands only -- per-group dispersion changes the effective
## weighting relative to lmer()'s pooled residual variance, so this is not
## expected to match tightly; it is here to catch gross regressions (e.g. a
## stray scalar-dispersion assumption silently using only the first group's
## value everywhere).
stopifnot(abs(post_mean_intercept - lmer_fe[["(Intercept)"]]) < 20)
stopifnot(abs(post_mean_days - lmer_fe[["Days"]]) < 5)

cat("\ntest_lmerb_fixed_vector_dispersion: OK\n")
