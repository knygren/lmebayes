# Manual validation: glmerb Poisson on airbnb_small -- ridge distance, DIC,
# and per-group summary.
#
# Model: reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood).
# Sampler: rwt set explicitly below (edit to test other values); otherwise
# glmerb() defaults (tv_tol, gap_tol, mode_gap_max).
# No set.seed(). Main draws n = 4000 (optional trailing arg sets n).
# progbar = TRUE only.
#
# Checks:
#   - the "[ridge] distance to the rwt = 0 posterior" message is emitted and
#     its Pareto k (when loo is installed) is below 0.7
#   - stored Dbar, Dthetabar, pD, DIC match a dpois() recomputation from
#     fit$groupef, and p_re < pD < J * p_re
#   - extractDIC() returns c(pD, DIC); summary() carries and prints the DIC
#
#   Rscript tests/manual/test_glmerb_airbnb_small_dic.R
#   Rscript tests/manual/test_glmerb_airbnb_small_dic.R 1000

args <- commandArgs(trailingOnly = TRUE)
n_draws <- if (length(args) >= 1L) as.integer(args[1L]) else 4000L
rwt <- 0.01

source("tests/manual/_load.R")
.manual_test_load()

data("airbnb_small", package = "bayesrules")
dat <- airbnb_small
dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
dat$walk_c   <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
dat <- dat[stats::complete.cases(dat[, c(
  "reviews", "rating_c", "walk_c", "neighborhood"
)]), ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))

form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)
ps <- Prior_Setup_GLMM(form, data = dat, family = poisson(), pop.pwt = 0.01)

cat("\n=== glmerb Poisson airbnb_small; n = ", n_draws, ", rwt = ", rwt,
    " ===\n\n", sep = "")
ridge_msg <- NA_character_
fit <- withCallingHandlers(
  glmerb(form, data = dat, family = poisson(),
         pfamily_list = pfamily_list(ps), n = n_draws, rwt = rwt,
         progbar = TRUE),
  message = function(m) {
    if (grepl("^\\[ridge\\]", conditionMessage(m))) {
      ridge_msg <<- trimws(conditionMessage(m))
    }
  }
)
stopifnot(inherits(fit, "glmerb"))

## --- ridge distance ---------------------------------------------------------

stopifnot(!is.na(ridge_msg))
cat("ridge: ", ridge_msg, "\n", sep = "")
k_txt <- regmatches(ridge_msg, regexpr("Pareto k [0-9.]+", ridge_msg))
if (length(k_txt)) {
  stopifnot(as.numeric(sub("Pareto k ", "", k_txt)) < 0.7)
}

## --- DIC vs dpois() recomputation -------------------------------------------

ms <- fit$model_setup
re <- ms$groupef.names
D  <- as.matrix(ms$D)[, re, drop = FALSE]
g  <- as.character(ms$group)
d  <- fit$groupef
Dbar <- 0
Dt   <- 0
for (l in unique(g)) {
  r  <- which(g == l)
  Dr <- D[r, , drop = FALSE]
  B  <- as.matrix(d[as.character(d[[ms$group_name]]) == l, re])
  ll <- apply(B, 1, function(b) sum(dpois(ms$y[r], exp(Dr %*% b), log = TRUE)))
  Dbar <- Dbar - 2 * mean(ll)
  Dt   <- Dt - 2 * sum(dpois(ms$y[r], exp(Dr %*% colMeans(B)), log = TRUE))
}
hand   <- c(Dbar = Dbar, Dthetabar = Dt, pD = Dbar - Dt, DIC = 2 * Dbar - Dt)
stored <- c(Dbar = fit$Dbar, Dthetabar = fit$Dthetabar, pD = fit$pD, DIC = fit$DIC)
stopifnot(length(stored) == 4L, all(is.finite(stored)))

err  <- max(abs(stored - hand))
J    <- length(unique(g))
p_re <- length(re)
cat(sprintf("DIC %.2f  pD %.2f  (between %d and %d)  max|stored - hand| %.2e\n",
            stored[["DIC"]], stored[["pD"]], p_re, J * p_re, err))
stopifnot(err < 1e-6 * max(1, abs(stored[["DIC"]])))
stopifnot(stored[["pD"]] > p_re, stored[["pD"]] < J * p_re)

e <- glmbayes::extractDIC(fit)
stopifnot(identical(names(e), c("pD", "DIC")),
          isTRUE(all.equal(unname(e), unname(stored[c("pD", "DIC")]))))

## --- summary ----------------------------------------------------------------

s <- summary(fit, groups = "Jefferson Park")
stopifnot(isTRUE(all.equal(unname(s$dic), unname(stored[c("Dbar", "pD", "DIC")]))))
out <- utils::capture.output(print(s))
stopifnot(any(grepl("^DIC:", out)))
print(s)

cat("\ntest_glmerb_airbnb_small_dic: OK\n")
