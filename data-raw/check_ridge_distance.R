## Check the in-package ridge distance report against the reweighting script.
##
## Fits a short Poisson ridge run with the installed lmebayesCore, captures
## the "[ridge] distance to the rwt = 0 posterior" message, and recomputes
## TV, ESS and Pareto k from the fit's group draws with ridge_gap() and
## reweight() from ridge_reweight_to_rwt0.R. Same h on the same draws, so the
## two should agree to rounding. Uses the installed package (no compilation).
##
## Run from C:/Rpackages/lmebayes:
##   Rscript data-raw/check_ridge_distance.R          (rwt = 0.25, 500 chains)
##   Rscript data-raw/check_ridge_distance.R 0.75 1000

suppressPackageStartupMessages(library(lmebayesCore))
args <- commandArgs(trailingOnly = TRUE)
rwt <- if (length(args) >= 1L) as.numeric(args[1L]) else 0.25
n_chains <- if (length(args) >= 2L) as.integer(args[2L]) else 500L

data(airbnb_small, package = "bayesrules")
dat <- airbnb_small
dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
dat$walk_c <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
dat <- dat[stats::complete.cases(dat[, c(
  "reviews", "rating_c", "walk_c", "neighborhood"
)]), ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))
form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

design <- model_setup(form, data = dat, family = poisson())
ps <- Prior_Setup_GLMM(form, data = dat, family = poisson(), pop.pwt = 0.01)
pf <- pfamily_list(ps)

ridge_msg <- NA_character_
fit <- withCallingHandlers(
  rglmerb(n = n_chains, design = design, pfamily_list = pf,
          family = poisson(), rwt = rwt, verbose = FALSE, progbar = TRUE),
  message = function(m) {
    if (grepl("^\\[ridge\\]", conditionMessage(m))) {
      ridge_msg <<- trimws(conditionMessage(m))
    }
  }
)

## ridge_gap(), group_draws() and reweight() from the reweighting script,
## with the globals they read
want <- c("group_draws", "ridge_gap", "reweight")
for (x in parse("data-raw/ridge_reweight_to_rwt0.R")) {
  if (is.call(x) && identical(x[[1]], as.name("<-")) &&
      as.character(x[[2]]) %in% want) eval(x)
}
re_cols <- c("(Intercept)", "rating_c")
y_all <- dat$reviews
X_all <- cbind(1, dat$rating_c)
g_all <- as.character(dat$neighborhood)

h <- ridge_gap(fit)
rw <- reweight(h, rwt, 0)
k_hat <- if (requireNamespace("loo", quietly = TRUE)) {
  suppressWarnings(loo::psis(rwt * h, r_eff = 1)$diagnostics$pareto_k)
} else NA_real_

cat("\n===== ridge distance: package message vs reweighting script =====\n")
cat("package: ", ridge_msg, "\n", sep = "")
cat(sprintf("script:  TV ~ %.3g (ESS %.0f of %d, Pareto k %.2f)\n",
            rw$tv, rw$ess, length(h), k_hat))
