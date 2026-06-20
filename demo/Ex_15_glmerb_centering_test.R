## Demo: glmerb Poisson centering test on airbnb_small subset
##
## Replica of tests/testthat/test-glmerb-example.R (removed from testthat because
## the full pilot + main run takes several minutes).  Exercises the two-stage
## sampler (pilot from ICM mode, main from pilot mean) and prints multivariate
## centering diagnostics: posterior mean vs pilot mean vs mode.
##
##   demo("Ex_15_glmerb_centering_test", package = "lmebayes")
##
## WHY A SUBSET:
##   The full airbnb_small has 17 neighborhood levels (J=17).  By a CLT-like
##   argument over J groups, the marginal posterior of the fixed-effect
##   hyperparameter gamma is approximately normal for large J, making the mode
##   nearly equal to the mean and rendering the pilot-vs-mode comparison
##   insensitive.  The commonly cited threshold is J ~ 30; at J=17 we are
##   already close to the normal regime.
##
##   To make the mode-vs-mean gap detectable we keep only the smallest
##   neighborhoods (5-20 observations, J=6, 72 rows total).
##   Two criteria amplify Poisson non-normality:
##     1. Few groups (J small) -> posterior of gamma far from normal.
##     2. Small counts per observation -> individual Poisson posteriors skewed;
##        ~43% of the review counts in these groups are single-digit.
##   Larger groups (>20 obs) and the very large ones (Logan Square n=330,
##   Rogers Park n=123) are excluded because for large n_j the Poisson
##   likelihood is well-approximated by a Gaussian, reducing skewness.
##   For Poisson (concave h): ICM mode < E[gamma|y] < gamma* (Banach fixed pt).

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

.print_glmerb_test_table <- function(x, digits = 6L) {
  op <- options(width = max(300L, getOption("width")))
  on.exit(options(op), add = TRUE)
  paste(
    capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE)),
    collapse = "\n"
  )
}

.print_glmerb_test_summary <- function(fit) {
  op <- options(width = max(300L, getOption("width")))
  on.exit(options(op), add = TRUE)
  paste(capture.output(print(summary(fit))), collapse = "\n")
}

data(airbnb_small, package = "bayesrules", envir = environment())
dat <- airbnb_small
dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
dat$walk_c   <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
dat <- dat[complete.cases(dat[, c("reviews", "rating_c", "walk_c",
                                  "neighborhood")]), ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))

grp_counts  <- table(dat$neighborhood)
keep_groups <- names(grp_counts[grp_counts >= 5L & grp_counts <= 20L])
dat <- dat[dat$neighborhood %in% keep_groups, ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))
stopifnot(nlevels(dat$neighborhood) == 6L)

form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)

fit <- glmerb(
  form,
  data         = dat,
  family       = poisson(),
  pfamily_list = pfamily_list(ps),
  n            = 10000L,
  gap_tol      = 0.0196,
  mode_gap_max = 1.0,
  seed         = 42L
)

re_names <- fit$model_setup$re_coef_names
stopifnot(identical(re_names, c("(Intercept)", "rating_c")))

n_draws <- nrow(fit$fixef[[re_names[1L]]])
stopifnot(identical(n_draws, 10000L))
stopifnot(!is.null(fit$fixef.init))
stopifnot(is.list(fit$convergence), is.finite(fit$convergence$m_convergence))
stopifnot(identical(fit$pilot_chisq$n_pilot, 10000L))
stopifnot(identical(fit$convergence$m_convergence_pilot, 9L))
stopifnot(identical(fit$convergence$mode_gap_max, 1.0))

X <- do.call(cbind, lapply(re_names, function(k) fit$fixef[[k]]))
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
colnames(X) <- cn
stopifnot(all(is.finite(X)))

beta_bar <- colMeans(X)
theta_pilot <- unlist(lapply(re_names, function(k) fit$fixef.init[[k]]))
theta_mode  <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
names(theta_pilot) <- cn
names(theta_mode)  <- cn

center_tab <- data.frame(
  parameter = cn,
  mode = unname(theta_mode),
  pilot_mean = unname(theta_pilot),
  main_mean = unname(beta_bar),
  stringsAsFactors = FALSE
)
rownames(center_tab) <- NULL
diff_tab <- data.frame(
  parameter = cn,
  pilot_minus_mode = unname(theta_pilot - theta_mode),
  main_minus_pilot = unname(beta_bar - theta_pilot),
  main_minus_mode = unname(beta_bar - theta_mode),
  stringsAsFactors = FALSE
)
rownames(diff_tab) <- NULL

sd_main <- apply(X, 2L, stats::sd)
se_main <- sd_main / sqrt(n_draws)
z_main_vs_pilot <- unname((beta_bar - theta_pilot) / se_main)
z_main_vs_mode  <- unname((beta_bar - theta_mode) / se_main)
uni_tab <- data.frame(
  parameter = cn,
  z_vs_pilot = z_main_vs_pilot,
  p_vs_pilot = 2 * stats::pnorm(abs(z_main_vs_pilot), lower.tail = FALSE),
  z_vs_mode = z_main_vs_mode,
  p_vs_mode = 2 * stats::pnorm(abs(z_main_vs_mode), lower.tail = FALSE),
  stringsAsFactors = FALSE
)
rownames(uni_tab) <- NULL

cat("glmerb centers table:\n")
cat(.print_glmerb_test_table(center_tab), "\n\n")
cat("glmerb differences table:\n")
cat(.print_glmerb_test_table(diff_tab), "\n\n")
cat("glmerb univariate z/p table:\n")
cat(.print_glmerb_test_table(uni_tab), "\n\n")

n_tot <- nrow(X)
p_tot <- ncol(X)
S <- stats::cov(X)
V <- S / n_tot
V_inv <- solve(V)

d_pilot <- beta_bar - theta_pilot
d_mode  <- beta_bar - theta_mode

Q_pilot <- as.numeric(t(d_pilot) %*% V_inv %*% d_pilot)
Q_mode  <- as.numeric(t(d_mode)  %*% V_inv %*% d_mode)
p_pilot <- stats::pchisq(Q_pilot, df = p_tot, lower.tail = FALSE)
p_mode  <- stats::pchisq(Q_mode,  df = p_tot, lower.tail = FALSE)

cat(sprintf(
  "glmerb centering test: m_convergence=%d, p(mean=pilot)=%.4g, p(mean=mode)=%.4g\n\n",
  fit$convergence$m_convergence, p_pilot, p_mode
))

stopifnot(is.finite(p_pilot), is.finite(p_mode))
if (!((p_pilot > 0.05) || (p_pilot >= p_mode))) {
  warning(sprintf(
    "Centering criterion not met: p_pilot=%.4g, p_mode=%.4g",
    p_pilot, p_mode
  ), call. = FALSE)
}

lmebayes:::print_coef_means(fit)
print(fit)
cat(.print_glmerb_test_summary(fit), "\n")
