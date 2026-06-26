## Demo: glmerb Poisson centering test on bayesrules::airbnb (full model)
##
## Replica of demo/Ex_15_glmerb_centering_test.R on the full airbnb dataset
## (not run in testthat because the pilot + main fit takes several minutes).  Exercises cost-optimal pilot sizing, the two-stage sampler
## (pilot from ICM mode, main from pilot mean), and multivariate centering
## diagnostics: posterior mean vs pilot mean vs mode.
##
##   demo("Ex_13_glmerb_Airbnb", package = "lmebayes")
##
## Model: listing-level rating and log-price with neighborhood random slopes;
## walkability and transit as level-2 covariates on the intercept; cross-level
## moderation via walk_c:rating_c and transit_c:log_price_c (same structure as
## demo/Ex_12_lmerb_BigWordClub.R on the Gaussian side).

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

data(airbnb, package = "bayesrules", envir = environment())

dat <- airbnb
dat$rating_c    <- dat$rating - mean(dat$rating)
dat$log_price_c <- scale(log(dat$price + 1))[, 1]
dat$walk_c      <- dat$walk_score - mean(dat$walk_score)
dat$transit_c   <- dat$transit_score - mean(dat$transit_score)
dat <- dat[complete.cases(dat[, c(
  "reviews", "rating", "rating_c", "price", "log_price_c",
  "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
)]), ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))

form_glmer <- reviews ~
  walk_c + transit_c +
  rating_c + log_price_c +
  walk_c:rating_c + transit_c:log_price_c +
  (1 + rating_c + log_price_c || neighborhood)

ps <- Prior_Setup_lmebayes(form_glmer, data = dat, family = poisson(), pwt = 0.01)

fit <- glmerb(
  form_glmer,
  data         = dat,
  family       = poisson(),
  pfamily_list = pfamily_list(ps),
  n            = 1000L,
  mode_gap_max = 1.0,
)

re_names <- fit$model_setup$re_coef_names
stopifnot(identical(re_names, c("(Intercept)", "rating_c", "log_price_c")))

n_draws <- nrow(fit$fixef[[re_names[1L]]])
stopifnot(identical(n_draws, 1000L))
stopifnot(!is.null(fit$fixef.init))
stopifnot(is.list(fit$convergence), is.finite(fit$convergence$m_convergence))
stopifnot(fit$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))
stopifnot(identical(fit$convergence$n_pilot_source, "cost"))
stopifnot(is.finite(fit$convergence$m_convergence_pilot))
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

cat("glmerb airbnb (full model) centers table:\n")
cat(.print_glmerb_test_table(center_tab), "\n\n")
cat("glmerb airbnb differences table:\n")
cat(.print_glmerb_test_table(diff_tab), "\n\n")
cat("glmerb airbnb univariate z/p table:\n")
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
  paste0(
    "glmerb airbnb centering test: n_pilot=%d, m_convergence=%d, ",
    "p(mean=pilot)=%.4g, p(mean=mode)=%.4g\n\n"
  ),
  fit$pilot_chisq$n_pilot,
  fit$convergence$m_convergence,
  p_pilot, p_mode
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
