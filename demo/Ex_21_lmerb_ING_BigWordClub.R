## Demo: lmerb() with ING Block~2 priors + pilot on the full big_word_club model
##
## Combines demo/Ex_12_lmerb_BigWordClub.R (school random intercept and slopes,
## cross-level moderation) with demo/Ex_20_lmerb_ING_pilot.R (dIndependent_Normal_
## Gamma dispersion, two-stage pilot/main sampling).
##
## With diag_sweeps = TRUE the fit uses the sweep-outer driver and prints one
## combined Block~2 chain-mean table per stage (pilot, then main), like
## print(fit$sweep_history$pilot) / $main.
##
##   demo("Ex_21_lmerb_ING_BigWordClub", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
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

design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data             = dat,
  pwt              = 0.01,
  pwt_dispersion   = 0.2
)
cat("\n=== Prior_Setup_lmebayes (ING calibration) ===\n\n")
print(ps)

pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

## ING fit with sweep diagnostics (diag_sweeps uses the R sweep-outer driver;
## slower than the default C++ engine but prints/plots per-inner-sweep tables).
## gap_tol = 0.05 => n_pilot from the Hotelling bound (~16 for this model).
## progbar defaults to FALSE when diag_sweeps = TRUE (tables use print()).
fit <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = ps$dispersion_ranef,
  n                = 10000L,
  gap_tol          = 0.05,
  mode_gap_max     = 1.0,
  diag_sweeps      = TRUE
)

stopifnot(isTRUE(fit$prior$any_non_normal))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(!is.null(fit$sweep_history))
stopifnot(!is.null(fit$sweep_history$pilot))
stopifnot(!is.null(fit$sweep_history$main))
stopifnot(fit$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))
stopifnot(is.finite(fit$pilot_chisq$p_value))

cat("\n=== Stored sweep history (last 3 sweeps per stage) ===\n\n")
print(fit$sweep_history$pilot, max_sweeps = 3000)
print(fit$sweep_history$main, max_sweeps = 3000)

cat("\n=== summary(fit) ===\n\n")
print(summary(fit))

re_names <- fit$model_setup$re_coef_names

## tau^2_k per RE component: sampled inside [disp_lower, disp_upper].
for (k in re_names) {
  pr_k <- pf[[k]]$prior_list
  t2   <- fit$fixef.dispersion[, k]
  stopifnot(
    all(is.finite(t2)), all(t2 > 0),
    all(t2 >= pr_k$disp_lower),
    all(t2 <= pr_k$disp_upper),
    stats::sd(t2) > 0
  )
  cat(sprintf(
    "\n%s tau^2: post mean = %.4f  [window (%.4f, %.4f); plugin disp_lower = %.4f]\n",
    k,
    fit$fixef.dispersion.mean[[k]],
    pr_k$disp_lower,
    pr_k$disp_upper,
    pr_k$disp_lower
  ))
}

cat(sprintf(
  "\nPilot vs plug-in start (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  fit$convergence$m_convergence_pilot
))

## Block~2 hyperparameters: prior mean, gamma @ lmer tau2, pilot init, MCMC mean.
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
beta_bar <- unlist(lapply(re_names, function(k) fit$fixef.means[[k]]))
theta_plug <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
theta_prior <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  unname(fit$prior$prior_list[[k]]$mu_fixef[nms])
}))
theta_pilot <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  unname(fit$fixef.init[[k]][nms])
}))
names(beta_bar) <- names(theta_plug) <- names(theta_prior) <- names(theta_pilot) <- cn

block2_cmp <- data.frame(
  prior_mean      = unname(theta_prior),
  gamma_lmer_tau2 = unname(theta_plug),
  pilot_mean      = unname(theta_pilot),
  mcmc_mean       = unname(beta_bar),
  row.names       = cn,
  check.names     = FALSE
)
cat("\n=== Block 2 hyperparameters (prior / plug-in / pilot / MCMC) ===\n\n")
print(round(block2_cmp, 4))

## Multivariate centering: main mean vs pilot start vs plug-in start.
X <- do.call(cbind, lapply(re_names, function(k) fit$fixef[[k]]))
colnames(X) <- cn
V <- stats::cov(X)
V_inv <- tryCatch(
  solve(V),
  error = function(e) solve(V + diag(1e-8 * mean(diag(V)), ncol(V)))
)
d_pilot <- beta_bar - theta_pilot
d_plug  <- beta_bar - theta_plug
Q_pilot <- as.numeric(t(d_pilot) %*% V_inv %*% d_pilot)
Q_plug  <- as.numeric(t(d_plug) %*% V_inv %*% d_plug)
p_pilot <- stats::pchisq(Q_pilot, df = ncol(X), lower.tail = FALSE)
p_plug  <- stats::pchisq(Q_plug, df = ncol(X), lower.tail = FALSE)

cat(sprintf(
  "\nOverall centering (chi-squared, p = %d hyperparameters): p(mean=pilot)=%.4g, p(mean=plug-in)=%.4g\n",
  ncol(X), p_pilot, p_plug
))

pilot_mean <- unname(theta_pilot)
post_mean  <- unname(beta_bar)
post_sd    <- unlist(lapply(re_names, function(k) apply(fit$fixef[[k]], 2L, sd)))
plug_in    <- unname(theta_plug)
names(post_sd) <- cn

n_main <- nrow(fit$fixef[[re_names[1L]]])
mc_se  <- post_sd / sqrt(n_main)

tab <- data.frame(
  pilot_mean  = round(pilot_mean, 4),
  post_mean   = round(post_mean, 4),
  difference  = round(post_mean - pilot_mean, 4),
  post_sd     = round(post_sd, 4),
  mc_se       = round(mc_se, 4),
  z_vs_pilot  = round((post_mean - pilot_mean) / mc_se, 2),
  plug_in     = round(plug_in, 4),
  row.names   = cn,
  check.names = FALSE
)

cat("\n=== Block 2: pilot mean vs posterior mean ===\n\n")
print(tab)

## Sweep-history diagnostics: cross-chain mean and SD vs inner sweep (pilot and main).
## Helps spot coefficients whose chains spread or drift across inner sweeps.
coef_focus <- list(
  c("(Intercept)", "(Intercept)"),
  c("(Intercept)", "private_school"),
  c("(Intercept)", "title1"),
  c("(Intercept)", "free_reduced_lunch"),
  c("distracted_ppvt", "(Intercept)"),
  c("distracted_a1", "(Intercept)"),
  c("distracted_a1", "free_reduced_lunch")
)

plot_sweep_history_diag <- function(
    hist,
    coef_focus,
    what = c("sd", "mean"),
    stage_label = hist$stage
) {
  what <- match.arg(what)
  sh_tab <- hist$table
  sh_sweeps <- subset(sh_tab, sweep > 0L)
  if (!nrow(sh_sweeps)) {
    warning("No sweep rows in sweep history for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }

  ylab <- if (what == "sd") "Cross-chain SD" else "Cross-chain mean"
  cat(sprintf(
    "\n=== %s sweep history plots (cross-chain %s) ===\n\n",
    stage_label, what
  ))

  plot_one <- function(re_comp, cov) {
    sub <- subset(
      sh_sweeps,
      re_component == re_comp & covariate == cov
    )
    if (!nrow(sub)) {
      warning("No sweep rows for ", re_comp, " | ", cov, call. = FALSE)
      return(invisible(NULL))
    }
    y <- if (what == "sd") sub$sd else sub$mean
    plot(
      sub$sweep, y,
      type = "b", pch = 16,
      xlab = "Inner sweep", ylab = ylab,
      main = paste(re_comp, cov, sep = " | ")
    )
    if (what == "mean") {
      mode_val <- subset(
        sh_tab,
        re_component == re_comp & covariate == cov & sweep == 0L
      )$mean
      abline(h = mode_val, lty = 2, col = "gray40")
    }
    invisible(sub)
  }

  op <- par(
    mfrow = c(length(coef_focus), 1L),
    mar = c(4, 4, 2.5, 1),
    oma = c(0, 0, 2, 0)
  )
  on.exit(par(op), add = TRUE)
  for (cc in coef_focus) {
    plot_one(cc[1L], cc[2L])
  }
  mtext(
    sprintf("%s Block 2 fixef: cross-chain %s by inner sweep", stage_label, what),
    outer = TRUE, line = 0.5, cex = 0.95
  )
  if (what == "mean") {
    mtext(
      "Dashed line = ICM mode (sweep 0)",
      outer = TRUE, line = -1.5, cex = 0.85
    )
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    sh_sweeps$coef <- interaction(
      sh_sweeps$re_component, sh_sweeps$covariate, sep = " | "
    )
    y_var <- if (what == "sd") "sd" else "mean"
    p <- ggplot2::ggplot(
      sh_sweeps,
      ggplot2::aes(sweep, .data[[y_var]], group = coef, colour = coef)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ coef, scales = "free_y") +
      ggplot2::labs(
        x = "Inner sweep",
        y = ylab,
        title = sprintf(
          "%s Block 2 fixef - cross-chain %s by sweep",
          stage_label, what
        )
      ) +
      ggplot2::theme(legend.position = "none")
    if (what == "mean") {
      mode_df <- subset(sh_tab, sweep == 0L)
      mode_df$coef <- interaction(
        mode_df$re_component, mode_df$covariate, sep = " | "
      )
      p <- p + ggplot2::geom_hline(
        ggplot2::aes(yintercept = mean, linetype = "ICM mode"),
        data = mode_df,
        colour = "gray40"
      ) +
        ggplot2::scale_linetype_manual(name = NULL, values = c("ICM mode" = "dashed"))
    }
    print(p)
  }
  invisible(hist)
}

for (st in list(fit$sweep_history$pilot, fit$sweep_history$main)) {
  plot_sweep_history_diag(st, coef_focus, what = "sd")
  plot_sweep_history_diag(st, coef_focus, what = "mean")
}
