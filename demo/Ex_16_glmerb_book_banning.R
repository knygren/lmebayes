## Demo: glmerb() binomial workflow on bayesrules::book_banning
##
## Bayes Rules! Ch. 18.2--18.4 hierarchical logistic regression of book-removal
## outcomes on challenge reasons with state random effects.  The textbook uses
##
##   removed ~ violent + antifamily + language + (1 | state)
##
## lmebayes requires climber-/challenge-level predictors to appear as population
## mean slopes with matching random effects (see demo/Ex_13_glmerb_Airbnb.R).
## This demo fits a minimal binomial glmerb model with violent coded as 0/1
## (violent_i) and a random slope on violent_i, and compares the classical
## glmer fit for the Ch. 18 random-intercept formula.
##
## After fitting:
##   print(fit, sweep_history = TRUE, max_sweeps = 5)
##   plot_sweep_history_diag(fit$sweep_history$main, coef_focus, what = "mean")
##
##   demo("Ex_16_glmerb_book_banning", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

data(book_banning, package = "bayesrules", envir = environment())

dat <- book_banning[, c(
  "state", "removed", "violent", "antifamily", "language"
)]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)

form_book <- removed ~ violent + antifamily + language + (1 | state)
form_glmerb <- removed_i ~ violent_i + (1 + violent_i || state)

ps <- Prior_Setup_GLMM(form_glmerb, data = dat, family = binomial(), pwt = 0.01)

print(ps)

fit <- glmerb(
  form_glmerb,
  data         = dat,
  family       = binomial(),
  pfamily_list = pfamily_list(ps),
  n            = 3000L,
  mode_gap_max = 1.0,
  progbar = TRUE
)

cat("m_convergence used:", fit$convergence$m_convergence, "\n")
cat(sprintf(
  "Pilot vs mode (chi-squared): p = %.4g (n_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot
))

cat("\n--- Ch. 18 reference glmer (random intercept; all three reasons) ---\n")
fit_book <- lme4::glmer(form_book, data = dat, family = binomial())
print(lme4::fixef(fit_book))

lmebayes:::print_coef_means(fit)
print(fit)
summary(fit)

coef_focus <- list(
  c("(Intercept)", "(Intercept)"),
  c("violent_i", "(Intercept)")
)

for (st in list(fit$sweep_history$pilot, fit$sweep_history$main)) {
  if (is.null(st)) next
  plot_sweep_history_diag(st, coef_focus, what = "sd")
}

## --- Temporary: cross-chain variance by inner sweep (coef_focus) ---------------
## Raw variance printout/plots, then Var / Var_main(final) ratios (Claim 3
## proxy in non-whitened space). Stacked pilot/main figures; pilot on top.

.ex16_var_by_coef_focus <- function(tab, coef_focus, min_sweep = 1L) {
  stats::setNames(
    lapply(coef_focus, function(cc) {
      re_comp <- as.character(cc[1L])
      cov <- as.character(cc[2L])
      sub <- tab[
        tab$re_component == re_comp &
          tab$covariate == cov &
          tab$sweep >= min_sweep,
        ,
        drop = FALSE
      ]
      sub <- sub[order(sub$sweep), , drop = FALSE]
      if (!nrow(sub)) {
        stop(
          "Missing variance rows for ", re_comp, " | ", cov,
          call. = FALSE
        )
      }
      stats::setNames(sub$sd^2, sub$sweep)
    }),
    vapply(coef_focus, function(cc) {
      paste(as.character(cc[1L]), as.character(cc[2L]), sep = " | ")
    }, character(1L))
  )
}

.ex16_cat_var_stage <- function(stage_label, var_list, value_label = "variance") {
  for (nm in names(var_list)) {
    v <- var_list[[nm]]
    cat(
      "  ", stage_label, " ", value_label, " (", nm, "): ",
      paste(names(v), format(v, digits = 3), sep = "=", collapse = ", "),
      "\n",
      sep = ""
    )
  }
}

.ex16_var_ref_main_final <- function(var_main) {
  vapply(
    var_main,
    function(v) {
      m_final <- max(as.integer(names(v)))
      ref <- as.numeric(v[as.character(m_final)])
      if (length(ref) != 1L || !is.finite(ref) || ref <= 0) {
        stop(
          "Need a positive main-final variance reference at sweep ", m_final,
          call. = FALSE
        )
      }
      ref
    },
    numeric(1L)
  )
}

.ex16_var_ratio <- function(var_list, ref) {
  stats::setNames(
    lapply(names(var_list), function(nm) {
      stats::setNames(var_list[[nm]] / ref[nm], names(var_list[[nm]]))
    }),
    names(var_list)
  )
}

.ex16_ylim_var_ratio <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- vals[is.finite(vals)]
  y_top <- max(1, vals, na.rm = TRUE)
  if (!is.finite(y_top) || y_top <= 0) {
    y_top <- 1
  }
  c(0, y_top * 1.05)
}

.ex16_caption_below <- function(text, line = 4.2) {
  if (!nzchar(text)) {
    return(invisible(NULL))
  }
  graphics::mtext(text, side = 1, line = line, cex = 0.85)
}

.ex16_legend_below <- function(labels, cols, pch = 16L) {
  if (length(labels) <= 1L) {
    return(invisible(NULL))
  }
  old_xpd <- graphics::par("xpd")
  on.exit(graphics::par(xpd = old_xpd), add = TRUE)
  graphics::par(xpd = TRUE)
  usr <- graphics::par("usr")
  graphics::legend(
    x = mean(usr[1:2]),
    y = usr[3] - 0.24 * diff(usr[3:4]),
    legend = labels,
    col = cols,
    pch = pch,
    horiz = TRUE,
    bty = "n",
    cex = 0.85,
    xjust = 0.5
  )
}

.ex16_plot_crosschain_var <- function(var_list, stage_label, xlab = "Inner sweep",
                                     ylim = NULL, n_chains = NULL,
                                     ylab = "Cross-chain variance",
                                     sub = "(cross-chain variance from sweep history)",
                                     ref_line = NULL) {
  if (!length(var_list)) {
    warning("No variance series to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  sweeps <- sort(unique(as.integer(unlist(lapply(var_list, names)))))
  if (!length(sweeps)) {
    warning("No sweep rows to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  if (is.null(ylim)) {
    y_top <- max(unlist(var_list), na.rm = TRUE)
    if (!is.finite(y_top) || y_top <= 0) {
      y_top <- 1
    }
    ylim <- c(0, y_top * 1.05)
  }
  if (!is.null(n_chains)) {
    sub <- paste0(sub, sprintf("; n = %d chains", as.integer(n_chains)))
  }
  cols <- grDevices::hcl.colors(length(var_list), palette = "Dark 3")
  y_first <- var_list[[1L]][match(sweeps, names(var_list[[1L]]))]
  graphics::plot(
    sweeps, y_first,
    type = "n",
    xlab = xlab,
    ylab = ylab,
    ylim = ylim,
    main = paste0(stage_label, ": ", ylab, " vs sweep")
  )
  graphics::grid()
  if (!is.null(ref_line)) {
    graphics::abline(h = ref_line, lty = 2, col = "gray40")
  }
  for (i in seq_along(var_list)) {
    v <- var_list[[i]]
    x <- as.integer(names(v))
    graphics::lines(x, v, type = "b", pch = 16, col = cols[i])
  }
  .ex16_caption_below(sub)
  .ex16_legend_below(names(var_list), cols)
  invisible(var_list)
}

if (!is.null(fit$sweep_history$main)) {
  re_names_var <- fit$sweep_history$main$re_names
  n_main_var <- if (!is.null(fit$fixef)) {
    nrow(as.matrix(fit$fixef[[re_names_var[1L]]]))
  } else {
    NULL
  }
  n_pilot_var <- fit$pilot_chisq$n_pilot
  if (is.null(n_pilot_var) && !is.null(fit$convergence)) {
    n_pilot_var <- fit$convergence$n_pilot
  }

  cat("\n--- Temporary cross-chain variance (coef_focus) ---\n")
  if (!is.null(n_pilot_var) && n_pilot_var > 0L && !is.null(n_main_var)) {
    cat(sprintf("  Pilot n = %d chains; Main n = %d chains\n", n_pilot_var, n_main_var))
  } else if (!is.null(n_main_var)) {
    cat(sprintf("  Main n = %d chains\n", n_main_var))
  }

  var_pilot <- NULL
  if (!is.null(fit$sweep_history$pilot)) {
    var_pilot <- .ex16_var_by_coef_focus(
      fit$sweep_history$pilot$table,
      coef_focus,
      min_sweep = 0L
    )
    .ex16_cat_var_stage("Pilot", var_pilot)
  }

  var_main <- .ex16_var_by_coef_focus(
    fit$sweep_history$main$table,
    coef_focus,
    min_sweep = 1L
  )
  .ex16_cat_var_stage("Main", var_main)

  if (!is.null(var_pilot)) {
    y_top <- max(c(unlist(var_pilot), unlist(var_main)), na.rm = TRUE)
    if (!is.finite(y_top) || y_top <= 0) {
      y_top <- 1
    }
    ylim_shared <- c(0, y_top * 1.05)
    graphics::par(mfrow = c(2L, 1L), mar = c(8, 4, 3, 1) + 0.1)
    .ex16_plot_crosschain_var(
      var_pilot,
      "Pilot stage",
      xlab = "Inner sweep (0 = common start)",
      ylim = ylim_shared,
      n_chains = n_pilot_var
    )
    .ex16_plot_crosschain_var(
      var_main,
      "Main stage",
      xlab = "Inner sweep",
      ylim = ylim_shared,
      n_chains = n_main_var
    )
    graphics::par(mfrow = c(1L, 1L), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    graphics::par(mar = c(8, 4, 3, 1) + 0.1)
    .ex16_plot_crosschain_var(
      var_main,
      "Main stage",
      xlab = "Inner sweep",
      n_chains = n_main_var
    )
    graphics::par(mar = c(5, 4, 4, 2) + 0.1)
  }

  ref_var <- .ex16_var_ref_main_final(var_main)
  m_final_var <- max(as.integer(names(var_main[[1L]])))
  cat(sprintf(
    "\n--- Temporary variance ratio: Var / Var_main(sweep %d) ---\n",
    m_final_var
  ))
  cat("  (non-whitened; fraction of main-final cross-chain variance, Claim 3)\n")
  for (nm in names(ref_var)) {
    cat(sprintf("  Reference variance (%s): %.4g\n", nm, ref_var[[nm]]))
  }

  var_ratio_pilot <- NULL
  if (!is.null(var_pilot)) {
    var_ratio_pilot <- .ex16_var_ratio(var_pilot, ref_var)
    .ex16_cat_var_stage("Pilot", var_ratio_pilot, value_label = "variance ratio")
  }
  var_ratio_main <- .ex16_var_ratio(var_main, ref_var)
  .ex16_cat_var_stage("Main", var_ratio_main, value_label = "variance ratio")

  ylim_ratio <- if (!is.null(var_ratio_pilot)) {
    .ex16_ylim_var_ratio(var_ratio_pilot, var_ratio_main)
  } else {
    .ex16_ylim_var_ratio(var_ratio_main)
  }
  sub_ratio <- paste0(
    "Var / Var_main(sweep ", m_final_var,
    "); 1 = main-final scale (may exceed 1 from sampling)"
  )
  if (!is.null(var_ratio_pilot)) {
    graphics::par(mfrow = c(2L, 1L), mar = c(8, 4, 3, 1) + 0.1)
    .ex16_plot_crosschain_var(
      var_ratio_pilot,
      "Pilot stage",
      xlab = "Inner sweep (0 = common start)",
      ylim = ylim_ratio,
      n_chains = n_pilot_var,
      ylab = "Var / Var_main(final)",
      sub = sub_ratio,
      ref_line = 1
    )
    .ex16_plot_crosschain_var(
      var_ratio_main,
      "Main stage",
      xlab = "Inner sweep",
      ylim = ylim_ratio,
      n_chains = n_main_var,
      ylab = "Var / Var_main(final)",
      sub = sub_ratio,
      ref_line = 1
    )
    graphics::par(mfrow = c(1L, 1L), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    graphics::par(mar = c(8, 4, 3, 1) + 0.1)
    .ex16_plot_crosschain_var(
      var_ratio_main,
      "Main stage",
      xlab = "Inner sweep",
      ylim = ylim_ratio,
      n_chains = n_main_var,
      ylab = "Var / Var_main(final)",
      sub = sub_ratio,
      ref_line = 1
    )
    graphics::par(mar = c(5, 4, 4, 2) + 0.1)
  }
}

## --- Temporary: proxy D_l = (E[gamma|start] - mu)^T Sigma^{-1} (.) -----------
## Uses main-final cross-chain mean and covariance (n stored draws) as proxies
## for mu_1 and Sigma_11.  Pilot and main each get one plot vs inner sweep.
## (Sweep-history tables supply cross-chain means; Sigma_hat from cov of fit$fixef.)

.ex16_stack_sweep_means <- function(tab, sweep, re_names, fixef_template) {
  v <- numeric(0L)
  for (k in re_names) {
    cn <- colnames(fixef_template[[k]])
    for (nm in cn) {
      hit <- tab$sweep == sweep &
        tab$re_component == k &
        tab$covariate == nm
      if (!any(hit)) {
        stop(
          "Missing sweep mean for ", k, " | ", nm, " at sweep ", sweep,
          call. = FALSE
        )
      }
      v <- c(v, tab$mean[hit][1L])
    }
  }
  v
}

.ex16_mahalanobis_sq <- function(gamma, mu_hat, Sigma_inv) {
  delta <- gamma - mu_hat
  as.numeric(t(delta) %*% Sigma_inv %*% delta)
}

.ex16_stack_fixef_list <- function(fixef_list, re_names, fixef_template) {
  v <- numeric(0L)
  for (k in re_names) {
    cn <- colnames(fixef_template[[k]])
    vec <- fixef_list[[k]]
    if (is.matrix(vec)) {
      vec <- vec[1L, , drop = TRUE]
    }
    for (nm in cn) {
      v <- c(v, as.numeric(vec[nm]))
    }
  }
  v
}

.ex16_proxy_D_by_sweep <- function(tab, re_names, fixef_template, mu_hat, Sigma_inv,
                                  min_sweep = 1L) {
  sweeps <- sort(unique(tab$sweep[tab$sweep >= min_sweep]))
  D <- vapply(
    sweeps,
    function(m) {
      gm <- .ex16_stack_sweep_means(tab, m, re_names, fixef_template)
      .ex16_mahalanobis_sq(gm, mu_hat, Sigma_inv)
    },
    numeric(1L)
  )
  stats::setNames(D, sweeps)
}

.ex16_plot_proxy_D <- function(D, stage_label, xlab = "Inner sweep", ylim = NULL,
                              n_chains = NULL) {
  sweeps <- as.integer(names(D))
  if (!length(sweeps)) {
    warning("No sweep rows to plot for stage ", stage_label, call. = FALSE)
    return(invisible(NULL))
  }
  if (is.null(ylim)) {
    y_top <- max(D, na.rm = TRUE)
    if (!is.finite(y_top) || y_top <= 0) {
      y_top <- 1
    }
    ylim <- c(0, y_top * 1.05)
  }
  sub <- "(ref: main-final mu, Sigma from fit$fixef)"
  if (!is.null(n_chains)) {
    sub <- paste0(sub, sprintf("; n = %d chains", as.integer(n_chains)))
  }
  graphics::plot(
    sweeps, D,
    type = "b", pch = 16,
    xlab = xlab,
    ylab = expression(hat(D)[l] ~ "(Mahalanobis"^2*")"),
    ylim = ylim,
    main = paste0(stage_label, ": proxy D[l] vs sweep")
  )
  graphics::grid()
  .ex16_caption_below(sub)
  invisible(D)
}

if (!is.null(fit$sweep_history$main) && !is.null(fit$fixef)) {
  re_names <- fit$sweep_history$main$re_names
  G_main <- do.call(
    cbind,
    lapply(re_names, function(k) as.matrix(fit$fixef[[k]]))
  )
  mu_hat <- colMeans(G_main)
  Sigma_inv <- solve(stats::cov(G_main))
  n_main <- nrow(G_main)
  n_pilot <- fit$pilot_chisq$n_pilot
  if (is.null(n_pilot) && !is.null(fit$convergence)) {
    n_pilot <- fit$convergence$n_pilot
  }
  gamma_mode <- .ex16_stack_fixef_list(fit$fixef.mode, re_names, fit$fixef)
  D_mode <- .ex16_mahalanobis_sq(gamma_mode, mu_hat, Sigma_inv)

  cat("\n--- Temporary proxy D_l (main-final mu, Sigma from fit$fixef) ---\n")
  cat(sprintf(
    "  q = %d; mu_hat and Sigma_hat from n = %d main chains\n",
    ncol(G_main), n_main
  ))
  if (!is.null(n_pilot) && n_pilot > 0L) {
    cat(sprintf("  Pilot n = %d chains; Main n = %d chains\n", n_pilot, n_main))
  } else {
    cat(sprintf("  Main n = %d chains\n", n_main))
  }
  cat(sprintf("  D at ICM mode: %.4g\n", D_mode))

  D_main_sweeps <- .ex16_proxy_D_by_sweep(
    fit$sweep_history$main$table,
    re_names,
    fit$fixef,
    mu_hat,
    Sigma_inv
  )
  D_main <- c(`-1` = D_mode, D_main_sweeps)
  if (!is.null(fit$fixef.init)) {
    gamma_pilot <- .ex16_stack_fixef_list(fit$fixef.init, re_names, fit$fixef)
    D_pilot_start <- .ex16_mahalanobis_sq(gamma_pilot, mu_hat, Sigma_inv)
    cat(sprintf("  D at pilot mean (main start): %.4g\n", D_pilot_start))
    D_main <- c(`-1` = D_mode, `0` = D_pilot_start, D_main_sweeps)
  }
  cat("  Main  D_l:", paste(names(D_main), format(D_main, digits = 3),
                            sep = "=", collapse = ", "), "\n")

  D_pilot <- NULL
  if (!is.null(fit$sweep_history$pilot)) {
    D_pilot <- .ex16_proxy_D_by_sweep(
      fit$sweep_history$pilot$table,
      re_names,
      fit$fixef,
      mu_hat,
      Sigma_inv,
      min_sweep = 0L
    )
    cat("  Pilot D_l:", paste(names(D_pilot), format(D_pilot, digits = 3),
                              sep = "=", collapse = ", "), "\n")
  }

  if (!is.null(D_pilot)) {
    y_top <- max(c(D_pilot, D_main), na.rm = TRUE)
    if (!is.finite(y_top) || y_top <= 0) {
      y_top <- 1
    }
    ylim_shared <- c(0, y_top * 1.05)
    graphics::par(mfrow = c(2L, 1L), mar = c(7, 4, 3, 1) + 0.1)
    .ex16_plot_proxy_D(
      D_pilot,
      "Pilot stage",
      xlab = "Inner sweep (0 = ICM mode)",
      ylim = ylim_shared,
      n_chains = n_pilot
    )
    .ex16_plot_proxy_D(
      D_main,
      "Main stage",
      xlab = paste(
        "Stage step (-1 = ICM mode, 0 = pilot mean,",
        "1+ = main inner sweep)"
      ),
      ylim = ylim_shared,
      n_chains = n_main
    )
    graphics::par(mfrow = c(1L, 1L), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    graphics::par(mar = c(7, 4, 3, 1) + 0.1)
    .ex16_plot_proxy_D(
      D_main,
      "Main stage",
      xlab = paste(
        "Stage step (-1 = ICM mode, 0 = pilot mean,",
        "1+ = main inner sweep)"
      ),
      n_chains = n_main
    )
    graphics::par(mar = c(5, 4, 4, 2) + 0.1)
  }
}

cat("\n=== Random effects: glmer reference vs glmerb chain mean ===\n\n")
lmebayes:::print_mer_bayes_re_compare(fit)
