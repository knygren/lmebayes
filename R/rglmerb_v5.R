#' Raw two-stage Gibbs sampler (v5 short-chain driver, sweep-outer loop)
#'
#' Pilot and main sampling call \code{\link{run_short_chains_v5}} (a single
#' \code{\link[glmbayesCore]{two_block_rNormal_reg_v5}} call per stage with
#' sweep-outer C++ loop order and per-sweep chain progress bars).
#'
#' Development driver with \code{fixef_temp} in C++.  \code{\link{glmerb}} calls
#' \code{rglmerb_v5} for Poisson/non-Gaussian sampling via \code{\link{rglmerb}}.
#'
#' @inheritParams rglmerb
#' @param n_pilot Optional integer pilot chain count. When \code{NULL}, derived
#'   from \code{tv_tol} / \code{mode_gap_max} (cost optimizer) or \code{gap_tol}.
#' @param m_convergence_pilot Optional integer inner sweeps for the pilot stage.
#' @return An object of class \code{c("rglmerb_v5", "list")} with staging fields
#'   compatible with \code{\link{rglmerb}}.
#' @seealso \code{\link{rglmerb}}, \code{\link{glmerb}},
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg_v5}}
#' @export
rglmerb_v5 <- function(
    n,
    design,
    prior,
    family              = poisson(),
    fixef_start         = NULL,
    m_convergence       = NULL,
    n_pilot             = NULL,
    m_convergence_pilot = NULL,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    gap_tol             = 0.0196,
    collect_block1      = TRUE,
    verbose             = TRUE,
    progbar             = FALSE
) {
  cl <- match.call()

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  if (!inherits(design, "model_setup")) {
    stop("'design' must be a model_setup object.", call. = FALSE)
  }

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  if (!is.null(m_convergence)) {
    if (!is.numeric(m_convergence) || length(m_convergence) != 1L ||
        !is.finite(m_convergence) || m_convergence < 1) {
      stop("'m_convergence' must be NULL or a single integer >= 1.", call. = FALSE)
    }
    m_convergence <- as.integer(m_convergence)
  }

  if (!is.null(m_convergence_pilot)) {
    if (!is.numeric(m_convergence_pilot) ||
        length(m_convergence_pilot) != 1L ||
        !is.finite(m_convergence_pilot) || m_convergence_pilot < 1) {
      stop("'m_convergence_pilot' must be NULL or a single integer >= 1.",
           call. = FALSE)
    }
    m_convergence_pilot <- as.integer(m_convergence_pilot)
  }

  n_pilot_int <- if (is.null(n_pilot) || identical(n_pilot, 0L) ||
                     identical(n_pilot, 0)) {
    0L
  } else {
    as.integer(n_pilot[1L])
  }
  n_pilot_source <- if (n_pilot_int > 0L) "explicit" else "none"

  if (!is.null(gap_tol) && identical(n_pilot_int, 0L) && is.null(n_pilot)) {
    if (!is.numeric(gap_tol) || length(gap_tol) != 1L ||
        !is.finite(gap_tol) || gap_tol <= 0) {
      stop("'gap_tol' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  run_pilot <- n_pilot_int > 0L

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  is_gaussian  <- identical(family$family, "gaussian")

  ranef_mode <- NULL
  if (is.null(fixef_start)) {
    pm          <- glmbayesCore::glmerb_posterior_mode(design, family, prior)
    fixef_start <- pm$fixef
    ranef_mode  <- pm$b_mean

    if (verbose) {
      fixef_glmer <- lapply(prior$prior_list, `[[`, "mu_fixef")
      names(fixef_glmer) <- re_names
      hdr <- sprintf("  %-18s  %-30s  %12s  %12s",
                     "RE component", "parameter", "glmer (start)", "post mode (ICM)")
      sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
      cat("--- glmerb: Block 2 fixed effects ---\n")
      cat(hdr, "\n")
      cat(sep, "\n")
      for (k in re_names) {
        nms_k   <- names(fixef_glmer[[k]])
        glmer_v <- fixef_glmer[[k]]
        pm_v    <- fixef_start[[k]]
        for (nm in nms_k) {
          cat(sprintf("  %-18s  %-30s  %12.4f  %12.4f\n",
                      k, nm, glmer_v[[nm]], pm_v[[nm]]))
        }
      }
      cat(sprintf("  (ICM converged: %s, %d iter, delta = %.2e)\n\n",
                  pm$converged, pm$iterations, pm$delta))
    }
  }

  .lmebayes_print_ranef_mode_reference(
    ranef_mode, re_names, group_levels, verbose
  )

  fixef_mode_ref <- fixef_start
  b_mode_ref     <- ranef_mode
  diag_sweeps    <- FALSE

  block1_prior <- .lmebayes_block1_prior_list(prior)

  if (is_gaussian) {
    rate <- glmbayesCore::two_block_rate_v2(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = block1_prior,
      pfamily_list      = prior$pfamily_list,
      family            = gaussian(),
      group_levels      = group_levels
    )
  } else {
    mode_w <- glmbayesCore::two_block_mode_weights(
      x            = design$Z,
      block        = design$groups,
      b_mode       = ranef_mode,
      family       = family,
      group_levels = group_levels
    )
    rate <- glmbayesCore::two_block_rate_v2(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = block1_prior,
      pfamily_list      = prior$pfamily_list,
      weights           = mode_w$weights,
      family            = family,
      group_levels      = group_levels
    )
  }

  m_min <- glmbayesCore::two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L

  p_dim            <- sum(vapply(fixef_start, length, integer(1L)))
  D_max            <- if (!is.null(mode_gap_max)) sqrt(p_dim) * mode_gap_max else 0
  m_pilot_from_gap <- NULL
  m_convergence_pilot_plan <- m_convergence_pilot

  if (is.null(m_convergence_pilot_plan) && !is.null(mode_gap_max) && !is.null(tv_tol)) {
    erf1_inv_tv      <- stats::qnorm((tv_tol + 1) / 2) / sqrt(2)
    c_tol            <- erf1_inv_tv * 2 * sqrt(2)
    m_pilot_from_gap <- if (D_max <= c_tol || rate$lambda_star <= 0) {
      m_min
    } else {
      as.integer(ceiling(log(D_max / c_tol) / log(1 / rate$lambda_star)))
    }
    m_convergence_pilot_plan <- max(m_min, m_pilot_from_gap)
  }

  if (identical(n_pilot_int, 0L) && is.null(n_pilot)) {
    if (!is.null(tv_tol) && !is.null(m_convergence_pilot_plan)) {
      pilot_cost_opt <- glmbayesCore::two_block_optimize_pilot_cost(
        n                   = n,
        rate                = rate,
        tv_tol              = tv_tol,
        m_convergence_pilot = m_convergence_pilot_plan,
        p                   = p_dim
      )
      n_pilot_int <- pilot_cost_opt$n_pilot_opt
      n_pilot_source <- "cost"
      if (is.null(m_convergence)) {
        m_convergence <- pilot_cost_opt$m_convergence_opt
      }
    } else if (!is.null(gap_tol)) {
      n_pilot_int <- as.integer(ceiling((stats::qnorm(0.975) / gap_tol)^2))
      n_pilot_source <- "gap_tol"
    }
  }

  run_pilot <- n_pilot_int > 0L

  if (run_pilot && is.null(m_convergence_pilot)) {
    m_convergence_pilot <- if (!is.null(m_convergence)) {
      m_convergence
    } else if (!is.null(m_convergence_pilot_plan)) {
      m_convergence_pilot_plan
    } else if (!is.null(tv_tol)) {
      m_min
    } else {
      10L
    }
  }

  if (is.null(m_convergence)) {
    m_convergence <- m_min
  } else if (m_convergence < m_min) {
    warning(
      "rglmerb_v5: m_convergence = ", m_convergence, " is below the derived ",
      "minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
      "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }

  calib_label <- if (is_gaussian) {
    "exact (Gaussian posterior)"
  } else {
    sprintf("approximate (local-Gaussian at mode, %s)", family$family)
  }
  if (isTRUE(prior$any_non_normal)) {
    calib_label <- paste0(calib_label, "; conservative: ING tau^2_k = disp_lower")
  }

  if (verbose) {
    cat(sprintf(
      "--- glmerb: convergence calibration [%s]:\n    lambda* = %.4f, tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n",
      calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
    if (run_pilot && !is.null(mode_gap_max) && !is.null(m_pilot_from_gap)) {
      cat(sprintf(
        "--- glmerb: pilot sweep calibration [mode_gap_max = %g SD/dim, p = %d, D_max = %.4f]:\n    m_min = %d, lambda* = %.4f => m_convergence_pilot = %d ---\n\n",
        mode_gap_max, p_dim, D_max, m_min, rate$lambda_star, m_convergence_pilot
      ))
    }
  }

  method_label <- if (is_gaussian) "exact" else "local_gaussian_mode"
  if (isTRUE(prior$any_non_normal)) {
    method_label <- paste0(method_label, "+disp_lower_bound")
  }

  convergence_info <- list(
    method              = method_label,
    tv_tol              = tv_tol,
    gap_tol             = gap_tol,
    n_pilot             = if (run_pilot) n_pilot_int else 0L,
    n_pilot_source      = n_pilot_source,
    lambda_star         = rate$lambda_star,
    eigenvalues         = rate$eigenvalues,
    m_min               = m_min,
    m_convergence       = m_convergence,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    mode_gap_max        = if (run_pilot) mode_gap_max else NULL,
    m_pilot_from_gap    = if (run_pilot) m_pilot_from_gap else NULL,
    draw_engine         = "two_block_rNormal_reg_v5"
  )

  pilot           <- NULL
  pilot_mode_test <- NULL
  fixef_main_start <- fixef_start

  if (run_pilot) {
    if (verbose) {
      cat(sprintf(
        "--- glmerb [v5 / two_block_rNormal_reg_v5]: pilot stage (%d independent chains from ICM mode; m_convergence_pilot = %d) ---\n\n",
        n_pilot_int, m_convergence_pilot
      ))
    }

    pilot <- run_short_chains_v5(
      n_chains       = n_pilot_int,
      start_fixef    = fixef_start,
      inner_sweeps   = m_convergence_pilot,
      design         = design,
      block1_prior   = block1_prior,
      pfamily_list   = prior$pfamily_list,
      family         = family,
      re_names       = re_names,
      group_levels   = group_levels,
      seed_offset    = 0L,
      collect_block1 = TRUE,
      progbar        = progbar,
      stage_label    = "pilot",
      diag_sweeps    = diag_sweeps,
      fixef_mode     = fixef_mode_ref,
      b_mode         = b_mode_ref
    )

    X_pilot <- do.call(cbind, lapply(re_names, function(k) pilot$fixef_draws[[k]]))
    pnames  <- unlist(lapply(re_names, function(k) {
      paste0(k, "::", colnames(pilot$fixef_draws[[k]]))
    }))
    colnames(X_pilot) <- pnames
    mu_pilot  <- colMeans(X_pilot)
    mode_vec  <- unlist(lapply(re_names, function(k) fixef_start[[k]]))
    names(mode_vec) <- pnames
    d_pm   <- mu_pilot - mode_vec
    S_pilot <- stats::cov(X_pilot)
    p_dim2  <- ncol(X_pilot)
    S_inv   <- tryCatch(
      solve(S_pilot),
      error = function(e) {
        ridge <- 1e-8 * mean(diag(S_pilot))
        solve(S_pilot + diag(ridge, p_dim2))
      }
    )
    Q_pm <- as.numeric(n_pilot_int * t(d_pm) %*% S_inv %*% d_pm)
    p_pm <- stats::pchisq(Q_pm, df = p_dim2, lower.tail = FALSE)
    pilot_mode_test <- list(
      Q       = Q_pm,
      df      = p_dim2,
      p_value = p_pm,
      n_pilot = n_pilot_int
    )
    if (verbose) {
      cat(sprintf(
        "--- glmerb: pilot vs mode chi-squared test: p = %.4g (df = %d, n_pilot = %d) ---\n\n",
        p_pm, p_dim2, n_pilot_int
      ))
    }

    fixef_main_start <- lapply(pilot$fixef_draws, colMeans)

    .lmebayes_print_fixef_init(
      fixef_main_start, re_names, verbose,
      header = "--- glmerb: main-stage fixef.init (pilot colMeans) ---"
    )

    n_grp        <- nlevels(design$groups)
    re_col_names <- colnames(design$Z)
    grp_col_name <- design$group_name
    n_eigs       <- length(rate$eigenvalues)

    lambda_star_vec  <- numeric(n_pilot_int)
    max_eigenvalues  <- rep(-Inf, n_eigs)
    rate_upper       <- NULL
    lambda_star_best <- -Inf
    i_max_rate       <- NA_integer_

    for (i in seq_len(n_pilot_int)) {
      rows_i   <- ((i - 1L) * n_grp + 1L):(i * n_grp)
      block_df <- pilot$coefficients[rows_i, , drop = FALSE]
      if (!is.null(grp_col_name) && grp_col_name %in% colnames(block_df)) {
        ord <- match(group_levels, block_df[[grp_col_name]])
        b_i <- as.matrix(block_df[ord, re_col_names, drop = FALSE])
      } else {
        b_i <- as.matrix(block_df[, re_col_names, drop = FALSE])
      }
      rownames(b_i) <- group_levels
      mode_w_i <- glmbayesCore::two_block_mode_weights(
        x            = design$Z,
        block        = design$groups,
        b_mode       = b_i,
        family       = family,
        group_levels = group_levels
      )
      rate_i <- glmbayesCore::two_block_rate_v2(
        x                 = design$Z,
        block             = design$groups,
        x_hyper           = design$X_hyper,
        prior_list_block1 = block1_prior,
        pfamily_list      = prior$pfamily_list,
        weights           = mode_w_i$weights,
        family            = family,
        group_levels      = group_levels
      )
      lambda_star_vec[i] <- rate_i$lambda_star
      max_eigenvalues    <- pmax(max_eigenvalues, rate_i$eigenvalues)
      if (rate_i$lambda_star > lambda_star_best) {
        lambda_star_best <- rate_i$lambda_star
        rate_upper       <- rate_i
        i_max_rate       <- i
      }
    }

    rate_upper_eig             <- rate_upper
    rate_upper_eig$eigenvalues <- max_eigenvalues
    rate_upper_eig$lambda_star <- max_eigenvalues[1L]

    m_min_upper <- glmbayesCore::two_block_l_for_tv(
      rate_upper_eig, tv_tol, method = "theorem3"
    ) + 1L

    m_convergence_upper <- m_min_upper
    if (m_convergence_upper > m_convergence) {
      m_convergence <- m_convergence_upper
    }

    convergence_info$lambda_star_upper <- rate_upper_eig$lambda_star
    convergence_info$eigenvalues_upper <- max_eigenvalues
    convergence_info$m_min_upper       <- m_min_upper
    convergence_info$i_max_rate        <- i_max_rate
    convergence_info$lambda_star_vec   <- lambda_star_vec
    convergence_info$m_convergence     <- m_convergence

    if (verbose) {
      .fmt_eigs <- function(ev) paste(sprintf("%.4f", ev), collapse = ", ")
      cat(sprintf(
        "--- glmerb: post-pilot convergence bounds (%d pilot draws) ---\n    ML estimate (local-Gaussian at mode):    lambda* = %.4f, m_min = %d, eigenvalues = [%s]\n    Pilot upper bound (per-eig max, #%d/%d): lambda* = %.4f, m_min = %d, eigenvalues = [%s]\n    => using m_convergence = %d ---\n\n",
        n_pilot_int,
        rate$lambda_star,          m_min,       .fmt_eigs(rate$eigenvalues),
        i_max_rate, n_pilot_int,
        rate_upper_eig$lambda_star, m_min_upper, .fmt_eigs(max_eigenvalues),
        m_convergence
      ))
      cat(sprintf(
        "--- glmerb [v5 / two_block_rNormal_reg_v5]: pilot complete; main stage (%d independent chains from pilot mean; m_convergence = %d) ---\n\n",
        n, m_convergence
      ))
    }
  } else {
    if (verbose) {
      cat(sprintf(
        "--- glmerb [v5 / two_block_rNormal_reg_v5]: main stage (%d independent chains from ICM mode; m_convergence = %d) ---\n\n",
        n, m_convergence
      ))
    }
  }

  sampler <- run_short_chains_v5(
    n_chains       = n,
    start_fixef    = fixef_main_start,
    inner_sweeps   = m_convergence,
    design         = design,
    block1_prior   = block1_prior,
    pfamily_list   = prior$pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    seed_offset    = if (run_pilot) n_pilot_int else 0L,
    collect_block1 = collect_block1,
    progbar        = progbar,
    stage_label    = "main",
    diag_sweeps    = diag_sweeps,
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref
  )

  structure(
    list(
      call                   = cl,
      fixef_draws            = sampler$fixef_draws,
      coefficients           = sampler$coefficients,
      dispersion_fixef_draws = sampler$dispersion_fixef_draws,
      iters_fixef_draws      = sampler$iters_fixef_draws,
      iters_ranef_draws      = sampler$iters_ranef_draws,
      mu_all_last            = sampler$mu_all_last,
      sweep_history          = sampler$sweep_history,
      coef.mode              = fixef_start,
      ranef.mode             = ranef_mode,
      fixef_main_start       = fixef_main_start,
      pilot                  = pilot,
      pilot_mode_test        = pilot_mode_test,
      m_convergence_used     = m_convergence,
      n_pilot                = if (run_pilot) n_pilot_int else 0L,
      convergence_info       = convergence_info,
      Prior                  = list(block1_prior = block1_prior,
                                    pfamily_list = prior$pfamily_list),
      design                 = design
    ),
    class = c("rglmerb_v5", "list")
  )
}
