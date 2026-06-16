#' Experimental \code{rGLMM}-based two-stage Gibbs sampler
#'
#' Batch pilot and main sampling via a single call to
#' \code{\link[glmbayesCore]{rGLMM}}.  Kept for development and comparison;
#' \code{\link{rglmerb}} uses independent short chains via
#' \code{\link{run_short_chains}} instead.
#'
#' @inheritParams rglmerb
#' @return An object of class \code{c("rglmerb_experimental", "list")} with the
#'   same components as \code{\link{rglmerb}}.
#' @seealso \code{\link{rglmerb}}, \code{\link{glmerb}}, \code{\link{rlmerb}},
#'   \code{\link[glmbayesCore]{rGLMM}}
#' @export
rglmerb_experimental <- function(
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
    seed                = NULL,
    collect_block1      = TRUE,
    verbose             = TRUE,
    progbar             = FALSE
) {
  cl <- match.call()

  # ---- argument validation --------------------------------------------------
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

  # ---- ICM posterior mode ---------------------------------------------------
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

  block1_prior <- .lmebayes_block1_prior_list(prior)

  # ---- convergence calibration ----------------------------------------------
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

  if (is.null(m_convergence)) {
    m_convergence <- m_min
  } else if (m_convergence < m_min) {
    warning(
      "rglmerb_experimental: m_convergence = ", m_convergence,
      " is below the derived minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
      "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }

  p_dim            <- sum(vapply(fixef_start, length, integer(1L)))
  D_max            <- if (!is.null(mode_gap_max)) sqrt(p_dim) * mode_gap_max else 0
  m_pilot_from_gap <- NULL

  if (run_pilot && is.null(m_convergence_pilot)) {
    erf1_inv_tv      <- stats::qnorm((tv_tol + 1) / 2) / sqrt(2)
    c_tol            <- erf1_inv_tv * 2 * sqrt(2)
    m_pilot_from_gap <- if (D_max <= c_tol || rate$lambda_star <= 0) {
      m_min
    } else {
      as.integer(ceiling(log(D_max / c_tol) / log(1 / rate$lambda_star)))
    }
    m_convergence_pilot <- max(m_min, m_pilot_from_gap)
  }

  calib_label <- if (is_gaussian) {
    "exact (Gaussian posterior)"
  } else {
    sprintf("approximate (local-Gaussian at mode, %s)", family$family)
  }
  if (prior$any_ing) {
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
  if (prior$any_ing) method_label <- paste0(method_label, "+disp_lower_bound")

  convergence_info <- list(
    method              = method_label,
    tv_tol              = tv_tol,
    lambda_star         = rate$lambda_star,
    eigenvalues         = rate$eigenvalues,
    m_min               = m_min,
    m_convergence       = m_convergence,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    mode_gap_max        = if (run_pilot) mode_gap_max else NULL,
    m_pilot_from_gap    = if (run_pilot) m_pilot_from_gap else NULL,
    draw_engine         = "rGLMM"
  )

  if (verbose) {
    if (run_pilot) {
      cat(sprintf(
        "--- glmerb: pilot + main via rGLMM (%d pilot chains, %d main chains; m_pilot = %d, m_main = %d) ---\n\n",
        n_pilot_int, n, m_convergence_pilot, m_convergence
      ))
    } else {
      cat(sprintf(
        "--- glmerb: main stage via rGLMM (%d chains from ICM mode; m_convergence = %d) ---\n\n",
        n, m_convergence
      ))
    }
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  res <- glmbayesCore::rGLMM(
    n                   = n,
    y                   = design$y,
    x                   = design$Z,
    block               = design$groups,
    x_hyper             = design$X_hyper,
    prior_list          = block1_prior,
    pfamily_list        = prior$pfamily_list,
    start               = fixef_start,
    family              = family,
    m_convergence       = m_convergence,
    re_coef_names       = re_names,
    group_levels        = group_levels,
    group_name          = design$group_name,
    n_pilot             = n_pilot_int,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    tv_tol              = if (run_pilot) tv_tol else NULL,
    mode_gap_max        = mode_gap_max,
    verbose             = FALSE,
    progbar             = isTRUE(progbar),
    stage_verbose       = verbose,
    rate_calibration    = if (run_pilot) {
      list(
        lambda_star = rate$lambda_star,
        eigenvalues = rate$eigenvalues,
        m_min       = m_min
      )
    } else {
      NULL
    }
  )

  m_convergence_used <- res$m_convergence
  convergence_info$m_convergence <- m_convergence_used
  fixef_main_start   <- res$fixef.init
  pilot_mode_test    <- res$pilot_chisq
  pilot              <- if (run_pilot) {
    .rglmerb_experimental_sampler_from_rGLMM(res$pilot, collect_block1 = TRUE)
  } else {
    NULL
  }

  if (run_pilot && !is.null(res$pilot_ub)) {
    ub <- res$pilot_ub
    convergence_info$lambda_star_upper <- ub$rate_upper$lambda_star
    convergence_info$eigenvalues_upper <- ub$max_eigenvalues
    convergence_info$m_min_upper       <- ub$m_min_upper
    convergence_info$i_max_rate        <- ub$i_max_rate
    convergence_info$lambda_star_vec   <- ub$lambda_star_vec
  }

  sampler <- .rglmerb_experimental_sampler_from_rGLMM(res, collect_block1 = collect_block1)

  structure(
    list(
      call                   = cl,
      fixef_draws            = sampler$fixef_draws,
      coefficients           = sampler$coefficients,
      dispersion_fixef_draws = sampler$dispersion_fixef_draws,
      iters_fixef_draws      = sampler$iters_fixef_draws,
      mu_all_last            = sampler$mu_all_last,
      coef.mode              = fixef_start,
      ranef.mode             = ranef_mode,
      fixef_main_start       = fixef_main_start,
      pilot                  = pilot,
      pilot_mode_test        = pilot_mode_test,
      m_convergence_used     = m_convergence_used,
      convergence_info       = convergence_info,
      Prior                  = list(block1_prior = block1_prior,
                                    pfamily_list = prior$pfamily_list),
      design                 = design
    ),
    class = c("rglmerb_experimental", "list")
  )
}

#' Map an \code{rGLMM} result to legacy sampler fields
#' @noRd
.rglmerb_experimental_sampler_from_rGLMM <- function(res, collect_block1 = TRUE) {
  list(
    fixef_draws            = res$fixef,
    dispersion_fixef_draws = res$fixef.dispersion,
    iters_fixef_draws      = res$fixef.iters,
    coefficients           = if (collect_block1) res$coefficients else NULL,
    mu_all_last            = res$fixef.mu
  )
}
