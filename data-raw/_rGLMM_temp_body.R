#' Replicate-chain Gibbs sampling for Bayesian GLMMs (v6 sweep-outer driver)
#'
#' Temporary low-level GLMM sampler in \pkg{lmebayes}, parallel to
#' \code{\link[glmbayesCore]{rGLMM}} but using
#' \code{\link[glmbayesCore]{run_sweep_outer_chains_v6}} instead of the C++ v2
#' engine.  See \code{\link{rglmerb}} for the \pkg{lmebayes} wrapper.
#'
#' @inheritParams rGLMM
#' @param m_convergence Inner Gibbs steps per main-stage stored draw. When
#'   \code{NULL} and \code{tv_tol} is set, derived from Theorem~3 at \code{start};
#'   when \code{NULL} and \code{tv_tol} is \code{NULL}, defaults to \code{10L}.
#' @param b_start Optional \code{J x p_re} Block~1 mode matrix for
#'   \code{two_block_mode_weights} and v6 batch init.  Required for
#'   non-Gaussian families.
#' @param collect_block1 Logical. Collect Block~1 \code{coefficients} from each
#'   chain (needed for post-pilot eigenvalue upper bounds).
#' @param any_ing Logical. Label convergence calibration as ING-conservative.
#' @return Object of class \code{c("rGLMM_temp", "list")} with \code{fixef.*}
#'   fields (as \code{rGLMM}), plus \code{draw_engine} and \code{convergence_info}.
#' @family simfuncs
#' @seealso \code{\link[glmbayesCore]{rGLMM}}, \code{\link{rglmerb}},
#'   \code{\link[glmbayesCore]{run_sweep_outer_chains_v6}}
#' @export
rGLMM_temp <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    start,
    offset              = NULL,
    weights             = 1,
    family              = gaussian(),
    m_convergence       = NULL,
    re_coef_names       = colnames(x),
    group_levels        = levels(block),
    group_name          = NULL,
    n_pilot             = 0L,
    m_convergence_pilot = NULL,
    tv_tol              = NULL,
    mode_gap_max        = 1.0,
    Gridtype            = 2,
    n_envopt            = NULL,
    use_parallel        = TRUE,
    use_opencl          = FALSE,
    verbose             = FALSE,
    progbar             = FALSE,
    stage_verbose       = FALSE,
    rate_calibration    = NULL,
    b_start             = NULL,
    collect_block1      = TRUE,
    any_ing             = FALSE
) {
  cl <- match.call()

  family <- glmbayesCore:::.two_block_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  n_pilot <- as.integer(n_pilot[1L])
  if (n_pilot < 0L) stop("'n_pilot' must be non-negative.", call. = FALSE)

  if (!is.null(m_convergence)) {
    m_convergence <- as.integer(m_convergence[1L])
    if (m_convergence < 1L) {
      stop("'m_convergence' must be at least 1.", call. = FALSE)
    }
  }

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be a single positive finite number.",
           call. = FALSE)
    }
  }

  run_pilot <- n_pilot > 0L
  run_ub    <- run_pilot && !is.null(tv_tol)

  if (!is.null(tv_tol)) {
    if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
        !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
      stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
    }
  }

  if (run_pilot && is.null(m_convergence_pilot)) {
    m_convergence_pilot <- if (!is.null(m_convergence)) {
      m_convergence
    } else if (!is.null(tv_tol)) {
      NULL
    } else {
      10L
    }
  } else if (run_pilot) {
    m_convergence_pilot <- as.integer(m_convergence_pilot[1L])
    if (m_convergence_pilot < 1L) {
      stop("'m_convergence_pilot' must be at least 1 when n_pilot > 0.",
           call. = FALSE)
    }
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) {
        paste0("RE", seq_len(ncol(x)))
      } else {
        cn
      }
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop("length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
         call. = FALSE)
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  pfamily_list <- glmbayesCore:::.two_block_validate_pfamily_list(
    pfamily_list, re_names, J = length(group_levels)
  )

  if (!is.list(start) || is.null(names(start))) {
    stop("'start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(start), re_names)) {
    stop("names(start) must match re_coef_names.", call. = FALSE)
  }
  start <- start[re_names]
  fixef_mode <- start

  if (!is_gaussian && is.null(b_start)) {
    stop("'b_start' is required for non-Gaussian families.", call. = FALSE)
  }

  glmbayesCore:::.two_block_validate_block1_prior(
    prior_list, family = family
  )

  ptypes <- vapply(pfamily_list, function(pf) pf$pfamily, character(1))
  names(ptypes) <- re_names

  design <- list(
    y             = y,
    Z             = x,
    groups        = block,
    X_hyper       = x_hyper,
    re_coef_names = re_names,
    group_name    = group_name
  )

  fixef_mode_ref <- fixef_mode
  b_mode_ref     <- b_start
  diag_sweeps    <- isTRUE(verbose) || isTRUE(stage_verbose)
  progbar_use    <- isTRUE(progbar) && !diag_sweeps

  rate <- .rGLMM_temp_rate_at_mode(
    design       = design,
    prior_list   = prior_list,
    pfamily_list = pfamily_list,
    family       = family,
    b_mode       = b_start,
    group_levels = group_levels,
    is_gaussian  = is_gaussian
  )

  m_min <- NULL
  if (!is.null(tv_tol)) {
    m_min <- glmbayesCore::two_block_l_for_tv(
      rate, tv_tol, method = "theorem3"
    ) + 1L
    if (is.null(m_convergence)) {
      m_convergence <- m_min
    } else if (m_convergence < m_min) {
      warning(
        "rGLMM_temp: m_convergence = ", m_convergence, " is below the derived ",
        "minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
        "; using m_min instead.",
        call. = FALSE
      )
      m_convergence <- m_min
    }
  } else if (is.null(m_convergence)) {
    m_convergence <- 10L
  }

  p_dim            <- sum(vapply(fixef_mode, length, integer(1L)))
  D_max            <- if (!is.null(mode_gap_max)) sqrt(p_dim) * mode_gap_max else 0
  m_pilot_from_gap <- NULL

  if (run_pilot && is.null(m_convergence_pilot) && !is.null(tv_tol)) {
    erf1_inv_tv <- stats::qnorm((tv_tol + 1) / 2) / sqrt(2)
    c_tol       <- erf1_inv_tv * 2 * sqrt(2)
    m_pilot_from_gap <- if (D_max <= c_tol || rate$lambda_star <= 0) {
      m_min
    } else {
      as.integer(ceiling(log(D_max / c_tol) / log(1 / rate$lambda_star)))
    }
    m_convergence_pilot <- max(m_min, m_pilot_from_gap)
  }

  if (is.null(rate_calibration) && !is.null(tv_tol)) {
    rate_calibration <- list(
      lambda_star = rate$lambda_star,
      eigenvalues = rate$eigenvalues,
      m_min       = m_min
    )
  }

  calib_label <- if (is_gaussian) {
    "exact (Gaussian posterior)"
  } else {
    sprintf("approximate (local-Gaussian at mode, %s)", family$family)
  }
  if (isTRUE(any_ing)) {
    calib_label <- paste0(calib_label, "; conservative: ING tau^2_k = disp_lower")
  }

  if (isTRUE(verbose) && !is.null(tv_tol)) {
    cat(sprintf(
      "--- rGLMM_temp: convergence calibration [%s]:\n    lambda* = %.4f, tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n",
      calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
    if (run_pilot && !is.null(mode_gap_max) && !is.null(m_pilot_from_gap)) {
      cat(sprintf(
        "--- rGLMM_temp: pilot sweep calibration [mode_gap_max = %g SD/dim, p = %d, D_max = %.4f]:\n    m_min = %d, lambda* = %.4f => m_convergence_pilot = %d ---\n\n",
        mode_gap_max, p_dim, D_max, m_min, rate$lambda_star, m_convergence_pilot
      ))
    }
  }

  method_label <- if (is_gaussian) "exact" else "local_gaussian_mode"
  if (isTRUE(any_ing)) method_label <- paste0(method_label, "+disp_lower_bound")

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
    draw_engine         = "run_sweep_outer_chains_v6"
  )

  m_convergence_used <- m_convergence
  fixef_init         <- fixef_mode
  pilot_res          <- NULL
  pilot_chisq        <- NULL
  pilot_ub           <- NULL

  run_sweep_stage <- function(n_chains, start_fixef, inner_sweeps, stage_label) {
    glmbayesCore::run_sweep_outer_chains_v6(
      n_chains       = n_chains,
      start_fixef    = start_fixef,
      inner_sweeps   = inner_sweeps,
      design         = design,
      block1_prior   = prior_list,
      pfamily_list   = pfamily_list,
      family         = family,
      re_names       = re_names,
      group_levels   = group_levels,
      collect_block1 = collect_block1,
      progbar        = progbar_use,
      stage_label    = stage_label,
      diag_sweeps    = diag_sweeps,
      fixef_mode     = fixef_mode_ref,
      b_mode         = b_mode_ref,
      b_start        = b_mode_ref,
      ptypes         = ptypes
    )
  }

  if (run_pilot) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "--- rGLMM_temp [sweep-outer]: pilot stage (%d chains; m_convergence_pilot = %d) ---\n\n",
        n_pilot, m_convergence_pilot
      ))
    }

    pilot_raw <- run_sweep_stage(
      n_chains     = n_pilot,
      start_fixef  = fixef_mode,
      inner_sweeps = m_convergence_pilot,
      stage_label  = "pilot"
    )

    pilot_chisq <- glmbayesCore:::.two_block_pilot_chisq_test(
      fixef_draws = pilot_raw$fixef_draws,
      re_names    = re_names,
      fixef_mode  = fixef_mode,
      n_pilot     = n_pilot
    )

    if (isTRUE(stage_verbose) || isTRUE(verbose)) {
      cat(sprintf(
        "--- rGLMM_temp: pilot vs mode chi-squared test: p = %.4g (df = %d, n_pilot = %d) ---\n\n",
        pilot_chisq$p_value, pilot_chisq$df, pilot_chisq$n_pilot
      ))
    }

    fixef_init <- glmbayesCore:::.two_block_fixef_colmeans(
      pilot_raw$fixef_draws, re_names, fixef_mode
    )

    if (run_ub) {
      pilot_ub <- .rGLMM_temp_pilot_ub_from_coefficients(
        pilot_coefficients = pilot_raw$coefficients,
        n_pilot            = n_pilot,
        re_names           = re_names,
        group_levels       = group_levels,
        group_name         = group_name,
        x                  = x,
        block              = block,
        x_hyper            = x_hyper,
        prior_list         = prior_list,
        pfamily_list       = pfamily_list,
        family             = family,
        tv_tol             = tv_tol
      )
      if (pilot_ub$m_min_upper > m_convergence_used) {
        m_convergence_used <- pilot_ub$m_min_upper
      }
      convergence_info$lambda_star_upper <- pilot_ub$rate_upper$lambda_star
      convergence_info$eigenvalues_upper <- pilot_ub$max_eigenvalues
      convergence_info$m_min_upper       <- pilot_ub$m_min_upper
      convergence_info$i_max_rate        <- pilot_ub$i_max_rate
      convergence_info$lambda_star_vec   <- pilot_ub$lambda_star_vec
      convergence_info$m_convergence     <- m_convergence_used
    }

    if (isTRUE(stage_verbose) && run_ub) {
      glmbayesCore:::.rGLMM_print_pilot_stage_diagnostics(
        n_pilot            = n_pilot,
        n_main             = n,
        pilot_chisq        = pilot_chisq,
        pilot_ub           = pilot_ub,
        rate_calibration   = rate_calibration,
        m_convergence_used = m_convergence_used
      )
    } else if (isTRUE(verbose)) {
      cat(sprintf(
        "--- rGLMM_temp [sweep-outer]: pilot complete; main stage (%d chains; m_convergence = %d) ---\n\n",
        n, m_convergence_used
      ))
    }

    pilot_res <- .rGLMM_temp_format_v6_out(
      v6_out       = pilot_raw,
      n            = n_pilot,
      re_names     = re_names,
      group_levels = group_levels,
      fixef_mode   = fixef_mode,
      fixef_init   = fixef_mode
    )
  } else if (isTRUE(verbose)) {
    cat(sprintf(
      "--- rGLMM_temp [sweep-outer]: main stage (%d chains; m_convergence = %d) ---\n\n",
      n, m_convergence_used
    ))
  }

  main_raw <- run_sweep_stage(
    n_chains     = n,
    start_fixef  = fixef_init,
    inner_sweeps = m_convergence_used,
    stage_label  = "main"
  )

  draw_engine_args <- list(
    n_chains       = n,
    start_fixef    = fixef_init,
    inner_sweeps   = m_convergence_used,
    design         = design,
    block1_prior   = prior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    stage_label    = "main",
    diag_sweeps    = diag_sweeps,
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes
  )

  main_res <- .rGLMM_temp_format_v6_out(
    v6_out       = main_raw,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )

  main_res$call                <- cl
  main_res$n_pilot             <- n_pilot
  main_res$m_convergence       <- m_convergence_used
  main_res$m_convergence_pilot <- if (run_pilot) m_convergence_pilot else NULL
  main_res$convergence_info    <- convergence_info
  main_res$draw_engine         <- "run_sweep_outer_chains_v6"
  main_res$draw_engine_call    <- quote(glmbayesCore::run_sweep_outer_chains_v6)
  main_res$draw_engine_args    <- draw_engine_args
  main_res$pfamily_list        <- pfamily_list
  main_res$family              <- family
  main_res$prior_list          <- prior_list

  if (run_pilot) {
    main_res$pilot       <- pilot_res
    main_res$pilot_chisq <- pilot_chisq
  }
  if (run_ub) {
    main_res$pilot_ub <- pilot_ub
    main_res$tv_tol   <- tv_tol
  }

  class(main_res) <- c("rGLMM_temp", "list")
  main_res
}

#' Local-Gaussian rate at the ICM mode
#' @noRd
.rGLMM_temp_rate_at_mode <- function(
    design,
    prior_list,
    pfamily_list,
    family,
    b_mode,
    group_levels,
    is_gaussian
) {
  if (is_gaussian) {
    glmbayesCore::two_block_rate_v2(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      family            = gaussian(),
      group_levels      = group_levels
    )
  } else {
    mode_w <- glmbayesCore::two_block_mode_weights(
      x            = design$Z,
      block        = design$groups,
      b_mode       = b_mode,
      family       = family,
      group_levels = group_levels
    )
    glmbayesCore::two_block_rate_v2(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      weights           = mode_w$weights,
      family            = family,
      group_levels      = group_levels
    )
  }
}

#' Format v6 batch output for staged \code{fixef.*} naming
#' @noRd
.rGLMM_temp_format_v6_out <- function(
    v6_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  x <- list(
    fixef_draws            = v6_out$fixef_draws,
    coefficients           = v6_out$coefficients,
    dispersion_fixef_draws = v6_out$dispersion_fixef_draws,
    iters_fixef_draws      = v6_out$iters_fixef_draws,
    mu_all_last            = v6_out$mu_all_last,
    re_coef_names          = re_names,
    group_levels           = group_levels,
    n                      = n
  )
  glmbayesCore:::.two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}

#' Post-pilot eigenvalue UB from v6 stacked \code{coefficients}
#' @noRd
.rGLMM_temp_pilot_ub_from_coefficients <- function(
    pilot_coefficients,
    n_pilot,
    re_names,
    group_levels,
    group_name,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    family,
    tv_tol
) {
  n_grp        <- length(group_levels)
  re_col_names <- re_names
  grp_col_name <- group_name

  lambda_star_vec  <- numeric(n_pilot)
  max_eigenvalues  <- NULL
  rate_upper       <- NULL
  lambda_star_best <- -Inf
  i_max_rate       <- NA_integer_

  for (i in seq_len(n_pilot)) {
    rows_i   <- ((i - 1L) * n_grp + 1L):(i * n_grp)
    block_df <- pilot_coefficients[rows_i, , drop = FALSE]
    if (!is.null(grp_col_name) && grp_col_name %in% colnames(block_df)) {
      ord <- match(group_levels, block_df[[grp_col_name]])
      b_i <- as.matrix(block_df[ord, re_col_names, drop = FALSE])
    } else {
      b_i <- as.matrix(block_df[, re_col_names, drop = FALSE])
    }
    rownames(b_i) <- group_levels
    mode_w_i <- glmbayesCore::two_block_mode_weights(
      x            = x,
      block        = block,
      b_mode       = b_i,
      family       = family,
      group_levels = group_levels
    )
    rate_i <- glmbayesCore::two_block_rate_v2(
      x                 = x,
      block             = block,
      x_hyper           = x_hyper,
      prior_list_block1 = prior_list,
      pfamily_list      = pfamily_list,
      weights           = mode_w_i$weights,
      family            = family,
      group_levels      = group_levels
    )
    if (is.null(max_eigenvalues)) {
      max_eigenvalues <- rep(-Inf, length(rate_i$eigenvalues))
    }
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

  list(
    rate_upper      = rate_upper_eig,
    m_min_upper     = m_min_upper,
    lambda_star_vec = lambda_star_vec,
    i_max_rate      = i_max_rate,
    max_eigenvalues = max_eigenvalues
  )
}