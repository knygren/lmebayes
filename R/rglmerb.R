#' Bayesian generalized linear mixed-effects model sampler (sweep-outer engine)
#'
#' Full sampling engine for non-Gaussian (and Gaussian) generalized linear
#' mixed models, parallel to \code{\link{rlmerb}} for the LMM path.
#' \code{rglmerb} wraps \code{\link{rGLMM}} after ICM mode computation
#' and \pkg{lmebayes} prior setup; see \code{\link{glmerb}} for the user API.
#'
#' @param n Integer. Number of independent chains in the main stage.
#' @param design A \code{\link{model_setup}} object.
#' @param prior A \code{lmebayes_prior_setup} object.
#' @param family A \code{\link[stats]{family}} object. Default \code{poisson()}.
#' @param fixef_start Optional named list of Block~2 starting vectors. When
#'   \code{NULL}, the ICM posterior mode is computed via
#'   \code{\link[glmbayesCore]{glmerb_posterior_mode}}.
#' @param m_convergence Optional integer inner Gibbs sweeps per main draw.
#' @param n_pilot Integer pilot chains; \code{0} or \code{NULL} skips pilot.
#' @param m_convergence_pilot Optional pilot inner sweeps.
#' @param tv_tol Total variation tolerance for convergence calibration.
#' @param mode_gap_max Pilot sweep calibration when \code{m_convergence_pilot}
#'   is \code{NULL}.
#' @param collect_block1 Collect Block~1 \code{coefficients} from main chains.
#' @param verbose Print stage headers and diagnostics.
#' @param progbar Progress bars when \code{verbose} is \code{FALSE}.
#' @return Object of class \code{c("rglmerb", "list")} with Block~2 fields in
#'   the \code{fixef.*} namespace (as \code{\link{rGLMM}}), plus
#'   \code{ranef.mode}, \code{Prior}, and \code{design}.
#' @seealso \code{\link{glmerb}}, \code{\link{rGLMM}}, \code{\link{rLMM}},
#'   \code{\link[glmbayesCore]{run_sweep_outer_chains_v6}}
#' @title The Bayesian Generalized Linear Mixed-Effects Model Distribution
#' @export
rglmerb <- function(
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

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

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

  block1_prior <- .lmebayes_block1_prior_list(prior)

  if (verbose && n_pilot_int > 0L) {
    cat(sprintf(
      "--- glmerb [rglmerb / sweep-outer]: pilot stage (%d independent chains from ICM mode) ---\n\n",
      n_pilot_int
    ))
  } else if (verbose) {
    cat(sprintf(
      "--- glmerb [rglmerb / sweep-outer]: main stage (%d independent chains from ICM mode) ---\n\n",
      n
    ))
  }

  out <- glmbayesCore::rGLMM(
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
    m_convergence_pilot = m_convergence_pilot,
    tv_tol              = tv_tol,
    mode_gap_max        = mode_gap_max,
    verbose             = verbose,
    progbar             = progbar,
    stage_verbose       = verbose,
    b_start             = ranef_mode,
    collect_block1      = collect_block1,
    any_ing             = isTRUE(prior$any_ing)
  )

  if (n_pilot_int > 0L) {
    .lmebayes_print_fixef_init(out$fixef.init, re_names, verbose)
  }

  out <- .lmebayes_add_fixef_summaries(out)
  out$call       <- cl
  out$ranef.mode <- ranef_mode
  out$convergence <- out$convergence_info
  out$Prior      <- list(
    block1_prior   = block1_prior,
    pfamily_list   = prior$pfamily_list
  )
  out$design     <- design

  class(out) <- c("rglmerb", "list")
  out
}
