#' Bayesian generalized linear mixed-effects model sampler
#'
#' Matrix-level sampler for \pkg{lmebayes} \code{model_setup} objects and prior
#' containers.  Routes by response family:
#' \itemize{
#'   \item \code{family = gaussian()} delegates to
#'     \code{\link[glmbayesCore]{rLMMNormal_reg}} or
#'     \code{\link[glmbayesCore]{rLMMindepNormalGamma_reg}} when \code{dispersion_ranef} is
#'     a \code{dGamma()} pfamily.
#'   \item Non-Gaussian families delegate to \code{\link[glmbayesCore]{rGLMM}}
#'     (sweep-outer engine with optional pilot/main staging).
#' }
#' See \code{\link{glmerb}} for the formula-level API.
#'
#' @param n Integer. Number of independent chains in the main stage.
#' @param design A \code{\link{model_setup}} object.
#' @param prior A \code{lmebayes_prior_setup} object.
#' @param family A \code{\link[stats]{family}} object. Default \code{poisson()}.
#' @param dispersion_ranef Observation-level measurement dispersion \eqn{\sigma^2}:
#'   required positive scalar for \code{family = gaussian()}, or a
#'   \code{\link{dGamma}()} pfamily with \code{Inv_Dispersion = TRUE}; must be
#'   \code{NULL} (default) for \code{poisson()} and \code{binomial()}.
#' @param fixef_start Optional named list of Block~2 starting vectors. When
#'   \code{NULL}, the ICM start is computed inside the Core engine
#'   (\code{\link[glmbayesCore]{lmerb_posterior_mean}} for Gaussian,
#'   \code{\link[glmbayesCore]{glmerb_posterior_mode}} otherwise).
#' @param m_convergence Optional integer inner Gibbs sweeps per main draw.
#' @param n_pilot Optional integer pilot chains (non-Gaussian only);
#'   \code{NULL} (default) uses cost-optimal count when \code{tv_tol} is set,
#'   else legacy \code{gap_tol}; \code{0L} skips pilot.
#' @param gap_tol Legacy mode--mean gap for deriving \code{n_pilot} when
#'   \code{n_pilot} is \code{NULL} and \code{tv_tol} is \code{NULL}. Set
#'   \code{NULL} to skip pilot unless \code{n_pilot} is explicit or
#'   \code{tv_tol} is set. Ignored for Gaussian.
#' @param m_convergence_pilot Optional pilot inner sweeps (non-Gaussian only).
#' @param tv_tol Total variation tolerance for convergence calibration.
#' @param mode_gap_max Pilot sweep calibration when \code{m_convergence_pilot}
#'   is \code{NULL} (non-Gaussian only).
#' @param collect_block1 Collect Block~1 \code{coefficients} from main chains
#'   (non-Gaussian only).
#' @param seed Optional RNG seed (Gaussian path only).
#' @param verbose Print stage headers and diagnostics.
#' @param progbar Progress bars when \code{verbose} is \code{FALSE}.
#' @return Object of class \code{c("rglmerb", "list")} with Block~2 fields in
#'   the \code{fixef.*} namespace, plus \code{ranef.mode}, \code{Prior},
#'   \code{design}, and \code{family}.  Non-Gaussian fits may include
#'   \code{n_pilot}, \code{pilot}, and \code{pilot_chisq}; Gaussian fits set
#'   \code{n_pilot = 0L} and omit pilot output.
#' @seealso \code{\link{glmerb}},
#'   \code{\link[glmbayesCore]{rLMMNormal_reg}},
#'   \code{\link[glmbayesCore]{rLMMindepNormalGamma_reg}},
#'   \code{\link[glmbayesCore]{rGLMM}}
#' @title The Bayesian Generalized Linear Mixed-Effects Model Distribution
#' @export
rglmerb <- function(
    n,
    design,
    prior,
    family              = poisson(),
    dispersion_ranef    = NULL,
    fixef_start         = NULL,
    m_convergence       = NULL,
    n_pilot             = NULL,
    gap_tol             = 0.0196,
    m_convergence_pilot = NULL,
    tv_tol              = 0.01,
    mode_gap_max        = 1.0,
    collect_block1      = TRUE,
    seed                = NULL,
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

  if (!inherits(family, "family") || is.null(family$family)) {
    stop("'family' must be a family object.", call. = FALSE)
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

  is_gaussian <- identical(family$family, "gaussian")

  disp_info <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = design,
    fn_name          = "rglmerb"
  )

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

  if (is_gaussian) {
    block1_prior <- .lmebayes_block1_prior_list(
      prior,
      dispersion_ranef = disp_info$dispersion_fix
    )

    out <- .lmebayes_run_lmm_engine(
      n             = n,
      design        = design,
      prior         = prior,
      disp_info     = disp_info,
      fixef_start   = fixef_start,
      m_convergence = m_convergence,
      tv_tol        = tv_tol,
      seed          = seed,
      progbar       = progbar,
      verbose       = verbose
    )

    if (is.null(fixef_start)) {
      .lmebayes_print_icm_fixef_table(
        prior_list = prior$prior_list,
        re_names   = re_names,
        fixef_icm  = out$fixef.mode,
        icm_info   = out$icm_info,
        ref_label  = "glmer (start)",
        icm_label  = "post mean (ICM)",
        header     = "--- glmerb: Block 2 fixed effects ---",
        verbose    = verbose
      )
    }

    out <- .lmebayes_add_fixef_summaries(out)
    out$call        <- cl
    out$convergence <- out$convergence_info
    out$Prior       <- list(
      block1_prior          = block1_prior,
      pfamily_list          = prior$pfamily_list,
      dispersion_ranef      = disp_info$dispersion_fix,
      dispersion_mode       = disp_info$mode,
      dispersion_pfamily    = disp_info$dispersion_pfamily,
      dispersion_prior_list = disp_info$dispersion_prior_list
    )
    out$design      <- design
    out$family      <- family
    out$n_pilot     <- 0L
    out$pilot_chisq <- NULL
    class(out)      <- c("rglmerb", "list")
    return(out)
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

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  block1_prior <- .lmebayes_block1_prior_list(prior, dispersion_ranef = NULL)

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
    n_pilot             = n_pilot,
    gap_tol             = gap_tol,
    m_convergence_pilot = m_convergence_pilot,
    tv_tol              = tv_tol,
    mode_gap_max        = mode_gap_max,
    verbose             = verbose,
    progbar             = progbar,
    stage_verbose       = verbose,
    b_start             = NULL,
    collect_block1      = collect_block1
  )

  if (is.null(fixef_start)) {
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$prior_list,
      re_names   = re_names,
      fixef_icm  = out$fixef.mode,
      icm_info   = out$icm_info,
      ref_label  = "glmer (start)",
      icm_label  = "post mode (ICM)",
      header     = "--- glmerb: Block 2 fixed effects ---",
      verbose    = verbose
    )
  }

  .lmebayes_print_ranef_mode_reference(
    out$ranef.mode, re_names, group_levels, verbose
  )

  if (!is.null(out$n_pilot) && out$n_pilot > 0L) {
    .lmebayes_print_fixef_init(out$fixef.init, re_names, verbose)
  }

  out <- .lmebayes_add_fixef_summaries(out)
  out$call        <- cl
  out$convergence <- out$convergence_info
  out$Prior       <- list(
    block1_prior          = block1_prior,
    pfamily_list          = prior$pfamily_list,
    dispersion_ranef      = disp_info$dispersion_fix,
    dispersion_mode       = disp_info$mode,
    dispersion_pfamily    = disp_info$dispersion_pfamily,
    dispersion_prior_list = disp_info$dispersion_prior_list
  )
  out$design      <- design
  out$family      <- family

  class(out) <- c("rglmerb", "list")
  out
}
