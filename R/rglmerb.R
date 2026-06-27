#' Bayesian generalized linear mixed-effects model sampler
#'
#' Matrix-level sampler for \pkg{lmebayes} \code{model_setup} objects and prior
#' containers.  Routes by response family:
#' \itemize{
#'   \item \code{family = gaussian()} delegates to
#'     \code{\link[glmbayesCore]{rLMMNormal_reg}} or
#'     \code{\link[glmbayesCore]{rLMMindepNormalGamma_reg}} when \code{dispersion_ranef} is
#'     a \code{dGamma()} pfamily.
#'   \item Non-Gaussian families delegate to \code{\link{rglmerb_v5}}
#'     (C++ sweep-outer \code{two_block_rNormal_reg_v5}).  The inactive R
#'     sweep-outer path is \code{\link{.rglmerb_v6_rGLMM}}.
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
#' @param gap_tol Legacy mode--mean gap for deriving the pilot chain count when
#'   \code{tv_tol} is \code{NULL}. Set \code{NULL} to skip pilot unless
#'   \code{tv_tol} is set. Ignored for Gaussian without ING Block~2 components.
#' @param tv_tol Total variation tolerance for convergence calibration.
#' @param mode_gap_max Pilot inner-sweep calibration (non-Gaussian and
#'   Gaussian+ING only).
#' @param collect_block1 Collect Block~1 \code{coefficients} from main chains
#'   (non-Gaussian only).
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
    gap_tol             = 0.0196,
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
      progbar       = progbar,
      verbose       = verbose,
      gap_tol       = gap_tol,
      mode_gap_max  = mode_gap_max
    )

    if (is.null(fixef_start)) {
      icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
      .lmebayes_print_icm_fixef_table(
        prior_list = prior$prior_list,
        re_names   = re_names,
        fixef_icm  = out$fixef.mode,
        icm_info   = out$icm_info,
        ref_label  = icm_lbl$ref_label,
        icm_label  = icm_lbl$icm_label,
        conv_label = icm_lbl$conv_label,
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

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  block1_prior <- .lmebayes_block1_prior_list(prior, dispersion_ranef = NULL)

  v5 <- rglmerb_v5(
    n                   = n,
    design              = design,
    prior               = prior,
    family              = family,
    fixef_start         = fixef_start,
    m_convergence       = m_convergence,
    gap_tol             = gap_tol,
    tv_tol              = tv_tol,
    mode_gap_max        = mode_gap_max,
    collect_block1      = collect_block1,
    verbose             = verbose,
    progbar             = progbar
  )

  out <- .lmebayes_pack_rglmerb_v5(v5, re_names, group_levels)
  out$icm_info <- NULL

  if (is.null(fixef_start)) {
    icm_lbl <- .lmebayes_block2_icm_labels(prior, family)
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$prior_list,
      re_names   = re_names,
      fixef_icm  = out$fixef.mode,
      icm_info   = out$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm_lbl$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- glmerb: Block 2 fixed effects ---",
      verbose    = verbose
    )
  }

  .lmebayes_print_ranef_mode_reference(
    out$ranef.mode, re_names, group_levels, verbose
  )

  if (!is.null(out$n_pilot) && out$n_pilot > 0L) {
    .lmebayes_print_fixef_init(
      out$fixef.init,
      re_names,
      verbose,
      header = "--- glmerb: main-stage fixef.init (pilot colMeans) ---"
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

  class(out) <- c("rglmerb", "list")
  out
}
