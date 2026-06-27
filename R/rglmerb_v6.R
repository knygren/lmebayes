#' Inactive R sweep-outer GLMM sampler (v6 driver via \code{rGLMM})
#'
#' Legacy development path: pilot and main stages call
#' \code{\link[glmbayesCore]{run_sweep_outer_chains_v6}} through
#' \code{\link[glmbayesCore]{rGLMM}}.  Retained for comparison and debugging;
#' \code{\link{rglmerb}} routes non-Gaussian sampling to
#' \code{\link{rglmerb_v5}} (C++ sweep-outer) instead.
#'
#' @inheritParams rglmerb
#' @return Object compatible with \code{\link{rglmerb}} staging fields.
#' @keywords internal
#' @noRd
.rglmerb_v6_rGLMM <- function(
    n,
    design,
    prior,
    family,
    dispersion_ranef,
    fixef_start,
    m_convergence,
    gap_tol,
    tv_tol,
    mode_gap_max,
    collect_block1,
    verbose,
    progbar,
    cl
) {
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

  disp_info <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = design,
    fn_name          = "rglmerb"
  )

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
    gap_tol             = gap_tol,
    tv_tol              = tv_tol,
    mode_gap_max        = mode_gap_max,
    verbose             = verbose,
    progbar             = progbar,
    stage_verbose       = verbose,
    b_start             = NULL,
    collect_block1      = collect_block1
  )

  out$call <- cl
  out
}
