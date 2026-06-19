#' Run independent short Gibbs chains via \code{two_block_rNormal_reg_v5}
#'
#' Same contract as \code{\link{run_short_chains_v4}}, but calls
#' \code{\link[glmbayesCore]{two_block_rNormal_reg_v5}} (v5 C++ driver with
#' sweep-outer loop order and per-sweep chain progress bars).
#' Used by \code{\link{rglmerb_v5}} and \code{\link{glmerb}}.
#'
#' @inheritParams run_short_chains
#' @return As \code{\link{run_short_chains}}.
#' @keywords internal
run_short_chains_v5 <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    pfamily_list,
    family,
    re_names,
    group_levels,
    seed = NULL,
    seed_offset = 0L,
    collect_block1 = TRUE,
    progbar = FALSE,
    stage_label = "",
    diag_sweeps = FALSE,
    fixef_mode = NULL,
    b_mode = NULL
) {
  out <- glmbayesCore::two_block_rNormal_reg_v5(
    n                 = n_chains,
    y                 = design$y,
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = block1_prior,
    pfamily_list      = pfamily_list,
    fixef_start       = start_fixef,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = design$group_name,
    family            = family,
    m_convergence     = inner_sweeps,
    collect_block1    = collect_block1,
    seed              = seed,
    seed_offset       = seed_offset,
    progbar           = progbar,
    stage_label       = stage_label,
    diag_sweeps       = diag_sweeps,
    fixef_mode        = fixef_mode,
    b_mode            = b_mode
  )

  list(
    fixef_draws            = out$fixef_draws,
    dispersion_fixef_draws = out$dispersion_fixef_draws,
    iters_fixef_draws      = out$iters_fixef_draws,
    coefficients           = out$coefficients,
    mu_all_last            = out$mu_all_last
  )
}
