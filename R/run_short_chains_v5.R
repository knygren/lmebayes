#' Run independent short Gibbs chains via \code{two_block_rNormal_reg_v5}
#'
#' Calls \code{\link[glmbayesCore]{two_block_rNormal_reg_v5}} (v5 C++ driver with
#' sweep-outer loop order and per-sweep chain progress bars).
#' Used by \code{\link{rglmerb_v5}} and \code{\link{glmerb}}.
#'
#' @param n_chains Number of independent short chains.
#' @param start_fixef Named list of Block~2 starting vectors (one per chain).
#' @param inner_sweeps Inner Gibbs sweeps per chain (\code{m_convergence}).
#' @param design A \code{\link{model_setup}} object.
#' @param block1_prior Block~1 prior list for \code{block_rNormalGLM}.
#' @param pfamily_list Named list of Block~2 \code{pfamily} objects.
#' @param family Response \code{\link[stats]{family}}.
#' @param re_names Random-effect coefficient names.
#' @param group_levels Factor levels for the grouping variable.
#' @param seed Base RNG seed for the v5 C++ driver (see
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg_v5}}).
#' @param seed_offset Passed to the v5 C++ driver (offsets per-chain reseeds).
#' @param collect_block1 If \code{TRUE}, row-bind Block~1 coefficient draws.
#' @param progbar Show per-sweep progress bars.
#' @param stage_label Stage label stored on \code{$sweep_history} (e.g. \code{"pilot"}).
#' @param diag_sweeps Unused; live sweep diagnostics are disabled on the v5 path.
#' @param fixef_mode ICM mode reference for \code{$sweep_history}.
#' @param b_mode Random-effects mode matrix for optional diagnostics.
#' @return A list with \code{fixef_draws}, \code{coefficients},
#'   \code{dispersion_fixef_draws}, \code{iters_fixef_draws},
#'   \code{iters_ranef_draws}, \code{mu_all_last}, and \code{sweep_history}.
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
    use_parallel      = FALSE,
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
    iters_ranef_draws      = out$iters_ranef_draws,
    coefficients           = out$coefficients,
    mu_all_last            = out$mu_all_last,
    sweep_history          = out$sweep_history
  )
}
