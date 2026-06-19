#' Run independent short Gibbs chains via sweep-outer R driver (v6)
#'
#' Same contract as \code{\link{run_short_chains_v4}}, but uses the pure-R
#' sweep-outer loop: for each inner sweep, all chains run Block~1, then all
#' chains run Block~2. Used by \code{\link{rglmerb_v6}} and \code{\link{glmerb}}.
#'
#' @inheritParams run_short_chains
#' @param stage_label Character label for diagnostics (e.g. \code{"pilot"}).
#' @param diag_sweeps If \code{TRUE}, print per-block means after each sweep.
#' @param fixef_mode ICM mode reference for fixef diagnostics.
#' @param b_mode ICM mode reference for random-effect diagnostics.
#' @param b_start Initial random-effect matrix for all chains (J x p_re).
#' @param ptypes Per-component pfamily names (optional; derived if \code{NULL}).
#' @return As \code{\link{run_short_chains}}.
#' @keywords internal
run_sweep_outer_chains_v6 <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    pfamily_list,
    family,
    re_names,
    group_levels,
    collect_block1 = TRUE,
    progbar        = FALSE,
    stage_label    = "",
    diag_sweeps    = FALSE,
    fixef_mode     = NULL,
    b_mode         = NULL,
    b_start        = NULL,
    ptypes         = NULL
) {
  if (is.null(ptypes)) {
    ptypes <- vapply(pfamily_list, function(pf) pf$pfamily, character(1))
    names(ptypes) <- re_names
  }

  tau2_start <- .lmebayes_tau2_start_from_pfamily(pfamily_list, re_names)
  if (is.null(b_start)) {
    if (is.null(b_mode)) {
      stop("'b_start' or 'b_mode' required for v6 batch init.", call. = FALSE)
    }
    b_start <- b_mode
  }

  batch <- .lmebayes_init_batch_state(
    n_chains     = n_chains,
    start_fixef  = start_fixef,
    b_start      = b_start,
    tau2_start   = tau2_start,
    re_names     = re_names,
    group_levels = group_levels
  )

  verbose_block_diag <- isTRUE(diag_sweeps)
  progbar_use <- isTRUE(progbar) && !verbose_block_diag

  for (m in seq_len(inner_sweeps)) {
    .lmebayes_print_sweep_boundary(
      stage_label  = stage_label,
      sweep        = m,
      inner_sweeps = inner_sweeps,
      phase        = "Block1",
      boundary     = "enter"
    )
    batch <- block1_all_chains(
      batch        = batch,
      design       = design,
      block1_prior = block1_prior,
      family       = family,
      ptypes       = ptypes,
      progbar      = progbar_use
    )
    .lmebayes_print_sweep_boundary(
      stage_label  = stage_label,
      sweep        = m,
      inner_sweeps = inner_sweeps,
      phase        = "Block1",
      boundary     = "exit"
    )
    if (verbose_block_diag) {
      .lmebayes_print_block_diag(
        stage_label  = stage_label,
        sweep        = m,
        inner_sweeps = inner_sweeps,
        phase        = "Block1",
        batch        = batch,
        fixef_mode   = fixef_mode,
        b_mode       = b_mode,
        re_names     = re_names,
        group_levels = group_levels
      )
    }

    .lmebayes_print_sweep_boundary(
      stage_label  = stage_label,
      sweep        = m,
      inner_sweeps = inner_sweeps,
      phase        = "Block2",
      boundary     = "enter"
    )
    batch <- block2_all_chains(
      batch        = batch,
      design       = design,
      pfamily_list = pfamily_list,
      ptypes       = ptypes,
      progbar      = progbar_use
    )
    .lmebayes_print_sweep_boundary(
      stage_label  = stage_label,
      sweep        = m,
      inner_sweeps = inner_sweeps,
      phase        = "Block2",
      boundary     = "exit"
    )
    if (verbose_block_diag) {
      .lmebayes_print_block_diag(
        stage_label  = stage_label,
        sweep        = m,
        inner_sweeps = inner_sweeps,
        phase        = "Block2",
        batch        = batch,
        fixef_mode   = fixef_mode,
        b_mode       = b_mode,
        re_names     = re_names,
        group_levels = group_levels
      )
    }
  }

  .lmebayes_pack_batch_draws(
    batch          = batch,
    design         = design,
    collect_block1 = collect_block1
  )
}
