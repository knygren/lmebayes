#' Run independent short Gibbs chains via \code{two_block_rNormal_reg_v2}
#'
#' Same contract as \code{\link{run_short_chains}}: an R loop over
#' \code{n_chains}, each iteration calling
#' \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}} with \code{n = 1L}
#' and \code{set.seed(seed + seed_offset + i)} per chain.  Entry point for
#' \code{\link{rglmerb_v2}} / \code{\link{glmerb}}.  For the C++ chain loop
#' (development v3), see \code{\link{run_short_chains_v3}}.
#'
#' @inheritParams run_short_chains
#' @return As \code{\link{run_short_chains}}.
#' @keywords internal
run_short_chains_v2 <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    pfamily_list,
    family,
    re_names,
    group_levels,
    seed_offset    = 0L,
    seed           = NULL,
    collect_block1 = TRUE,
    progbar        = FALSE
) {
  fe_draws <- lapply(start_fixef, function(beta0) {
    mat <- matrix(NA_real_, nrow = n_chains, ncol = length(beta0))
    colnames(mat) <- names(beta0)
    mat
  })
  names(fe_draws) <- names(start_fixef)

  tau2_draws_local  <- matrix(NA_real_, nrow = n_chains, ncol = length(re_names))
  colnames(tau2_draws_local) <- re_names
  iters_draws_local <- matrix(NA_real_, nrow = n_chains, ncol = length(re_names))
  colnames(iters_draws_local) <- re_names

  coef_rows <- if (collect_block1) vector("list", n_chains) else NULL
  mu_last   <- NULL

  show_chain_bar <- isTRUE(progbar) && n_chains > 1L
  inner_progbar  <- isTRUE(progbar) && n_chains == 1L

  for (i in seq_len(n_chains)) {
    if (show_chain_bar) {
      .lmebayes_progress_bar(i, n_chains)
    }
    seed_i <- if (!is.null(seed)) as.integer(seed + seed_offset + i) else NULL
    out_i  <- glmbayesCore::two_block_rNormal_reg_v2(
      n                 = 1L,
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
      seed              = seed_i,
      progbar           = inner_progbar
    )
    for (k in re_names) {
      fe_draws[[k]][i, ] <- out_i$fixef_draws[[k]][1L, ]
    }
    tau2_draws_local[i, ]  <- out_i$dispersion_fixef_draws[1L, re_names]
    iters_draws_local[i, ] <- out_i$iters_fixef_draws[1L, re_names]
    if (collect_block1) {
      coef_rows[[i]] <- out_i$coefficients
    }
    mu_last <- out_i$mu_all_last
  }

  if (show_chain_bar) {
    .lmebayes_progress_bar_finish()
  }

  coefficients_local <- if (collect_block1) {
    out <- do.call(rbind, coef_rows)
    rownames(out) <- NULL
    out
  } else {
    NULL
  }

  list(
    fixef_draws            = fe_draws,
    dispersion_fixef_draws = tau2_draws_local,
    iters_fixef_draws      = iters_draws_local,
    coefficients           = coefficients_local,
    mu_all_last            = mu_last
  )
}
