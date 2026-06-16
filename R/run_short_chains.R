#' Run independent short Gibbs chains for the two-block mixed model
#'
#' Executes \code{n_chains} independent two-block Gibbs chains, each
#' initialised at \code{start_fixef} and run for \code{inner_sweeps} sweeps.
#' Each chain stores exactly one draw (\code{n = 1L}), so the result is a
#' collection of \code{n_chains} approximately independent draws from the
#' target posterior (for large enough \code{inner_sweeps}).  Used by
#' \code{\link{glmerb}} for both the pilot and main sampling stages.
#'
#' @param n_chains Integer. Number of independent chains to run.
#' @param start_fixef Named list of starting hyper-parameter vectors, one
#'   named numeric vector per RE component (matches \code{re_names}).
#' @param inner_sweeps Integer. Number of inner Gibbs sweeps per chain
#'   (\code{m_convergence} argument of
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}).
#' @param design A \code{\link{model_setup}} object supplying \code{y},
#'   \code{Z}, \code{groups}, \code{X_hyper}, and \code{group_name}.
#' @param block1_prior Block 1 prior list as returned by
#'   \code{.lmebayes_block1_prior_list()}.
#' @param pfamily_list Named list of \code{\link[glmbayesCore]{pfamily}}
#'   objects (one per RE component), as supplied to \code{\link{glmerb}}.
#' @param family A \code{\link[stats]{family}} object for the response model.
#' @param re_names Character vector of RE coefficient names
#'   (\code{design$re_coef_names}).
#' @param group_levels Character vector of group levels
#'   (\code{levels(design$groups)}).
#' @param seed_offset Integer added to \code{seed} for each chain
#'   (\code{seed + seed_offset + i}).  Default \code{0L}.
#' @param seed Base RNG seed, or \code{NULL} for no seeding.  Default
#'   \code{NULL}.
#' @param collect_block1 Logical. If \code{TRUE}, collect and rbind the
#'   Block 1 (\code{coefficients}) matrix from every chain.  Default
#'   \code{TRUE}.
#' @param progbar Logical. When \code{TRUE} and \code{n_chains > 1}, show a
#'   text progress bar over independent chains. When \code{n_chains == 1},
#'   passed through to \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}
#'   (inner Gibbs sweeps). Default \code{FALSE}.
#' @return A list with components:
#'   \describe{
#'     \item{\code{fixef_draws}}{Named list of \code{n_chains x q_k} matrices,
#'       one per RE component.}
#'     \item{\code{dispersion_fixef_draws}}{\code{n_chains x p_re} matrix of
#'       Block 2 dispersion (\eqn{\tau^2}) draws.}
#'     \item{\code{iters_fixef_draws}}{\code{n_chains x p_re} matrix of
#'       envelope iteration counts.}
#'     \item{\code{coefficients}}{If \code{collect_block1 = TRUE}: an
#'       \code{(n_chains * J) x p_re} matrix of Block 1 endpoint draws with
#'       row names stripped; \code{NULL} otherwise.}
#'     \item{\code{mu_all_last}}{Per-observation fitted means vector from the
#'       final chain.}
#'   }
#' @keywords internal
run_short_chains <- function(
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

#' Text progress bar matching \pkg{glmbayesCore} C++ style
#' @param current Completed step (1-based, up to \code{total}).
#' @param total Total number of steps.
#' @noRd
.lmebayes_progress_bar <- function(current, total) {
  if (total <= 0L) {
    return(invisible())
  }
  totaldotz <- 40L
  fraction  <- current / total
  dotz      <- round(fraction * totaldotz)
  cat("\r", strrep(" ", 80L), "\r", sep = "")
  cat(sprintf("%3.0f%% [", fraction * 100), sep = "")
  cat(paste0(rep("=", dotz), collapse = ""))
  cat(paste0(rep(" ", totaldotz - dotz), collapse = ""))
  cat("]", sep = "")
  utils::flush.console()
}

#' Finish a progress bar started by \code{.lmebayes_progress_bar}
#' @noRd
.lmebayes_progress_bar_finish <- function() {
  cat("\n")
}
