#' Bayesian linear mixed-effects model sampler (two-block Gibbs engine)
#'
#' Full sampling engine for Gaussian linear mixed models, parallel to
#' \code{\link[glmbayes]{rlmb}} in \pkg{glmbayes} and \code{\link{rglmerb}}
#' / \code{\link{glmerb}} in \pkg{lmebayes}.  Takes structured \code{design}
#' and \code{prior} objects, computes the ICM posterior mean internally, performs TV-calibrated
#' convergence calibration, and runs the two-block Gibbs sampler.
#'
#' \code{rlmerb} is called internally by \code{\link{lmerb}} after
#' \code{\link{model_setup}} and prior construction are complete.  It can also
#' be called directly in simulation or Gibbs-sampling workflows where formula
#' parsing and model-fit overhead are unnecessary.  It encapsulates:
#' \enumerate{
#'   \item ICM posterior mean (\code{\link[glmbayesCore]{lmerb_posterior_mean}})
#'   \item Block 1 prior construction (\code{.lmebayes_block1_prior_list})
#'   \item TV-calibrated convergence constants
#'         (\code{\link[glmbayesCore]{two_block_rate_v2}},
#'          \code{\link[glmbayesCore]{two_block_l_for_tv}})
#'   \item The two-block Gibbs sampler
#'         (\code{\link[glmbayesCore]{two_block_rNormal_reg_v2}})
#' }
#'
#' @param n Integer. Number of stored draws (each draw is one full pass through
#'   \code{m_convergence} inner Gibbs sweeps).
#' @param design A \code{\link{model_setup}} object as returned by
#'   \code{\link{model_setup}}, supplying \code{y}, \code{Z}, \code{groups},
#'   \code{X_hyper}, \code{group_name}, and \code{re_coef_names}.
#' @param prior A \code{lmebayes_prior_setup} object as returned by
#'   \code{\link{.lmebayes_priors_from_pfamily_list}}.
#' @param fixef_start Optional named list of starting hyper-parameter vectors
#'   (one per RE component).  When \code{NULL} (default), the ICM posterior
#'   mean is computed internally via
#'   \code{\link[glmbayesCore]{lmerb_posterior_mean}}.
#' @param m_convergence Optional integer. Number of inner Gibbs sweeps per
#'   stored draw.  When \code{NULL} (default), derived from \code{tv_tol} via
#'   Theorem 3 (Nygren 2020) and floored at the derived \code{m_min}.  A
#'   user-supplied value is floored at \code{m_min} with a warning if it had
#'   to be raised.
#' @param tv_tol Single numeric in \code{(0, 1)}. Total variation tolerance
#'   used for convergence calibration.  Default \code{0.01}.
#' @param seed Optional integer RNG seed.  Default \code{NULL}.
#' @param progbar Logical. Show a text progress bar during sampling.
#'   Default \code{TRUE}.
#' @param verbose Logical. Print the lmer-vs-ICM table and the convergence
#'   calibration line.  Default \code{TRUE}.
#' @return An object of class \code{c("rlmerb", "list")} with Block~2 fields in
#'   the \code{fixef.*} namespace (as \code{\link{rGLMM}}):
#'   \code{fixef}, \code{fixef.mode}, \code{fixef.init}, \code{fixef.means},
#'   \code{fixef.dispersion}, \code{fixef.dispersion.mean}, \code{fixef.iters},
#'   \code{fixef.iters.mean}, \code{fixef.mu}; Block~1 draws in
#'   \code{coefficients}; \code{ranef.mode}; \code{m_convergence};
#'   \code{convergence}; \code{Prior}; \code{design}.
#' @seealso \code{\link{lmerb}}, \code{\link{glmerb}}, \code{\link{rglmerb}},
#'   \code{\link[glmbayes]{rlmb}},
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}
#' @title The Bayesian Linear Mixed-Effects Model Distribution
#' @export
rlmerb <- function(
    n,
    design,
    prior,
    fixef_start   = NULL,
    m_convergence = NULL,
    tv_tol        = 0.01,
    seed          = NULL,
    progbar       = TRUE,
    verbose       = TRUE
) {
  cl <- match.call()

  # ---- argument validation --------------------------------------------------
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
      stop("'m_convergence' must be NULL or a single integer >= 1.",
           call. = FALSE)
    }
    m_convergence <- as.integer(m_convergence)
  }

  # ---- helpers --------------------------------------------------------------
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

  # ---- ICM posterior mean ---------------------------------------------------
  ranef_mode <- NULL
  if (is.null(fixef_start)) {
    pm          <- glmbayesCore::lmerb_posterior_mean(design, prior)
    fixef_start <- pm$fixef
    ranef_mode  <- pm$b_mean

    if (verbose) {
      fixef_lmer <- lapply(prior$prior_list, `[[`, "mu_fixef")
      names(fixef_lmer) <- re_names
      hdr <- sprintf("  %-18s  %-30s  %12s  %12s",
                     "RE component", "parameter", "lmer (start)", "post mean (ICM)")
      sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
      cat("--- lmerb: Block 2 fixed effects ---\n")
      cat(hdr, "\n")
      cat(sep, "\n")
      for (k in re_names) {
        nms_k  <- names(fixef_lmer[[k]])
        lmer_v <- fixef_lmer[[k]]
        pm_v   <- fixef_start[[k]]
        for (nm in nms_k) {
          cat(sprintf("  %-18s  %-30s  %12.4f  %12.4f\n",
                      k, nm, lmer_v[[nm]], pm_v[[nm]]))
        }
      }
      cat(sprintf("  (ICM converged: %s, %d iter, delta = %.2e)\n\n",
                  pm$converged, pm$iterations, pm$delta))
    }
  }

  # ---- block 1 prior --------------------------------------------------------
  block1_prior <- .lmebayes_block1_prior_list(prior)

  # ---- convergence calibration ----------------------------------------------
  rate <- glmbayesCore::two_block_rate_v2(
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = block1_prior,
    pfamily_list      = prior$pfamily_list,
    family            = gaussian(),
    group_levels      = group_levels
  )
  m_min <- glmbayesCore::two_block_l_for_tv(
    rate, tv_tol, method = "theorem3"
  ) + 1L

  if (is.null(m_convergence)) {
    m_convergence <- m_min
  } else if (m_convergence < m_min) {
    warning(
      "rlmerb: m_convergence = ", m_convergence, " is below the derived ",
      "minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
      "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }

  calib_label <- if (prior$any_ing) {
    "conservative: ING tau^2_k = disp_lower"
  } else {
    "exact"
  }
  if (verbose) {
    cat(sprintf(
      "--- lmerb: convergence calibration [%s]: lambda* = %.4f, tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n",
      calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
  }

  convergence_info <- list(
    method        = if (prior$any_ing) "disp_lower_bound" else "exact",
    tv_tol        = tv_tol,
    lambda_star   = rate$lambda_star,
    eigenvalues   = rate$eigenvalues,
    m_min         = m_min,
    m_convergence = m_convergence
  )

  # ---- sampling -------------------------------------------------------------
  out <- glmbayesCore::two_block_rNormal_reg_v2(
    n                 = n,
    y                 = design$y,
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = block1_prior,
    pfamily_list      = prior$pfamily_list,
    fixef_start       = fixef_start,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = design$group_name,
    family            = gaussian(),
    m_convergence     = m_convergence,
    seed              = seed,
    progbar           = progbar
  )

  staged <- .lmebayes_stage_v2_fixef(
    out          = out,
    fixef_mode   = fixef_start,
    fixef_init   = fixef_start,
    re_names     = re_names,
    group_levels = group_levels,
    n            = n
  )
  staged <- .lmebayes_add_fixef_summaries(staged)
  staged$call        <- cl
  staged$ranef.mode  <- ranef_mode
  staged$convergence <- convergence_info
  staged$Prior       <- list(
    block1_prior = block1_prior,
    pfamily_list = prior$pfamily_list
  )
  staged$design      <- design

  class(staged) <- c("rlmerb", "list")
  staged
}
