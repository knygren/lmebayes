#' Bayesian linear mixed-effects model sampler (two-block Gibbs engine)

#'

#' Full sampling engine for Gaussian linear mixed models, parallel to

#' \code{\link[glmbayes]{rlmb}} in \pkg{glmbayes} and \code{\link{rglmerb}}

#' / \code{\link{glmerb}} in \pkg{lmebayes}.  Takes structured \code{design}

#' and \code{prior} objects, computes the ICM posterior mean internally, and

#' delegates replicate-chain sampling to
#' \code{\link[glmbayesCore]{rLMMNormal_reg}} or
#' \code{\link[glmbayesCore]{rLMMindepNormalGamma_reg}} when \code{dispersion_ranef} is a
#' \code{dGamma()} pfamily.

#'

#' \code{rlmerb} is called internally by \code{\link{lmerb}} after

#' \code{\link{model_setup}} and prior construction are complete.  It can also

#' be called directly in simulation or Gibbs-sampling workflows where formula

#' parsing and model-fit overhead are unnecessary.

#'

#' @param n Integer. Number of stored draws (each draw is one full pass through

#'   \code{m_convergence} inner Gibbs sweeps).

#' @param design A \code{\link{model_setup}} object as returned by

#'   \code{\link{model_setup}}, supplying \code{y}, \code{Z}, \code{groups},

#'   \code{X_hyper}, \code{group_name}, and \code{re_coef_names}.

#' @param prior Normalized prior container with \code{Sigma_ranef}, \code{prior_list},
#'   and related Block~2 fields (as built internally by \code{\link{lmerb}} from
#'   \code{pfamily_list} and \code{dispersion_ranef}).

#' @param dispersion_ranef Required observation-level dispersion: a positive

#'   scalar \eqn{\sigma^2} (fixed) or a \code{\link{dGamma}()} pfamily with

#'   \code{Inv_Dispersion = TRUE} for a Gamma prior on \eqn{\sigma^2}.

#'   Typically \code{Prior_Setup_lmebayes(...)$dispersion_ranef} for a fixed plug-in.

#' @param fixef_start Optional named list of starting hyper-parameter vectors

#'   (one per RE component).  When \code{NULL} (default), the ICM posterior

#'   mean is computed inside the Core engine
#'   (\code{\link[glmbayesCore]{rLMMNormal_reg}} or
#'   \code{\link[glmbayesCore]{rLMMindepNormalGamma_reg}}).

#' @param m_convergence Optional integer. Number of inner Gibbs sweeps per

#'   stored draw.  When \code{NULL} (default), derived from \code{tv_tol} via

#'   Theorem 3 (Nygren 2020) and floored at the derived \code{m_min}.  A

#'   user-supplied value is floored at \code{m_min} with a warning if it had

#'   to be raised.

#' @param tv_tol Single numeric in \code{(0, 1)}. Total variation tolerance

#'   used for convergence calibration.  Default \code{0.01}.

#' @param gap_tol Legacy mode--mean gap tolerance for the pilot stage when
#'   any Block~2 component uses \code{dIndependent_Normal_Gamma} and
#'   \code{tv_tol} is \code{NULL}.  Ignored for all-\code{dNormal} models.

#' @param mode_gap_max Pilot inner-sweep calibration for ING Block~2 models
#'   (default \code{1.0}).  Ignored for all-\code{dNormal} models.

#' @param progbar Logical. Show a text progress bar during sampling.

#'   Default \code{TRUE}.

#' @param verbose Logical. Print the lmer-vs-ICM table and the convergence

#'   calibration line.  Default \code{TRUE}.

#' @param print_icm_table Logical. When \code{FALSE}, skip the reference-vs-ICM

#'   table (e.g. when \code{\link{rglmerb}} prints glmer-labelled output).

#'   The convergence calibration line from the Core engine still follows

#'   \code{verbose}.  Default \code{TRUE}.

#' @param diag_sweeps Temporary diagnostic flag for ING Block~2 models with a
#'   pilot stage.  When \code{TRUE}, print one combined Block~2 chain-mean table
#'   per stage when each stage finishes (as \code{print()} on
#'   \code{$sweep_history}) and attach sweep history on the fit.
#'   the result.  Default \code{FALSE}.

#' @return An object of class \code{c("rlmerb", "list")} with Block~2 fields in

#'   the \code{fixef.*} namespace (Core LMM engines):

#'   \code{fixef}, \code{fixef.mode}, \code{fixef.init}, \code{fixef.means},

#'   \code{fixef.dispersion}, \code{fixef.dispersion.mean}, \code{fixef.iters},

#'   \code{fixef.iters.mean}, \code{fixef.mu}; Block~1 draws in

#'   \code{coefficients}; \code{ranef.mode}; \code{m_convergence};

#'   \code{convergence}; \code{Prior}; \code{design}.

#' @seealso \code{\link{lmerb}},
#'   \code{\link[glmbayesCore]{rLMMNormal_reg}},
#'   \code{\link[glmbayesCore]{rLMMindepNormalGamma_reg}}, \code{\link{glmerb}},

#'   \code{\link{rglmerb}}, \code{\link[glmbayes]{rlmb}}

#' @title The Bayesian Linear Mixed-Effects Model Distribution

#' @export

rlmerb <- function(

    n,

    design,

    prior,

    dispersion_ranef,

    fixef_start   = NULL,

    m_convergence = NULL,

    tv_tol        = 0.01,

    progbar         = TRUE,

    verbose         = TRUE,

    print_icm_table = TRUE,

    gap_tol             = 0.0196,

    mode_gap_max        = 1.0,

    diag_sweeps         = FALSE

) {

  cl <- match.call()



  if (length(n) > 1L) n <- length(n)

  n <- as.integer(n[1L])

  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)



  if (!inherits(design, "model_setup")) {

    stop("'design' must be a model_setup object.", call. = FALSE)

  }

  if (missing(dispersion_ranef)) {
    stop(
      "'dispersion_ranef' is required for rlmerb(). Typically ",
      "Prior_Setup_lmebayes(...)$dispersion_ranef.",
      call. = FALSE
    )
  }

  disp_info <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = gaussian(),
    design           = design,
    fn_name          = "rlmerb"
  )

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

  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  block1_prior <- .lmebayes_block1_prior_list(
    prior,
    dispersion_ranef = disp_info$dispersion_fix
  )

  out <- .lmebayes_run_lmm_engine(
    n               = n,
    design          = design,
    prior           = prior,
    disp_info       = disp_info,
    fixef_start     = fixef_start,
    m_convergence   = m_convergence,
    tv_tol          = tv_tol,
    progbar         = progbar,
    verbose         = verbose,
    gap_tol             = gap_tol,
    mode_gap_max        = mode_gap_max,
    diag_sweeps         = diag_sweeps
  )



  if (is.null(fixef_start) && isTRUE(print_icm_table)) {
    icm_lbl <- .lmebayes_block2_icm_labels(prior, gaussian())
    .lmebayes_print_icm_fixef_table(
      prior_list = prior$prior_list,
      re_names   = re_names,
      fixef_icm  = out$fixef.mode,
      icm_info   = out$icm_info,
      ref_label  = icm_lbl$ref_label,
      icm_label  = icm_lbl$icm_label,
      conv_label = icm_lbl$conv_label,
      header     = "--- lmerb: Block 2 fixed effects ---",
      verbose    = verbose
    )
  }



  out <- .lmebayes_add_fixef_summaries(out)

  out$call       <- cl

  out$convergence <- out$convergence_info

  out$Prior      <- list(

    block1_prior         = block1_prior,

    pfamily_list         = prior$pfamily_list,

    dispersion_ranef     = disp_info$dispersion_fix,

    dispersion_mode      = disp_info$mode,

    dispersion_pfamily   = disp_info$dispersion_pfamily,

    dispersion_prior_list = disp_info$dispersion_prior_list

  )

  out$design     <- design

  if (!is.null(out$n_pilot) && out$n_pilot > 0L) {
    .lmebayes_print_fixef_init(
      out$fixef.init,
      re_names,
      verbose,
      header = "--- lmerb: main-stage fixef.init (pilot colMeans) ---"
    )
  }

  class(out) <- c("rlmerb", "list")

  out

}

