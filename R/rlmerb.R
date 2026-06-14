#' Raw two-block Gibbs sampler for Bayesian linear mixed models
#'
#' Matrix-level sampling engine for the two-block Gibbs sampler, parallel to
#' \code{\link[glmbayes]{rlmb}} in \pkg{glmbayes}.  Takes raw numeric inputs
#' (no formula, no data frame, no \code{\link[lme4]{lmer}} call) and returns a
#' lean \code{"rlmerb"} object containing the Block 1 and Block 2 posterior
#' draws.
#'
#' \code{rlmerb} is a thin wrapper over
#' \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}.  It is called
#' internally by \code{\link{lmerb}} after \code{\link{model_setup}} and prior
#' construction are complete.  It can also be called directly in simulation or
#' Gibbs-sampling workflows where formula parsing and model-fit overhead are
#' unnecessary.
#'
#' @param n Integer. Number of stored draws (each draw is one full pass through
#'   \code{m_convergence} inner Gibbs sweeps).
#' @param y Numeric response vector of length \code{nrow(Z)}.
#' @param Z Level-1 random-effects design matrix (\code{l2 x p_re}).  Passed
#'   as \code{x} to \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}.
#' @param groups Factor of length \code{nrow(Z)} giving the group membership
#'   for each observation.
#' @param X_hyper Named list of group-level design matrices (\code{J x q_k}),
#'   one per RE coefficient.  Names must match \code{re_names}.
#' @param block1_prior Block 1 prior list as returned by
#'   \code{.lmebayes_block1_prior_list()}: contains \code{P} or \code{Sigma}
#'   and (for Gaussian) \code{dispersion}.
#' @param pfamily_list Named list of \code{\link[glmbayesCore]{pfamily}}
#'   objects (one per RE component), specifying Block 2 hyperpriors.
#' @param fixef_start Named list of starting hyper-parameter vectors, one per
#'   RE component (Block 2 initialisation).
#' @param family A \code{\link[stats]{family}} object for the response model.
#'   Default \code{gaussian()}.
#' @param re_names Character vector of RE coefficient names.  Default
#'   \code{names(pfamily_list)}.
#' @param group_levels Character vector defining the row order of Block 1
#'   draws (group levels).  Default \code{levels(groups)}.
#' @param group_name Optional character(1): name of the grouping variable,
#'   stored as a column label in \code{coefficients}.  Default \code{NULL}.
#' @param m_convergence Integer. Number of inner Gibbs sweeps per stored draw.
#'   Default \code{10L}.
#' @param seed Optional integer RNG seed.  Default \code{NULL}.
#' @param progbar Logical. Show a text progress bar.  Default \code{TRUE}.
#' @return An object of class \code{"rlmerb"}: a list with components:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{fixef_draws}}{Named list of \code{n x q_k} matrices of
#'       Block 2 (hyper-parameter) draws, one per RE component.}
#'     \item{\code{coefficients}}{Matrix or data frame of Block 1 (random-effect)
#'       endpoint draws as returned by \code{two_block_rNormal_reg_v2}.}
#'     \item{\code{dispersion_fixef_draws}}{\code{n x p_re} matrix of
#'       \eqn{\tau^2_k} draws per stored draw.}
#'     \item{\code{iters_fixef_draws}}{\code{n x p_re} matrix of total Block 2
#'       envelope candidates per stored draw.}
#'     \item{\code{mu_all_last}}{Per-observation fitted means at the final
#'       Gibbs state.}
#'     \item{\code{coef.mode}}{\code{fixef_start} echoed: the Block 2 starting
#'       point (typically the ICM posterior mean from \code{lmerb}).}
#'     \item{\code{Prior}}{List with \code{block1_prior} and \code{pfamily_list}
#'       echoed.}
#'     \item{\code{y}}{Response vector echoed.}
#'     \item{\code{x}}{\code{Z} echoed (mirrors \code{rlmb$x}).}
#'     \item{\code{groups}}{Grouping factor echoed.}
#'   }
#' @seealso \code{\link{lmerb}}, \code{\link[glmbayes]{rlmb}},
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}
#' @export
rlmerb <- function(
    n,
    y,
    Z,
    groups,
    X_hyper,
    block1_prior,
    pfamily_list,
    fixef_start,
    family        = gaussian(),
    re_names      = names(pfamily_list),
    group_levels  = levels(groups),
    group_name    = NULL,
    m_convergence = 10L,
    seed          = NULL,
    progbar       = TRUE
) {
  cl <- match.call()

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  m_convergence <- as.integer(m_convergence[1L])
  if (m_convergence < 1L) {
    stop("'m_convergence' must be at least 1.", call. = FALSE)
  }

  out <- glmbayesCore::two_block_rNormal_reg_v2(
    n                 = n,
    y                 = y,
    x                 = Z,
    block             = groups,
    x_hyper           = X_hyper,
    prior_list_block1 = block1_prior,
    pfamily_list      = pfamily_list,
    fixef_start       = fixef_start,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = group_name,
    family            = family,
    m_convergence     = m_convergence,
    seed              = seed,
    progbar           = progbar
  )

  structure(
    list(
      call                   = cl,
      fixef_draws            = out$fixef_draws,
      coefficients           = out$coefficients,
      dispersion_fixef_draws = out$dispersion_fixef_draws,
      iters_fixef_draws      = out$iters_fixef_draws,
      mu_all_last            = out$mu_all_last,
      coef.mode              = fixef_start,
      Prior                  = list(block1_prior = block1_prior,
                                    pfamily_list = pfamily_list),
      y                      = y,
      x                      = Z,
      groups                 = groups
    ),
    class = c("rlmerb", "list")
  )
}
