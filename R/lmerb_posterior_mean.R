#' Joint posterior mean of the two-block Gaussian model
#'
#' Finds the joint posterior mean (= joint mode, since the posterior is exactly
#' multivariate normal when variance components are fixed) of the two-block
#' model fit by \code{\link{lmerb}}, using an \emph{iterated conditional means}
#' (ICM) algorithm.
#'
#' @details
#' \strong{Algorithm.}
#' For any jointly Gaussian distribution the conditional mean of each block is
#' an affine function of the other block's value.  ICM alternates between the
#' two closed-form conditional mean updates:
#'
#' \describe{
#'   \item{Block 1 mean}{For each group \eqn{j}:
#'     \deqn{
#'       E[b_j \mid \gamma] =
#'       \bigl(Z_j^\top Z_j / \sigma^2 + P_b\bigr)^{-1}
#'       \bigl(Z_j^\top y_j / \sigma^2 + P_b \,\mu_j(\gamma)\bigr)
#'     }
#'     where \eqn{P_b = \Sigma_b^{-1}} and
#'     \eqn{\mu_j(\gamma)} is the Block 2 prior mean from
#'     \code{\link{build_mu_all}}.
#'     The quantities \eqn{Z_j^\top Z_j / \sigma^2} and
#'     \eqn{Z_j^\top y_j / \sigma^2} are constant across iterations and are
#'     pre-computed once.
#'   }
#'   \item{Block 2 mean}{For each RE component \eqn{k}:
#'     \deqn{
#'       E[\gamma_k \mid b_k] =
#'       \bigl(X_k^\top X_k / \tau^2_k + P_{\gamma_k}\bigr)^{-1}
#'       \bigl(X_k^\top b_k / \tau^2_k + P_{\gamma_k} \mu_{\gamma_k}\bigr)
#'     }
#'     where \eqn{b_k} is the \eqn{k}-th column of the current Block 1 mean
#'     matrix, \eqn{X_k = } \code{design$X_hyper[[k]]},
#'     \eqn{\tau^2_k = } \code{dispersion_fixef},
#'     and \eqn{P_{\gamma_k} = \Sigma_{\gamma_k}^{-1}}.
#'   }
#' }
#'
#' \strong{Convergence.}
#' Because the iteration is a linear contraction, it converges to the unique
#' fixed point (the joint posterior mean) with a geometric rate equal to
#' \eqn{\lambda^*}, the same eigenvalue that bounds Gibbs mixing in Nygren
#' (2020).  Convergence is declared when
#' \eqn{\max_k \|\gamma_k^{(\text{new})} - \gamma_k\|_\infty < \delta}.
#'
#' \strong{Use in \code{lmerb}.}
#' The joint posterior mean is the optimal starting point for the inner Gibbs
#' loop: initialising at the mean minimises the TV distance achievable in
#' \code{m_convergence} steps.  Replace the placeholder \code{fixef_start}
#' in \code{\link{lmerb}} with \code{lmerb_posterior_mean(...)$fixef} once
#' this function has been validated on the target data.
#'
#' @param design Object of class \code{"model_setup"} from
#'   \code{\link{model_setup}}.
#' @param measurement_prior_list Object of class \code{"lmebayes_prior_setup"}
#'   from \code{\link{Prior_Setup_lmebayes}}.
#' @param tol Convergence tolerance on the \eqn{\ell_\infty} change in
#'   \code{fixef} between successive iterations.  Default \code{1e-10}.
#' @param maxit Maximum number of ICM iterations.  Default \code{200L}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{fixef}}{Named list of posterior mean vectors for the
#'       level-2 fixed effects \eqn{\gamma_k}, one entry per
#'       \code{design$re_coef_names}.  Same structure and names as the
#'       \code{fixef} argument of \code{\link{lmerb}}.  Use directly as
#'       \code{fixef_start} in \code{\link{lmerb}}.}
#'     \item{\code{b_mean}}{\eqn{J \times p_{\mathrm{re}}} numeric matrix of
#'       posterior mean random effects.  Rows are group levels
#'       (\code{levels(design$groups)}); columns are
#'       \code{design$re_coef_names}.}
#'     \item{\code{converged}}{Logical; \code{TRUE} if \code{tol} was reached
#'       before \code{maxit} iterations.}
#'     \item{\code{iterations}}{Integer; number of ICM iterations performed.}
#'     \item{\code{delta}}{Numeric; final \eqn{\ell_\infty} change in
#'       \code{fixef} at the last iteration.}
#'   }
#'
#' @references
#' Nygren, K. (2020). \emph{On the total variation distance between multivariate
#' normal densities with applications to two-block Gibbs samplers.}
#' Unpublished manuscript.
#'
#' @seealso \code{\link{lmerb}}, \code{\link{build_mu_all}},
#'   \code{\link{Prior_Setup_lmebayes}}
#' @export
lmerb_posterior_mean <- function(design,
                                 measurement_prior_list,
                                 tol   = 1e-10,
                                 maxit = 200L) {

  if (!inherits(design, "model_setup")) {
    stop("'design' must be a model_setup object.", call. = FALSE)
  }
  if (!inherits(measurement_prior_list, "lmebayes_prior_setup")) {
    stop(
      "'measurement_prior_list' must be an lmebayes_prior_setup object.",
      call. = FALSE
    )
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  J            <- length(group_levels)
  p_re         <- length(re_names)
  g_chr        <- as.character(design$groups)

  sigma2 <- measurement_prior_list$dispersion_ranef
  P_b    <- solve(measurement_prior_list$Sigma_ranef)   # p_re x p_re

  # Per-RE Block 2 quantities (constant across ICM iterations)
  P_gamma  <- stats::setNames(
    lapply(re_names, function(k) {
      solve(measurement_prior_list$prior_list[[k]]$Sigma_fixef)
    }),
    re_names
  )
  mu_gamma <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$mu_fixef),
    re_names
  )
  tau2 <- stats::setNames(
    lapply(re_names, function(k) measurement_prior_list$prior_list[[k]]$dispersion_fixef),
    re_names
  )

  # Pre-compute per-group Z_j'Z_j/sigma2 and Z_j'y_j/sigma2.
  # These depend only on the data and are constant across all ICM iterations.
  ZtZ_scaled <- vector("list", J)
  Zty_scaled <- vector("list", J)
  names(ZtZ_scaled) <- names(Zty_scaled) <- group_levels

  for (lev in group_levels) {
    rows <- which(g_chr == lev)
    Z_j  <- design$Z[rows, , drop = FALSE]
    y_j  <- design$y[rows]
    ZtZ_scaled[[lev]] <- crossprod(Z_j) / sigma2         # p_re x p_re
    Zty_scaled[[lev]] <- crossprod(Z_j, y_j) / sigma2    # p_re vector
  }

  # Initialise fixef to lmer mu_fixef (best available proxy for posterior mean)
  fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
  names(fixef) <- re_names

  # b_mean: J x p_re, rows = group levels, cols = RE names
  b_mean <- matrix(
    0.0, nrow = J, ncol = p_re,
    dimnames = list(group_levels, re_names)
  )

  converged <- FALSE
  delta     <- NA_real_

  for (iter in seq_len(maxit)) {

    # -- Block 1 mean: E[b_j | fixef] ------------------------------------
    mu_all <- as.matrix(build_mu_all(design, fixef)$mu_all)   # p_re x J

    for (jj in seq_len(J)) {
      lev      <- group_levels[jj]
      mu_j     <- mu_all[, jj]                               # p_re vector
      post_P_j <- ZtZ_scaled[[lev]] + P_b                   # p_re x p_re
      post_v_j <- Zty_scaled[[lev]] + P_b %*% mu_j          # p_re vector
      b_mean[jj, ] <- solve(post_P_j, post_v_j)
    }

    # -- Block 2 mean: E[gamma_k | b_mean_k] ----------------------------
    fixef_new <- vector("list", p_re)
    names(fixef_new) <- re_names

    for (k in re_names) {
      X_k      <- design$X_hyper[[k]]                       # J x q_k
      b_k      <- b_mean[, k]                               # J vector
      tau2_k   <- tau2[[k]]
      P_gam_k  <- P_gamma[[k]]                              # q_k x q_k
      mu_gam_k <- mu_gamma[[k]]                             # q_k vector

      post_P_k  <- crossprod(X_k) / tau2_k + P_gam_k       # q_k x q_k
      post_v_k  <- crossprod(X_k, b_k) / tau2_k +
                   P_gam_k %*% mu_gam_k                     # q_k vector
      gam_k <- as.vector(solve(post_P_k, post_v_k))
      names(gam_k) <- colnames(X_k)
      fixef_new[[k]] <- gam_k
    }

    # -- Convergence check -----------------------------------------------
    delta <- max(vapply(re_names, function(k) {
      max(abs(fixef_new[[k]] - fixef[[k]]))
    }, numeric(1L)))

    fixef <- fixef_new

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      "lmerb_posterior_mean() did not converge in ", maxit, " iterations ",
      "(final delta = ", signif(delta, 3L), "). ",
      "Consider increasing 'maxit' or checking model identifiability.",
      call. = FALSE
    )
  }

  list(
    fixef      = fixef,
    b_mean     = b_mean,
    converged  = converged,
    iterations = iter,
    delta      = delta
  )
}
