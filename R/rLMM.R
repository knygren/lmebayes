#' Replicate-chain Gibbs sampling for Bayesian LMMs (v2 two-block driver)
#'
#' Low-level matrix-level LMM replicate-chain sampler using
#' \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}.  For formula-level
#' fitting with ICM posterior mean and \pkg{lmebayes} priors, see
#' \code{\link{rlmerb}} / \code{\link{lmerb}}.
#'
#' @param n Number of stored draws. If \code{length(n) > 1}, the length is
#'   taken to be the number required.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices \code{X_k}
#'   (\code{J x q_k}), one per column of \code{x}.
#' @param prior_list Prior for Block~1: \code{P} or \code{Sigma},
#'   \code{dispersion} (required for Gaussian), optional \code{ddef}.
#' @param pfamily_list Named list of \code{pfamily} objects for Block~2.
#' @param start Named list of Block~2 hyper-parameter vectors at which each
#'   replicate chain is initialised (typically the ICM posterior mean).
#' @param m_convergence Inner Gibbs steps per stored draw. When \code{NULL},
#'   derived from Theorem~3 at \code{start} using \code{tv_tol}.
#' @param tv_tol Total variation tolerance in \code{(0, 1)} for convergence
#'   calibration.
#' @param re_coef_names Character vector naming columns of \code{x}.
#' @param group_levels Character vector defining row order of Block~1 draws.
#' @param group_name Name for the grouping column in \code{coefficients}.
#' @param seed Optional RNG seed passed to the sampler.
#' @param progbar Logical; show a text progress bar during sampling.
#' @param verbose Print the convergence calibration line.
#' @param any_ing Logical. Label convergence calibration as ING-conservative
#'   when Block~2 uses \code{dIndependent_Normal_Gamma} components.
#' @return Object of class \code{c("rLMM", "list")} with Block~2 fields in
#'   the \code{fixef.*} namespace, Block~1 draws in \code{coefficients},
#'   plus \code{convergence_info}, \code{m_convergence}, and
#'   \code{draw_engine}.
#' @seealso \code{\link{rlmerb}}, \code{\link{rGLMM}},
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}}
#' @export
rLMM <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list,
    pfamily_list,
    start,
    m_convergence = NULL,
    tv_tol          = 0.01,
    re_coef_names   = colnames(x),
    group_levels    = levels(block),
    group_name      = NULL,
    seed            = NULL,
    progbar         = TRUE,
    verbose         = FALSE,
    any_ing         = FALSE
) {
  cl <- match.call()

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

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

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) {
        paste0("RE", seq_len(ncol(x)))
      } else {
        cn
      }
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop("length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
         call. = FALSE)
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  if (!is.list(start) || is.null(names(start))) {
    stop("'start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(start), re_names)) {
    stop("names(start) must match re_coef_names.", call. = FALSE)
  }
  start <- start[re_names]
  fixef_mode <- start

  if (!is.list(prior_list)) {
    stop("'prior_list' must be a list (Block~1 prior).", call. = FALSE)
  }
  if (!is.list(pfamily_list) || is.null(names(pfamily_list))) {
    stop("'pfamily_list' must be a named list.", call. = FALSE)
  }

  rate <- glmbayesCore::two_block_rate_v2(
    x                 = x,
    block             = block,
    x_hyper           = x_hyper,
    prior_list_block1 = prior_list,
    pfamily_list      = pfamily_list,
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
      "rLMM: m_convergence = ", m_convergence, " is below the derived ",
      "minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
      "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }

  calib_label <- if (isTRUE(any_ing)) {
    "conservative: ING tau^2_k = disp_lower"
  } else {
    "exact"
  }
  if (isTRUE(verbose)) {
    cat(sprintf(
      "--- rLMM: convergence calibration [%s]: lambda* = %.4f, tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n",
      calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
  }

  convergence_info <- list(
    method        = if (isTRUE(any_ing)) "disp_lower_bound" else "exact",
    tv_tol        = tv_tol,
    lambda_star   = rate$lambda_star,
    eigenvalues   = rate$eigenvalues,
    m_min         = m_min,
    m_convergence = m_convergence,
    draw_engine   = "two_block_rNormal_reg_v2"
  )

  out <- glmbayesCore::two_block_rNormal_reg_v2(
    n                 = n,
    y                 = y,
    x                 = x,
    block             = block,
    x_hyper           = x_hyper,
    prior_list_block1 = prior_list,
    pfamily_list      = pfamily_list,
    fixef_start       = start,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = group_name,
    family            = gaussian(),
    m_convergence     = m_convergence,
    seed              = seed,
    progbar           = progbar
  )

  staged <- .lmebayes_stage_v2_fixef(
    out          = out,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_mode,
    re_names     = re_names,
    group_levels = group_levels,
    n            = n
  )

  staged$call             <- cl
  staged$m_convergence    <- m_convergence
  staged$convergence_info <- convergence_info
  staged$draw_engine      <- "two_block_rNormal_reg_v2"
  staged$pfamily_list     <- pfamily_list
  staged$prior_list       <- prior_list
  staged$family           <- gaussian()

  class(staged) <- c("rLMM", "list")
  staged
}
