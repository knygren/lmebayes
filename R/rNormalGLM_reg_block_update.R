#' One Gibbs block update via \code{rNormalGLM_reg_block}
#'
#' Draw a single blockwise GLM posterior sample (\code{n = 1}) and return
#' latent-level draws (e.g. updated \code{theta}) for two-block Gibbs samplers
#' where each observation or group is its own block with scalar (or vector)
#' coefficients per block.
#'
#' This is a thin wrapper around \code{\link{rNormalGLM_reg_block}} that always
#' uses \code{n = 1}. When \code{prior_lists} and \code{prior_list} are both
#' omitted, per-block \code{dNormal} priors are built from \code{mu_all} and
#' \code{sigma_theta_sq} (scalar intercept per block).
#'
#' @param mu_all Length-\code{k} vector of prior means per block (linear predictors
#'   from the population block). Required unless \code{prior_lists} or
#'   \code{prior_list} is supplied.
#' @param sigma_theta_sq Shared prior variance for scalar blocks when building
#'   \code{prior_lists} from \code{mu_all}.
#' @param y,x,block,family Passed to \code{\link{rNormalGLM_reg_block}}.
#' @param prior_list,prior_lists Optional; passed through when supplied.
#' @param offset,weights Passed to \code{\link{rNormalGLM_reg_block}}.
#' @inheritParams rNormalGLM_reg_block
#' @param seed Optional; passed to \code{\link{set.seed}} before sampling.
#' @param theta_coef_col Column index of \code{coefficients} to return as
#'   \code{theta} (default \code{1} for scalar intercept blocks).
#' @return A list with:
#'   \describe{
#'     \item{theta}{Vector of length \code{nrow(x)} (one draw per row/block).}
#'     \item{coefficients,coef.mode}{Matrices from \code{rNormalGLM_reg_block}.}
#'     \item{rNormalGLM_reg_block}{Full block sampler output.}
#'   }
#' @seealso \code{\link{rNormalGLM_reg_block}}
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("glmbayes", quietly = TRUE)) {
#'   n <- 5L
#'   y <- rpois(n, 2)
#'   x <- matrix(1, n, 1)
#'   mu <- rep(0, n)
#'   rNormalGLM_reg_block_update(
#'     mu_all = mu, sigma_theta_sq = 1,
#'     y = y, x = x, block = seq_len(n),
#'     family = poisson(), use_parallel = FALSE
#'   )
#' }
#' }
rNormalGLM_reg_block_update <- function(mu_all,
                                        sigma_theta_sq = NULL,
                                        y,
                                        x,
                                        block,
                                        family = poisson(),
                                        prior_list = NULL,
                                        prior_lists = NULL,
                                        offset = NULL,
                                        weights = 1,
                                        Gridtype = 2L,
                                        n_envopt = 1L,
                                        use_parallel = TRUE,
                                        use_opencl = FALSE,
                                        verbose = FALSE,
                                        progbar = FALSE,
                                        seed = NULL,
                                        theta_coef_col = 1L) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(prior_lists) && is.null(prior_list)) {
    if (missing(mu_all)) {
      stop("Provide 'mu_all' (with 'sigma_theta_sq'), 'prior_lists', or 'prior_list'.",
           call. = FALSE)
    }
    if (is.null(sigma_theta_sq)) {
      stop("'sigma_theta_sq' is required when building per-block priors from 'mu_all'.",
           call. = FALSE)
    }
    prior_lists <- lapply(mu_all, function(m) {
      list(
        mu = m,
        Sigma = matrix(sigma_theta_sq, 1, 1),
        dispersion = 1,
        ddef = FALSE
      )
    })
  }

  out <- rNormalGLM_reg_block(
    n = 1L,
    y = y,
    x = x,
    block = block,
    prior_list = prior_list,
    prior_lists = prior_lists,
    offset = offset,
    weights = weights,
    family = family,
    Gridtype = as.integer(Gridtype),
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose,
    progbar = progbar
  )

  theta <- as.vector(out$coefficients[, theta_coef_col, drop = TRUE])

  list(
    theta = theta,
    coefficients = out$coefficients,
    coef.mode = out$coef.mode,
    rNormalGLM_reg_block = out
  )
}
