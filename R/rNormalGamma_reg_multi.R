#' Normal--Gamma Regression with Multiple Responses
#' @family simfuncs
#' @description
#' Normal--Gamma regression with a **matrix** response: run
#' \code{\link{rNormalGamma_reg}} once per column of \code{y}, sharing the
#' same design matrix \code{x} (as \code{lm} with \code{cbind} responses).
#'
#' Same argument list as \code{\link{rNormalGamma_reg}} except \code{prior_list}:
#' a **list** of length \code{ncol(y)}; element \code{j} is the
#' \code{prior_list} passed to \code{\link{rNormalGamma_reg}} for column
#' \code{j} of \code{y} (each with \code{mu}, \code{Sigma} or \code{P},
#' \code{shape}, and \code{rate}).
#'
#' @inheritParams rNormalGamma_reg
#' @param prior_list List of length \code{ncol(y)} of per-column prior lists (see Details).
#' @return A named list of class \code{"mrglmb"}.  Each element is the
#'   \code{\link{rNormalGamma_reg}} output (class \code{"rglmb"}) for the
#'   corresponding column of \code{y}, named by \code{colnames(y)}.  Metadata
#'   are stored as attributes; access column \code{j} with \code{out[[j]]}.
#' @seealso \code{\link{rNormalGamma_reg}}, \code{\link{rindepNormalGamma_reg_multi}}
#' @example inst/examples/Ex_rNormalGamma_reg_multi.R
#' @export
rNormalGamma_reg_multi <- function(n,
                                   y,
                                   x,
                                   prior_list,
                                   offset = NULL,
                                   weights = 1,
                                   family = gaussian(),
                                   Gridtype = 2,
                                   n_envopt = NULL,
                                   use_parallel = TRUE,
                                   use_opencl = FALSE,
                                   verbose = FALSE,
                                   progbar = TRUE) {
  call <- match.call()
  inp <- .mrglmb_check_inputs(y, x, prior_list)
  prior_lists <- .mrglmb_normalize_prior_lists(
    prior_list, inp$l1, inp$p, .validate_normal_gamma_prior_list
  )
  n_draw <- .mrglmb_n_draw(n)

  block_results <- vector("list", inp$l1)
  for (j in seq_len(inp$l1)) {
    block_results[[j]] <- rNormalGamma_reg(
      n = n_draw,
      y = inp$y_mat[, j],
      x = inp$x,
      prior_list = prior_lists[[j]],
      offset = offset,
      weights = weights,
      family = family,
      Gridtype = Gridtype,
      n_envopt = n_envopt,
      use_parallel = use_parallel,
      use_opencl = use_opencl,
      verbose = verbose,
      progbar = progbar && (j == 1L)
    )
  }

  .mrglmb_assemble(
    block_results,
    inp$coef_names,
    call,
    inp$y_mat,
    inp$x,
    inp$l1,
    inp$p,
    prior_lists,
    inp$pred_names
  )
}
