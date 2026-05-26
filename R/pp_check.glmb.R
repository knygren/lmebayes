#' @title Posterior predictive checks (\code{pp_check}) for \code{glmb} fits
#'
#' @description
#' Thin wrapper around [bayesplot::pp_check()] for objects of class \code{glmb}:
#' builds posterior predictive replications with \code{\link[=simulate.glmb]{simulate.glmb}()}
#' and passes
#' \code{y} and \code{yrep} to the default \code{pp_check} method in \pkg{bayesplot}.
#'
#' @param object A fitted model of class \code{glmb}.
#' @param ndraws Number of posterior draws to retain for \code{yrep} (first rows
#'   of \code{predict(object, type = "response")}). The default is \code{100}.
#'   Use \code{NULL} for all draws (can be slow or memory-heavy for large fits).
#' @param fun Passed to [bayesplot::pp_check()] as the \code{fun} argument: a
#'   `ppc_*` plotting function or its name without the \code{ppc_} prefix (see
#'   [bayesplot::PPC-overview]).
#' @param ... Additional arguments passed to [bayesplot::pp_check()].
#'
#' @details
#' For binomial models with a two-column matrix response \code{cbind(successes,
#' failures)}, both \code{y} and \code{yrep} are expressed on the **proportion**
#' scale so they match \code{\link[=simulate.glmb]{simulate.glmb}()}.
#'
#' Call this function as `bayesplot::pp_check(object, ...)` (or attach
#' \pkg{bayesplot} first) so the S3 generic is visible.
#'
#' @return
#' A \pkg{ggplot2} object (from the \pkg{bayesplot} \code{ppc_*} function used).
#'
#' @seealso \code{\link[=simulate.glmb]{simulate.glmb}},
#' \code{\link[=predict.glmb]{predict.glmb}},
#' \code{\link[bayesplot:pp_check]{bayesplot::pp_check}}
#' @importFrom bayesplot pp_check
#' @method pp_check glmb
#' @export
pp_check.glmb <- function(object, ndraws = 100L, fun = "dens_overlay", ...) {
  y <- .glmb_pp_y(object)
  pred <- predict(object, type = "response")
  n <- nrow(pred)
  if (!is.null(ndraws)) {
    nd <- max(1L, min(as.integer(ndraws), n))
    pred <- pred[seq_len(nd), , drop = FALSE]
  }
  yrep <- stats::simulate(
    object,
    pred = pred,
    prior.weights = object$prior.weights
  )
  pp_fun <- if (is.function(fun)) {
    fun
  } else {
    nm <- as.character(fun)[1L]
    if (!startsWith(nm, "ppc_")) nm <- paste0("ppc_", nm)
    get(nm, envir = asNamespace("bayesplot"))
  }
  bayesplot::pp_check(y, yrep, fun = pp_fun, ...)
}

## Prepare observed outcome on the same scale as simulate.glmb() for posterior
## predictive checks.
## @param object A fitted `glmb` object.
## @return Numeric vector (length = number of observations).
## @noRd
.glmb_pp_y <- function(object) {
  y <- object$y
  fam <- object$family$family
  if (fam %in% c("binomial", "quasibinomial") && is.matrix(y)) {
    rs <- rowSums(y)
    rs[rs == 0] <- NA_real_
    return(y[, 1L] / rs)
  }
  if (is.matrix(y) && ncol(y) == 1L) {
    return(as.numeric(y))
  }
  as.numeric(y)
}
