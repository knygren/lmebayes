#' Multi-response Bayesian regression and simulation
#'
#' @description
#' These functions run the corresponding single-response sampler once per column
#' of a matrix \code{y}, sharing the same design matrix \code{x} (as \code{lm}
#' with \code{cbind} responses). Each returns a named list of class
#' \code{"mrglmb"}; element \code{j} is the fit for column \code{j} of
#' \code{y}. Use \code{\link{summary.mrglmb}} for column-wise summaries.
#'
#' @details
#' \describe{
#'   \item{\code{multi_rlmb}}{
#'     Same arguments as \code{\link[glmbayes]{rlmb}} except \code{pfamily} is replaced by
#'     \code{pfamily_list} (length \code{ncol(y)} of \code{pfamily} objects).
#'     Each element is class \code{"rlmb"} (and \code{"rglmb"}).
#'   }
#'   \item{\code{multi_rNormal_reg}}{
#'     Same arguments as \code{\link[glmbayesCore]{rNormal_reg}} except \code{prior_list} is a
#'     list of per-column prior lists (\code{mu}, \code{Sigma} or \code{P}, optional
#'     \code{dispersion}).
#'   }
#'   \item{\code{multi_rNormalGamma_reg}}{
#'     Same arguments as \code{\link[glmbayesCore]{rNormalGamma_reg}} except \code{prior_list} is a
#'     list of per-column prior lists (\code{mu}, \code{Sigma} or \code{P},
#'     \code{shape}, \code{rate}).
#'   }
#'   \item{\code{multi_rindepNormalGamma_reg}}{
#'     Same arguments as \code{\link[glmbayesCore]{rindepNormalGamma_reg}} except \code{prior_list}
#'     is a list of per-column prior lists (\code{mu}, \code{Sigma}, \code{shape},
#'     \code{rate}, optional dispersion bounds).
#'   }
#'   \item{\code{\link[glmbayes]{multi_prior_setup}}}{
#'     In \pkg{glmbayes}: same arguments as \code{\link[glmbayes]{Prior_Setup}}, but the
#'     formula left-hand side may be several responses (\code{cbind(...)}). Returns a
#'     named list of \code{"PriorSetup"} objects (one per response column).
#'   }
#' }
#'
#' @return
#' A named list of class \code{"mrglmb"}. Metadata (\code{call}, \code{y},
#' \code{x}, \code{l1}, \code{p}, \code{coef_names}, \code{pred_names}) are
#' attributes; per-column priors are in \code{attr(..., "prior_lists")} or
#' \code{attr(..., "pfamily_lists")} for \code{multi_rlmb}.
#'
#' @seealso
#' \code{\link{summary.mrglmb}}, \code{\link[glmbayes]{lmb}} (\code{cbind} responses), \code{\link[glmbayes]{Prior_Setup}},
#' \code{\link[glmbayes]{rlmb}}, \code{\link[glmbayesCore]{rNormal_reg}}, \code{\link[glmbayesCore]{rNormalGamma_reg}},
#' \code{\link[glmbayesCore]{rindepNormalGamma_reg}}
#'
#' @name multi_rlmb
#' @aliases multi_rlmb multi_rNormalGamma_reg multi_rNormal_reg
#'   multi_rindepNormalGamma_reg
#' @example inst/examples/Ex_multi_rlmb.R
NULL

#' @describeIn multi_rlmb Gaussian \code{\link[glmbayes]{rlmb}} simulation with multiple responses.
#' @inheritParams glmbayesCore::rlmb
#' @param pfamily_list List of length \code{ncol(y)} of \code{pfamily} objects.
#' @family modelfuns
#' @export
multi_rlmb <- function(n = 1,
                       y,
                       x,
                       pfamily_list,
                       offset = NULL,
                       weights = NULL,
                       Gridtype = 2,
                       n_envopt = NULL,
                       use_parallel = TRUE,
                       use_opencl = FALSE,
                       verbose = FALSE,
                       progbar = FALSE) {
  call <- match.call()
  inp <- .mrglmb_check_inputs(y, x, pfamily_list, spec_name = "pfamily_list")
  pfamily_lists <- .mrglmb_normalize_pfamily_lists(
    pfamily_list, inp$l1, inp$p, .validate_pfamily_for_rlmb
  )
  n_draw <- .mrglmb_n_draw(n)

  block_results <- vector("list", inp$l1)
  for (j in seq_len(inp$l1)) {
    block_results[[j]] <- rlmb(
      n = n_draw,
      y = inp$y_mat[, j],
      x = inp$x,
      pfamily = pfamily_lists[[j]],
      offset = offset,
      weights = weights,
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
    prior_lists = NULL,
    inp$pred_names,
    pfamily_lists = pfamily_lists
  )
}

#' @describeIn multi_rlmb Normal-prior regression with multiple responses.
#' @inheritParams glmbayesCore::rNormal_reg
#' @param prior_list List of length \code{ncol(y)} of per-column prior lists.
#' @family simfuncs
#' @export
multi_rNormal_reg <- function(n,
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
    prior_list, inp$l1, inp$p, .validate_normal_prior_list
  )
  n_draw <- .mrglmb_n_draw(n)

  block_results <- vector("list", inp$l1)
  for (j in seq_len(inp$l1)) {
    block_results[[j]] <- rNormal_reg(
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

#' @describeIn multi_rlmb Normal-prior regression with multiple responses and
#'   per-column design matrices.
#'
#' Extension of \code{multi_rNormal_reg} where \code{x} accepts either a
#' single shared design matrix (identical behaviour to \code{multi_rNormal_reg})
#' or a \strong{list of matrices} with one entry per response column.  The list
#' path is the key new capability: it allows each column of \code{y} to have a
#' different number of predictors, which is required for Block 2 of the
#' two-block Gibbs sampler in \code{\link{lmerb}}, where the hyper design
#' matrices \code{design$X_hyper[[k]]} can have differing column dimensions
#' across random-effect components \eqn{k}.
#'
#' @inheritParams glmbayesCore::rNormal_reg
#' @param n Number of draws to request from \code{\link[glmbayesCore]{rNormal_reg}}
#'   for each column of \code{y}.  For a single Gibbs step inside
#'   \code{\link{lmerb}} set \code{n = 1L}; the first row of each
#'   \code{$coefficients} matrix is then the draw.
#' @param y Numeric matrix (or any object coercible via \code{as.matrix}) with
#'   one column per random-effect component \eqn{k}.  In \code{\link{lmerb}}
#'   Block 2 this is the \eqn{J \times p_{\mathrm{re}}} matrix of current
#'   Block 1 draws \eqn{b_j}, where rows are groups and columns correspond to
#'   RE coefficients (\code{design$re_coef_names}).
#' @param x Either
#'   \itemize{
#'     \item a numeric matrix (shared across all response columns) — in this
#'       case the function behaves identically to \code{multi_rNormal_reg} and
#'       returns an \code{"mrglmb"} object; or
#'     \item a \strong{list} of numeric matrices, one per column of \code{y},
#'       where \code{x[[k]]} is the \eqn{J \times q_k} design matrix for RE
#'       component \eqn{k}.  Supply \code{design$X_hyper} from
#'       \code{\link{model_setup}} for the Block 2 path in \code{\link{lmerb}}.
#'       Predictor dimensions \eqn{q_k} may differ across columns.
#'   }
#' @param prior_list
#'   \itemize{
#'     \item \emph{Shared-x path} (matrix \code{x}): a list of length
#'       \code{ncol(y)}, each element a prior list with components \code{mu}
#'       (length \eqn{p}), \code{Sigma} or \code{P} (\eqn{p \times p}), and
#'       optionally \code{dispersion}.
#'     \item \emph{List-x path}: same structure, but for column \eqn{k} the
#'       \code{mu} and \code{Sigma}/\code{P} must conform to
#'       \code{ncol(x[[k]])} = \eqn{q_k}.  For \code{\link{lmerb}} Block 2,
#'       supply the renamed fields from
#'       \code{\link{Prior_Setup_lmebayes}}\code{$prior_list[[k]]}: use
#'       \code{mu = mu_fixef}, \code{Sigma = Sigma_fixef}, and
#'       \code{dispersion = dispersion_fixef}.
#'   }
#' @return
#'   \describe{
#'     \item{Shared-x path (matrix \code{x})}{An \code{"mrglmb"} S3 object
#'       identical in structure to the return value of
#'       \code{\link{multi_rNormal_reg}}.}
#'     \item{List-x path (list \code{x})}{A plain named list of length
#'       \code{ncol(y)}, named by \code{colnames(y)} (or \code{"Y1"},
#'       \code{"Y2"}, \ldots{} if \code{y} has no column names).  Each element
#'       is the \code{\link[glmbayesCore]{rNormal_reg}} output for column
#'       \eqn{k}; access the single Gibbs draw via
#'       \code{result[[k]]$coefficients[1L, ]}.}
#'   }
#' @family simfuncs
#' @export
multi_rNormal_reg_v2 <- function(n,
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
  call  <- match.call()
  n_draw <- .mrglmb_n_draw(n)

  x_is_list <- is.list(x) && !is.data.frame(x)

  if (x_is_list) {
    # ------------------------------------------------------------------
    # Per-column design path: x is a list of matrices, p_k may vary.
    # Returns a plain named list (Option R2) -- not an mrglmb.
    # ------------------------------------------------------------------
    y_mat <- as.matrix(y)
    l1    <- ncol(y_mat)
    n_obs <- nrow(y_mat)

    if (l1 < 1L) stop("y must have at least one column.", call. = FALSE)

    coef_names <- colnames(y_mat)
    if (is.null(coef_names) || length(coef_names) != l1) {
      coef_names <- paste0("Y", seq_len(l1))
    }

    if (length(x) != l1) {
      stop(
        "When x is a list, length(x) must equal ncol(y) = ", l1, ".",
        call. = FALSE
      )
    }
    x_list <- lapply(x, as.matrix)
    for (j in seq_len(l1)) {
      if (nrow(x_list[[j]]) != n_obs) {
        stop(
          "nrow(x[[", j, "]]) (", nrow(x_list[[j]]),
          ") must equal nrow(y) (", n_obs, ").",
          call. = FALSE
        )
      }
    }
    p_vec <- vapply(x_list, ncol, integer(1L))

    if (!is.list(prior_list)) {
      stop(
        "prior_list must be a list of length ncol(y) = ", l1, ".",
        call. = FALSE
      )
    }
    if (!is.null(prior_list$mu) || !is.null(prior_list$Sigma)) {
      stop(
        "prior_list must be a list of prior_list objects (one per column ",
        "of y), not a single prior_list with components mu and Sigma.",
        call. = FALSE
      )
    }
    if (length(prior_list) != l1) {
      stop(
        "length(prior_list) must equal ncol(y) = ", l1, ".",
        call. = FALSE
      )
    }
    prior_lists <- lapply(seq_len(l1), function(j) {
      .validate_normal_prior_list(prior_list[[j]], j = j, p = p_vec[j])
    })

    block_results <- vector("list", l1)
    names(block_results) <- coef_names
    for (j in seq_len(l1)) {
      block_results[[j]] <- rNormal_reg(
        n            = n_draw,
        y            = y_mat[, j],
        x            = x_list[[j]],
        prior_list   = prior_lists[[j]],
        offset       = offset,
        weights      = weights,
        family       = family,
        Gridtype     = Gridtype,
        n_envopt     = n_envopt,
        use_parallel = use_parallel,
        use_opencl   = use_opencl,
        verbose      = verbose,
        progbar      = progbar && (j == 1L)
      )
    }
    block_results

  } else {
    # ------------------------------------------------------------------
    # Shared-design path: x is a single matrix. Identical to
    # multi_rNormal_reg; returns an mrglmb object.
    # ------------------------------------------------------------------
    inp         <- .mrglmb_check_inputs(y, x, prior_list)
    prior_lists <- .mrglmb_normalize_prior_lists(
      prior_list, inp$l1, inp$p, .validate_normal_prior_list
    )

    block_results <- vector("list", inp$l1)
    for (j in seq_len(inp$l1)) {
      block_results[[j]] <- rNormal_reg(
        n            = n_draw,
        y            = inp$y_mat[, j],
        x            = inp$x,
        prior_list   = prior_lists[[j]],
        offset       = offset,
        weights      = weights,
        family       = family,
        Gridtype     = Gridtype,
        n_envopt     = n_envopt,
        use_parallel = use_parallel,
        use_opencl   = use_opencl,
        verbose      = verbose,
        progbar      = progbar && (j == 1L)
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
}

#' @describeIn multi_rlmb Normal--Gamma regression with multiple responses.
#' @inheritParams glmbayesCore::rNormalGamma_reg
#' @param prior_list List of length \code{ncol(y)} of per-column prior lists.
#' @family simfuncs
#' @export
multi_rNormalGamma_reg <- function(n,
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

#' @describeIn multi_rlmb Independent Normal--Gamma regression with multiple responses.
#' @inheritParams glmbayesCore::rindepNormalGamma_reg
#' @param prior_list List of length \code{ncol(y)} of per-column prior lists.
#' @family simfuncs
#' @export
multi_rindepNormalGamma_reg <- function(n,
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
    prior_list, inp$l1, inp$p, .validate_rindep_prior_list
  )
  n_draw <- .mrglmb_n_draw(n)

  block_results <- vector("list", inp$l1)
  for (j in seq_len(inp$l1)) {
    block_results[[j]] <- rindepNormalGamma_reg(
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

# ---- helpers (multi-response) -----------------------------------------------

#' @keywords internal
.mrglmb_check_inputs <- function(y, x, spec_list, spec_name = "prior_list") {
  if (missing(spec_list)) {
    stop("'", spec_name, "' is required.", call. = FALSE)
  }
  y_mat <- as.matrix(y)
  x <- as.matrix(x)
  l1 <- ncol(y_mat)
  if (l1 < 1L) {
    stop("y must have at least one column.", call. = FALSE)
  }
  p <- ncol(x)
  if (p < 1L) {
    stop("x must have at least one column.", call. = FALSE)
  }
  if (nrow(x) != nrow(y_mat)) {
    stop("nrow(x) must equal nrow(y).", call. = FALSE)
  }
  coef_names <- colnames(y_mat)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }
  pred_names <- colnames(x)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }
  list(
    y_mat = y_mat,
    x = x,
    l1 = l1,
    p = p,
    coef_names = coef_names,
    pred_names = pred_names
  )
}

#' @keywords internal
.mrglmb_n_draw <- function(n) {
  n_draw <- if (length(n) > 1L) length(n) else as.integer(n)
  if (!is.finite(n_draw) || n_draw < 1L) {
    stop(
      "'n' must be a positive scalar or a vector whose length defines the number of draws.",
      call. = FALSE
    )
  }
  n_draw
}

#' @keywords internal
.mrglmb_normalize_prior_lists <- function(prior_list, l1, p, validate_fn) {
  if (!is.list(prior_list)) {
    stop(
      "prior_list must be a list of length ncol(y) of per-column prior lists.",
      call. = FALSE
    )
  }
  if (!is.null(prior_list$mu) || !is.null(prior_list$Sigma)) {
    stop(
      "prior_list must be a list of prior_list objects (one per column of y), ",
      "not a single prior_list with components mu and Sigma.",
      call. = FALSE
    )
  }
  if (length(prior_list) != l1) {
    stop("length(prior_list) must equal ncol(y) = ", l1, ".", call. = FALSE)
  }
  lapply(seq_len(l1), function(j) {
    validate_fn(prior_list[[j]], j = j, p = p)
  })
}

#' @keywords internal
.mrglmb_assemble <- function(block_results,
                             coef_names,
                             call,
                             y_mat,
                             x,
                             l1,
                             p,
                             prior_lists,
                             pred_names,
                             pfamily_lists = NULL) {
  outlist <- setNames(block_results, coef_names)
  attr(outlist, "call")       <- call
  attr(outlist, "y")          <- y_mat
  attr(outlist, "x")          <- x
  attr(outlist, "l1")         <- l1
  attr(outlist, "p")          <- p
  attr(outlist, "coef_names") <- coef_names
  attr(outlist, "pred_names") <- pred_names
  if (!is.null(prior_lists)) {
    attr(outlist, "prior_lists") <- prior_lists
  }
  if (!is.null(pfamily_lists)) {
    attr(outlist, "pfamily_lists") <- pfamily_lists
  }
  class(outlist) <- "mrglmb"
  outlist
}

#' @keywords internal
.mrglmb_normalize_pfamily_lists <- function(pfamily_list, l1, p, validate_fn) {
  if (!is.list(pfamily_list)) {
    stop(
      "pfamily_list must be a list of length ncol(y) of per-column pfamily objects.",
      call. = FALSE
    )
  }
  if (!is.null(pfamily_list$pfamily) && !is.null(pfamily_list$prior_list)) {
    stop(
      "pfamily_list must be a list of pfamily objects (one per column of y), ",
      "not a single pfamily with components pfamily and prior_list.",
      call. = FALSE
    )
  }
  if (length(pfamily_list) != l1) {
    stop("length(pfamily_list) must equal ncol(y) = ", l1, ".", call. = FALSE)
  }
  lapply(seq_len(l1), function(j) {
    validate_fn(pfamily_list[[j]], j = j, p = p)
  })
}

#' @keywords internal
.validate_pfamily_for_rlmb <- function(pl, j, p) {
  if (!inherits(pl, "pfamily")) {
    stop("pfamily_list[[", j, "]] must inherit from class \"pfamily\".", call. = FALSE)
  }
  if (is.null(pl$pfamily) || is.null(pl$prior_list) || is.null(pl$simfun)) {
    stop(
      "pfamily_list[[", j, "]] must contain 'pfamily', 'prior_list', and 'simfun'.",
      call. = FALSE
    )
  }
  mu <- pl$prior_list$mu
  if (is.null(mu)) {
    stop("pfamily_list[[", j, "]]$prior_list must contain 'mu'.", call. = FALSE)
  }
  mu <- as.numeric(mu)
  if (length(mu) != p) {
    stop(
      "pfamily_list[[", j, "]]$prior_list$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }
  pl
}

#' @keywords internal
.validate_rindep_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma)) {
    stop("prior_list[[", j, "]] must contain 'Sigma'.", call. = FALSE)
  }
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop("prior_list[[", j, "]] must contain 'shape' and 'rate'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  S <- as.matrix(pl$Sigma)
  if (nrow(S) != p || ncol(S) != p) {
    stop(
      "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
      call. = FALSE
    )
  }
  .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))

  shape <- as.numeric(pl$shape)
  rate <- as.numeric(pl$rate)
  if (length(shape) != 1L || !is.finite(shape)) {
    stop("prior_list[[", j, "]]$shape must be a finite scalar.", call. = FALSE)
  }
  if (length(rate) != 1L || !is.finite(rate)) {
    stop("prior_list[[", j, "]]$rate must be a finite scalar.", call. = FALSE)
  }

  out <- list(mu = mu, Sigma = S, shape = shape, rate = rate)
  if (!is.null(pl$max_disp_perc)) {
    out$max_disp_perc <- as.numeric(pl$max_disp_perc)
  }
  if (!is.null(pl$disp_lower)) {
    out$disp_lower <- pl$disp_lower
  }
  if (!is.null(pl$disp_upper)) {
    out$disp_upper <- pl$disp_upper
  }
  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  out
}

#' @keywords internal
.validate_normal_gamma_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma) && is.null(pl$P)) {
    stop("prior_list[[", j, "]] must contain 'Sigma' or 'P'.", call. = FALSE)
  }
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop("prior_list[[", j, "]] must contain 'shape' and 'rate'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  out <- list(mu = mu)
  if (!is.null(pl$Sigma)) {
    S <- as.matrix(pl$Sigma)
    if (nrow(S) != p || ncol(S) != p) {
      stop(
        "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))
    out$Sigma <- S
  }
  if (!is.null(pl$P)) {
    P <- as.matrix(pl$P)
    if (nrow(P) != p || ncol(P) != p) {
      stop(
        "prior_list[[", j, "]]$P must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(P, label = paste0("prior_list[[", j, "]]$P"))
    out$P <- P
  }

  shape <- as.numeric(pl$shape)
  rate <- as.numeric(pl$rate)
  if (length(shape) != 1L || !is.finite(shape)) {
    stop("prior_list[[", j, "]]$shape must be a finite scalar.", call. = FALSE)
  }
  if (length(rate) != 1L || !is.finite(rate)) {
    stop("prior_list[[", j, "]]$rate must be a finite scalar.", call. = FALSE)
  }
  out$shape <- shape
  out$rate <- rate

  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  if (!is.null(pl$max_disp_perc)) {
    out$max_disp_perc <- as.numeric(pl$max_disp_perc)
  }
  if (!is.null(pl$disp_lower)) {
    out$disp_lower <- pl$disp_lower
  }
  if (!is.null(pl$disp_upper)) {
    out$disp_upper <- pl$disp_upper
  }
  if (!is.null(pl$Precision)) {
    out$Precision <- pl$Precision
  }
  out
}

#' @keywords internal
.validate_normal_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma) && is.null(pl$P)) {
    stop("prior_list[[", j, "]] must contain 'Sigma' or 'P'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  out <- list(mu = mu)
  if (!is.null(pl$Sigma)) {
    S <- as.matrix(pl$Sigma)
    if (nrow(S) != p || ncol(S) != p) {
      stop(
        "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))
    out$Sigma <- S
  }
  if (!is.null(pl$P)) {
    P <- as.matrix(pl$P)
    if (nrow(P) != p || ncol(P) != p) {
      stop(
        "prior_list[[", j, "]]$P must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(P, label = paste0("prior_list[[", j, "]]$P"))
    out$P <- P
  }
  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  if (!is.null(pl$shape)) {
    out$shape <- pl$shape
  }
  if (!is.null(pl$rate)) {
    out$rate <- pl$rate
  }
  if (!is.null(pl$ddef)) {
    out$ddef <- pl$ddef
  }
  out
}

#' @keywords internal
.check_symmetric_pd <- function(M, label) {
  if (!isSymmetric(M)) {
    stop(label, " must be symmetric.", call. = FALSE)
  }
  tol <- 1e-6
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (!all(ev >= -tol * abs(ev[1L]))) {
    stop(label, " is not positive definite.", call. = FALSE)
  }
  invisible(TRUE)
}
