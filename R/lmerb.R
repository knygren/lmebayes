#' Bayesian linear mixed model fit (draft)
#'
#' Entry point for \pkg{lmebayes} models with an \code{lmer}-like interface,
#' analogous to \code{\link{lmb}} and \code{\link{glmb}} for fixed-effects models.
#'
#' Current implementation performs Block 1 sampling: \code{n} iid draws of
#' group-level random effects \eqn{b_j} conditional on the supplied measurement
#' priors and iter-0 hyper means from \code{measurement_prior_list}.
#' A full two-block Gibbs sampler is not yet implemented.
#'
#' @param formula Mixed-model formula. Must match
#'   \code{measurement_prior_list$formula} when both are supplied.
#' @param measurement_prior_list Required object from
#'   \code{\link{Prior_Setup_lmebayes}}. Supplies \code{design},
#'   \code{dispersion_ranef}, \code{Sigma_ranef}, and iter-0 \code{fixef}
#'   means via \code{prior_list}. Not created inside \code{lmerb}; call
#'   \code{\link{Prior_Setup_lmebayes}} explicitly first.
#' @param n Number of iid draws per group (default \code{1000L}, as in \code{\link{lmb}}).
#' @param data Optional data frame; must be consistent with
#'   \code{measurement_prior_list$design}.
#' @param REML Logical; reserved for \code{lmer}-like API compatibility (ignored).
#' @param control \code{\link[lme4]{lmerControl}} settings; reserved (ignored).
#' @param start Optional starting values; reserved (ignored).
#' @param verbose Verbosity flag; reserved (ignored).
#' @param subset Optional subset; reserved (ignored).
#' @param weights Optional weights; reserved (ignored).
#' @param na.action Missing-data handler; reserved (ignored).
#' @param offset Optional offset; reserved (ignored).
#' @param contrasts Optional contrasts; reserved (ignored).
#' @param devFunOnly If \code{TRUE}, return deviance function only; reserved (ignored).
#' @param fixef Optional named list of hyper-parameter vectors (Block 2 state).
#'   When \code{NULL} (default), iter-0 means are taken from
#'   \code{measurement_prior_list$prior_list$mu_fixef}.
#' @param seed Optional; sets the RNG seed before sampling.
#' @param ... Reserved for future use.
#' @return Object of class \code{"lmerb"}: a list of three components:
#'   \describe{
#'     \item{\code{model_setup}}{The \code{\link{model_setup}} object (from
#'       \code{measurement_prior_list$design}).}
#'     \item{\code{lmer}}{The \code{\link[lme4]{lmer}} fit on full-rank groups
#'       (from \code{measurement_prior_list$fit_fr}).}
#'     \item{\code{coefficients}}{\code{data.frame} with \code{n * J} rows and
#'       \code{2 + p_re} columns: \code{draw}, the grouping-factor column
#'       (\code{design$group_name}), and one column per random-effect variable
#'       (\code{design$re_coef_names}).}
#'   }
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{build_mu_all}},
#'   \code{\link[glmbayesCore]{block_rNormalReg}},
#'   \code{\link{lmb}}, \code{\link{glmb}}
#' @export
lmerb <- function(
    formula,
    measurement_prior_list,
    n = 1000L,
    data = NULL,
    REML = TRUE,
    control = lme4::lmerControl(),
    start = NULL,
    verbose = 0L,
    subset,
    weights,
    na.action,
    offset,
    contrasts = NULL,
    devFunOnly = FALSE,
    fixef = NULL,
    seed = NULL,
    ...
) {
  if (missing(measurement_prior_list) || is.null(measurement_prior_list)) {
    stop(
      "'measurement_prior_list' is required. ",
      "Build it with Prior_Setup_lmebayes() and pass the result to lmerb().",
      call. = FALSE
    )
  }
  if (!inherits(measurement_prior_list, "lmebayes_prior_setup")) {
    stop(
      "'measurement_prior_list' must be an object from Prior_Setup_lmebayes().",
      call. = FALSE
    )
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }

  if (missing(formula)) {
    formula <- measurement_prior_list$formula
  } else if (!identical(deparse(formula), deparse(measurement_prior_list$formula))) {
    stop(
      "'formula' does not match measurement_prior_list$formula.",
      call. = FALSE
    )
  }

  design <- measurement_prior_list$design
  if (!inherits(design, "model_setup")) {
    stop("measurement_prior_list$design must be a model_setup object.", call. = FALSE)
  }

  if (is.null(fixef)) {
    fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
    names(fixef) <- design$re_coef_names
  }

  mu_all <- build_mu_all(design, fixef)$mu_all

  Sigma_ranef <- measurement_prior_list$Sigma_ranef
  if (is.null(Sigma_ranef)) {
    stop("measurement_prior_list must contain 'Sigma_ranef'.", call. = FALSE)
  }
  P <- solve(Sigma_ranef)

  dispersion <- measurement_prior_list$dispersion_ranef
  if (is.null(dispersion)) {
    stop("measurement_prior_list must contain 'dispersion_ranef'.", call. = FALSE)
  }

  prior_list <- list(
    mu         = mu_all,
    P          = P,
    dispersion = dispersion,
    ddef       = FALSE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  block_args <- list(
    n          = 1L,
    y          = design$y,
    x          = design$Z,
    block      = design$groups,
    prior_list = prior_list
  )

  grp_col  <- design$group_name
  re_names <- design$re_coef_names
  coef_cols <- c("draw", grp_col, re_names)

  draw_rows <- vector("list", n)

  for (i in seq_len(n)) {
    block_i <- do.call(block_rNormalReg, block_args)
    b_i <- block_i$coefficients
    if (is.null(rownames(b_i))) {
      rownames(b_i) <- block_i$block_info$ids
    }
    colnames(b_i) <- re_names

    J_i <- nrow(b_i)
    draw_df <- data.frame(
      draw = rep(i, J_i),
      stringsAsFactors = FALSE
    )
    draw_df[[grp_col]] <- rownames(b_i)
    for (nm in re_names) {
      draw_df[[nm]] <- b_i[, nm]
    }
    draw_rows[[i]] <- draw_df
  }

  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      model_setup  = design,
      lmer         = measurement_prior_list$fit_fr,
      coefficients = coefficients
    ),
    class = c("lmerb", "list")
  )
}
