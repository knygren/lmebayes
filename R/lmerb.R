#' Bayesian linear mixed model fit (draft)
#'
#' Entry point for \pkg{lmebayes} models with an \code{lmer}-like interface,
#' analogous to \code{\link[glmbayes]{lmb}} and \code{\link[glmbayes]{glmb}} for
#' fixed-effects models. Initial implementation delegates to
#' \code{\link{model_setup}}; sampling and prior integration will be added in
#' later versions.
#'
#' @param formula Mixed-model formula (single grouping factor; same constraints
#'   as \code{\link{model_setup}}).
#' @param data Optional data frame.
#' @param measurement_prior_list Optional prior specification for the
#'   measurement model (Block 1: residual variance and random-effect structure).
#'   Reserved for future use; currently ignored.
#' @param REML Logical; passed to \code{\link[lme4]{lmer}} inside
#'   \code{\link{model_setup}}.
#' @param control \code{\link[lme4]{lmerControl}} settings.
#' @param start Optional starting values for the inner optimization.
#' @param verbose Passed to \code{\link[lme4]{lmer}}.
#' @param subset,weights,na.action,offset,contrasts Passed to
#'   \code{\link[lme4]{lmer}}.
#' @param devFunOnly If \code{TRUE}, return the deviance function only.
#' @param ... Passed to \code{\link{model_setup}}.
#' @return Object of class \code{"model_setup"} (for now): design matrices,
#'   \code{lmer} fits, variance components, and rank diagnostics from
#'   \code{\link{model_setup}}.
#' @seealso \code{\link{model_setup}}, \code{\link{Prior_Setup_lmebayes}},
#'   \code{\link{build_mu_all}},,
#'   \code{\link[glmbayes]{lmb}}, \code{\link[glmbayes]{glmb}}
#' @export
lmerb <- function(
    formula,
    data = NULL,
    measurement_prior_list = NULL,
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
    ...
) {
  model_setup(
    formula = formula,
    data = data,
    REML = REML,
    control = control,
    start = start,
    verbose = verbose,
    subset = subset,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    devFunOnly = devFunOnly,
    ...
  )
}
