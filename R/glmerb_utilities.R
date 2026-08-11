#' Human-readable sampler-engine label for \code{print()}/\code{summary()} headers
#'
#' \code{sim_method_used} is \code{"DEFAULT"} exactly when the exact-iid
#' engine ran (\code{\link[lmebayesCore]{rLMMNormal_reg_known_vcov_iid}} /
#' \code{\link[lmebayesCore]{rLMMNormal_joint_iid}}), and
#' \code{"TWO_BLOCK_GIBBS"} whenever the two-block Gibbs engine ran instead
#' (either because \code{sim_method = "TWO_BLOCK_GIBBS"} was requested, or
#' because the model has no iid engine at all -- e.g. any
#' \code{dIndependent_Normal_Gamma} Block~2 component). \code{NULL}/missing
#' (e.g. \code{simulate = FALSE}, no sampler ran) also reports
#' \code{"two-block Gibbs"}, matching pre-\code{sim_method} behavior.
#' @noRd
.lmerb_engine_label <- function(sim_method_used) {
  if (identical(sim_method_used, "DEFAULT")) "exact iid" else "two-block Gibbs"
}

#' Reference fit backing calibration/summary tables for an \code{lmerb}/
#' \code{glmerb} object.
#'
#' Prefers \code{object$glmmTMB} (the per-group-dispersion \code{glmmTMB}
#' reference; only non-\code{NULL} when \code{dispformula} requested
#' per-group dispersion) when present, since \code{object$lmer}/
#' \code{object$glmer} is always the plain pooled-dispersion fit and would
#' misreport per-group residual variance / fixed-effect SEs for those
#' models. Falls back to \code{object$lmer}/\code{object$glmer} otherwise.
#' @noRd
.lmerb_reference_fit <- function(object) {
  if (!is.null(object$glmmTMB)) {
    return(object$glmmTMB)
  }
  if (inherits(object, "glmerb")) {
    if (is.null(object$glmer)) {
      stop("glmerb fit is missing component 'glmer'.", call. = FALSE)
    }
    return(object$glmer)
  }
  if (is.null(object$lmer)) {
    stop("lmerb fit is missing component 'lmer'.", call. = FALSE)
  }
  object$lmer
}


#' @noRd
.lmebayes_stage_v2_fixef <- function(
    out,
    fixef_mode,
    fixef_init,
    re_names,
    group_levels,
    n
) {
  x <- list(
    fixef_draws            = out$fixef_draws,
    coefficients           = out$coefficients,
    dispersion_fixef_draws = out$dispersion_fixef_draws,
    iters_fixef_draws      = out$iters_fixef_draws,
    iters_ranef_draws      = out$iters_ranef_draws,
    mu_all_last            = out$mu_all_last,
    re_coef_names          = re_names,
    group_levels           = group_levels,
    n                      = n
  )
  lmebayesCore:::.two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}

