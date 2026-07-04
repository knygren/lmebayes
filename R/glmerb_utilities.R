#' @noRd
.lmerb_reference_fit <- function(object) {
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
  glmbayesCore:::.two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}
