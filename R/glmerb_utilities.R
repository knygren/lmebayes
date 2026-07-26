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

#' ICM at fixed variance-component plug-ins (\code{simulate = FALSE} path)
#' @noRd
.lmebayes_icm_at_fixed_vc <- function(design, prior, family) {
  measurement_prior_list <- list(
    Sigma_ranef      = prior$Sigma_ranef,
    prior_list       = prior$prior_list,
    dispersion_ranef = prior$dispersion_ranef
  )
  re_names <- design$re_coef_names
  is_gaussian <- identical(family$family, "gaussian")
  if (is_gaussian) {
    pm <- lmebayesCore::lmerb_posterior_mean(
      design                 = design,
      measurement_prior_list = measurement_prior_list
    )
    icm_label <- "ICM mean"
  } else {
    pm <- lmebayesCore::glmerb_posterior_mode(
      design                 = design,
      family                 = family,
      measurement_prior_list = measurement_prior_list
    )
    icm_label <- "ICM mode"
  }
  fixef_init <- lapply(prior$prior_list, `[[`, "mu_fixef")
  names(fixef_init) <- re_names
  tau2_mode <- stats::setNames(
    vapply(re_names, function(k) prior$prior_list[[k]]$dispersion_fixef, numeric(1)),
    re_names
  )
  list(
    fixef      = pm$fixef,
    fixef_init = fixef_init,
    b_mean     = pm$b_mean,
    icm_info   = list(
      converged  = pm$converged,
      iterations = pm$iterations,
      delta      = pm$delta
    ),
    icm_label  = icm_label,
    joint_mode = FALSE,
    sigma2     = prior$dispersion_ranef,
    tau2       = tau2_mode
  )
}

#' Print Block~2 ICM table for \code{simulate = FALSE}
#' @noRd
.lmebayes_print_icm_simulate_false <- function(prior, re_names, icm, header) {
  lmebayesCore:::.lmebayes_print_icm_fixef_table(
    prior_list = prior$prior_list,
    re_names   = re_names,
    fixef_icm  = icm$fixef,
    icm_info   = icm$icm_info,
    ref_label  = "prior mean",
    icm_label  = icm$icm_label,
    conv_label = "ICM",
    header     = header,
    verbose    = TRUE
  )
}
