## Internal helpers for the `dispformula` argument on lmerb() / glmerb().
## `dispformula` gates which `dispersion_ranef` shape (and therefore which
## lmebayesCore sampler route) is accepted; it does not touch the embedded
## lme4::lmer/glmer reference fit (x$lmer / x$glmer stay unchanged in every
## case). When per-group dispersion is requested it additionally fits a
## glmmTMB::glmmTMB() diagnostic reference model, stored separately.

#' Validate `dispformula` against the grouping factor, family, and the
#' resolved `dispersion_ranef` mode.
#'
#' `dispformula = ~1` (pooled) accepts `disp_mode` \code{"none"}, \code{"fixed"},
#' or \code{"gamma"} (no dispersion parameter, a fixed scalar, or a single
#' pooled \code{dGamma()}). `dispformula = ~<group_name>` (per-group) requires
#' `disp_mode` to be `"gamma_list"` (a named list from
#' \code{\link[lmebayesCore:dGamma_list.lmebayes_prior_setup]{dGamma_list}()})
#' or `"fixed_vector"` (a plain named numeric vector, one known fixed
#' dispersion per group). No other `dispformula` or mode combination is
#' accepted.
#'
#' @param dispformula One-sided formula: `~1` or `~<group_name>`.
#' @param group_name Character scalar, the random-effects grouping factor
#'   (`design$group_name`).
#' @param family A \code{\link[stats]{family}} object.
#' @param disp_mode `"none"`, `"fixed"`, `"gamma"`, `"gamma_list"`, or
#'   `"fixed_vector"` (from `prior$dispersion_mode`, already resolved by
#'   `lmebayesCore:::.lmebayes_resolve_dispersion_ranef()` inside
#'   `lmebayesCore::priors_from_pfamily_list()`).
#' @return `"pooled"` or `"group"`.
#' @noRd
.lmebayes_validate_dispformula <- function(dispformula, group_name, family, disp_mode) {
  if (!inherits(dispformula, "formula") || length(dispformula) != 2L) {
    stop(
      "'dispformula' must be a one-sided formula, either ~1 (pooled) or ~",
      group_name, " (per-group).",
      call. = FALSE
    )
  }
  vars <- all.vars(dispformula)
  if (length(vars) == 0L) {
    kind <- "pooled"
  } else if (length(vars) == 1L && identical(vars, group_name)) {
    kind <- "group"
  } else {
    stop(
      "'dispformula' must be ~1 (pooled) or ~", group_name,
      " (per-group, matching the random-effects grouping factor); got ",
      deparse(dispformula), ".",
      call. = FALSE
    )
  }

  if (identical(kind, "group")) {
    has_dispersion <- family$family %in%
      c("gaussian", "Gamma", "quasipoisson", "quasibinomial")
    if (!has_dispersion) {
      stop(
        "'dispformula' must be ~1 for family = ", family$family,
        "() (no observation-level dispersion parameter).",
        call. = FALSE
      )
    }
    if (!disp_mode %in% c("gamma_list", "fixed_vector")) {
      stop(
        "dispformula = ~", group_name, " requires 'dispersion_ranef' to be ",
        "a dGamma_list(...) result or a named numeric vector (one value per ",
        group_name, " level); got dispersion_ranef mode = '", disp_mode,
        "'. Use dispformula = ~1 for a fixed scalar or a pooled dGamma().",
        call. = FALSE
      )
    }
  } else if (disp_mode %in% c("gamma_list", "fixed_vector")) {
    stop(
      "'dispersion_ranef' is per-group (", disp_mode, "); this requires ",
      "dispformula = ~", group_name, ", not dispformula = ~1.",
      call. = FALSE
    )
  }

  kind
}

#' Fit a \code{glmmTMB::glmmTMB()} reference model with per-group residual
#' dispersion, for diagnostics only.
#'
#' \code{x$lmer}/\code{x$glmer} (\pkg{lme4}) stay the primary embedded
#' reference fit in every case; this is stored separately as
#' \code{x$dispersion_fit} and is only fit when \code{dispformula} requests
#' per-group dispersion (see \code{\link{.lmebayes_validate_dispformula}}).
#'
#' Thin wrapper around lmebayesCore's
#' \code{.lmebayes_fit_glmmtmb_reference()} (the same helper
#' \code{Prior_Setup_lmebayes()} uses to calibrate priors when
#' \code{dispformula} requests per-group dispersion), so \code{lmerb()} only
#' re-fits glmmTMB here when the caller's \code{dispersion_ranef} was not
#' already produced by
#' \code{dGamma_list(Prior_Setup_lmebayes(..., dispformula = dispformula))}
#' (which carries the fit forward as an attribute; see \code{lmerb()}).
#'
#' @noRd
.lmebayes_fit_glmmtmb_dispersion <- function(
    formula, data, family, dispformula, REML, mer_optional_args = list(), ...
) {
  do.call(
    lmebayesCore:::.lmebayes_fit_glmmtmb_reference,
    c(
      list(
        formula     = formula,
        data        = data,
        family      = family,
        dispformula = dispformula,
        REML        = isTRUE(REML)
      ),
      mer_optional_args,
      list(...)
    )
  )
}
