## Internal helpers for the `dispformula` argument on lmerb() / glmerb().
## `dispformula` gates which `dispersion_ranef` shape (and therefore which
## glmbayesCore sampler route) is accepted; it does not touch the embedded
## lme4::lmer/glmer reference fit (x$lmer / x$glmer stay unchanged in every
## case). When per-group dispersion is requested it additionally fits a
## glmmTMB::glmmTMB() diagnostic reference model, stored separately.

#' Validate `dispformula` against the grouping factor, family, and the
#' resolved `dispersion_ranef` mode.
#'
#' `dispformula = ~1` (pooled) accepts `disp_mode` \code{"none"}, \code{"fixed"},
#' or \code{"gamma"} (no dispersion parameter, a fixed scalar, or a single
#' pooled \code{dGamma()}). `dispformula = ~<group_name>` (per-group) requires
#' `disp_mode == "gamma_list"` (a named list from
#' \code{\link[glmbayesCore:dGamma_list.lmebayes_prior_setup]{dGamma_list}()}).
#' No other `dispformula` or mode combination is accepted.
#'
#' @param dispformula One-sided formula: `~1` or `~<group_name>`.
#' @param group_name Character scalar, the random-effects grouping factor
#'   (`design$group_name`).
#' @param family A \code{\link[stats]{family}} object.
#' @param disp_mode `"none"`, `"fixed"`, `"gamma"`, or `"gamma_list"` (from
#'   `prior$dispersion_mode`, already resolved by
#'   `glmbayesCore:::.lmebayes_resolve_dispersion_ranef()` inside
#'   `.lmebayes_priors_from_pfamily_list()`).
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
    if (!identical(disp_mode, "gamma_list")) {
      stop(
        "dispformula = ~", group_name, " requires 'dispersion_ranef' to be ",
        "a dGamma_list(...) result (one dGamma() per ", group_name,
        " level); got dispersion_ranef mode = '", disp_mode, "'. Use ",
        "dispformula = ~1 for a fixed scalar or a pooled dGamma().",
        call. = FALSE
      )
    }
  } else if (identical(disp_mode, "gamma_list")) {
    stop(
      "'dispersion_ranef' is a dGamma_list(...) (per-group); this requires ",
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
#' @noRd
.lmebayes_fit_glmmtmb_dispersion <- function(
    formula, data, family, dispformula, REML, mer_optional_args = list(), ...
) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for dispformula = ", deparse(dispformula),
      " (per-group residual dispersion reference fit). Install it with ",
      "install.packages(\"glmmTMB\").",
      call. = FALSE
    )
  }
  tmb_args <- c(
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
  do.call(glmmTMB::glmmTMB, tmb_args)
}
