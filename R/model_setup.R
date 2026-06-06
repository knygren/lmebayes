#' Bayesian mixed model setup (single-factor \code{lmer} gate)
#'
#' Wrapper around \code{\link[lme4]{lmer}} for models with exactly one
#' grouping factor. Design matrices come from \code{formula} (including
#' cross-level RE moderation terms). Variance components use \code{vcov_formula}
#' (defaults to \code{\link{lmerb_default_vcov_formula}}): level-2 fixed only,
#' same \code{||} or \code{|} random structure as \code{formula}, without
#' cross-level fixed interactions (so RE moderation is not double-coded).
#'
#' @param formula Mixed-model formula for design extraction and \code{lmer_fit}
#'   (fixed effects / hyper calibration).
#' @param vcov_formula Optional formula for the \code{lmer} fit used to extract
#'   \code{vcov_re}. Defaults to \code{\link{lmerb_default_vcov_formula}(formula)}.
#' @param data Optional data frame.
#' @param REML Logical; passed to \code{\link[lme4]{lmer}}.
#' @param control \code{\link[lme4]{lmerControl}} settings.
#' @param start Optional starting values for the inner optimization.
#' @param verbose Passed to \code{\link[lme4]{lmer}}.
#' @param subset,weights,na.action,offset,contrasts Passed to
#'   \code{\link[lme4]{lmer}}.
#' @param devFunOnly If \code{TRUE}, return the deviance function only.
#' @param ... Passed to design extraction and \code{\link[lme4]{lmer}}.
#' @return Object of class \code{"model_setup"}: \code{y}, \code{Z},
#'   \code{groups}, \code{X_hyper}, \code{formula}, \code{vcov_formula},
#'   \code{lmer_fit} (full formula), \code{lmer_vcov_fit}, \code{varcorr},
#'   \code{vcov_re}, and \code{residual_var}.
#' @seealso \code{\link{extract_re_hyper_matrices}},
#'   \code{\link{lmerb_default_vcov_formula}},
#'   \code{\link{extract_lmer_variance_components}}
#' @export
model_setup <- function(
    formula,
    data = NULL,
    vcov_formula = NULL,
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
  design <- extract_re_hyper_matrices(formula = formula, data = data, ...)
  design$formula <- formula

  if (is.null(vcov_formula)) {
    vcov_formula <- lmerb_default_vcov_formula(
      formula = formula,
      data = data,
      ...
    )
  }
  design$vcov_formula <- vcov_formula

  lmer_args <- list(
    data = data,
    REML = REML,
    control = control,
    verbose = verbose,
    devFunOnly = devFunOnly
  )
  if (!missing(start) && !is.null(start)) {
    lmer_args$start <- start
  }
  if (!missing(subset)) {
    lmer_args$subset <- subset
  }
  if (!missing(weights)) {
    lmer_args$weights <- weights
  }
  if (!missing(na.action)) {
    lmer_args$na.action <- na.action
  }
  if (!missing(offset)) {
    lmer_args$offset <- offset
  }
  if (!missing(contrasts)) {
    lmer_args$contrasts <- contrasts
  }

  fit_full <- do.call(lme4::lmer, c(lmer_args, list(formula = formula), list(...)))
  fit_vcov <- do.call(
    lme4::lmer,
    c(lmer_args, list(formula = vcov_formula), list(...))
  )

  if (lme4::isSingular(fit_vcov)) {
    message(
      "lmer vcov fit is singular -- check VarCorr; ",
      "RE variances may be on boundary."
    )
  }

  vc <- extract_lmer_variance_components(fit_vcov, design$re_coef_names)
  design$lmer_fit <- fit_full
  design$lmer_vcov_fit <- fit_vcov
  design$varcorr <- vc$varcorr
  design$vcov_re <- vc$vcov_re
  design$residual_var <- vc$residual_var

  design
}

#' @export
print.model_setup <- function(x, ...) {
  cat("model setup (group =", x$group_name, ")\n")
  cat("Random coefficients:", paste(x$re_coef_names, collapse = ", "), "\n")
  cat("y: length ", length(x$y), "\n", sep = "")
  cat("Z:", paste(dim(x$Z), collapse = " x "),
      " (columns: ", paste(colnames(x$Z), collapse = ", "), ")\n", sep = "")
  cat("groups: length ", length(x$groups), ", ",
      nlevels(x$groups), " levels\n", sep = "")

  for (coef in x$re_coef_names) {
    Xj <- x$X_hyper[[coef]]
    cat("\nX_hyper[[\"", coef, "\"]] ", paste(dim(Xj), collapse = " x "),
        " ~ 1", sep = "")
    other <- setdiff(colnames(Xj), "(Intercept)")
    if (length(other) > 0L) {
      cat(" + ", paste(other, collapse = " + "), sep = "")
    }
    cat("\n  columns: ", paste(colnames(Xj), collapse = ", "), "\n", sep = "")
  }

  if (nrow(x$re_slope_moderation) > 0L) {
    cat("\nCross-level RE moderation:\n")
    print(x$re_slope_moderation, row.names = FALSE)
  }

  if (!is.null(x$vcov_formula)) {
    cat("\nvcov_formula:", deparse1(x$vcov_formula), "\n")
  }

  if (!is.null(x$vcov_re)) {
    cat("\nVariance components (lmer_vcov_fit):\n")
    print(x$vcov_re)
    cat("residual_var:", x$residual_var, "\n")
  }

  invisible(x)
}
