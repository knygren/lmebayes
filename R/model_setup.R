#' Bayesian mixed model setup (single-factor \code{lmer} gate)
#'
#' Wrapper around \code{\link[lme4]{lmer}} for models with exactly one
#' grouping factor. Design matrices come from \code{formula} (including
#' cross-level RE moderation terms). Variance components use \code{vcov_formula}
#' (defaults to \code{\link{lmerb_default_vcov_formula}}): level-2 fixed only,
#' same \code{||} or \code{|} random structure as \code{formula}, without
#' cross-level fixed interactions (so RE moderation is not double-coded).
#'
#' @details
#' The example uses \code{big_word_club} from the Suggested package
#' \pkg{bayesrules} (see \code{?bayesrules::big_word_club}).
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
#'   \code{vcov_re}, \code{residual_var}, and \code{re_rank} (named logical
#'   vector: \code{TRUE} if \code{Z_j} is full column rank for that group).
#' @seealso \code{\link{extract_re_hyper_matrices}},
#'   \code{\link{lmerb_default_vcov_formula}},
#'   \code{\link{extract_lmer_variance_components}}
#' @examplesIf requireNamespace("bayesrules", quietly = TRUE)
#' @example inst/examples/Ex_model_setup_big_word_club.R
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
  cl <- match.call()
  design <- extract_re_hyper_matrices(formula = formula, data = data, ...)
  design$call    <- cl
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

  # Per-group rank check: is Z_j full column rank for each factor level?
  p_re  <- ncol(design$Z)
  g_chr <- as.character(design$groups)
  design$re_rank <- vapply(
    levels(design$groups),
    function(lev) {
      rows <- which(g_chr == lev)
      Z_j  <- design$Z[rows, , drop = FALSE]
      nrow(Z_j) >= p_re &&
        Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
    },
    logical(1L)
  )

  # Hyper-design rank check: for each RE coefficient, is the level-2 design
  # matrix X_hyper[[nm]] full column rank when restricted to the full-rank
  # groups?  Rank-deficient groups contribute a zero BLUP for the missing
  # slope and are excluded here so the check reflects only groups that
  # actually supply information about each RE.
  full_rank_levs <- names(design$re_rank)[design$re_rank]
  design$hyper_rank <- vapply(
    design$re_coef_names,
    function(nm) {
      Xh <- design$X_hyper[[nm]][full_rank_levs, , drop = FALSE]
      p  <- ncol(Xh)
      nrow(Xh) >= p && Matrix::rankMatrix(Xh, method = "qr")[1L] == p
    },
    logical(1L)
  )

  # Convenience summaries:
  #   hyper_deficient : named logical, TRUE = that RE's hyper-matrix is
  #                     rank-deficient (inverse of hyper_rank)
  #   rank_ok         : scalar TRUE only when every Z_j AND every hyper-matrix
  #                     is full-rank -- a quick go/no-go indicator
  design$hyper_deficient <- !design$hyper_rank

  # rank_ok reflects only the hyper-design matrices (level-2 estimability):
  # TRUE  = all X_hyper are full-rank after restricting to full-rank groups
  #         => the random-effects model can be estimated
  # FALSE = at least one X_hyper is rank-deficient => hyper parameters are
  #         not identified; Z_j rank deficiency is reported separately above
  design$rank_ok <- all(design$hyper_rank)

  design
}

#' @export
print.model_setup <- function(x, ...) {

  resp     <- deparse(x$formula[[2L]])
  re_names <- x$re_coef_names
  grp      <- x$group_name
  n_obs    <- length(x$y)
  n_lev    <- nlevels(x$groups)

  # ---- Call ------------------------------------------------------------------
  if (!is.null(x$call)) {
    cat("Call:\n  ", deparse1(x$call), "\n\n", sep = "")
  }

  # ---- Section 1: Measurement Model -----------------------------------------
  cat("--- Measurement Model ---\n")
  cat(sprintf("  %s ~ %s\n\n", resp, paste(re_names, collapse = " + ")))
  cat(sprintf("  Observations : %d\n", n_obs))
  cat(sprintf("  RE predictors: %d\n", length(re_names)))
  cat(sprintf("  Group        : %s  [%d levels]\n", grp, n_lev))
  if (!is.null(x$re_rank)) {
    n_full <- sum(x$re_rank)
    cat(sprintf("  Full-rank Z_j: %d of %d groups\n", n_full, n_lev))
    if (n_full < n_lev) {
      deficient <- names(x$re_rank)[!x$re_rank]
      shown     <- deficient[seq_len(min(10L, length(deficient)))]
      suffix    <- if (length(deficient) > 10L)
        sprintf(", ... (%d more)", length(deficient) - 10L) else ""
      cat(sprintf("    rank-deficient: %s%s\n",
                  paste(shown, collapse = ", "), suffix))
    }
  }
  cat("\n")

  # ---- Section 2: Random Effects Model --------------------------------------
  cat("--- Random Effects Model ---\n")

  w <- max(nchar(re_names))

  for (nm in re_names) {
    Xj    <- x$X_hyper[[nm]]
    other <- setdiff(colnames(Xj), "(Intercept)")

    hyper_rhs <- if (length(other) == 0L) "1" else paste(c("1", other), collapse = " + ")

    cat(sprintf("  %-*s ~ %s\n", w, nm, hyper_rhs))
  }
  cat("\n")

  # ---- Section 3: Hyper-design rank (full-rank groups only) -----------------
  if (!is.null(x$hyper_rank) && !is.null(x$re_rank)) {
    n_full_groups <- sum(x$re_rank)
    cat("--- Random Effects Model: Hyper-Design Rank ---\n")
    cat(sprintf("  (Restricted to %d full-rank %s)\n\n", n_full_groups, grp))
    deficient_nms <- character(0)
    for (nm in re_names) {
      Xh      <- x$X_hyper[[nm]]
      p_hyper <- ncol(Xh)
      is_fr   <- if (nm %in% names(x$hyper_rank)) x$hyper_rank[[nm]] else NA
      status  <- if (isTRUE(is_fr)) "full-rank" else if (isFALSE(is_fr)) "RANK-DEFICIENT" else "unknown"
      cat(sprintf("  %-*s  groups=%-3d  predictors=%-2d  %s\n",
                  w, nm, n_full_groups, p_hyper, status))
      if (isFALSE(is_fr)) deficient_nms <- c(deficient_nms, nm)
    }
    # Per-RE deficient flags
    cat("\n")
    flag_strs <- ifelse(x$hyper_deficient[re_names], "TRUE (deficient)", "FALSE")
    cat("  Rank-deficient flags:\n")
    for (nm in re_names) {
      cat(sprintf("    %-*s  %s\n", w, nm, flag_strs[nm]))
    }

    # Overall indicator
    ok_label <- if (isTRUE(x$rank_ok)) "TRUE  -- model rank looks OK" else
                  "FALSE -- rank issues detected (see above)"
    cat(sprintf("\n  rank_ok: %s\n", ok_label))

    if (length(deficient_nms) > 0L) {
      cat("\n")
      for (nm in deficient_nms) {
        cat(sprintf(
          "  NOTE: X_hyper for '%s' is rank-deficient after restricting to\n",
          nm))
        cat(sprintf(
          "  %d full-rank %s. Consider removing predictors or merging\n",
          n_full_groups, grp))
        cat("  factor levels.\n")
      }
    }
    cat("\n")
  }

  invisible(x)
}
