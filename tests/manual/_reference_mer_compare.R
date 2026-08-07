#' Print lme4 and glmmTMB reference fits for manual smoke-test diagnostics.
#'
#' Intended for scripts under \code{tests/manual/} and \code{data-raw/} only.
#' Mirrors \code{Prior_Setup_GLMM()} routing: \code{dispformula = ~1}
#' shows \code{lme4::lmer} only (the calibration reference); per-group
#' \code{dispformula} adds the matching \code{glmmTMB} fit used as
#' \code{fit_ref}.
#'
#' @param form Mixed-model formula.
#' @param dat Data frame.
#' @param group_name Grouping variable name (default \code{"school_id"}).
#' @param dispformula One-sided formula, \code{~1} (pooled) or
#'   \code{~group_name} (per-group residual variance). Must match the
#'   \code{dispformula} passed to \code{Prior_Setup_GLMM()} /
#'   \code{lmerb()} in the calling script.
#' @return Invisibly a list with \code{lmer} and optional \code{glmmTMB}.
#' @keywords internal
.print_reference_mer_compare <- function(
    form,
    dat,
    group_name = "school_id",
    dispformula = ~1
) {
  disp_kind <- .reference_mer_compare_dispformula_kind(dispformula, group_name)

  cat("\n=== lme4::lmer reference (REML) ===\n\n")
  fit_lmer <- lme4::lmer(form, data = dat, REML = TRUE)
  print(summary(fit_lmer))
  cat("\nfixef(lmer):\n")
  print(lme4::fixef(fit_lmer))
  cat("\nVarCorr(lmer):\n")
  print(lme4::VarCorr(fit_lmer), comp = "Std.Dev.")
  cat("\nranef(lmer):\n")
  print(lme4::ranef(fit_lmer))
  cat("\ncoef(lmer):\n")
  print(stats::coef(fit_lmer))
  cat("\nResidual SD (lmer): ", stats::sigma(fit_lmer), "\n", sep = "")

  out <- list(lmer = fit_lmer)

  if (identical(disp_kind, "pooled")) {
    cat(
      "\n(dispformula = ~1: Prior_Setup_GLMM() uses lme4 only; ",
      "skipping glmmTMB reference fit.)\n",
      sep = ""
    )
    return(invisible(out))
  }

  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    cat("\n(glmmTMB not installed — skipping glmmTMB reference fit.)\n")
    return(invisible(out))
  }

  print_glmmtmb <- function(fit, label) {
    cat("\n=== glmmTMB reference: ", label, " ===\n\n", sep = "")
    conv <- fit$fit$convergence
    pd   <- fit$sdr$pdHess
    cat(
      "convergence = ", conv,
      if (!is.null(pd)) paste0(", pdHess = ", pd) else "",
      "\n",
      sep = ""
    )
    if (!identical(as.numeric(conv), 0)) {
      cat("  ** optimizer did not report convergence (code != 0) **\n")
    }
    if (!is.null(pd) && !isTRUE(pd)) {
      cat("  ** Hessian is not positive-definite **\n")
    }
    print(summary(fit))
    cat("\nfixef(glmmTMB, cond):\n")
    print(glmmTMB::fixef(fit)$cond)
    cat("\nVarCorr(glmmTMB, cond):\n")
    vc <- glmmTMB::VarCorr(fit)
    if (is.list(vc) && !is.null(vc$cond)) {
      vc <- vc$cond
    }
    print(vc, comp = "Std.Dev.")
    cat("\nranef(glmmTMB, cond):\n")
    print(glmmTMB::ranef(fit)$cond)
    disp <- tryCatch(
      as.numeric(stats::predict(fit, type = "disp")),
      error = function(e) NULL
    )
    if (!is.null(disp)) {
      cat("\nResidual dispersion (predict type = 'disp'):\n")
      if (length(disp) == 1L) {
        cat("  pooled sigma^2 = ", disp, ", sigma = ", sqrt(disp), "\n", sep = "")
      } else if (group_name %in% names(dat)) {
        by_grp <- stats::aggregate(
          disp,
          by = list(group = dat[[group_name]]),
          FUN = mean
        )
        names(by_grp)[2L] <- "sigma2"
        by_grp$sigma <- sqrt(by_grp$sigma2)
        print(by_grp[order(by_grp$group), ], row.names = FALSE)
      } else {
        cat(
          "  ", length(disp), " values; mean sigma^2 = ", mean(disp),
          ", sigma = ", sqrt(mean(disp)), "\n",
          sep = ""
        )
      }
    }
  }

  fit_tmb <- glmmTMB::glmmTMB(
    formula     = form,
    data        = dat,
    family      = gaussian(),
    dispformula = dispformula,
    REML        = TRUE
  )
  print_glmmtmb(
    fit_tmb,
    paste0("dispformula = ", deparse(dispformula), " (Prior_Setup fit_ref)")
  )
  out$glmmTMB <- fit_tmb

  invisible(out)
}

## Classify dispformula as "pooled" (~1) or "group" (~<group_name>).
.reference_mer_compare_dispformula_kind <- function(dispformula, group_name) {
  if (!inherits(dispformula, "formula") || length(dispformula) != 2L) {
    stop(
      "'dispformula' must be a one-sided formula, either ~1 or ~",
      group_name,
      ".",
      call. = FALSE
    )
  }
  vars <- all.vars(dispformula)
  if (length(vars) == 0L) {
    return("pooled")
  }
  if (length(vars) == 1L && identical(vars, group_name)) {
    return("group")
  }
  stop(
    "'dispformula' must be ~1 (pooled) or ~", group_name,
    "; got ", deparse(dispformula), ".",
    call. = FALSE
  )
}
