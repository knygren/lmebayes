#' Summarize per-group measurement dispersion (sigma^2)
#'
#' @description
#' Per-group observation-level \eqn{\sigma^2} summaries when
#' \code{dispersion_ranef} is a named list of \code{\link{dGamma}()} priors
#' (\code{dispersion_mode == "gamma_list"}).  Layout mirrors the Block~2
#' \eqn{\tau^2_k} tables in \code{\link{summary.lmerb}}: a prior reference
#' table (with a common \code{lmer}/\code{glmer} residual column), an
#' optional posterior overview from stored draws (including per-group
#' \code{Cand/draw}, envelope candidates per inner sweep), and per-group
#' distribution percentiles of the \eqn{\sigma^2} draws (1\%/2.5\%/5\%/
#' median/95\%/97.5\%/99\%) -- comparing these against \code{disp_lower}/
#' \code{disp_upper} in the prior table shows whether the simulation tends
#' to pile up against the truncation bounds. The prior table also reports
#' each group's implied effective prior sample size (\code{n_prior},
#' back-computed from \code{shape_ING} under the \code{Prior_Setup()}
#' default \code{k = 1} calibration) and observed data sample size
#' (\code{n_data}); the overview table reports their sum (\code{n_total}).
#'
#' @param object An \code{lmerb} or \code{glmerb} fit.
#' @param type \code{"both"} (default), \code{"prior"}, or \code{"overview"}.
#' @param digits Number of significant digits for printing.
#' @param \ldots Ignored.
#' @return An object of class \code{"summary.sigma2.lmerb"} with components
#'   \code{prior}, \code{overview}, \code{percentiles} (all when applicable),
#'   \code{group_name}, \code{n_groups}, \code{simulated}, and
#'   \code{mer_label}.  When \code{window_diagnostics} were attached by
#'   \code{\link[glmbayesCore]{dGamma_list}()}, the prior table includes
#'   \code{blup_infl}, \code{R_lo}, \code{R_hi}, and \code{asymmetric_window}.
#' @seealso \code{\link{summary.lmerb}}, \code{\link{lmerb}}.
#' @export
summary_sigma2 <- function(object, ...) {
  UseMethod("summary_sigma2")
}

#' @rdname summary_sigma2
#' @export
#' @method summary_sigma2 lmerb
summary_sigma2.lmerb <- function(
    object,
    type = c("both", "prior", "overview"),
    digits = max(3L, getOption("digits") - 3L),
    ...
) {
  type <- match.arg(type)
  .lmerb_summary_sigma2_build(object, type, digits)
}

#' @rdname summary_sigma2
#' @export
#' @method summary_sigma2 glmerb
summary_sigma2.glmerb <- function(
    object,
    type = c("both", "prior", "overview"),
    digits = max(3L, getOption("digits") - 3L),
    ...
) {
  type <- match.arg(type)
  .lmerb_summary_sigma2_build(object, type, digits)
}

#' @keywords internal
#' @noRd
.lmerb_summary_sigma2_build <- function(object, type, digits) {
  if (!inherits(object, c("lmerb", "glmerb"))) {
    stop("object must be an lmerb or glmerb fit.", call. = FALSE)
  }
  if (!identical(object$prior$dispersion_mode, "gamma_list")) {
    stop(
      "summary_sigma2() applies only when dispersion_ranef is a per-group ",
      "list of dGamma() priors (dispersion_mode = \"gamma_list\").",
      call. = FALSE
    )
  }
  re_names  <- object$model_setup$re_coef_names
  simulated <- !is.null(object$coefficients)
  n_draws   <- if (simulated && length(re_names)) {
    nrow(object$fixef[[re_names[1L]]])
  } else {
    NULL
  }
  mer_label <- if (inherits(object, "glmerb")) "glmer" else "lmer"

  prior <- .lmerb_sigma2_gamma_list_prior_overview(object)
  if (is.null(prior) || nrow(prior) == 0L) {
    stop(
      "No per-group sigma^2 prior information found on object$prior.",
      call. = FALSE
    )
  }
  overview <- .lmerb_sigma2_gamma_list_posterior_overview(
    object, simulated = simulated, n_draws = n_draws
  )
  percentiles <- .lmerb_sigma2_gamma_list_percentiles_overview(
    object, simulated = simulated
  )

  structure(
    list(
      call        = object$call,
      formula     = object$formula,
      type        = type,
      digits      = digits,
      group_name  = object$model_setup$group_name,
      n_groups    = nlevels(object$model_setup$groups),
      simulated   = simulated,
      n           = n_draws,
      mer_label   = mer_label,
      prior       = if (type %in% c("both", "prior")) prior else NULL,
      overview    = if (type %in% c("both", "overview")) overview else NULL,
      percentiles = if (type %in% c("both", "overview")) percentiles else NULL
    ),
    class = "summary.sigma2.lmerb"
  )
}

#' @rdname summary_sigma2
#' @param x An object of class \code{"summary.sigma2.lmerb"}.
#' @export
#' @method print summary.sigma2.lmerb
print.summary.sigma2.lmerb <- function(x, digits = NULL, ...) {
  if (is.null(digits)) {
    digits <- x$digits
  }
  if (is.null(digits)) {
    digits <- max(3L, getOption("digits") - 3L)
  }
  mer_label <- if (!is.null(x$mer_label)) x$mer_label else "lmer"

  cat("Measurement dispersion (sigma^2, per group)\n")
  if (!is.null(x$formula)) {
    cat("Formula:", deparse1(x$formula), "\n")
  }
  cat(sprintf(
    "Groups: %s (%d level(s))",
    x$group_name, x$n_groups
  ))
  if (isTRUE(x$simulated)) {
    cat(sprintf("  [%d draws]\n\n", x$n))
  } else {
    cat("  [ICM only; no stored draws]\n\n")
  }

  if (!is.null(x$prior) && nrow(x$prior) > 0L) {
    cat(sprintf("Prior and %s reference:\n\n", mer_label))
    .lmerb_print_summary_table(x$prior, digits = digits)
  }

  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview:\n")
    stats::printCoefmat(
      round(x$overview, digits = digits),
      digits = digits,
      quote = FALSE
    )
    if (!isTRUE(x$simulated)) {
      cat("\n  (Run with simulate = TRUE for MCMC means, SDs, and tail probabilities.)\n")
    }
  }

  if (!is.null(x$percentiles) && nrow(x$percentiles) > 0L) {
    cat("\nDistribution percentiles (sigma^2):\n\n")
    print(round(x$percentiles, digits = digits), right = TRUE)
  }

  invisible(x)
}

#' @rdname summary_sigma2
#' @param x An object of class \code{"summary.sigma2.lmerb"}.
#' @param row.names Ignored.
#' @export
#' @method as.data.frame summary.sigma2.lmerb
as.data.frame.summary.sigma2.lmerb <- function(x, row.names = NULL, ...) {
  if (is.null(x$prior)) {
    out <- x$overview
  } else if (is.null(x$overview)) {
    out <- x$prior
  } else {
    pri <- x$prior
    ovw <- x$overview
    ovw_cols <- setdiff(colnames(ovw), intersect(colnames(pri), colnames(ovw)))
    out <- cbind(pri, ovw[, ovw_cols, drop = FALSE])
  }
  if (!is.null(out)) {
    out[[x$group_name]] <- rownames(out)
    out <- out[, c(x$group_name, setdiff(names(out), x$group_name)), drop = FALSE]
    rownames(out) <- NULL
  }
  out
}
