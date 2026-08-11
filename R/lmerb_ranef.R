#' Random effects and coefficient extractors for \code{lmerb} / \code{glmerb}
#'
#' @name lmerb_ranef
#' @keywords internal
NULL

## Block 1 / Block 2 containers were renamed coefficients -> groupef,
## ranef.mode -> groupef.mode, fixef.* -> popef.*.  Fits produced before the
## rename may still be in a user's workspace, so read either spelling.

#' @keywords internal
#' @noRd
.lmerb_groupef_draws <- function(object) {
  if (!is.null(object$groupef)) object$groupef else object$coefficients
}

#' @keywords internal
#' @noRd
.lmerb_groupef_mode <- function(object) {
  if (!is.null(object$groupef.mode)) object$groupef.mode else object$ranef.mode
}

#' @keywords internal
#' @noRd
.lmerb_popef_mode <- function(object) {
  if (!is.null(object$popef.mode)) object$popef.mode else object$fixef.mode
}

#' @keywords internal
#' @noRd
.lmerb_popef_means <- function(object) {
  if (!is.null(object$popef.means)) object$popef.means else object$fixef.means
}

#' @keywords internal
#' @noRd
.lmerb_popef_draws <- function(object) {
  if (!is.null(object$popef)) object$popef else object$fixef
}

## Per-group population expectation W_j %*% gamma, as a p_re x J matrix.
## Legacy fits stored this directly as fixef.mu; current fits do not, so it is
## rebuilt from the level-2 design and the population mode.
#' @keywords internal
#' @noRd
.lmerb_popef_mu <- function(object) {
  if (!is.null(object$fixef.mu)) return(as.matrix(object$fixef.mu))

  design <- object$model_setup
  gamma  <- .lmerb_popef_mode(object)
  if (is.null(design) || is.null(gamma)) return(NULL)

  re_names <- design$groupef.names
  grp_levs <- levels(design$group)
  mu <- matrix(
    NA_real_,
    nrow = length(re_names), ncol = length(grp_levs),
    dimnames = list(re_names, grp_levs)
  )
  for (k in re_names) {
    mu[k, ] <- as.vector(design$W[[k]] %*% as.numeric(gamma[[k]]))
  }
  mu
}

#' Extract random effects from Bayesian mixed model fits
#'
#' @description
#' Methods mirroring \code{\link[lme4]{ranef}} for \code{\link{lmerb}} and
#' \code{\link{glmerb}} fits.  Values summarize Block~1 random-effects draws
#' \eqn{b} stored in \code{fit$groupef} (posterior chain mean by default,
#' or the exact ICM \code{fit$groupef.mode}).
#'
#' @param object An \code{lmerb} or \code{glmerb} object.
#' @param condVar Logical; if \code{TRUE}, attach a \code{postVar} attribute
#'   (posterior variance from stored draws, diagonal by RE term).
#' @param drop Logical; if \code{TRUE}, convert single-column grouping
#'   \code{data.frame}s to named vectors (as in \pkg{lme4}).
#' @param type \code{"mean"} (default when draws are stored) or \code{"mode"}
#'   (ICM / \code{groupef.mode}).  When \code{simulate = FALSE}, \code{"mean"}
#'   falls back to \code{"mode"} with a warning.
#' @param \ldots Ignored.
#' @return An object of class \code{"ranef.lmerb"}: a named list with one
#'   \code{data.frame} per grouping factor (rows = levels, columns = RE terms).
#'   Attributes include \code{type}, \code{n_draws}, and optionally
#'   \code{postVar}.
#' @seealso \code{\link[lme4]{ranef}}, \code{\link{coef.lmerb}},
#'   \code{\link{fixef.lmerb}}, \code{\link{lmerb}}.
#' @importFrom lme4 fixef ranef VarCorr
#' @importFrom stats coef
#' @export
#' @method ranef lmerb
ranef.lmerb <- function(
    object,
    condVar = FALSE,
    drop = FALSE,
    type = c("mean", "mode"),
    ...
) {
  type <- match.arg(type)
  .lmerb_ranef_extract(object, condVar, drop, type)
}

#' @rdname ranef.lmerb
#' @export
#' @method ranef glmerb
ranef.glmerb <- function(
    object,
    condVar = FALSE,
    drop = FALSE,
    type = c("mean", "mode"),
    ...
) {
  type <- match.arg(type)
  .lmerb_ranef_extract(object, condVar, drop, type)
}

#' @keywords internal
#' @noRd
.lmerb_ranef_extract <- function(object, condVar, drop, type) {
  if (!inherits(object, c("lmerb", "glmerb"))) {
    stop("object must be an lmerb or glmerb fit.", call. = FALSE)
  }
  re_names <- object$model_setup$groupef.names
  grp_col  <- object$model_setup$group_name
  if (!length(re_names)) {
    stop("object has no random-effects components.", call. = FALSE)
  }

  n_draws <- NULL
  draws <- .lmerb_groupef_draws(object)
  if (type == "mean") {
    if (is.null(draws)) {
      warning(
        "No stored draws (simulate = FALSE); using groupef.mode instead.",
        call. = FALSE
      )
      type <- "mode"
    } else {
      mat <- .bayes_ranef_chain_means(object)
      n_draws <- length(unique(draws[["draw"]]))
    }
  }
  if (type == "mode") {
    mat <- .lmerb_groupef_mode(object)
    if (is.null(mat)) {
      stop("object$groupef.mode is NULL.", call. = FALSE)
    }
  }

  df <- as.data.frame(mat, check.names = FALSE)
  if (isTRUE(drop) && ncol(df) == 1L) {
    df <- setNames(df[[1L]], rownames(df))
  }
  out <- stats::setNames(list(df), grp_col)

  postVar <- NULL
  if (isTRUE(condVar)) {
    if (is.null(draws)) {
      warning(
        "condVar requested but no stored draws; postVar omitted.",
        call. = FALSE
      )
    } else {
      postVar <- stats::setNames(
        list(.lmerb_ranef_postVar(object)),
        grp_col
      )
    }
  }

  structure(
    out,
    class    = "ranef.lmerb",
    postVar  = postVar,
    type     = type,
    n_draws  = n_draws,
    fit      = object
  )
}

#' @keywords internal
#' @noRd
.lmerb_ranef_postVar <- function(object) {
  re_names <- object$model_setup$groupef.names
  grp_col  <- object$model_setup$group_name
  coef_df  <- .lmerb_groupef_draws(object)
  grp_levs <- levels(factor(coef_df[[grp_col]]))
  p_re     <- length(re_names)
  J        <- length(grp_levs)
  pv <- array(
    0,
    dim = c(p_re, p_re, J),
    dimnames = list(re_names, re_names, grp_levs)
  )
  for (lev in grp_levs) {
    idx <- coef_df[[grp_col]] == lev
    for (k in seq_len(p_re)) {
      v <- stats::var(coef_df[idx, re_names[k]])
      pv[k, k, lev] <- if (is.finite(v)) v else NA_real_
    }
  }
  pv
}

#' Coerce \code{ranef.lmerb} objects to long \code{data.frame} format
#'
#' @param x An object of class \code{"ranef.lmerb"}.
#' @param \ldots Ignored.
#' @return A \code{data.frame} with columns \code{grpvar}, \code{term},
#'   \code{grp}, \code{condval}, and \code{condsd} (when \code{postVar} is
#'   available on \code{x}).
#' @seealso \code{\link[lme4]{ranef}}.
#' @export
#' @method as.data.frame ranef.lmerb
as.data.frame.ranef.lmerb <- function(x, ...) {
  pv <- attr(x, "postVar")
  pieces <- lapply(names(x), function(gv) {
    elem <- x[[gv]]
    if (is.data.frame(elem)) {
      levs <- rownames(elem)
      terms <- names(elem)
      do.call(rbind, lapply(terms, function(trm) {
        condsd <- if (!is.null(pv) && !is.null(pv[[gv]])) {
          sqrt(pv[[gv]][trm, trm, levs])
        } else {
          rep(NA_real_, length(levs))
        }
        data.frame(
          grpvar  = gv,
          term    = trm,
          grp     = levs,
          condval = unname(elem[[trm]]),
          condsd  = unname(condsd),
          row.names = NULL,
          check.names = FALSE
        )
      }))
    } else {
      data.frame(
        grpvar  = gv,
        term    = names(elem),
        grp     = names(elem),
        condval = unname(elem),
        condsd  = NA_real_,
        row.names = NULL,
        check.names = FALSE
      )
    }
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  class(out) <- c("ranef.lmerb", "data.frame")
  out
}

#' Print random effects in \pkg{lme4} layout
#'
#' @param x An object of class \code{"ranef.lmerb"}.
#' @param \ldots Ignored.
#' @return \code{x}, invisibly.
#' @export
#' @method print ranef.lmerb
print.ranef.lmerb <- function(x, ...) {
  ## as.data.frame() keeps the class, so the long format lands here too.
  if (is.data.frame(x)) {
    y <- x
    class(y) <- "data.frame"
    print(y, ...)
    return(invisible(x))
  }
  for (gv in names(x)) {
    cat("$", gv, "\n", sep = "")
    print(x[[gv]])
  }
  invisible(x)
}

#' Summarize \code{ranef.lmerb} objects
#'
#' @param object An object of class \code{"ranef.lmerb"}.
#' @param groups Optional character vector of grouping levels to highlight.
#' @param digits Number of digits for printing.
#' @param \ldots Ignored.
#' @return A \code{data.frame} (invisibly) with long-format random effects.
#' @export
#' @method summary ranef.lmerb
summary.ranef.lmerb <- function(object, groups = NULL, digits = 4L, ...) {
  tab <- as.data.frame(object)
  if (!is.null(groups)) {
    groups <- as.character(groups)
    tab <- tab[tab$grp %in% groups, , drop = FALSE]
  }
  num_cols <- vapply(tab, is.numeric, logical(1L))
  tab[num_cols] <- lapply(tab[num_cols], round, digits = digits)
  print(tab, row.names = FALSE)
  invisible(tab)
}

#' Extract combined coefficients (fixed + random) by group
#'
#' @description
#' Mirrors \code{\link[lme4]{coef.merMod}}: Block~1 stored coefficients
#' (\code{ranef.lmerb}) are already on the per-group coefficient scale used
#' internally by the sampler (not lme4-style deviations).
#'
#' @inheritParams ranef.lmerb
#' @param \ldots Ignored.
#' @return A named list (one element per grouping factor) of \code{data.frame}s
#'   with the same layout as \code{coef(merMod)}.
#' @seealso \code{\link[lme4]{coef.merMod}}, \code{\link{ranef.lmerb}}.
#' @export
#' @method coef lmerb
coef.lmerb <- function(object, type = c("mean", "mode"), ...) {
  type <- match.arg(type)
  .lmerb_coef_extract(object, type)
}

#' @rdname coef.lmerb
#' @export
#' @method coef glmerb
coef.glmerb <- function(object, type = c("mean", "mode"), ...) {
  type <- match.arg(type)
  .lmerb_coef_extract(object, type)
}

#' @keywords internal
#' @noRd
.lmerb_coef_extract <- function(object, type) {
  mer <- if (inherits(object, "lmerb")) object$lmer else object$glmer
  if (is.null(mer)) {
    stop("object has no reference lmer/glmer fit.", call. = FALSE)
  }
  re <- lme4::ranef(object, condVar = FALSE, drop = FALSE, type = type)
  grp_col <- names(re)[1L]
  out_df <- re[[1L]]
  if (!is.data.frame(out_df)) {
    out_df <- data.frame(`(Intercept)` = out_df, check.names = FALSE)
    rownames(out_df) <- names(re[[1L]])
  }
  re_names <- colnames(out_df)
  ## Block~1 coefficients stored in fit$groupef / groupef.mode are
  ## already full per-group random coefficients (same scale as coef(mer)),
  ## not lme4-style deviations; do not add fixef(mer) again.
  structure(stats::setNames(list(out_df), grp_col), class = "coef.lmerb")
}

#' Extract fixed effects from \code{lmerb} / \code{glmerb} fits
#'
#' @description
#' By default returns \code{\link[lme4]{fixef}} from the reference
#' \code{lmer}/\code{glmer} fit (population-level regression coefficients).
#' Use \code{type = "hyper"} or \code{"hyper.mean"} for Block~2 hyperparameter
#' summaries stored on the Bayesian fit.
#'
#' @param object An \code{lmerb} or \code{glmerb} object.
#' @param type \code{"mer"} (default), \code{"hyper"} (\code{popef.mode}), or
#'   \code{"hyper.mean"} (\code{popef.means}; requires stored draws).
#' @param \ldots Ignored.
#' @return A named numeric vector.
#' @seealso \code{\link[lme4]{fixef}}, \code{\link{coef.lmerb}}.
#' @export
#' @method fixef lmerb
fixef.lmerb <- function(object, type = c("mer", "hyper", "hyper.mean"), ...) {
  type <- match.arg(type)
  .lmerb_fixef_extract(object, type)
}

#' @rdname fixef.lmerb
#' @export
#' @method fixef glmerb
fixef.glmerb <- function(object, type = c("mer", "hyper", "hyper.mean"), ...) {
  type <- match.arg(type)
  .lmerb_fixef_extract(object, type)
}

#' @keywords internal
#' @noRd
.lmerb_fixef_extract <- function(object, type) {
  switch(
    type,
    mer = {
      mer <- if (inherits(object, "lmerb")) object$lmer else object$glmer
      if (is.null(mer)) {
        stop("object has no reference lmer/glmer fit.", call. = FALSE)
      }
      lme4::fixef(mer)
    },
    hyper = {
      out <- .lmerb_popef_mode(object)
      if (is.null(out)) {
        stop("object$popef.mode is NULL.", call. = FALSE)
      }
      unlist(out, use.names = TRUE)
    },
    hyper.mean = {
      out <- .lmerb_popef_means(object)
      if (is.null(out)) {
        stop(
          "object$popef.means is NULL (run with simulate = TRUE).",
          call. = FALSE
        )
      }
      unlist(out, use.names = TRUE)
    }
  )
}
