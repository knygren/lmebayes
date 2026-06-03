#' Summarize and print bglmb fits
#'
#' @description
#' Methods for row-block \code{\link{glmb}} fits (class \code{"bglmb"}).
#' \code{summary.bglmb} applies \code{\link{summary.glmb}} to each block;
#' printing follows \code{\link[glmbayes]{summary.mlmb}} / \code{\link{summary.blmb}} with
#' per-block sections.
#'
#' @param object An object of class \code{"bglmb"} from \code{\link{block_glmb}}.
#' @param x An object of class \code{"summary.bglmb"}.
#' @param digits Number of significant digits for printing.
#' @param \ldots Passed to \code{\link{summary.glmb}} or print methods.
#' @return \code{summary.bglmb} returns a named list of \code{"summary.glmb"}
#'   objects with class \code{"summary.bglmb"}.
#' @seealso \code{\link{block_glmb}}, \code{\link{block_lmb}},
#'   \code{\link[glmbayes]{summary.mlmb}}, \code{\link{summary.blmb}},
#'   \code{\link{summary.glmb}}
#' @name summary.bglmb
#' @aliases summary.bglmb print.bglmb print.summary.bglmb
NULL

#' @rdname summary.bglmb
#' @export
#' @method summary bglmb
summary.bglmb <- function(object, ...) {
  res <- lapply(object, function(fit) {
    s <- summary(fit, ...)
    s$call <- fit$call
    s
  })
  names(res) <- names(object)
  attr(res, "bglmb_call") <- attr(object, "call")
  attr(res, "family") <- attr(object, "family")
  attr(res, "coef_means") <- .blmb_coef_means_matrix(object)
  attr(res, "dic_table") <- .blmb_dic_table(object)
  class(res) <- "summary.bglmb"
  res
}

#' @rdname summary.bglmb
#' @export
#' @method print bglmb
print.bglmb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cl <- attr(x, "call")
  if (is.null(cl) && length(x) >= 1L) {
    cl <- x[[1L]]$call
  }
  cat("\nCall:\n")
  if (!is.null(cl)) {
    if (is.call(cl)) {
      cat(paste(deparse(cl, width.cutoff = 500L), collapse = "\n"), "\n")
    } else {
      print(cl)
    }
  } else {
    cat("  (not recorded)\n")
  }

  fam <- attr(x, "family")
  if (!is.null(fam) && !is.null(fam$family)) {
    cat("Family:", fam$family)
    if (!is.null(fam$link) && nzchar(fam$link)) {
      cat(" (link = ", fam$link, ")", sep = "")
    }
    cat("\n")
  }

  cm <- .blmb_coef_means_matrix(x)
  if (!is.null(cm) && length(cm)) {
    cat("\nPosterior mean coefficients (rows = blocks):\n")
    print.default(
      format(cm, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  } else {
    cat("\nNo coefficients\n")
  }

  dic_tab <- .blmb_dic_table(x)
  if (!is.null(dic_tab) && nrow(dic_tab) >= 1L) {
    rownames(dic_tab) <- names(x)
    cat("\nBayesian fit (per block; independent BY model):\n")
    print.default(
      format(dic_tab, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
    cat(
      "Sum DIC:",
      format(sum(dic_tab[, "DIC"]), digits = digits),
      "  Sum pD:",
      format(sum(dic_tab[, "pD"]), digits = digits),
      "\n",
      sep = ""
    )
  }

  cat("\n")
  invisible(x)
}

#' @rdname summary.bglmb
#' @export
#' @method print summary.bglmb
print.summary.bglmb <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...
) {
  cl <- attr(x, "bglmb_call")
  if (!is.null(cl)) {
    cat("\nCall:\n")
    if (is.call(cl)) {
      cat(paste(deparse(cl, width.cutoff = 500L), collapse = "\n"), "\n")
    } else {
      print(cl)
    }
  }

  fam <- attr(x, "family")
  if (!is.null(fam) && !is.null(fam$family)) {
    cat("Family:", fam$family)
    if (!is.null(fam$link) && nzchar(fam$link)) {
      cat(" (link = ", fam$link, ")", sep = "")
    }
    cat("\n")
  }

  cm <- attr(x, "coef_means")
  if (!is.null(cm) && length(cm)) {
    cat("\nPosterior mean coefficients (rows = blocks):\n")
    print.default(
      format(cm, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  }

  dic_tab <- attr(x, "dic_table")
  if (!is.null(dic_tab) && nrow(dic_tab) >= 1L) {
    cat("\nBayesian fit (per block):\n")
    print.default(
      format(dic_tab, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
    cat(
      "Sum DIC:",
      format(sum(dic_tab[, "DIC"]), digits = digits),
      "  Sum pD:",
      format(sum(dic_tab[, "pD"]), digits = digits),
      "\n",
      sep = ""
    )
  }

  for (nm in names(x)) {
    cat("\nBlock", nm, ":\n")
    print(x[[nm]], digits = digits, ...)
  }
  invisible(x)
}
