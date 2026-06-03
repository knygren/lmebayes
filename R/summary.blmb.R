#' Summarize and print blmb fits
#'
#' @description
#' Methods for row-block \code{\link[glmbayes]{lmb}} fits (class \code{"blmb"}).
#'
#' @param object An object of class \code{"blmb"} from \code{\link{block_lmb}}.
#' @param x An object of class \code{"summary.blmb"}.
#' @param digits Number of significant digits for printing.
#' @param \ldots Passed to \code{\link{summary.glmb}} or print methods.
#' @name summary.blmb
#' @aliases summary.blmb print.blmb print.summary.blmb
NULL

#' @rdname summary.blmb
#' @export
#' @method summary blmb
summary.blmb <- function(object, ...) {
  res <- lapply(object, function(fit) {
    s <- summary(fit, ...)
    s$call <- fit$call
    s
  })
  names(res) <- names(object)
  attr(res, "blmb_call") <- attr(object, "call")
  attr(res, "coef_means") <- .blmb_coef_means_matrix(object)
  attr(res, "dic_table") <- .blmb_dic_table(object)
  class(res) <- "summary.blmb"
  res
}

#' @rdname summary.blmb
#' @export
#' @method print blmb
print.blmb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
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

#' @rdname summary.blmb
#' @export
#' @method print summary.blmb
print.summary.blmb <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...
) {
  cl <- attr(x, "blmb_call")
  if (!is.null(cl)) {
    cat("\nCall:\n")
    if (is.call(cl)) {
      cat(paste(deparse(cl, width.cutoff = 500L), collapse = "\n"), "\n")
    } else {
      print(cl)
    }
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
