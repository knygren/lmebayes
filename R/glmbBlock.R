#' @title Fitting Blocked Bayesian Generalized Linear Models
#'
#' @description
#' Fits one \code{\link[glmbayes]{glmb}} per observation block. Same row partition as
#' \code{\link{lmbBlock}}, but supports GLM \code{\link{family}} objects.
#' Counterpart to \code{\link{lmbBlock}}; see \code{\link{summary.bglmb}} for
#' the summary method and \code{\link{rNormalGLMBlock}} for Gibbs sampling.
#'
#' @inheritParams lmbBlock
#' @name glmbBlock
#' @family modelfuns
NULL

#' @describeIn glmbBlock \code{\link[glmbayes]{glmb}} fit per row block.
#' @param pfamily Recycled to all blocks, or use \code{pfamily_list} of length \code{k}.
#' @param pfamily_list Optional list of \code{pfamily} objects, one per block.
#' @inheritParams glmbayes::glmb
#' @inheritParams lmbBlock
#' @param model,x,y For \code{glmbBlock}, logicals passed to each block
#'   \code{\link[glmbayes]{glmb}} fit (see \code{\link[stats]{glm}}). For
#'   \code{print}, \code{x} is the \code{"bglmb"} object.
#' @return A named list of class \code{"bglmb"} (list of \code{"glmb"} fits from **glmbayes**).
#' @param digits Number of significant digits to use when printing.
#' @aliases print.bglmb
#' @export
glmbBlock <- function(
    formula,
    block,
    family = gaussian(),
    pfamily = NULL,
    pfamily_list = NULL,
    n = 1000,
    data,
    subset,
    weights,
    na.action,
    offset,
    contrasts = NULL,
    model = TRUE,
    x = FALSE,
    y = TRUE,
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    ...
) {
  mc <- match.call(expand.dots = FALSE)
  call <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }

  if (is.null(pfamily_list) && is.null(pfamily)) {
    stop("Provide 'pfamily' or 'pfamily_list'.", call. = FALSE)
  }
  if (!is.null(pfamily_list) && !is.null(pfamily)) {
    stop("Specify only one of 'pfamily' or 'pfamily_list'.", call. = FALSE)
  }

  meta <- .blmb_formula_block_meta(
    formula = formula,
    block = block,
    data = data,
    subset = if (!missing(subset)) subset else NULL,
    weights = if (!missing(weights)) weights else NULL,
    na.action = if (!missing(na.action)) na.action else NULL,
    offset = if (!missing(offset)) offset else NULL,
    contrasts = if (!missing(contrasts)) contrasts else NULL
  )
  k <- meta$block_info$k
  p <- meta$p

  if (is.null(pfamily_list)) {
    pfamily_list <- rep(list(pfamily), k)
  }
  pfamily_lists <- .mrglmb_normalize_pfamily_lists(
    pfamily_list,
    k,
    p,
    .validate_pfamily_for_rlmb
  )

  glmb_args <- list(
    n = n,
    family = family,
    data = data,
    contrasts = contrasts,
    model = model,
    x = x,
    y = y,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose
  )
  if (!missing(subset)) glmb_args$subset <- subset
  if (!missing(weights)) glmb_args$weights <- weights
  if (!missing(na.action)) glmb_args$na.action <- na.action
  if (!missing(offset)) glmb_args$offset <- offset
  if (length(list(...))) {
    glmb_args <- c(glmb_args, list(...))
  }

  block_results <- vector("list", k)
  for (b in seq_len(k)) {
    rows_b <- .blmb_rows_to_data_subset(
      meta$block_info$rows[[b]], meta$mf, data
    )
    fit_b <- do.call(
      glmb,
      c(
        list(
          formula = formula,
          pfamily = pfamily_lists[[b]],
          subset = rows_b
        ),
        glmb_args
      )
    )
    fit_b$call <- .blmb_glmb_display_call(mc, formula, rows_b)
    block_results[[b]] <- fit_b
  }

  .bglmb_assemble(
    block_results = block_results,
    block_ids = meta$block_info$ids,
    call = call,
    formula = formula,
    block = block,
    block_info = meta$block_info,
    family = family,
    p = p,
    pred_names = meta$pred_names,
    pfamily_lists = pfamily_lists
  )
}

#' @keywords internal
.blmb_glmb_display_call <- function(mc_block, formula, rows_b) {
  cl <- call("glmb", formula = formula, subset = rows_b)
  pass <- c(
    "n", "family", "data", "weights", "na.action", "offset",
    "use_parallel", "use_opencl", "verbose", "Gridtype", "n_envopt"
  )
  for (nm in pass) {
    if (nm %in% names(mc_block)) {
      cl[[nm]] <- mc_block[[nm]]
    }
  }
  cl
}

#' @keywords internal
.bglmb_assemble <- function(
    block_results,
    block_ids,
    call,
    formula,
    block,
    block_info,
    family,
    p,
    pred_names,
    pfamily_lists = NULL
) {
  outlist <- setNames(block_results, block_ids)
  attr(outlist, "call") <- call
  attr(outlist, "formula") <- formula
  attr(outlist, "block") <- block
  attr(outlist, "block_info") <- block_info
  attr(outlist, "family") <- family
  attr(outlist, "k") <- block_info$k
  attr(outlist, "p") <- p
  attr(outlist, "pred_names") <- pred_names
  if (!is.null(pfamily_lists)) {
    attr(outlist, "pfamily_lists") <- pfamily_lists
  }
  class(outlist) <- "bglmb"
  outlist
}

#' @rdname glmbBlock
#' @method print bglmb
#' @export
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
