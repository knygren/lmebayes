#' Prior setup for row-block \code{\link[glmbayes]{lmb}} / \code{\link{glmbBlock}}
#'
#' Runs \code{\link[glmbayesCore]{Prior_Setup}} on each block subset of the data.
#'
#' @param formula A \code{\link{formula}} with a single response.
#' @param block Block partition: \code{factor} or vector of length \code{nrow(data)}
#'   (after \code{model.frame}), a column name in \code{data}, \code{l2_blocks}
#'   counts, or a list of row index vectors (see \code{\link[glmbayesCore]{normalize_block}}).
#' @inheritParams glmbayesCore::Prior_Setup
#' @return A named list of class \code{"block_PriorSetup"}. Each element is a
#'   \code{\link[glmbayesCore]{Prior_Setup}} result for one block.
#' @family prior
#' @seealso \code{\link{lmbBlock}}, \code{\link[glmbayesCore]{multi_prior_setup}},
#'   \code{\link{Prior_Setup_lmebayes}}, \code{\link[glmbayesCore]{normalize_block}}
#' @export
Prior_SetupBlock <- function(
    formula,
    block,
    family = gaussian(),
    data = NULL,
    weights = NULL,
    subset = NULL,
    na.action = na.fail,
    offset = NULL,
    contrasts = NULL,
    pwt = NULL,
    pwt_default_low = 0.01,
    pwt_default_high = 0.05,
    n_prior = NULL,
    sd = NULL,
    dispersion = NULL,
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model"),
    mu = NULL,
    k = 1,
    ...
) {
  call <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  fam_ok <- family$family %in% c("gaussian", "poisson")
  if (is.null(family$family) || !fam_ok) {
    stop(
      "Prior_SetupBlock() supports family = gaussian() or poisson() only.",
      call. = FALSE
    )
  }
  if (missing(data)) {
    data <- environment(formula)
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

  ps_args <- list(
    family = family,
    data = data,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    pwt = pwt,
    pwt_default_low = pwt_default_low,
    pwt_default_high = pwt_default_high,
    n_prior = n_prior,
    sd = sd,
    dispersion = dispersion,
    intercept_source = intercept_source,
    effects_source = effects_source,
    mu = mu,
    k = k
  )

  setups <- vector("list", meta$block_info$k)
  for (b in seq_len(meta$block_info$k)) {
    rows_b <- .blmb_rows_to_data_subset(
      meta$block_info$rows[[b]], meta$mf, data
    )
    setups[[b]] <- do.call(
      Prior_Setup,
      c(
        list(formula = formula, subset = rows_b),
        ps_args,
        list(...)
      )
    )
  }
  names(setups) <- meta$block_info$ids

  attr(setups, "call") <- call
  attr(setups, "formula") <- formula
  attr(setups, "block") <- block
  attr(setups, "block_info") <- meta$block_info
  class(setups) <- c("block_PriorSetup", "list")
  setups
}
