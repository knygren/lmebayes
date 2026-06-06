#' Row-block (BY-style) Bayesian linear models
#'
#' @description
#' Fits one \code{\link[glmbayes]{lmb}} per observation block (SAS \code{BY}-style split on
#' rows), sharing the same formula on each subset. Contrast with
#' \code{\link[glmbayes]{lmb}} on a \code{cbind(...)} response (several response columns) and
#' \code{\link{block_rNormalGLM}} (Gibbs conditional draws, matrix API).
#'
#' @param block Block partition: \code{factor} or vector of length \code{nrow(data)}
#'   (after \code{model.frame}), a column name in \code{data}, \code{l2_blocks}
#'   counts, or a list of row index vectors (see \code{\link{normalize_block}}).
#' @name block_lmb
#' @family modelfuns
NULL

#' Prior setup for row-block \code{\link[glmbayes]{lmb}} / \code{block_glmb}
#'
#' Runs \code{\link[glmbayes]{Prior_Setup}} on each block subset of the data.
#'
#' @param formula A \code{\link{formula}} with a single response.
#' @param block Block partition: \code{factor} or vector of length \code{nrow(data)}
#'   (after \code{model.frame}), a column name in \code{data}, \code{l2_blocks}
#'   counts, or a list of row index vectors (see \code{\link{normalize_block}}).
#' @inheritParams glmbayes::Prior_Setup
#' @return A named list of class \code{"block_PriorSetup"}. Each element is a
#'   \code{\link[glmbayes]{Prior_Setup}} result for one block.
#' @family prior
#' @seealso \code{\link{block_lmb}}, \code{\link[glmbayes]{multi_prior_setup}},
#'   \code{\link{normalize_block}}
#' @export
block_prior_setup <- function(
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
      "block_prior_setup() supports family = gaussian() or poisson() only.",
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
      glmbayes::Prior_Setup,
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

#' @describeIn block_lmb Gaussian \code{\link[glmbayes]{lmb}} fit per row block.
#' @param pfamily A single \code{\link[glmbayes]{pfamily}} recycled to every block, or
#'   use \code{pfamily_list} of length \code{k} (number of blocks).
#' @param pfamily_list Optional list of \code{pfamily} objects, one per block.
#' @inheritParams glmbayes::lmb
#' @return A named list of class \code{"blmb"} (list of \code{"lmb"} fits).
#' @example inst/examples/Ex_block_lmb.R
#' @export
block_lmb <- function(
    formula,
    block,
    pfamily = NULL,
    pfamily_list = NULL,
    n = 1000,
    data,
    subset,
    weights,
    na.action,
    method = "qr",
    model = TRUE,
    x = TRUE,
    y = TRUE,
    qr = TRUE,
    singular.ok = TRUE,
    contrasts = NULL,
    offset,
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

  lmb_args <- list(
    n = n,
    data = data,
    method = method,
    model = model,
    x = x,
    y = y,
    qr = qr,
    singular.ok = singular.ok,
    contrasts = contrasts,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose
  )
  if (!missing(subset)) lmb_args$subset <- subset
  if (!missing(weights)) lmb_args$weights <- weights
  if (!missing(na.action)) lmb_args$na.action <- na.action
  if (!missing(offset)) lmb_args$offset <- offset
  if (length(list(...))) {
    lmb_args <- c(lmb_args, list(...))
  }

  block_results <- vector("list", k)
  for (b in seq_len(k)) {
    rows_b <- .blmb_rows_to_data_subset(
      meta$block_info$rows[[b]], meta$mf, data
    )
    fit_b <- do.call(
      lmb,
      c(
        list(
          formula = formula,
          pfamily = pfamily_lists[[b]],
          subset = rows_b
        ),
        lmb_args
      )
    )
    fit_b$call <- .blmb_lmb_display_call(mc, formula, rows_b)
    block_results[[b]] <- fit_b
  }

  .blmb_assemble(
    block_results = block_results,
    block_ids = meta$block_info$ids,
    call = call,
    formula = formula,
    block = block,
    block_info = meta$block_info,
    p = p,
    pred_names = meta$pred_names,
    pfamily_lists = pfamily_lists
  )
}

#' @keywords internal
.blmb_formula_block_meta <- function(
    formula,
    block,
    data,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL
) {
  mf_args <- list(
    formula = formula,
    data = data,
    drop.unused.levels = TRUE
  )
  if (!is.null(subset)) mf_args$subset <- subset
  if (!is.null(weights)) mf_args$weights <- weights
  if (!is.null(na.action)) mf_args$na.action <- na.action
  if (!is.null(offset)) mf_args$offset <- offset
  if (!is.null(contrasts)) mf_args$contrasts <- contrasts
  mf <- do.call(stats::model.frame, mf_args)

  l2 <- nrow(mf)
  block_vec <- .blmb_resolve_block(block, data, mf, l2)
  block_info <- normalize_block(block_vec, l2)

  mt <- attr(mf, "terms")
  x_mat <- stats::model.matrix(mt, mf, contrasts)
  p <- ncol(x_mat)
  pred_names <- colnames(x_mat)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }

  list(
    mf = mf,
    block_info = block_info,
    p = p,
    pred_names = pred_names
  )
}

#' @keywords internal
.blmb_rows_to_data_subset <- function(rows_mf, mf, data) {
  rn <- rownames(mf)
  if (is.null(rn)) {
    return(rows_mf)
  }
  if (!is.null(rownames(data))) {
    out <- match(rn[rows_mf], rownames(data))
    if (anyNA(out)) {
      stop("Could not map model.frame rows to rownames(data).", call. = FALSE)
    }
    return(out)
  }
  as.integer(rn[rows_mf])
}

#' @keywords internal
.blmb_resolve_block <- function(block, data, mf, l2) {
  if (is.list(block)) {
    return(block)
  }
  if (is.character(block) && length(block) == 1L && block %in% names(data)) {
    rn_mf <- rownames(mf)
    col <- data[[block]]
    if (length(col) == l2 && is.null(rn_mf)) {
      return(col)
    }
    if (!is.null(rn_mf)) {
      if (!is.null(rownames(data))) {
        return(col[match(rn_mf, rownames(data))])
      }
      return(col[as.integer(rn_mf)])
    }
  }
  block <- as.vector(block)
  if (length(block) == l2) {
    return(block)
  }
  stop(
    "'block' must have length nrow(model.frame), be a list of row indices, ",
    "l2_blocks counts, or a single column name in 'data'.",
    call. = FALSE
  )
}

#' @keywords internal
.blmb_lmb_display_call <- function(mc_block, formula, rows_b) {
  cl <- call("lmb", formula = formula, subset = rows_b)
  pass <- c(
    "n", "data", "weights", "na.action", "offset",
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
.blmb_assemble <- function(
    block_results,
    block_ids,
    call,
    formula,
    block,
    block_info,
    p,
    pred_names,
    pfamily_lists = NULL
) {
  outlist <- setNames(block_results, block_ids)
  attr(outlist, "call") <- call
  attr(outlist, "formula") <- formula
  attr(outlist, "block") <- block
  attr(outlist, "block_info") <- block_info
  attr(outlist, "k") <- block_info$k
  attr(outlist, "p") <- p
  attr(outlist, "pred_names") <- pred_names
  if (!is.null(pfamily_lists)) {
    attr(outlist, "pfamily_lists") <- pfamily_lists
  }
  class(outlist) <- "blmb"
  outlist
}

#' Per-block design-matrix rank (row-block formulas)
#'
#' @description
#' For each block in \code{block}, builds the \code{\link{model.matrix}} on that
#' block's \code{model.frame} rows and marks blocks with \code{qr(x)$rank == ncol(x)}.
#' Saturated blocks (\code{n = p}) need \code{dispersion} in \code{\link[glmbayes]{Prior_Setup}}.
#'
#' @param formula Model formula shared across blocks.
#' @inheritParams block_lmb
#' @return A list with \code{keep} and \code{drop} (block id character vectors),
#'   \code{table} (data frame with \code{id}, \code{n}, \code{rank}, \code{p},
#'   \code{full_rank}), and \code{block_info}.
#' @keywords internal
.blmb_blocks_full_rank <- function(
    formula,
    block,
    data,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL
) {
  meta <- .blmb_formula_block_meta(
    formula = formula,
    block = block,
    data = data,
    subset = subset,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts
  )
  mt <- attr(meta$mf, "terms")
  k <- meta$block_info$k
  tab <- data.frame(
    id = character(k),
    n = integer(k),
    rank = integer(k),
    p = integer(k),
    full_rank = logical(k),
    stringsAsFactors = FALSE
  )
  for (b in seq_len(k)) {
    rows <- meta$block_info$rows[[b]]
    mf_b <- meta$mf[rows, , drop = FALSE]
    x <- stats::model.matrix(mt, mf_b, contrasts)
    p <- ncol(x)
    rk <- qr(x)$rank
    n_b <- nrow(x)
    tab[b, ] <- list(
      id = meta$block_info$ids[b],
      n = n_b,
      rank = rk,
      p = p,
      full_rank = rk == p
    )
  }
  keep <- tab$id[tab$full_rank]
  drop <- tab$id[!tab$full_rank]
  list(keep = keep, drop = drop, table = tab, block_info = meta$block_info)
}

#' Per-block design-matrix rank from a pre-formed matrix (matrix interface)
#'
#' @description
#' Matrix-interface analogue of \code{.blmb_blocks_full_rank}.  Accepts a
#' pre-formed design matrix \code{x} and a block specification (factor, integer
#' vector, or list of row-index vectors) rather than a formula and data frame.
#' Used internally by \code{block_check_identifiability_xy} and by
#' \code{\link{block_rNormalGLM}}.
#'
#' @param x Numeric matrix \code{(l2 x l1)}: the full design matrix.
#' @param block_info Block partition as returned by \code{normalize_block()}.
#' @return A list with \code{keep}, \code{drop}, \code{table}.
#' @keywords internal
.blmb_blocks_full_rank_xy <- function(x, block_info) {
  k  <- block_info$k
  l1 <- ncol(x)
  tab <- data.frame(
    id        = character(k),
    n         = integer(k),
    rank      = integer(k),
    p         = integer(k),
    full_rank = logical(k),
    stringsAsFactors = FALSE
  )
  for (b in seq_len(k)) {
    rows <- block_info$rows[[b]]
    x_b  <- x[rows, , drop = FALSE]
    rk   <- qr(x_b)$rank
    n_b  <- nrow(x_b)
    tab[b, ] <- list(
      id        = block_info$ids[b],
      n         = n_b,
      rank      = rk,
      p         = l1,
      full_rank = rk == l1
    )
  }
  keep <- tab$id[tab$full_rank]
  drop <- tab$id[!tab$full_rank]
  list(keep = keep, drop = drop, table = tab)
}

#' Check block-level and hyper-level identifiability (matrix interface)
#'
#' @description
#' Matrix-interface variant of \code{\link{block_check_identifiability}} for
#' use when the design matrix \code{x} and response \code{y} are already formed
#' (e.g. inside a Gibbs loop setup).  Applies the same two-level algorithm:
#'
#' \enumerate{
#'   \item \strong{Level 1:} full column rank of \code{x[rows_b, ]} per block.
#'   \item \strong{Level 2:} full column rank of \code{X_nbhd} restricted to
#'     Level-1-identified blocks.
#' }
#'
#' See \code{\link{block_check_identifiability}} and
#' \code{inst/BLOCK_GIBBS_ERGODICITY.md} for background.
#'
#' @param x Numeric design matrix \code{(l2 x l1)}.
#' @param block Block specification: factor, integer/character vector, or list
#'   of row-index vectors.  Passed to \code{normalize_block()}.
#' @param X_nbhd Optional \code{(k x q)} numeric matrix of group-level
#'   covariates (one row per block, in block-id order or with matching
#'   \code{rownames}).  \code{NULL} assumes an intercept-only hyper design.
#' @param on_failure One of \code{"warn"} (default) or \code{"stop"}.
#' @return Invisibly, the same list structure as
#'   \code{\link{block_check_identifiability}}.
#' @seealso \code{\link{block_check_identifiability}},
#'   \code{\link{block_rNormalGLM}}
#' @export
block_check_identifiability_xy <- function(
    x,
    block,
    X_nbhd     = NULL,
    on_failure = c("warn", "stop")
) {
  on_failure <- match.arg(on_failure)
  x <- as.matrix(x)
  l2 <- nrow(x)
  block_info <- normalize_block(block, l2)
  k  <- block_info$k

  ri   <- .blmb_blocks_full_rank_xy(x, block_info)
  keep <- ri$keep
  drop <- ri$drop
  tab  <- ri$table

  if (length(drop)) {
    message(
      "block_check_identifiability_xy: Level 1 -- ",
      length(drop), " of ", k, " block(s) are rank-deficient:\n  ",
      paste(drop, collapse = ", ")
    )
  } else {
    message("block_check_identifiability_xy: Level 1 -- all ", k, " blocks are full rank.")
  }

  if (is.null(X_nbhd)) {
    X_nbhd <- matrix(1, nrow = k, ncol = 1L,
                     dimnames = list(tab$id, "(Intercept)"))
  }
  if (!is.matrix(X_nbhd) || nrow(X_nbhd) != k) {
    stop("X_nbhd must be a matrix with one row per block (", k, " rows).",
         call. = FALSE)
  }
  if (!is.null(rownames(X_nbhd))) {
    idx <- match(tab$id, rownames(X_nbhd))
    if (anyNA(idx)) {
      stop("rownames(X_nbhd) do not match all block ids.", call. = FALSE)
    }
    X_nbhd <- X_nbhd[idx, , drop = FALSE]
  }

  q <- ncol(X_nbhd)
  if (length(keep) == 0L) {
    l2_rank <- 0L
    l2_ok   <- FALSE
  } else {
    keep_idx <- which(tab$id %in% keep)
    X_sub    <- X_nbhd[keep_idx, , drop = FALSE]
    l2_rank  <- qr(X_sub)$rank
    l2_ok    <- l2_rank == q
  }

  action <- if (l2_ok) "proceed" else on_failure
  if (!l2_ok) {
    msg <- paste0(
      "block_check_identifiability_xy: Level 2 FAILED -- rank of X_nbhd ",
      "restricted to Level-1 blocks is ", l2_rank, " (need ", q, "). ",
      "mu is not identified; the Gibbs chain will be null recurrent in ",
      q - l2_rank, " direction(s)."
    )
    if (on_failure == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  } else {
    message(
      "block_check_identifiability_xy: Level 2 -- mu identified (rank ",
      l2_rank, " == q = ", q, " from ", length(keep), " data-bearing blocks)."
    )
  }

  invisible(list(
    level1_table = tab,
    level1_keep  = keep,
    level1_drop  = drop,
    level2_rank  = l2_rank,
    level2_q     = q,
    level2_ok    = l2_ok,
    action       = action
  ))
}

#' Check block-level and hyper-level identifiability for a coupled block Gibbs model
#'
#' @description
#' Implements a two-level identifiability preflight for block Gibbs samplers:
#'
#' \enumerate{
#'   \item \strong{Level 1 (per block):} identify blocks whose data design matrix
#'     \eqn{X_b} is full column rank.  Rank-deficient blocks contribute zero
#'     likelihood precision and cannot identify their coefficient vector from the
#'     data alone.
#'   \item \strong{Level 2 (hyper):} among the level-1 identified blocks, check
#'     that the hyper design \code{X_nbhd} (neighborhood/group-level covariates
#'     for the prior \eqn{\beta_b \sim N(X_{nbhd,b}\,\mu,\,\Sigma)}) has full
#'     column rank.  If this fails, the population parameter \eqn{\mu} is not
#'     identified by the data-bearing groups even in the limit of an improper
#'     prior on \eqn{\mu}.
#' }
#'
#' When both levels pass, the chain is geometrically ergodic (normal fixed-variance
#' case) and well-behaved (Poisson/GLM case) even with a near-flat prior on
#' \eqn{\mu}.  Non-identified blocks from Level 1 may be retained in the chain
#' as "prior-draw" groups; they do not disrupt ergodicity once Level 2 holds.
#' See \code{inst/BLOCK_GIBBS_ERGODICITY.md} for derivations.
#'
#' For BY-style independent fits (\code{\link{block_lmb}}, \code{\link{block_glmb}})
#' only Level 1 applies; non-identified blocks should be dropped before fitting.
#'
#' @param formula Model formula (data-level, shared across blocks).
#' @param block Block specification: factor, column name, or list of row indices.
#'   See \code{\link{block_lmb}}.
#' @param data A data frame.
#' @param X_nbhd Optional numeric matrix with one row per block (in block-id order)
#'   and \eqn{q} columns of group-level covariates.  If \code{NULL}, an
#'   intercept-only hyper design (\code{matrix(1, k, 1)}) is assumed, and Level 2
#'   reduces to: at least one Level-1-identified block exists.
#' @param subset,weights,na.action,offset,contrasts Passed to
#'   \code{\link[stats]{model.frame}}.
#' @param on_failure One of \code{"warn"} (default) or \code{"stop"}.  Controls
#'   whether a Level-2 failure emits a warning or stops with an error.
#' @return A list (invisibly) with components:
#'   \describe{
#'     \item{level1_table}{Data frame: \code{id}, \code{n}, \code{rank}, \code{p},
#'       \code{full_rank} for every block.}
#'     \item{level1_keep}{Character: block ids passing Level 1.}
#'     \item{level1_drop}{Character: rank-deficient block ids.}
#'     \item{level2_rank}{Integer: rank of \code{X_nbhd[level1_keep, ]}.}
#'     \item{level2_q}{Integer: number of columns of \code{X_nbhd} (target rank).}
#'     \item{level2_ok}{Logical: Level 2 satisfied?}
#'     \item{action}{\code{"proceed"} or \code{"warn"} or \code{"stop"}.}
#'   }
#' @seealso \code{\link{block_lmb}}, \code{\link{block_prior_setup}},
#'   \code{\link{block_glmb}}
#' @export
block_check_identifiability <- function(
    formula,
    block,
    data,
    X_nbhd = NULL,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL,
    on_failure = c("warn", "stop")
) {
  on_failure <- match.arg(on_failure)
  if (missing(data)) data <- environment(formula)

  ri <- .blmb_blocks_full_rank(
    formula    = formula,
    block      = block,
    data       = data,
    subset     = subset,
    weights    = weights,
    na.action  = na.action,
    offset     = offset,
    contrasts  = contrasts
  )

  keep <- ri$keep
  drop <- ri$drop
  tab  <- ri$table
  k    <- nrow(tab)

  if (length(drop)) {
    message(
      "block_check_identifiability: Level 1 -- ",
      length(drop), " of ", k, " block(s) are rank-deficient and cannot",
      " identify their coefficients from the data:\n  ",
      paste(drop, collapse = ", ")
    )
  } else {
    message("block_check_identifiability: Level 1 -- all ", k, " blocks are full rank.")
  }

  if (is.null(X_nbhd)) {
    X_nbhd <- matrix(1, nrow = k, ncol = 1L,
                     dimnames = list(tab$id, "(Intercept)"))
  }

  if (!is.matrix(X_nbhd) || nrow(X_nbhd) != k) {
    stop(
      "X_nbhd must be a matrix with one row per block (", k, " rows).",
      call. = FALSE
    )
  }
  if (!is.null(rownames(X_nbhd))) {
    X_nbhd <- X_nbhd[match(tab$id, rownames(X_nbhd)), , drop = FALSE]
    if (anyNA(rownames(X_nbhd))) {
      stop(
        "rownames(X_nbhd) do not match all block ids. ",
        "Supply X_nbhd with rownames matching block ids, or unnamed.",
        call. = FALSE
      )
    }
  }

  q <- ncol(X_nbhd)

  if (length(keep) == 0L) {
    l2_rank <- 0L
    l2_ok   <- FALSE
  } else {
    keep_idx  <- which(tab$id %in% keep)
    X_sub     <- X_nbhd[keep_idx, , drop = FALSE]
    l2_rank   <- qr(X_sub)$rank
    l2_ok     <- l2_rank == q
  }

  action <- if (l2_ok) "proceed" else on_failure

  if (!l2_ok) {
    msg <- paste0(
      "block_check_identifiability: Level 2 FAILED -- ",
      "rank of X_nbhd restricted to Level-1-identified blocks is ",
      l2_rank, " (need ", q, "). ",
      "The population parameter mu is not identified by the data-bearing blocks; ",
      "the coupled Gibbs chain will be null recurrent in ",
      q - l2_rank, " direction(s) of mu even with a proper prior."
    )
    if (on_failure == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  } else {
    message(
      "block_check_identifiability: Level 2 -- mu identified (rank ",
      l2_rank, " == q = ", q, " from ", length(keep), " data-bearing blocks). ",
      "Proceed with full model."
    )
  }

  invisible(list(
    level1_table = tab,
    level1_keep  = keep,
    level1_drop  = drop,
    level2_rank  = l2_rank,
    level2_q     = q,
    level2_ok    = l2_ok,
    action       = action
  ))
}

#' @keywords internal
.blmb_coef_means_matrix <- function(object) {
  nm <- names(object)
  if (length(nm) < 1L) {
    return(NULL)
  }
  cm <- do.call(rbind, lapply(object, function(fit) fit$coef.means))
  pred <- names(object[[1L]]$coef.means)
  if (is.null(pred)) {
    pred <- colnames(object[[1L]]$x)
  }
  if (!is.null(pred) && ncol(cm) == length(pred)) {
    colnames(cm) <- pred
  }
  rownames(cm) <- nm
  cm
}

#' @keywords internal
.blmb_dic_table <- function(object) {
  nm <- names(object)
  if (length(nm) < 1L) {
    return(NULL)
  }
  pD <- vapply(object, function(fit) fit$pD, numeric(1))
  dic <- vapply(object, function(fit) fit$DIC, numeric(1))
  cbind(pD = pD, DIC = dic)
}
