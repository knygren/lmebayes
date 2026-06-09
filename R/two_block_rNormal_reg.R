#' Two-block Gaussian Gibbs sampler for hierarchical regression
#'
#' Runs the coupled Block~1 / Block~2 Gibbs sampler used by
#' \code{\link{lmerb}}: Block~1 draws group-level random effects \eqn{b_j}
#' via \code{\link[glmbayesCore]{block_rNormalReg}}; Block~2 updates
#' hyper means \eqn{\gamma_k} via \code{\link{multi_rNormal_reg}}.
#'
#' @details
#' With \code{sampling = "replicate"} (default), each of the \code{n} stored
#' draws is obtained by resetting to \code{fixef_start}, running
#' \code{m_convergence} Gibbs steps, and retaining only the final
#' \eqn{(\gamma, b)} state.  This matches the near-iid sampling strategy
#' described in \code{\link{lmerb}} (Nygren 2020 TV bound).
#'
#' @param n Number of stored draws.
#' @param y Response vector of length \code{nrow(x)}.
#' @param x Level-1 design matrix \code{Z} (\code{l2 x p_re}).
#' @param block Grouping factor or block partition of length \code{l2}.
#' @param x_hyper Named list of group-level design matrices \code{X_k}
#'   (\code{J x q_k}), one per column of \code{x}.  Supply
#'   \code{design$X_hyper} from \code{\link{model_setup}}.
#' @param prior_list_block1 Prior for Block~1: must contain \code{P} or
#'   \code{Sigma}, \code{dispersion}, and optionally \code{ddef}.  The
#'   \code{mu} field (\code{p_re x J}) is updated internally from
#'   \code{fixef} at each Gibbs step.
#' @param prior_list_block2 Named list of Block~2 prior lists (one per RE
#'   column), each with \code{mu}, \code{Sigma} or \code{P}, and
#'   \code{dispersion}.  Passed to \code{\link{multi_rNormal_reg}}.
#' @param fixef_start Named list of hyper-parameter vectors \eqn{\gamma_k}
#'   at which each inner chain is initialised.
#' @param re_coef_names Character vector naming columns of \code{x}.  Defaults
#'   to \code{colnames(x)}.
#' @param group_levels Character vector defining row order of Block~1 draws
#'   and column order of \code{mu_all}.  Defaults to \code{levels(block)}.
#' @param group_name Name for the grouping column in \code{coefficients}.
#'   When \code{NULL}, taken from \code{deparse(substitute(block))} or
#'   \code{"group"}.
#' @param m_convergence Number of inner Gibbs steps per stored draw.
#' @param sampling Sampling scheme.  Currently only \code{"replicate"} is
#'   implemented.
#' @param seed Optional RNG seed set before sampling.
#' @param progbar Logical; show a text progress bar over the \code{n} draws.
#' @return Object of class \code{"two_block_rNormal_reg"}: a list with
#'   \code{fixef_draws}, \code{coefficients}, \code{fixef_last},
#'   \code{b_last}, \code{mu_all_last}, and metadata (\code{n},
#'   \code{m_convergence}, \code{sampling}, \code{fixef_start},
#'   \code{re_coef_names}, \code{group_levels}, \code{group_name},
#'   \code{call}).
#' @family simfuncs
#' @seealso \code{\link{lmerb}}, \code{\link{build_mu_all}},
#'   \code{\link[glmbayesCore]{block_rNormalReg}},
#'   \code{\link{multi_rNormal_reg}}
#' @export
two_block_rNormal_reg <- function(
    n,
    y,
    x,
    block,
    x_hyper,
    prior_list_block1,
    prior_list_block2,
    fixef_start,
    re_coef_names = colnames(x),
    group_levels = levels(block),
    group_name = NULL,
    m_convergence = 10L,
    sampling = c("replicate", "chain"),
    seed = NULL,
    progbar = TRUE) {

  cl <- match.call()
  sampling <- match.arg(sampling)
  if (!identical(sampling, "replicate")) {
    stop("Only sampling = \"replicate\" is implemented.", call. = FALSE)
  }

  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  m_convergence <- as.integer(m_convergence[1L])
  if (m_convergence < 1L) {
    stop("'m_convergence' must be at least 1.", call. = FALSE)
  }

  y <- as.vector(y)
  x <- as.matrix(x)
  l2 <- nrow(x)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }

  if (is.null(re_coef_names) || length(re_coef_names) != ncol(x)) {
    re_coef_names <- if (ncol(x) >= 1L) {
      cn <- colnames(x)
      if (is.null(cn) || length(cn) != ncol(x)) paste0("RE", seq_len(ncol(x))) else cn
    } else {
      stop("'x' must have at least one column.", call. = FALSE)
    }
  }
  colnames(x) <- re_coef_names
  re_names <- re_coef_names

  group_levels <- as.character(group_levels)
  if (length(group_levels) < 1L) {
    stop("'group_levels' must contain at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    group_name <- tryCatch(
      deparse(substitute(block))[1L],
      error = function(e) "group"
    )
    if (!nzchar(group_name)) group_name <- "group"
  }

  if (!is.list(x_hyper) || is.data.frame(x_hyper)) {
    stop("'x_hyper' must be a list of design matrices.", call. = FALSE)
  }
  if (length(x_hyper) != length(re_names)) {
    stop(
      "length(x_hyper) must equal ncol(x) = ", length(re_names), ".",
      call. = FALSE
    )
  }
  if (!setequal(names(x_hyper), re_names)) {
    x_hyper <- x_hyper[re_names]
  }

  if (!is.list(prior_list_block2)) {
    stop("'prior_list_block2' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(prior_list_block2), re_names)) {
    stop(
      "names(prior_list_block2) must match re_coef_names.",
      call. = FALSE
    )
  }

  if (!is.list(fixef_start) || is.null(names(fixef_start))) {
    stop("'fixef_start' must be a named list.", call. = FALSE)
  }
  if (!setequal(names(fixef_start), re_names)) {
    stop("names(fixef_start) must match re_coef_names.", call. = FALSE)
  }
  fixef_start <- fixef_start[re_names]

  if (!is.list(prior_list_block1)) {
    stop("'prior_list_block1' must be a list.", call. = FALSE)
  }
  P <- prior_list_block1$P
  if (is.null(P)) {
    Sigma_b <- prior_list_block1$Sigma
    if (is.null(Sigma_b)) {
      stop("prior_list_block1 must contain 'P' or 'Sigma'.", call. = FALSE)
    }
    P <- solve(Sigma_b)
  }
  dispersion <- prior_list_block1$dispersion
  if (is.null(dispersion)) {
    stop("prior_list_block1 must contain 'dispersion'.", call. = FALSE)
  }
  ddef <- if (is.null(prior_list_block1$ddef)) FALSE else prior_list_block1$ddef

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Block 1 argument template; prior_list$mu is updated at each Gibbs step.
  fixef <- fixef_start
  mu_all <- .two_block_mu_all(fixef, x_hyper, re_names, group_levels)
  block1_args <- list(
    n          = 1L,
    y          = y,
    x          = x,
    block      = block,
    prior_list = list(
      mu         = mu_all,
      P          = P,
      dispersion = dispersion,
      ddef       = ddef
    )
  )

  coef_cols <- c("draw", group_name, re_names)
  draw_rows <- vector("list", n)

  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      q_k <- length(fixef_start[[k]])
      matrix(NA_real_, nrow = n, ncol = q_k,
             dimnames = list(NULL, names(fixef_start[[k]])))
    }),
    re_names
  )

  if (isTRUE(progbar)) {
    pb <- utils::txtProgressBar(min = 0L, max = n, style = 3L)
    on.exit(close(pb), add = TRUE)
  }

  b_i <- NULL

  for (i in seq_len(n)) {
    if (isTRUE(progbar)) {
      utils::setTxtProgressBar(pb, i)
    }

    fixef <- fixef_start

    for (m in seq_len(m_convergence)) {

      mu_all <- .two_block_mu_all(fixef, x_hyper, re_names, group_levels)
      block1_args$prior_list$mu <- mu_all
      block_i <- do.call(block_rNormalReg, block1_args)
      b_i <- block_i$coefficients
      if (is.null(rownames(b_i))) {
        rownames(b_i) <- block_i$block_info$ids
      }
      colnames(b_i) <- re_names

      fixef_draw <- multi_rNormal_reg(
        n          = 1L,
        y          = b_i,
        x          = x_hyper,
        prior_list = prior_list_block2,
        progbar    = FALSE
      )
      fixef <- stats::setNames(
        lapply(re_names, function(k) fixef_draw[[k]]$coefficients[1L, ]),
        re_names
      )
    }

    for (k in re_names) {
      fixef_draws[[k]][i, ] <- fixef[[k]]
    }

    J_i <- nrow(b_i)
    draw_df <- data.frame(
      draw = rep(i, J_i),
      stringsAsFactors = FALSE
    )
    draw_df[[group_name]] <- rownames(b_i)
    for (nm in re_names) {
      draw_df[[nm]] <- b_i[, nm]
    }
    draw_rows[[i]] <- draw_df
  }

  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      fixef_draws   = fixef_draws,
      coefficients  = coefficients,
      fixef_last    = fixef,
      b_last        = b_i,
      mu_all_last   = mu_all,
      n             = n,
      m_convergence = m_convergence,
      sampling      = sampling,
      fixef_start   = fixef_start,
      re_coef_names = re_names,
      group_levels  = group_levels,
      group_name    = group_name,
      call          = cl
    ),
    class = "two_block_rNormal_reg"
  )
}

#' @noRd
.two_block_mu_all <- function(fixef, x_hyper, re_names, group_levels) {
  p_re <- length(re_names)
  J    <- length(group_levels)
  mu_all <- matrix(NA_real_, nrow = p_re, ncol = J,
                   dimnames = list(re_names, group_levels))
  for (i in seq_len(p_re)) {
    k       <- re_names[i]
    gamma_k <- fixef[[k]]
    X_k     <- as.matrix(x_hyper[[k]])
    rn      <- rownames(X_k)
    if (is.null(rn)) {
      if (nrow(X_k) != J) {
        stop(
          "nrow(x_hyper[[", k, "]]) must equal length(group_levels).",
          call. = FALSE
        )
      }
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[j, , drop = TRUE] * gamma_k)
      }
    } else {
      for (j in seq_len(J)) {
        mu_all[i, j] <- sum(X_k[group_levels[j], , drop = TRUE] * gamma_k)
      }
    }
  }
  mu_all
}
