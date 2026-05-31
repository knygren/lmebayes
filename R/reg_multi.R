#' @keywords internal
.mrglmb_check_inputs <- function(y, x, prior_list) {
  if (missing(prior_list)) {
    stop("'prior_list' is required.", call. = FALSE)
  }
  y_mat <- as.matrix(y)
  x <- as.matrix(x)
  l1 <- ncol(y_mat)
  if (l1 < 1L) {
    stop("y must have at least one column.", call. = FALSE)
  }
  p <- ncol(x)
  if (p < 1L) {
    stop("x must have at least one column.", call. = FALSE)
  }
  if (nrow(x) != nrow(y_mat)) {
    stop("nrow(x) must equal nrow(y).", call. = FALSE)
  }
  coef_names <- colnames(y_mat)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }
  pred_names <- colnames(x)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }
  list(
    y_mat = y_mat,
    x = x,
    l1 = l1,
    p = p,
    coef_names = coef_names,
    pred_names = pred_names
  )
}

#' @keywords internal
.mrglmb_n_draw <- function(n) {
  n_draw <- if (length(n) > 1L) length(n) else as.integer(n)
  if (!is.finite(n_draw) || n_draw < 1L) {
    stop(
      "'n' must be a positive scalar or a vector whose length defines the number of draws.",
      call. = FALSE
    )
  }
  n_draw
}

#' @keywords internal
.mrglmb_normalize_prior_lists <- function(prior_list, l1, p, validate_fn) {
  if (!is.list(prior_list)) {
    stop(
      "prior_list must be a list of length ncol(y) of per-column prior lists.",
      call. = FALSE
    )
  }
  if (!is.null(prior_list$mu) || !is.null(prior_list$Sigma)) {
    stop(
      "prior_list must be a list of prior_list objects (one per column of y), ",
      "not a single prior_list with components mu and Sigma.",
      call. = FALSE
    )
  }
  if (length(prior_list) != l1) {
    stop("length(prior_list) must equal ncol(y) = ", l1, ".", call. = FALSE)
  }
  lapply(seq_len(l1), function(j) {
    validate_fn(prior_list[[j]], j = j, p = p)
  })
}

#' @keywords internal
.mrglmb_assemble <- function(block_results,
                             coef_names,
                             call,
                             y_mat,
                             x,
                             l1,
                             p,
                             prior_lists,
                             pred_names) {
  outlist <- setNames(block_results, coef_names)
  attr(outlist, "call")        <- call
  attr(outlist, "y")           <- y_mat
  attr(outlist, "x")           <- x
  attr(outlist, "l1")          <- l1
  attr(outlist, "p")           <- p
  attr(outlist, "prior_lists") <- prior_lists
  attr(outlist, "coef_names")  <- coef_names
  attr(outlist, "pred_names")  <- pred_names
  class(outlist) <- "mrglmb"
  outlist
}

#' @keywords internal
.validate_rindep_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma)) {
    stop("prior_list[[", j, "]] must contain 'Sigma'.", call. = FALSE)
  }
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop("prior_list[[", j, "]] must contain 'shape' and 'rate'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  S <- as.matrix(pl$Sigma)
  if (nrow(S) != p || ncol(S) != p) {
    stop(
      "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
      call. = FALSE
    )
  }
  .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))

  shape <- as.numeric(pl$shape)
  rate <- as.numeric(pl$rate)
  if (length(shape) != 1L || !is.finite(shape)) {
    stop("prior_list[[", j, "]]$shape must be a finite scalar.", call. = FALSE)
  }
  if (length(rate) != 1L || !is.finite(rate)) {
    stop("prior_list[[", j, "]]$rate must be a finite scalar.", call. = FALSE)
  }

  out <- list(mu = mu, Sigma = S, shape = shape, rate = rate)
  if (!is.null(pl$max_disp_perc)) {
    out$max_disp_perc <- as.numeric(pl$max_disp_perc)
  }
  if (!is.null(pl$disp_lower)) {
    out$disp_lower <- pl$disp_lower
  }
  if (!is.null(pl$disp_upper)) {
    out$disp_upper <- pl$disp_upper
  }
  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  out
}

#' @keywords internal
.validate_normal_gamma_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma) && is.null(pl$P)) {
    stop("prior_list[[", j, "]] must contain 'Sigma' or 'P'.", call. = FALSE)
  }
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop("prior_list[[", j, "]] must contain 'shape' and 'rate'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  out <- list(mu = mu)
  if (!is.null(pl$Sigma)) {
    S <- as.matrix(pl$Sigma)
    if (nrow(S) != p || ncol(S) != p) {
      stop(
        "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))
    out$Sigma <- S
  }
  if (!is.null(pl$P)) {
    P <- as.matrix(pl$P)
    if (nrow(P) != p || ncol(P) != p) {
      stop(
        "prior_list[[", j, "]]$P must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(P, label = paste0("prior_list[[", j, "]]$P"))
    out$P <- P
  }

  shape <- as.numeric(pl$shape)
  rate <- as.numeric(pl$rate)
  if (length(shape) != 1L || !is.finite(shape)) {
    stop("prior_list[[", j, "]]$shape must be a finite scalar.", call. = FALSE)
  }
  if (length(rate) != 1L || !is.finite(rate)) {
    stop("prior_list[[", j, "]]$rate must be a finite scalar.", call. = FALSE)
  }
  out$shape <- shape
  out$rate <- rate

  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  if (!is.null(pl$max_disp_perc)) {
    out$max_disp_perc <- as.numeric(pl$max_disp_perc)
  }
  if (!is.null(pl$disp_lower)) {
    out$disp_lower <- pl$disp_lower
  }
  if (!is.null(pl$disp_upper)) {
    out$disp_upper <- pl$disp_upper
  }
  if (!is.null(pl$Precision)) {
    out$Precision <- pl$Precision
  }
  out
}

#' @keywords internal
.validate_normal_prior_list <- function(pl, j, p) {
  if (!is.list(pl)) {
    stop("prior_list[[", j, "]] must be a list.", call. = FALSE)
  }
  if (is.null(pl$mu)) {
    stop("prior_list[[", j, "]] must contain 'mu'.", call. = FALSE)
  }
  if (is.null(pl$Sigma) && is.null(pl$P)) {
    stop("prior_list[[", j, "]] must contain 'Sigma' or 'P'.", call. = FALSE)
  }

  mu <- as.numeric(pl$mu)
  if (length(mu) != p) {
    stop(
      "prior_list[[", j, "]]$mu must have length ncol(x) = ", p, ".",
      call. = FALSE
    )
  }

  out <- list(mu = mu)
  if (!is.null(pl$Sigma)) {
    S <- as.matrix(pl$Sigma)
    if (nrow(S) != p || ncol(S) != p) {
      stop(
        "prior_list[[", j, "]]$Sigma must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(S, label = paste0("prior_list[[", j, "]]$Sigma"))
    out$Sigma <- S
  }
  if (!is.null(pl$P)) {
    P <- as.matrix(pl$P)
    if (nrow(P) != p || ncol(P) != p) {
      stop(
        "prior_list[[", j, "]]$P must be ", p, " x ", p, ".",
        call. = FALSE
      )
    }
    .check_symmetric_pd(P, label = paste0("prior_list[[", j, "]]$P"))
    out$P <- P
  }
  if (!is.null(pl$dispersion)) {
    out$dispersion <- pl$dispersion
  }
  if (!is.null(pl$shape)) {
    out$shape <- pl$shape
  }
  if (!is.null(pl$rate)) {
    out$rate <- pl$rate
  }
  if (!is.null(pl$ddef)) {
    out$ddef <- pl$ddef
  }
  out
}

#' @keywords internal
.check_symmetric_pd <- function(M, label) {
  if (!isSymmetric(M)) {
    stop(label, " must be symmetric.", call. = FALSE)
  }
  tol <- 1e-6
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (!all(ev >= -tol * abs(ev[1L]))) {
    stop(label, " is not positive definite.", call. = FALSE)
  }
  invisible(TRUE)
}
