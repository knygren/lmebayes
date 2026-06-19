#' Initialize batch state for sweep-outer R Gibbs driver (v6)
#' @noRd
.lmebayes_init_batch_state <- function(
    n_chains,
    start_fixef,
    b_start,
    tau2_start,
    re_names,
    group_levels
) {
  p_re <- length(re_names)
  J    <- length(group_levels)

  fixef <- lapply(start_fixef, function(beta0) {
    mat <- matrix(
      beta0,
      nrow = n_chains,
      ncol = length(beta0),
      byrow = TRUE,
      dimnames = list(NULL, names(beta0))
    )
    mat
  })
  names(fixef) <- re_names

  tau2 <- matrix(
    tau2_start,
    nrow = n_chains,
    ncol = p_re,
    byrow = TRUE,
    dimnames = list(NULL, re_names)
  )

  b <- array(b_start, dim = c(J, p_re, n_chains),
             dimnames = list(group_levels, re_names, NULL))

  iters <- matrix(0, nrow = n_chains, ncol = p_re)
  colnames(iters) <- re_names

  list(
    n     = n_chains,
    fixef = fixef,
    tau2  = tau2,
    b     = b,
    iters = iters,
    re_names     = re_names,
    group_levels = group_levels
  )
}

#' Extract chain-i fixef list from batch state
#' @noRd
.lmebayes_batch_fixef_chain <- function(batch, i) {
  lapply(batch$fixef, function(mat) {
    v <- mat[i, , drop = TRUE]
    stats::setNames(as.numeric(v), colnames(mat))
  })
}

#' Refresh Block 1 prior precision for ING components (mirrors C++ twoBlockGibbs)
#' @noRd
.lmebayes_block1_prior_with_tau2 <- function(
    base_prior,
    tau2_vec,
    ptypes,
    re_names,
    mu_all
) {
  out <- list(
    mu         = mu_all,
    dispersion = base_prior$dispersion,
    ddef       = base_prior$ddef
  )
  if (isTRUE(base_prior$ddef)) {
    out$P <- base_prior$P
    return(out)
  }
  if (!any(ptypes == "dIndependent_Normal_Gamma")) {
    out$P <- base_prior$P
    return(out)
  }
  P1 <- base_prior$P
  for (k in seq_along(re_names)) {
    if (ptypes[[k]] != "dIndependent_Normal_Gamma") next
    P1[k, ] <- 0
    P1[, k] <- 0
    P1[k, k] <- 1 / tau2_vec[k]
  }
  out$P <- P1
  out
}

#' One-chain Block 1 update (writes batch$b[,,i])
#' @noRd
.lmebayes_block1_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    family,
    ptypes,
    is_gaussian,
    sweep_m = NULL,
    seed = NULL,
    seed_offset = 0L
) {
  if (!is.null(seed) && identical(sweep_m, 1L)) {
    set.seed(as.integer(seed + seed_offset + i + 1L))
  }

  fixef_i <- .lmebayes_batch_fixef_chain(batch, i)
  mu_all  <- as.matrix(glmbayesCore::build_mu_all(
    design, fixef_i, batch$group_levels
  )$mu_all)

  tau2_i <- batch$tau2[i, ]
  pl1    <- .lmebayes_block1_prior_with_tau2(
    block1_prior, tau2_i, ptypes, batch$re_names, mu_all
  )

  wt <- rep(1, length(design$y))
  if (!is_gaussian) {
    b_i <- batch$b[, , i, drop = FALSE]
    b_i <- matrix(b_i, nrow = nrow(b_i), ncol = ncol(b_i),
                  dimnames = dimnames(b_i)[1:2])
    mode_w <- glmbayesCore::two_block_mode_weights(
      x            = design$Z,
      block        = design$groups,
      b_mode       = b_i,
      family       = family,
      group_levels = batch$group_levels
    )
    wt <- mode_w$weights
  }

  if (is_gaussian) {
    block_out <- glmbayesCore::block_rNormalReg(
      n          = 1L,
      y          = design$y,
      x          = design$Z,
      block      = design$groups,
      prior_list = pl1,
      weights    = wt
    )
  } else {
    block_out <- glmbayesCore::block_rNormalGLM(
      n            = 1L,
      y            = design$y,
      x            = design$Z,
      block        = design$groups,
      prior_list   = pl1,
      weights      = wt,
      family       = family,
      use_parallel = FALSE,
      verbose      = FALSE,
      progbar      = FALSE
    )
  }

  b_draw <- block_out$coefficients
  rn <- rownames(b_draw)
  if (!is.null(rn)) {
    ord <- match(batch$group_levels, rn)
    if (any(is.na(ord))) {
      stop("Block 1 group ids do not match group_levels.", call. = FALSE)
    }
    b_draw <- b_draw[ord, , drop = FALSE]
  }
  batch$b[, , i] <- b_draw
  batch
}

#' Align b vector to X_hyper row order (group levels)
#' @noRd
.lmebayes_align_b_to_xhyper <- function(b_vec, X_k, group_levels) {
  rn <- rownames(X_k)
  if (is.null(rn)) {
    if (length(b_vec) != nrow(X_k)) {
      stop(
        "length(b) (", length(b_vec), ") must equal nrow(X_hyper) (",
        nrow(X_k), ") when X_hyper has no rownames.",
        call. = FALSE
      )
    }
    return(b_vec)
  }
  if (!is.null(names(b_vec))) {
    miss <- setdiff(rn, names(b_vec))
    if (length(miss) > 0L) {
      stop(
        "Group level(s) missing from b: ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    return(unname(b_vec[rn]))
  }
  if (length(b_vec) != length(group_levels) ||
      length(b_vec) != nrow(X_k)) {
    stop(
      "b and X_hyper row counts disagree (b: ", length(b_vec),
      ", X_hyper: ", nrow(X_k), ", group_levels: ", length(group_levels), ").",
      call. = FALSE
    )
  }
  names(b_vec) <- group_levels
  miss <- setdiff(rn, group_levels)
  if (length(miss) > 0L) {
    stop(
      "X_hyper rownames do not match group_levels; missing in groups: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  b_vec[rn]
}

#' One-chain Block 2 update (writes batch$fixef, batch$tau2, batch$iters)
#' @noRd
.lmebayes_block2_one_chain <- function(
    batch,
    i,
    design,
    pfamily_list,
    ptypes
) {
  b_i <- batch$b[, , i, drop = FALSE]
  b_i <- matrix(b_i, nrow = nrow(b_i), ncol = ncol(b_i),
                dimnames = dimnames(b_i)[1:2])

  for (k in batch$re_names) {
    X_k <- as.matrix(design$X_hyper[[k]])
    y_k <- .lmebayes_align_b_to_xhyper(
      b_vec        = b_i[, k],
      X_k          = X_k,
      group_levels = batch$group_levels
    )
    pf  <- pfamily_list[[k]]

    fit_k <- glmbayesCore::rglmb(
      n       = 1L,
      y       = y_k,
      x       = X_k,
      family  = stats::gaussian(),
      pfamily = pf,
      verbose = FALSE
    )

    cn <- colnames(batch$fixef[[k]])
    coef_k <- fit_k$coef.mode
    if (!is.null(names(coef_k))) {
      batch$fixef[[k]][i, names(coef_k)] <- coef_k
    } else {
      batch$fixef[[k]][i, ] <- coef_k
    }
    if (ptypes[[k]] == "dIndependent_Normal_Gamma") {
      batch$tau2[i, k] <- fit_k$dispersion[1L]
      it_k <- if (!is.null(fit_k$iters)) fit_k$iters[1L, 1L] else 1L
      batch$iters[i, k] <- batch$iters[i, k] + it_k
    } else {
      batch$iters[i, k] <- batch$iters[i, k] + 1L
    }
  }
  batch
}

#' Block 1 batch: update random effects for all chains
#' @noRd
block1_all_chains <- function(
    batch,
    design,
    block1_prior,
    family,
    ptypes,
    sweep_m = NULL,
    seed = NULL,
    seed_offset = 0L,
    progbar = FALSE
) {
  is_gaussian <- identical(family$family, "gaussian")
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L

  for (i in seq_len(n)) {
    if (show_bar) .lmebayes_progress_bar(i, n)
    batch <- .lmebayes_block1_one_chain(
      batch        = batch,
      i            = i,
      design       = design,
      block1_prior = block1_prior,
      family       = family,
      ptypes       = ptypes,
      is_gaussian  = is_gaussian,
      sweep_m      = sweep_m,
      seed         = seed,
      seed_offset  = seed_offset
    )
  }
  if (show_bar) .lmebayes_progress_bar_finish()
  batch
}

#' Block 2 batch: update fixed effects / tau2 for all chains
#' @noRd
block2_all_chains <- function(
    batch,
    design,
    pfamily_list,
    ptypes,
    progbar = FALSE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L

  for (i in seq_len(n)) {
    if (show_bar) .lmebayes_progress_bar(i, n)
    batch <- .lmebayes_block2_one_chain(
      batch        = batch,
      i            = i,
      design       = design,
      pfamily_list = pfamily_list,
      ptypes       = ptypes
    )
  }
  if (show_bar) .lmebayes_progress_bar_finish()
  batch
}

#' Starting tau2 vector from pfamily_list (plug-in dispersions)
#' @noRd
.lmebayes_tau2_start_from_pfamily <- function(pfamily_list, re_names) {
  vapply(re_names, function(k) {
    pl <- pfamily_list[[k]]$prior_list
    pf <- pfamily_list[[k]]$pfamily
    if (pf == "dNormal") {
      pl$dispersion
    } else if (pf == "dIndependent_Normal_Gamma") {
      pl$disp_lower
    } else {
      stop("Unsupported pfamily: ", pf, call. = FALSE)
    }
  }, numeric(1))
}

#' Pack batch state into run_short_chains contract
#' @noRd
.lmebayes_pack_batch_draws <- function(
    batch,
    design,
    collect_block1 = TRUE
) {
  re_names     <- batch$re_names
  group_levels <- batch$group_levels
  n            <- batch$n
  J            <- length(group_levels)

  fixef_draws <- lapply(batch$fixef, function(mat) {
    out <- mat
    rownames(out) <- NULL
    out
  })
  names(fixef_draws) <- re_names

  if (isTRUE(collect_block1)) {
    grp_col <- design$group_name
    if (is.null(grp_col) || !nzchar(grp_col)) grp_col <- "group"
    coef_rows <- vector("list", n)
    for (i in seq_len(n)) {
      draw_df <- data.frame(
        draw = rep(i, J),
        stringsAsFactors = FALSE
      )
      draw_df[[grp_col]] <- group_levels
      for (k in re_names) {
        draw_df[[k]] <- batch$b[, k, i]
      }
      coef_rows[[i]] <- draw_df
    }
    coefficients <- do.call(rbind, coef_rows)
    rownames(coefficients) <- NULL
  } else {
    coefficients <- NULL
  }

  fixef_mean <- lapply(batch$fixef, colMeans)
  mu_all_last <- as.matrix(glmbayesCore::build_mu_all(
    design, fixef_mean, group_levels
  )$mu_all)

  list(
    fixef_draws            = fixef_draws,
    dispersion_fixef_draws = batch$tau2,
    iters_fixef_draws      = batch$iters,
    coefficients           = coefficients,
    mu_all_last            = mu_all_last
  )
}
