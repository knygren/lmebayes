# Pure-R batch driver for two-block Gibbs sampling (v6 sweep-outer layout).
#
# Two-block Gibbs alternates:
#   Block 1 — given hyperparameters (fixef gamma, tau2), sample random effects b
#             per neighborhood via block_rNormalGLM / block_rNormalReg.
#   Block 2 — given b, treat each RE column as a Gaussian pseudo-response and
#             update hyperparameters via rglmb (one call per RE component).
#
# v6 runs ALL chains through Block 1, then ALL chains through Block 2, for each
# inner sweep (sweep-outer). Block 1 is split into prep (mu_all + prior_list per
# chain) and draw (block_rNormalGLM / block_rNormalReg per chain); both phases
# are embarrassingly parallel over chains (optional n_cores on Unix/macOS).


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
  # J = number of grouping levels (e.g. neighborhoods); p_re = number of RE
  # components (e.g. intercept + rating_c slopes each get their own column).
  p_re <- length(re_names)
  J    <- length(group_levels)

  # fixef[[k]]: n_chains x p_fixef matrix of Block-2 coefficients (gamma_k)
  # for RE component k. Every chain starts at the same vector (ICM mode or pilot
  # mean), replicated by row — matrix(beta0, byrow = TRUE) fills each row.
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

  # tau2: n_chains x p_re matrix of ING dispersion parameters (one per RE
  # component when pfamily is dIndependent_Normal_Gamma). Same plug-in start
  # for every chain.
  tau2 <- matrix(
    tau2_start,
    nrow = n_chains,
    ncol = p_re,
    byrow = TRUE,
    dimnames = list(NULL, re_names)
  )

  # b: J x p_re x n_chains array of random effects. b[j, k, i] is the draw for
  # group j, RE component k, chain i. Initialized from b_mode (ICM posterior
  # mode of b) broadcast to all chains.
  b <- array(b_start, dim = c(J, p_re, n_chains),
             dimnames = list(group_levels, re_names, NULL))

  # iters: cumulative inner Gibbs iteration counts from Block-2 rglmb calls
  # (used for diagnostics; not part of the statistical model).
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
  # Convert row i of each fixef matrix back to a named vector suitable for
  # build_mu_all() and other glmbayesCore helpers expecting a list of gammas.
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
  # Block 1 prior is multivariate normal on b_j with mean mu_all[j, ] and
  # precision P. base_prior comes from Prior_Setup_lmebayes (template structure).
  out <- list(
    mu         = mu_all,              # J x p_re: prior mean per group (from gamma)
    dispersion = base_prior$dispersion,
    ddef       = base_prior$ddef
  )

  # ddef = TRUE: prior precision is data-defined (fixed); do not touch P.
  if (isTRUE(base_prior$ddef)) {
    out$P <- base_prior$P
    return(out)
  }

  # No ING components: P is fixed (e.g. all dNormal priors on hyperparameters).
  if (!any(ptypes == "dIndependent_Normal_Gamma")) {
    out$P <- base_prior$P
    return(out)
  }

  # ING: Block 2 just drew tau2_k. Refresh the k-th diagonal block of P so the
  # conditional prior on b is N(mu, diag(1/tau2)) across RE components — the
  # standard conjugate Normal-Gamma coupling between hyper variance and RE scale.
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

#' One-chain Block 1 prep: fixef -> mu_all -> prior_list (no sampling)
#' @noRd
.lmebayes_block1_prep_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    ptypes
) {
  fixef_i <- .lmebayes_batch_fixef_chain(batch, i)
  mu_all  <- as.matrix(glmbayesCore::build_mu_all(
    design, fixef_i, batch$group_levels
  )$mu_all)
  tau2_i  <- batch$tau2[i, ]
  prior_list <- .lmebayes_block1_prior_with_tau2(
    block1_prior, tau2_i, ptypes, batch$re_names, mu_all
  )
  list(mu_all = mu_all, prior_list = prior_list)
}

#' All-chain Block 1 prep: mu_all and prior_list for every chain
#'
#' Embarrassingly parallel over chain index (optional \code{n_cores}).
#' @noRd
block1_prep_all_chains <- function(
    batch,
    design,
    block1_prior,
    ptypes,
    n_cores = NULL,
    progbar = FALSE
) {
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)

  prep_i <- function(i) {
    if (show_bar) .lmebayes_progress_bar(i, n)
    .lmebayes_block1_prep_one_chain(
      batch        = batch,
      i            = i,
      design       = design,
      block1_prior = block1_prior,
      ptypes       = ptypes
    )
  }

  prep_list <- .lmebayes_lapply_chains(n, prep_i, n_cores = n_cores)
  if (show_bar) .lmebayes_progress_bar_finish()

  structure(
    list(
      mu_all      = lapply(prep_list, `[[`, "mu_all"),
      prior_lists = lapply(prep_list, `[[`, "prior_list")
    ),
    class = "lmebayes_block1_prep"
  )
}

#' One-chain Block 1 draw given a prepared prior_list
#' @noRd
.lmebayes_block1_draw_one_chain <- function(
    prior_list,
    design,
    family,
    is_gaussian,
    group_levels
) {
  if (is_gaussian) {
    block_out <- glmbayesCore::block_rNormalReg(
      n          = 1L,
      y          = design$y,
      x          = design$Z,
      block      = design$groups,
      prior_list = prior_list
    )
  } else {
    block_out <- glmbayesCore::block_rNormalGLM(
      n            = 1L,
      y            = design$y,
      x            = design$Z,
      block        = design$groups,
      prior_list   = prior_list,
      family       = family,
      use_parallel = FALSE,
      verbose      = FALSE,
      progbar      = FALSE
    )
  }

  b_draw <- block_out$coefficients
  rn <- rownames(b_draw)
  if (!is.null(rn)) {
    ord <- match(group_levels, rn)
    if (any(is.na(ord))) {
      stop("Block 1 group ids do not match group_levels.", call. = FALSE)
    }
    b_draw <- b_draw[ord, , drop = FALSE]
  }
  b_draw
}

#' All-chain Block 1 draw from prepared prior_lists
#'
#' Embarrassingly parallel over chain index (optional \code{n_cores}).
#' Updates \code{batch$b} in place.
#' @noRd
block1_draw_all_chains <- function(
    batch,
    prep,
    design,
    family,
    n_cores = NULL,
    progbar = FALSE
) {
  is_gaussian <- identical(family$family, "gaussian")
  n <- batch$n
  show_bar <- isTRUE(progbar) && n > 1L &&
    (is.null(n_cores) || as.integer(n_cores[1L]) < 2L)
  prior_lists <- prep$prior_lists
  if (length(prior_lists) != n) {
    stop("length(prep$prior_lists) must equal batch$n.", call. = FALSE)
  }

  draw_i <- function(i) {
    if (show_bar) .lmebayes_progress_bar(i, n)
    .lmebayes_block1_draw_one_chain(
      prior_list   = prior_lists[[i]],
      design       = design,
      family       = family,
      is_gaussian  = is_gaussian,
      group_levels = batch$group_levels
    )
  }

  b_draws <- .lmebayes_lapply_chains(n, draw_i, n_cores = n_cores)
  if (show_bar) .lmebayes_progress_bar_finish()

  for (i in seq_len(n)) {
    batch$b[, , i] <- b_draws[[i]]
  }
  batch
}

#' Apply FUN to each chain index, optionally in parallel (Unix/macOS only)
#' @noRd
.lmebayes_lapply_chains <- function(n, FUN, n_cores = NULL) {
  idx <- seq_len(n)
  if (is.null(n_cores)) {
    return(lapply(idx, FUN))
  }
  n_cores <- as.integer(n_cores[1L])
  if (!is.finite(n_cores) || n_cores < 2L) {
    return(lapply(idx, FUN))
  }
  n_cores <- min(n_cores, n)
  if (.Platform$OS.type == "windows") {
    warning(
      "Chain-parallel Block 1 (n_cores > 1) is not supported on Windows; ",
      "using sequential lapply.",
      call. = FALSE
    )
    return(lapply(idx, FUN))
  }
  parallel::mclapply(idx, FUN, mc.cores = n_cores)
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
    is_gaussian
) {
  prep <- .lmebayes_block1_prep_one_chain(
    batch        = batch,
    i            = i,
    design       = design,
    block1_prior = block1_prior,
    ptypes       = ptypes
  )
  b_draw <- .lmebayes_block1_draw_one_chain(
    prior_list   = prep$prior_list,
    design       = design,
    family       = family,
    is_gaussian  = is_gaussian,
    group_levels = batch$group_levels
  )
  batch$b[, , i] <- b_draw
  batch
}

#' Align b vector to X_hyper row order (group levels)
#' @noRd
.lmebayes_align_b_to_xhyper <- function(b_vec, X_k, group_levels) {
  # Block 2 regresses b_k on X_hyper for RE component k. Row order of X_k must
  # match the pseudo-response y_k = b_i[, k]. This helper handles three cases:
  # named b + rownamed X, named b only, or positional alignment via group_levels.
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
  # Current RE draw for chain i (just updated in Block 1 of this sweep).
  b_i <- batch$b[, , i, drop = FALSE]
  b_i <- matrix(b_i, nrow = nrow(b_i), ncol = ncol(b_i),
                dimnames = dimnames(b_i)[1:2])

  # One rglmb call per RE component: Gaussian GLM with response = b_k and
  # design X_hyper. This is the conditional draw of gamma_k (and tau2_k if ING).
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
      family  = stats::gaussian(),  # pseudo-likelihood; real prior is pfamily
      pfamily = pf,
      verbose = FALSE
    )

    # Store gamma draw (coef.mode is the single stored draw when n = 1).
    cn <- colnames(batch$fixef[[k]])
    coef_k <- fit_k$coef.mode
    if (!is.null(names(coef_k))) {
      batch$fixef[[k]][i, names(coef_k)] <- coef_k
    } else {
      batch$fixef[[k]][i, ] <- coef_k
    }

    # ING: also store dispersion (tau2) and accumulate inner iteration count.
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
    n_cores = NULL,
    progbar = FALSE
) {
  n <- batch$n
  .lmebayes_print_block1_phase("prep", "enter", n)
  # Phase 1 — prep (embarrassingly parallel): mu_all + prior_list per chain.
  prep <- block1_prep_all_chains(
    batch        = batch,
    design       = design,
    block1_prior = block1_prior,
    ptypes       = ptypes,
    n_cores      = n_cores,
    progbar      = FALSE
  )
  .lmebayes_print_block1_phase("prep", "exit", n)
  .lmebayes_print_block1_phase("draw", "enter", n)
  # Phase 2 — draw (embarrassingly parallel): block_rNormalReg / block_rNormalGLM.
  batch <- block1_draw_all_chains(
    batch   = batch,
    prep    = prep,
    design  = design,
    family  = family,
    n_cores = n_cores,
    progbar = progbar
  )
  .lmebayes_print_block1_phase("draw", "exit", n)
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
  # Sweep-outer step: every chain gets one Block-2 pass (rglmb per RE component).
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
  # Initial tau2 for each RE component before any Gibbs sweeps:
  # dNormal -> fixed prior dispersion; ING -> lower bound disp_lower (plug-in).
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
  # Convert internal batch arrays to the list structure expected by rglmerb_v6 /
  # run_short_chains (fixef_draws, coefficients long form, etc.).
  re_names     <- batch$re_names
  group_levels <- batch$group_levels
  n            <- batch$n
  J            <- length(group_levels)

  # fixef_draws[[k]]: n_chains x p_fixef matrix (one row per stored draw).
  fixef_draws <- lapply(batch$fixef, function(mat) {
    out <- mat
    rownames(out) <- NULL
    out
  })
  names(fixef_draws) <- re_names

  if (isTRUE(collect_block1)) {
    # Long-format RE table: n * J rows, one row per (chain, group) with b columns.
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

  # Posterior mean of gamma across chains, for mu_all_last diagnostic output.
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
