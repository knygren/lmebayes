#' Print ICM random-effects mode (all group x RE coefficients)
#' @noRd
.lmebayes_print_ranef_mode_reference <- function(
    ranef_mode,
    re_names,
    group_levels,
    verbose
) {
  if (!isTRUE(verbose) || is.null(ranef_mode)) {
    return(invisible(NULL))
  }
  if (is.null(colnames(ranef_mode))) {
    colnames(ranef_mode) <- re_names
  }
  if (is.null(rownames(ranef_mode))) {
    rownames(ranef_mode) <- group_levels
  }
  cat("--- glmerb: ICM mode reference (Block 1 random effects) ---\n")
  cat(sprintf(
    "  %d group x RE coefficients (%d groups x %d RE cols):\n",
    length(ranef_mode), nrow(ranef_mode), ncol(ranef_mode)
  ))
  for (g in rownames(ranef_mode)) {
    for (k in colnames(ranef_mode)) {
      cat(sprintf("    %s::%-18s  %12.4f\n", g, k, ranef_mode[g, k]))
    }
  }
  cat("\n")
  invisible(NULL)
}

#' Print fixef_main_start (pilot colMeans) before main stage
#' @noRd
.lmebayes_print_fixef_main_start <- function(
    fixef_main_start,
    re_names,
    verbose
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat("--- glmerb: main-stage fixef_start (pilot colMeans) ---\n")
  for (k in re_names) {
    for (nm in names(fixef_main_start[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_main_start[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}

#' Print Block 1 prep/draw sub-phase boundary with wall-clock timestamp
#' @noRd
.lmebayes_print_block1_phase <- function(phase, boundary, n_chains) {
  phase <- as.character(phase)[1L]
  boundary <- as.character(boundary)[1L]
  action <- if (identical(boundary, "enter")) "Entering" else "Exiting"
  phase_label <- if (identical(phase, "prep")) {
    "Block1 prep (mu_all + prior_list)"
  } else if (identical(phase, "draw")) {
    "Block1 draw (block_rNormalGLM / block_rNormalReg)"
  } else {
    phase
  }
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf(
    "[Block1] %s %s, n=%d chains (%s)\n",
    action, phase_label, as.integer(n_chains)[1L], ts
  ))
  utils::flush.console()
  invisible(NULL)
}

#' Print sweep/block enter or exit line with wall-clock timestamp
#' @noRd
.lmebayes_print_sweep_boundary <- function(
    stage_label,
    sweep,
    inner_sweeps,
    phase,
    boundary
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) stage_label <- "stage"
  phase_label <- if (identical(phase, "Block1")) {
    "random effects update"
  } else {
    "fixed effects update"
  }
  action <- if (identical(boundary, "enter")) "Entering" else "Exiting"
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf(
    "[%s] Sweep %d / %d: %s %s (%s)\n",
    stage_label, sweep, inner_sweeps, action, phase_label, ts
  ))
  utils::flush.console()
  invisible(NULL)
}

#' Print per-sweep block diagnostics (fixef + all J x p_re b means vs ICM mode)
#' @noRd
.lmebayes_print_block_diag <- function(
    stage_label,
    sweep,
    inner_sweeps,
    phase,
    batch,
    fixef_mode,
    b_mode,
    re_names,
    group_levels
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) stage_label <- "stage"
  cat(sprintf(
    "--- glmerb [%s sweep %d / %d after %s, n=%d] ---\n",
    stage_label, sweep, inner_sweeps, phase, batch$n
  ))

  if (identical(phase, "Block2") && !is.null(fixef_mode)) {
    cat("  fixef (chain colMeans vs ICM mode):\n")
    for (k in re_names) {
      fe_mean <- colMeans(batch$fixef[[k]])
      fe_mode <- fixef_mode[[k]]
      for (nm in names(fe_mean)) {
        cat(sprintf(
          "    %-18s  %-30s  mean %12.4f  mode %12.4f  delta %+.4f\n",
          k, nm, fe_mean[[nm]], fe_mode[[nm]], fe_mean[[nm]] - fe_mode[[nm]]
        ))
      }
    }
  } else if (identical(phase, "Block1")) {
    cat("  fixef (unchanged this phase; current chain colMeans vs ICM mode):\n")
    if (!is.null(fixef_mode)) {
      for (k in re_names) {
        fe_mean <- colMeans(batch$fixef[[k]])
        fe_mode <- fixef_mode[[k]]
        for (nm in names(fe_mean)) {
          cat(sprintf(
            "    %-18s  %-30s  mean %12.4f  mode %12.4f\n",
            k, nm, fe_mean[[nm]], fe_mode[[nm]]
          ))
        }
      }
    }
  }

  if (!is.null(b_mode)) {
    b_mean <- apply(batch$b, c(1, 2), mean)
    if (is.null(colnames(b_mean))) colnames(b_mean) <- re_names
    if (is.null(rownames(b_mean))) rownames(b_mean) <- group_levels
    cat(sprintf(
      "  b mean vs mode (%d group x RE coefficients):\n",
      length(b_mean)
    ))
    for (g in rownames(b_mean)) {
      for (k in colnames(b_mean)) {
        bm <- b_mean[g, k]
        md <- b_mode[g, k]
        cat(sprintf(
          "    %s::%-18s  mean %12.4f  mode %12.4f  delta %+.4f\n",
          g, k, bm, md, bm - md
        ))
      }
    }
  }
  cat("\n")
  invisible(NULL)
}
