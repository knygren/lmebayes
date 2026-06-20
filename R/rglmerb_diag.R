#' Print Block~2 reference vs ICM fixed effects
#' @noRd
.lmebayes_print_icm_fixef_table <- function(
    prior_list,
    re_names,
    fixef_icm,
    icm_info,
    ref_label,
    icm_label,
    header,
    verbose
) {
  if (!isTRUE(verbose) || is.null(fixef_icm)) {
    return(invisible(NULL))
  }
  fixef_ref <- lapply(prior_list, `[[`, "mu_fixef")
  names(fixef_ref) <- re_names
  hdr <- sprintf("  %-18s  %-30s  %12s  %12s",
                 "RE component", "parameter", ref_label, icm_label)
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(header, "\n")
  cat(hdr, "\n")
  cat(sep, "\n")
  for (k in re_names) {
    nms_k  <- names(fixef_ref[[k]])
    ref_v  <- fixef_ref[[k]]
    icm_v  <- fixef_icm[[k]]
    for (nm in nms_k) {
      cat(sprintf("  %-18s  %-30s  %12.4f  %12.4f\n",
                  k, nm, ref_v[[nm]], icm_v[[nm]]))
    }
  }
  if (!is.null(icm_info)) {
    cat(sprintf("  (ICM converged: %s, %d iter, delta = %.2e)\n\n",
                icm_info$converged, icm_info$iterations, icm_info$delta))
  } else {
    cat("\n")
  }
  invisible(NULL)
}

#' Print ICM random-effects mode (all group x RE coefficients)
#' @noRd
.lmebayes_print_ranef_mode_reference <- function(
    ranef_mode,
    re_names,
    group_levels,
    verbose
) {
  invisible(NULL)

  # Disabled: long per-group Block 1 mode listing at fit startup.
  # if (!isTRUE(verbose) || is.null(ranef_mode)) {
  #   return(invisible(NULL))
  # }
  # if (is.null(colnames(ranef_mode))) {
  #   colnames(ranef_mode) <- re_names
  # }
  # if (is.null(rownames(ranef_mode))) {
  #   rownames(ranef_mode) <- group_levels
  # }
  # cat("--- glmerb: ICM mode reference (Block 1 random effects) ---\n")
  # cat(sprintf(
  #   "  %d group x RE coefficients (%d groups x %d RE cols):\n",
  #   length(ranef_mode), nrow(ranef_mode), ncol(ranef_mode)
  # ))
  # for (g in rownames(ranef_mode)) {
  #   for (k in colnames(ranef_mode)) {
  #     cat(sprintf("    %s::%-18s  %12.4f\n", g, k, ranef_mode[g, k]))
  #   }
  # }
  # cat("\n")
  # invisible(NULL)
}

#' Print \code{fixef.init} (pilot colMeans) before main stage
#' @noRd
.lmebayes_print_fixef_init <- function(
    fixef_init,
    re_names,
    verbose
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat("--- glmerb: main-stage fixef.init (pilot colMeans) ---\n")
  for (k in re_names) {
    for (nm in names(fixef_init[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_init[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}
