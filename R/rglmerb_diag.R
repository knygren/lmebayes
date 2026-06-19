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
