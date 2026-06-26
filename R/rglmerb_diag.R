#' Column labels for Block~2 prior vs ICM comparison tables
#' @noRd
.lmebayes_block2_icm_labels <- function(prior, family = gaussian()) {
  any_ing <- isTRUE(prior$any_non_normal)
  is_gauss <- is.null(family) || identical(family$family, "gaussian")
  ref_label <- "prior mean"
  if (any_ing) {
    icm_label   <- "gamma @ lmer tau2"
    icm_verbose <- "Block 2 start at lmer tau^2 plug-in"
    conv_label  <- "Plug-in fixed point"
  } else if (is_gauss) {
    icm_label   <- "ICM mean"
    icm_verbose <- "ICM posterior mean"
    conv_label  <- "ICM"
  } else {
    icm_label   <- "ICM mode"
    icm_verbose <- "ICM posterior mode"
    conv_label  <- "ICM"
  }
  list(
    ref_label   = ref_label,
    icm_label   = icm_label,
    icm_verbose = icm_verbose,
    conv_label  = conv_label
  )
}

#' Print Block~2 reference vs ICM fixed effects
#' @noRd
.lmebayes_print_icm_fixef_table <- function(
    prior_list,
    re_names,
    fixef_icm,
    icm_info,
    ref_label,
    icm_label,
    conv_label = "ICM",
    header,
    verbose
) {
  if (!isTRUE(verbose) || is.null(fixef_icm)) {
    return(invisible(NULL))
  }
  fixef_ref <- lapply(prior_list, `[[`, "mu_fixef")
  names(fixef_ref) <- re_names
  hdr <- sprintf("  %-18s  %-30s  %14s  %18s",
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
      cat(sprintf("  %-18s  %-30s  %14.4f  %18.4f\n",
                  k, nm, ref_v[[nm]], icm_v[[nm]]))
    }
  }
  if (!is.null(icm_info)) {
    cat(sprintf("  (%s converged: %s, %d iter, delta = %.2e)\n\n",
                conv_label,
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
    verbose,
    header = "--- main-stage fixef.init (pilot colMeans) ---"
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(header, "\n")
  for (k in re_names) {
    for (nm in names(fixef_init[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_init[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}
