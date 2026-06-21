#' Check for a single grouping factor in a mixed model formula
#'
#' @param formula A model formula with at most one \code{|} grouping factor.
#' @param data Data frame passed to \code{\link[lme4]{lFormula}}.
#' @param ... Passed to \code{\link[lme4]{lFormula}}.
#' @return Logical: \code{TRUE} if exactly one grouping factor is present.
#' @keywords internal
is_single_factor_model <- function(formula, data, ...) {
  
  # 1. Use the reformulas package to extract random effects bars
  random_bars <- reformulas::findbars(formula)
  
  # If there are no random effects at all, it's a fixed-effects model (0 factors)
  if (is.null(random_bars) || length(random_bars) == 0) {
    return(FALSE)
  }
  
  # 2. Parse the formula with lme4's formula module to evaluate actual grouping execution
  parsed_formula <- tryCatch({
    lme4::lFormula(formula = formula, data = data, ...)
  }, error = function(e) {
    # If it fails to parse structural components, return FALSE
    return(NULL)
  })
  
  if (is.null(parsed_formula)) {
    return(FALSE)
  }
  
  # 3. Extract the names of the grouping factors from flist
  group_names <- names(parsed_formula$reTrms$flist)
  
  # 4. Check if the length of unique grouping variables is EXACTLY 1
  if (length(group_names) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



#' Check whether a formula has no random-effects terms
#'
#' @param formula Model formula.
#' @param data Optional data frame.
#' @param ... Passed to \code{\link[lme4]{lFormula}} when \code{data} is given.
#' @return Logical: \code{TRUE} if no \code{|} random-effects terms are present.
#' @keywords internal
is_fixed_effects_only <- function(formula, data = NULL, ...) {
  
  # Use the dedicated reformulas package to extract random effects bars
  # This avoids the lme4 deprecation warning entirely
  random_bars <- reformulas::findbars(formula)
  
  # Check if the list of random terms is completely empty
  if (is.null(random_bars) || length(random_bars) == 0) {
    return(TRUE)
  }
  
  # Double-check by attempting a lean formula parse via lme4
  parsed_formula <- tryCatch({
    lme4::lFormula(formula = formula, data = data, ...)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(parsed_formula) || length(parsed_formula$reTrms$flist) == 0) {
    return(TRUE)
  }
  
  return(FALSE)
}

#' Validate uncorrelated (diagonal) random effects
#'
#' \pkg{lmebayes} treats \code{Sigma_ranef} as diagonal. Multi-coefficient
#' random terms must use \code{||}; a single random intercept may use
#' \code{(1 | group)} (\code{(1 || group)} is not supported by \code{lme4}).
#'
#' @param formula Mixed-model formula.
#' @param data Optional data frame for \code{\link[lme4]{lFormula}}.
#' @param ... Passed to \code{\link[lme4]{lFormula}}.
#' @return \code{formula} invisibly.
#' @keywords internal
.lmebayes_validate_uncorrelated_re_formula <- function(formula, data = NULL, ...) {
  if (is_fixed_effects_only(formula, data = data, ...)) {
    return(invisible(formula))
  }

  f_chr <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")
  if (grepl("\\(\\s*1\\s*\\|\\|", f_chr)) {
    stop(
      "Intercept-only '(1 || group)' is not supported by lme4. ",
      "Use '(1 | group)' for a random intercept only, or ",
      "'(1 + slope || group)' when adding uncorrelated random slopes.",
      call. = FALSE
    )
  }

  parsed <- tryCatch(
    lme4::lFormula(formula = formula, data = data, ...),
    error = function(e) {
      stop(
        "Could not parse random-effects formula: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  cnms <- parsed$reTrms$cnms
  bad <- which(vapply(cnms, length, integer(1L)) > 1L)
  if (length(bad)) {
    stop(
      "Correlated random effects are not supported (off-diagonal RE covariance). ",
      "Use '||' for uncorrelated terms, e.g. ",
      "(1 + x || group) instead of (1 + x | group).",
      call. = FALSE
    )
  }

  invisible(formula)
}



.lme4_Z_random_column_map <- function(reTrms) {
  Ztlist <- reTrms$Ztlist
  cnms <- reTrms$cnms
  if (length(cnms) != length(Ztlist)) {
    stop("reTrms$cnms and reTrms$Ztlist have different lengths.", call. = FALSE)
  }
  term_names <- names(Ztlist)
  pieces <- vector("list", length(Ztlist))
  col_idx <- 0L
  for (i in seq_along(Ztlist)) {
    Zk <- Ztlist[[i]]
    grp <- sub("^.*\\| (.*)$", "\\1", term_names[i])
    coef_nm <- cnms[[i]]
    levels <- rownames(Zk)
    if (is.null(levels)) {
      levels <- as.character(seq_len(nrow(Zk)))
    }
    n_k <- length(levels)
    pieces[[i]] <- data.frame(
      column_index = (col_idx + 1L):(col_idx + n_k),
      label = paste(grp, coef_nm, levels, sep = "|"),
      group = grp,
      coef = coef_nm,
      level = levels,
      term = term_names[i],
      stringsAsFactors = FALSE
    )
    col_idx <- col_idx + n_k
  }
  do.call(rbind, pieces)
}

.lme4_Z_random_colnames <- function(reTrms) {
  .lme4_Z_random_column_map(reTrms)$label
}

.lme4_Z_random_rownames <- function(reTrms, n_obs) {
  Ztlist <- reTrms$Ztlist
  if (length(Ztlist) < 1L) {
    return(as.character(seq_len(n_obs)))
  }
  rn <- colnames(Ztlist[[1L]])
  if (is.null(rn) || length(rn) != n_obs) {
    return(as.character(seq_len(n_obs)))
  }
  rn
}

.lme4_Z_random_row_map <- function(reTrms, n_obs) {
  rn <- .lme4_Z_random_rownames(reTrms, n_obs)
  data.frame(
    row_index = seq_along(rn),
    observation_row = rn,
    stringsAsFactors = FALSE
  )
}

.lme4_label_Z_random_sparse <- function(Z_random_sparse, reTrms) {
  rn <- .lme4_Z_random_rownames(reTrms, nrow(Z_random_sparse))
  cn <- .lme4_Z_random_colnames(reTrms)
  if (ncol(Z_random_sparse) != length(cn)) {
    stop(
      "Z_random_sparse column count (", ncol(Z_random_sparse),
      ") does not match labelled RE columns (", length(cn), ").",
      call. = FALSE
    )
  }
  dimnames(Z_random_sparse) <- list(rn, cn)
  Z_random_sparse
}

#' Print row and column labels for an lme4 random-effects design matrix
#'
#' Sparse \code{Matrix} objects suppress most dimnames in \code{print()}.
#' Prefer \code{View(lme4_comps$Z_random_column_map)} and
#' \code{View(lme4_comps$Z_random_row_map)} from
#' \code{\link{get_lme4_components}}; use this helper for a quick console summary.
#'
#' @param matrices_list List from \code{\link{get_lme4_components}}.
#' @param max_rows Maximum observation row labels to print.
#' @param max_cols Maximum random-effect column labels to print.
#' @return Invisibly, a list with \code{rownames} and \code{colnames} vectors.
#' @keywords internal
show_lme4_Z_random <- function(
    matrices_list,
    max_rows = 10L,
    max_cols = 15L
) {
  Z <- matrices_list$Z_random_sparse
  if (is.null(Z)) {
    stop("matrices_list has no Z_random_sparse component.", call. = FALSE)
  }
  rn <- matrices_list$Z_random_rownames
  cn <- matrices_list$Z_random_colnames
  if (is.null(rn)) rn <- rownames(Z)
  if (is.null(cn)) cn <- colnames(Z)
  cat("Z_random_sparse:", nrow(Z), "x", ncol(Z), "\n")
  cat("\nRow names (observation / model.frame rows), n =", length(rn), ":\n")
  if (is.null(rn)) {
    cat("  (none)\n")
  } else {
    print(rn[seq_len(min(max_rows, length(rn)))])
    if (length(rn) > max_rows) {
      cat("  ...", length(rn) - max_rows, "more rows\n")
    }
  }
  cat("\nColumn names (group|coef|level), n =", length(cn), ":\n")
  if (is.null(cn)) {
    cat("  (none)\n")
  } else {
    print(cn[seq_len(min(max_cols, length(cn)))])
    if (length(cn) > max_cols) {
      cat("  ...", length(cn) - max_cols, "more columns\n")
    }
  }
  invisible(list(rownames = rn, colnames = cn))
}

#' Extract fixed and random design matrices from an lme4 formula
#'
#' Parses \code{formula} with \code{\link[lme4]{lFormula}} and returns the
#' fixed-effects matrix \code{X}, sparse random-effects matrix \code{Z}, and
#' grouping factors. Row and column label tables are returned explicitly so
#' they are available even when sparse-matrix \code{print()} hides dimnames.
#'
#' @param formula Mixed-model formula understood by \code{lme4}.
#' @param data Data frame.
#' @param ... Passed to \code{\link[lme4]{lFormula}}.
#' @return A list with:
#'   \describe{
#'     \item{X_fixed}{Fixed-effects model matrix.}
#'     \item{Z_random_sparse}{Random-effects model matrix (observations x RE).}
#'     \item{Z_random_rownames}{Character vector of observation/model.frame rows.}
#'     \item{Z_random_colnames}{Character vector \code{group|coef|level}.}
#'     \item{Z_random_row_map}{Data frame mapping row index to observation row.}
#'     \item{Z_random_column_map}{Data frame mapping column index to group/coef/level.}
#'     \item{groups}{Named list of grouping factors (same as \code{lme4} \code{flist}).}
#'     \item{group_names}{Names of grouping factors.}
#'     \item{reTrms}{Raw \code{lme4} random-terms structure.}
#'   }
#' @details
#' Sparse \code{Z_random_sparse} hides most dimnames in \code{print()}. In
#' interactive R, inspect labels with:
#' \preformatted{
#' View(lme4_comps$Z_random_column_map)
#' View(lme4_comps$Z_random_row_map)
#' }
#' @seealso \code{\link{show_lme4_Z_random}}
#' @keywords internal
get_lme4_components <- function(formula, data, ...) {
  # Parse formula into modular lme4 elements
  parsed_formula <- lme4::lFormula(formula = formula, data = data, ...)
  
  # lme4 fixed-effects matrix (contains population terms and Level-2 covariates)
  X_fixed <- parsed_formula$X
  
  # lme4 random-effects sparse design matrix (contains Level-1 repeated measures)
  dev_function <- do.call(lme4::mkLmerDevfun, parsed_formula)
  Z_random_sparse <- Matrix::t(environment(dev_function)$pp$Zt)
  reTrms <- parsed_formula$reTrms
  Z_random_sparse <- .lme4_label_Z_random_sparse(Z_random_sparse, reTrms)
  
  # Factor grouping list
  flist <- reTrms$flist
  group_names <- names(flist)
  
  Z_random_row_map <- .lme4_Z_random_row_map(reTrms, nrow(Z_random_sparse))
  Z_random_column_map <- .lme4_Z_random_column_map(reTrms)

  list(
    X_fixed               = X_fixed,
    Z_random_sparse       = Z_random_sparse,
    Z_random_rownames     = Z_random_row_map$observation_row,
    Z_random_colnames     = Z_random_column_map$label,
    Z_random_row_map      = Z_random_row_map,
    Z_random_column_map   = Z_random_column_map,
    groups                = flist,
    group_names           = group_names,
    reTrms                = reTrms
  )
}




#' Classify \code{lme4} fixed columns as level-1 vs level-2
#'
#' Level-2 columns are constant within each level of \code{group_factor}
#' (including \code{(Intercept)}). Level-1 columns vary within at least one
#' group level.
#'
#' @param X_fixed Fixed-effects model matrix from \code{\link[lme4]{lFormula}}.
#' @param group_factor Grouping factor aligned with rows of \code{X_fixed}.
#' @return A list with \code{level1_cols}, \code{level2_cols}, and named
#'   logical \code{level2_flags}.
#' @keywords internal
classify_lme4_fixed_columns <- function(X_fixed, group_factor) {
  x_colnames <- colnames(X_fixed)
  if (is.null(x_colnames) || length(x_colnames) == 0L) {
    return(list(
      level1_cols = character(0),
      level2_cols = character(0),
      level2_flags = logical(0)
    ))
  }

  level2_flags <- logical(length(x_colnames))
  names(level2_flags) <- x_colnames
  
  for (col in x_colnames) {
    if (col == "(Intercept)") {
      level2_flags[col] <- TRUE
      next
    }
    unique_counts_per_group <- tapply(
      X_fixed[, col], group_factor, function(vec) length(unique(vec))
    )
    level2_flags[col] <- all(unique_counts_per_group == 1L)
  }

  list(
    level1_cols = x_colnames[!level2_flags],
    level2_cols = x_colnames[level2_flags],
    level2_flags = level2_flags
  )
}

#' Classify level-1 fixed interactions as cross-level RE-slope moderation
#'
#' Identifies \code{lme4} fixed columns that vary within groups because they
#' are products of a group-constant predictor and a random-slope coefficient
#' (e.g. \code{private_school:age_c} when \code{age_c} is a random slope).
#' Such terms are allowed in \code{\link{model_setup}} but are not returned in
#' \code{X_fixed}; they are encoded in block priors via
#' \code{mu_b[random_slope] = X_hyper[s, ] \%*\% gamma}.
#'
#' A main fixed effect whose name matches a random slope name (e.g. \code{age_c}
#' alongside \code{(1 + age_c || group)}) is treated as the \emph{population
#' mean slope} (gamma_10) and is also allowed. These are returned in
#' \code{slope_mean_cols}; they do not change the \code{X_hyper} structure but
#' the corresponding \code{fixef()} coefficient is the hyper-model intercept for
#' that slope.
#'
#' @param level1_cols Character vector of level-1 fixed column names.
#' @param level2_cols Character vector of group-constant fixed column names.
#' @param re_coef_names Character vector of random-effects coefficient names
#'   from \code{reTrms$cnms}.
#' @return List with \code{re_slope_moderation} (data frame with columns
#'   \code{interaction_col}, \code{moderator}, \code{random_slope}),
#'   \code{slope_mean_cols} (main fixed effects that are random slope names,
#'   i.e. population mean slope terms), and
#'   \code{disallowed_level1_cols} (level-1 fixed not explained as either).
#' @keywords internal
classify_crosslevel_re_moderation <- function(
    level1_cols,
    level2_cols,
    re_coef_names
) {
  empty_mod <- data.frame(
    interaction_col = character(0),
    moderator = character(0),
    random_slope = character(0),
    stringsAsFactors = FALSE
  )
  if (length(level1_cols) == 0L) {
    return(list(
      re_slope_moderation = empty_mod,
      slope_mean_cols = character(0),
      disallowed_level1_cols = character(0)
    ))
  }

  level2_preds <- setdiff(level2_cols, "(Intercept)")
  re_slopes <- setdiff(unique(as.character(re_coef_names)), "(Intercept)")
  moderation <- vector("list", length(level1_cols))
  slope_means <- character(0)
  disallowed <- character(0)
  n_mod <- 0L

  for (col in level1_cols) {
    # Population mean slope: main effect matches a random slope name
    if (col %in% re_slopes) {
      slope_means <- c(slope_means, col)
      next
    }
    if (!grepl(":", col, fixed = TRUE)) {
      disallowed <- c(disallowed, col)
      next
    }
    parts <- strsplit(col, ":", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      disallowed <- c(disallowed, col)
      next
    }
    p1 <- parts[1L]
    p2 <- parts[2L]
    if (p1 %in% level2_preds && p2 %in% re_slopes) {
      n_mod <- n_mod + 1L
      moderation[[n_mod]] <- data.frame(
        interaction_col = col,
        moderator = p1,
        random_slope = p2,
        stringsAsFactors = FALSE
      )
    } else if (p2 %in% level2_preds && p1 %in% re_slopes) {
      n_mod <- n_mod + 1L
      moderation[[n_mod]] <- data.frame(
        interaction_col = col,
        moderator = p2,
        random_slope = p1,
        stringsAsFactors = FALSE
      )
    } else {
      disallowed <- c(disallowed, col)
    }
  }

  list(
    re_slope_moderation = if (n_mod > 0L) {
      do.call(rbind, moderation[seq_len(n_mod)])
    } else {
      empty_mod
    },
    slope_mean_cols = slope_means,
    disallowed_level1_cols = disallowed
  )
}

#' Hyper covariate matrices per random coefficient
#'
#' For a single-factor \code{lmer} formula, returns one group-level design
#' matrix per random coefficient (intercept and random slopes). These are the
#' \code{X_nbhd} matrices for a coupled block Gibbs hierarchy:
#' \deqn{\beta_{b,j} \mid \gamma_j \sim N(X_{\mathrm{hyper},j}[b,]\gamma_j, \tau_j^2)}
#'
#' Rules (from the \code{lmer} formula):
#' \itemize{
#'   \item \strong{Random intercept:} all group-constant fixed terms from the
#'     formula, with explicit intercept column: e.g.
#'     \eqn{\sim} \code{1 + private_school + title1} (columns
#'     \code{(Intercept)}, \code{private_school}, \code{title1}).
#'   \item \strong{Random slope with cross-level moderation}
#'     (\code{moderator:random_slope} in fixed): e.g.
#'     \code{female:private_school} \eqn{\Rightarrow}
#'     \eqn{\sim} \code{1 + private_school} for the \code{female} slope.
#'   \item \strong{Other random slopes:} \eqn{\sim} \code{1} (column of ones).
#' }
#'
#' @param formula Mixed-model formula understood by \code{lme4}.
#' @param data Data frame.
#' @param ... Passed to \code{\link[lme4]{lFormula}}.
#' @return Object of class \code{"model_setup"} with:
#'   \describe{
#'     \item{\code{group_name}}{Grouping factor name.}
#'     \item{\code{groups}}{Factor of length \code{nrow(Z)} for block subsetting.}
#'     \item{\code{re_coef_names}}{Random coefficient names from \code{lme4}.}
#'     \item{\code{y}}{Response vector, length \code{nrow(Z)} (aligned with
#'       \code{Z} and \code{groups}).}
#'     \item{\code{Z}}{Random-effects model matrix (\code{n_obs} x \code{p_re}):
#'       per-observation loadings on the within-group random coefficient
#'       vector (columns \code{re_coef_names}).}
#'     \item{\code{X_hyper}}{Named list of matrices (one row per group level),
#'       keyed by \code{re_coef_names}.}
#'     \item{\code{re_slope_moderation}}{Cross-level moderation metadata.}
#'   }
#' @seealso \code{\link{model_setup}}, \code{\link{extract_lme4_fixed_group_matrix}},
#'   \code{\link{extract_re_Z_obs}}
#' @keywords internal
extract_re_hyper_matrices <- function(formula, data = NULL, ...) {
  if (!is_single_factor_model(formula, data = data, ...)) {
    stop("extract_re_hyper_matrices() requires exactly one grouping factor.",
         call. = FALSE)
  }
  .lmebayes_validate_uncorrelated_re_formula(formula, data = data, ...)

  parsed <- lme4::lFormula(formula = formula, data = data, ...)
  group_name <- names(parsed$reTrms$flist)[1L]
  group_factor <- parsed$reTrms$flist[[1L]]
  re_coef_names <- unique(unlist(parsed$reTrms$cnms, use.names = FALSE))

  re_slope_moderation <- data.frame(
    interaction_col = character(0),
    moderator = character(0),
    random_slope = character(0),
    stringsAsFactors = FALSE
  )

  X_fixed <- parsed$X
  if (ncol(X_fixed) > 0L) {
    fixed_cls <- classify_lme4_fixed_columns(X_fixed, group_factor)
    if (length(fixed_cls$level1_cols) > 0L) {
      cross_cls <- classify_crosslevel_re_moderation(
        level1_cols = fixed_cls$level1_cols,
        level2_cols = fixed_cls$level2_cols,
        re_coef_names = re_coef_names
      )
      re_slope_moderation <- cross_cls$re_slope_moderation
      # slope_mean_cols (e.g. age_c alongside (1 + age_c || group)) are the
      # population mean slopes gamma_10; they are level-2 parameters estimated
      # as fixed effects and are allowed without changing X_hyper structure.
      if (length(cross_cls$disallowed_level1_cols) > 0L) {
        stop(
          "Fixed effects must be constant within ", group_name,
          " (level-2), a population mean slope (matching a random slope name),",
          " or cross-level RE moderation (moderator:random_slope); ",
          "not allowed: ",
          paste(cross_cls$disallowed_level1_cols, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  comps <- get_lme4_components(formula = formula, data = data, ...)
  X_group <- extract_lme4_fixed_group_matrix(comps, group_name)
  k <- nrow(X_group)

  moderators <- stats::setNames(
    re_slope_moderation$moderator,
    re_slope_moderation$random_slope
  )

  X_hyper <- vector("list", length(re_coef_names))
  names(X_hyper) <- re_coef_names

  for (coef in re_coef_names) {
    if (coef == "(Intercept)") {
      X_hyper[[coef]] <- X_group
    } else if (coef %in% names(moderators)) {
      mod <- moderators[[coef]]
      if (!mod %in% colnames(X_group)) {
        stop(
          "Moderator '", mod, "' for random slope '", coef,
          "' is not in the group-constant fixed matrix.",
          call. = FALSE
        )
      }
      X_hyper[[coef]] <- X_group[, c("(Intercept)", mod), drop = FALSE]
    } else {
      X_hyper[[coef]] <- matrix(
        1,
        nrow = k,
        ncol = 1L,
        dimnames = list(rownames(X_group), "(Intercept)")
      )
    }
  }

  Z <- extract_re_Z_obs(
    matrices_list = comps,
    group_name = group_name,
    re_coef_names = re_coef_names
  )

  y <- stats::model.response(parsed$fr)
  if (is.matrix(y) && ncol(y) == 1L) {
    y <- y[, 1L]
  }
  if (length(y) != nrow(Z)) {
    stop(
      "Length of y (", length(y), ") does not match nrow(Z) (", nrow(Z), ").",
      call. = FALSE
    )
  }

  structure(
    list(
      group_name = group_name,
      groups = group_factor,
      re_coef_names = re_coef_names,
      y = y,
      Z = Z,
      X_hyper = X_hyper,
      re_slope_moderation = re_slope_moderation
    ),
    class = "model_setup"
  )
}

#' Build default \code{lmer} formula for RE variance calibration
#'
#' Cross-level RE moderation (e.g. \code{female:free_reduced_lunch}) belongs
#' in the Bayesian hyper prior, not as a student-varying fixed effect in the
#' \code{lmer} fit used to calibrate random-effect variances. This helper
#' returns a formula with:
#' \itemize{
#'   \item group-constant fixed terms only (level-2 hyper covariates);
#'   \item the same uncorrelated random-effects structure as \code{formula}
#'     (\code{||} notation), without cross-level fixed interactions.
#' }
#'
#' @inheritParams extract_re_hyper_matrices
#' @return Formula suitable for \code{\link[lme4]{lmer}} variance calibration.
#' @keywords internal
lmerb_default_vcov_formula <- function(formula, data = NULL, ...) {
  if (!is_single_factor_model(formula, data = data, ...)) {
    stop("lmerb_default_vcov_formula() requires exactly one grouping factor.",
         call. = FALSE)
  }

  parsed <- lme4::lFormula(formula = formula, data = data, ...)
  group_name <- names(parsed$reTrms$flist)[1L]
  group_factor <- parsed$reTrms$flist[[1L]]
  re_coef_names <- unique(unlist(parsed$reTrms$cnms, use.names = FALSE))

  fixed_cls <- classify_lme4_fixed_columns(parsed$X, group_factor)
  cross_cls <- classify_crosslevel_re_moderation(
    level1_cols = fixed_cls$level1_cols,
    level2_cols = fixed_cls$level2_cols,
    re_coef_names = re_coef_names
  )

  fixed_form <- reformulas::nobars(formula)
  term_labels <- attr(terms(fixed_form), "term.labels")
  # Drop only the cross-level interaction terms; keep slope_mean_cols (gamma_10)
  drop_terms <- cross_cls$re_slope_moderation$interaction_col
  terms_keep <- setdiff(term_labels, drop_terms)

  resp <- all.vars(formula)[1L]
  fixed_rhs <- if (length(terms_keep) == 0L) {
    "1"
  } else {
    paste(terms_keep, collapse = " + ")
  }

  # Build random term from lFormula cnms (findbars() splits || into separate bars).
  re_terms <- ifelse(re_coef_names == "(Intercept)", "1", re_coef_names)
  if (length(re_coef_names) == 1L && identical(re_coef_names, "(Intercept)")) {
    random_chr <- paste0("(1 | ", group_name, ")")
  } else {
    random_chr <- paste0(
      "(", paste(re_terms, collapse = " + "), " || ", group_name, ")"
    )
  }

  stats::as.formula(paste(resp, "~", fixed_rhs, "+", random_chr))
}

#' Extract variance components from an \code{lmer} fit
#'
#' Maps \code{lme4::VarCorr()} output to a named vector of random-effect
#' variances (aligned with \code{re_coef_names}) plus residual variance.
#'
#' @param fit Object of class \code{"lmerMod"} from \code{\link[lme4]{lmer}}.
#' @param re_coef_names Random coefficient names (as in \code{reTrms$cnms}).
#' @return List with \code{varcorr}, \code{vcov_re}, and \code{residual_var}.
#' @keywords internal
extract_lmer_variance_components <- function(fit, re_coef_names) {
  if (!inherits(fit, "lmerMod")) {
    stop("fit must be an lmerMod object.", call. = FALSE)
  }
  extract_mer_variance_components(fit, re_coef_names)
}

#' Extract variance components from an \code{lmer} or \code{glmer} fit
#'
#' @param fit Object of class \code{"merMod"}.
#' @param re_coef_names Random coefficient names (as in \code{reTrms$cnms}).
#' @return List with \code{varcorr}, \code{vcov_re}, and \code{residual_var}.
#' @keywords internal
extract_mer_variance_components <- function(fit, re_coef_names) {
  if (!inherits(fit, "merMod")) {
    stop("fit must be an merMod object.", call. = FALSE)
  }

  vc_df <- as.data.frame(lme4::VarCorr(fit), stringsAsFactors = FALSE)
  vcov_re <- stats::setNames(
    rep(NA_real_, length(re_coef_names)),
    re_coef_names
  )

  for (nm in re_coef_names) {
    hit <- vc_df$var1 == nm & (vc_df$var2 == nm | is.na(vc_df$var2))
    hit[is.na(hit)] <- FALSE
    if (any(hit)) {
      vcov_re[[nm]] <- vc_df$vcov[which(hit)[1L]]
    }
  }

  if (anyNA(vcov_re)) {
    missing_coefs <- names(vcov_re)[is.na(vcov_re)]
    stop(
      "Could not find variance components for: ",
      paste(missing_coefs, collapse = ", "),
      call. = FALSE
    )
  }

  residual_var <- vc_df$vcov[vc_df$grp == "Residual"]
  if (length(residual_var) >= 1L) {
    residual_var <- unname(residual_var[1L])
  } else if (inherits(fit, "lmerMod")) {
    residual_var <- lme4::getME(fit, "sigma")^2
  } else {
    residual_var <- NA_real_
  }

  list(
    varcorr = vc_df,
    vcov_re = vcov_re,
    residual_var = unname(residual_var)
  )
}

#' Reference \code{lmer}/\code{glmer} fit embedded in an \code{lmerb}/\code{glmerb} object
#' @keywords internal
.lmerb_reference_fit <- function(object) {
  if (inherits(object, "glmerb")) {
    if (is.null(object$glmer)) {
      stop("glmerb fit is missing component 'glmer'.", call. = FALSE)
    }
    return(object$glmer)
  }
  if (is.null(object$lmer)) {
    stop("lmerb fit is missing component 'lmer'.", call. = FALSE)
  }
  object$lmer
}

#' Validate observation-level \code{dispersion_ranef} for an \code{lmerb}/\code{glmerb} family
#'
#' Returns a scalar plug-in \eqn{\sigma^2} or \code{NULL}.  For routing that
#' distinguishes fixed vs \code{dGamma()} priors, use
#' \code{\link{.lmebayes_resolve_dispersion_ranef}}.
#' @keywords internal
.lmebayes_validate_dispersion_ranef <- function(
    dispersion_ranef,
    family,
    fn_name = "lmerb"
) {
  resolved <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = NULL,
    fn_name          = fn_name
  )
  resolved$dispersion_fix
}

#' Resolve observation-level \code{dispersion_ranef} (fixed scalar or \code{dGamma()})
#' @return List with \code{mode} (\code{"none"}, \code{"fixed"}, or \code{"gamma"}),
#'   \code{dispersion_fix} (plug-in \eqn{\sigma^2}), \code{dispersion_prior_list}
#'   (\code{dGamma()} \code{prior_list} or \code{NULL}), and \code{dispersion_pfamily}.
#' @keywords internal
.lmebayes_resolve_dispersion_ranef <- function(
    dispersion_ranef,
    family,
    design = NULL,
    fn_name = "lmerb"
) {
  has_dispersion <- family$family %in%
    c("gaussian", "Gamma", "quasipoisson", "quasibinomial")

  if (!has_dispersion) {
    if (!is.null(dispersion_ranef)) {
      stop(
        "'dispersion_ranef' must be NULL for family = ", family$family,
        "() (no observation-level dispersion).",
        call. = FALSE
      )
    }
    return(list(
      mode                   = "none",
      dispersion_fix         = NULL,
      dispersion_prior_list  = NULL,
      dispersion_pfamily     = NULL
    ))
  }

  if (inherits(dispersion_ranef, "pfamily")) {
    if (!identical(dispersion_ranef$pfamily, "dGamma")) {
      stop(
        fn_name, "(): 'dispersion_ranef' pfamily must be dGamma(); got ",
        dispersion_ranef$pfamily, ". RE priors belong in 'pfamily_list'.",
        call. = FALSE
      )
    }
    pl <- dispersion_ranef$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): dGamma() observation-dispersion prior requires ",
        "Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    if (is.null(design) || is.null(design$residual_var)) {
      stop(
        fn_name, "(): a model_setup with residual_var is required for ",
        "dGamma() dispersion_ranef (plug-in sigma^2).",
        call. = FALSE
      )
    }
    return(list(
      mode                  = "gamma",
      dispersion_fix        = as.numeric(design$residual_var),
      dispersion_prior_list = pl,
      dispersion_pfamily    = dispersion_ranef
    ))
  }

  if (is.null(dispersion_ranef) || !is.numeric(dispersion_ranef) ||
      length(dispersion_ranef) != 1L || !is.finite(dispersion_ranef) ||
      dispersion_ranef <= 0) {
    stop(
      "'dispersion_ranef' must be a single positive number or a dGamma() ",
      "pfamily for family = ", family$family, "().",
      call. = FALSE
    )
  }
  list(
    mode                  = "fixed",
    dispersion_fix        = as.numeric(dispersion_ranef),
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' Call \code{rLMMNormal_reg} or \code{rLMMindepNormalGamma_reg} from matrix-level inputs
#' @keywords internal
.lmebayes_run_lmm_engine <- function(
    n,
    design,
    prior,
    disp_info,
    fixef_start   = NULL,
    m_convergence = NULL,
    tv_tol        = 0.01,
    seed          = NULL,
    progbar       = TRUE,
    verbose       = FALSE
) {
  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)
  P            <- solve(prior$Sigma_ranef)
  common_args  <- list(
    n             = n,
    y             = design$y,
    x             = design$Z,
    block         = design$groups,
    x_hyper       = design$X_hyper,
    P             = P,
    pfamily_list  = prior$pfamily_list,
    start         = fixef_start,
    m_convergence = m_convergence,
    tv_tol        = tv_tol,
    re_coef_names = re_names,
    group_levels  = group_levels,
    group_name    = design$group_name,
    seed          = seed,
    progbar       = progbar,
    verbose       = verbose
  )
  if (identical(disp_info$mode, "gamma")) {
    do.call(
      glmbayesCore::rLMMindepNormalGamma_reg,
      c(
        common_args,
        list(
          prior_list     = disp_info$dispersion_prior_list,
          dispersion_fix = disp_info$dispersion_fix
        )
      )
    )
  } else {
    do.call(
      glmbayesCore::rLMMNormal_reg,
      c(
        common_args,
        list(prior_list = list(dispersion = disp_info$dispersion_fix))
      )
    )
  }
}

#' Build Block~1 prior list from a normalized prior container
#' @param measurement_prior_list Prior container with \code{Sigma_ranef}.
#' @param dispersion_ranef Optional observation-level dispersion override
#'   (\eqn{\sigma^2}); when supplied, used instead of
#'   \code{measurement_prior_list$dispersion_ranef}.
#' @keywords internal
.lmebayes_block1_prior_list <- function(
    measurement_prior_list,
    dispersion_ranef = NULL
) {
  if (is.null(measurement_prior_list$Sigma_ranef)) {
    stop("measurement_prior_list must contain 'Sigma_ranef'.", call. = FALSE)
  }
  P <- solve(measurement_prior_list$Sigma_ranef)
  dispersion <- if (!is.null(dispersion_ranef)) {
    dispersion_ranef
  } else {
    measurement_prior_list$dispersion_ranef
  }
  if (is.null(dispersion)) {
    list(P = P, ddef = TRUE)
  } else {
    list(P = P, dispersion = dispersion, ddef = FALSE)
  }
}

#' Stage v2 sampler output to the \code{fixef.*} namespace
#' @keywords internal
.lmebayes_stage_v2_fixef <- function(
    out,
    fixef_mode,
    fixef_init,
    re_names,
    group_levels,
    n
) {
  x <- list(
    fixef_draws            = out$fixef_draws,
    coefficients           = out$coefficients,
    dispersion_fixef_draws = out$dispersion_fixef_draws,
    iters_fixef_draws      = out$iters_fixef_draws,
    mu_all_last            = out$mu_all_last,
    re_coef_names          = re_names,
    group_levels           = group_levels,
    n                      = n
  )
  glmbayesCore:::.two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}

#' Add \code{fixef.means} and related summary fields to a staged sampler object
#' @keywords internal
.lmebayes_add_fixef_summaries <- function(x) {
  if (!is.null(x$fixef)) {
    x$fixef.means <- lapply(x$fixef, colMeans)
  }
  if (!is.null(x$fixef.dispersion)) {
    x$fixef.dispersion.mean <- colMeans(x$fixef.dispersion)
  }
  if (!is.null(x$fixef.iters) && !is.null(x$m_convergence)) {
    x$fixef.iters.mean <- colMeans(x$fixef.iters) / x$m_convergence
  }
  x
}

#' Normalize a pfamily list + dispersion into the internal prior container
#'
#' Validates the \code{pfamily_list} / \code{dispersion_ranef} arguments of
#' \code{\link{lmerb}} and \code{\link{glmerb}} against the model design and
#' converts them into the internal prior container consumed downstream
#' (\code{Sigma_ranef}, \code{dispersion_ranef}, per-component
#' \code{prior_list} with \code{mu_fixef} / \code{Sigma_fixef} /
#' \code{dispersion_fixef}).
#'
#' The Block~1 random-effect covariance is reconstructed from the Block~2
#' dispersions: \code{Sigma_ranef = diag(tau^2_k)} where \code{tau^2_k} is the
#' \code{dispersion} of component \code{k}'s \code{dNormal} pfamily.
#'
#' \code{dIndependent_Normal_Gamma} components must carry \emph{both}
#' truncation bounds: \code{disp_lower} doubles as the plug-in
#' \eqn{\tau^2_k} for the eigenvalue / TV calibration (smaller \eqn{\tau^2}
#' increases the coupling between blocks and hence the contraction rate
#' \eqn{\lambda^*}, so the lower bound yields a conservative eigenvalue and
#' sweep count), and supplying \code{disp_upper} as well fixes the
#' \eqn{\tau^2_k} truncation window \code{[disp_lower, disp_upper]} across
#' all inner Gibbs sweeps of \code{two_block_rNormal_reg_v2}, making the
#' calibration valid over the chain's entire dispersion support.
#'
#' @param pfamily_list Named list of \code{"pfamily"} objects, one per
#'   random-effect coefficient (e.g. from
#'   \code{\link{pfamily_list.lmebayes_prior_setup}}).  Names must match
#'   \code{design$re_coef_names} (any order).
#' @param dispersion_ranef Observation-level dispersion for the measurement
#'   model.  Required positive scalar for \code{gaussian()}; must be
#'   \code{NULL} for families without a dispersion parameter.
#' @param design A \code{\link{model_setup}} object.
#' @param family A \code{\link[stats]{family}} object.
#' @param fn_name Calling function name used in error messages.
#' @return List with \code{pfamily_list} (reordered), \code{dispersion_ranef},
#'   \code{Sigma_ranef}, \code{prior_list}, \code{ptypes} (per-component
#'   constructor names), and \code{any_non_normal}.
#' @keywords internal
.lmebayes_priors_from_pfamily_list <- function(pfamily_list,
                                               dispersion_ranef,
                                               design,
                                               family,
                                               fn_name = "lmerb") {

  re_names <- design$re_coef_names
  p_re     <- length(re_names)

  ## --- dispersion_ranef (Block 1 measurement dispersion) -------------------
  disp_res <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = design,
    fn_name          = fn_name
  )
  dispersion_ranef <- disp_res$dispersion_fix

  ## --- pfamily_list ---------------------------------------------------------
  if (!is.list(pfamily_list) || length(pfamily_list) != p_re) {
    stop(
      "'pfamily_list' must be a list with one pfamily per random-effect ",
      "component (", p_re, " expected: ", paste(re_names, collapse = ", "),
      "). Build it with pfamily_list(Prior_Setup_lmebayes(...)).",
      call. = FALSE
    )
  }
  if (is.null(names(pfamily_list)) || !setequal(names(pfamily_list), re_names)) {
    stop(
      "Names of 'pfamily_list' must match the random-effect coefficient ",
      "names: ", paste(re_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  pfamily_list <- pfamily_list[re_names]

  prior_list <- stats::setNames(vector("list", p_re), re_names)
  tau2   <- stats::setNames(numeric(p_re), re_names)
  ptypes <- stats::setNames(character(p_re), re_names)

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily")) {
      stop("pfamily_list[[\"", k, "\"]] must be a pfamily object.",
           call. = FALSE)
    }
    if (!pf$pfamily %in% c("dNormal", "dIndependent_Normal_Gamma")) {
      stop(
        fn_name, "() supports only dNormal and dIndependent_Normal_Gamma ",
        "pfamilies in 'pfamily_list'; component \"", k, "\" is ",
        pf$pfamily, ".",
        call. = FALSE
      )
    }
    ptypes[[k]] <- pf$pfamily

    par_names <- colnames(design$X_hyper[[k]])
    q_k <- length(par_names)

    mu_k <- as.numeric(pf$prior_list$mu)
    if (length(mu_k) != q_k) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu has length ",
        length(mu_k), " but the hyper design has ", q_k, " column(s): ",
        paste(par_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pf$prior_list$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k, "\"]] (",
          paste(mu_nms, collapse = ", "), ") do not match the hyper design ",
          "columns (", paste(par_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pf$prior_list$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pf$prior_list$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    ## Keep the pfamily object itself aligned with the hyper-design column
    ## order: it is passed straight to the v2 sampler as the Block 2 source
    ## of truth, so its mu/Sigma must match x_hyper[[k]].
    pfamily_list[[k]]$prior_list$mu <-
      matrix(mu_k, ncol = 1L, dimnames = list(par_names, NULL))
    pfamily_list[[k]]$prior_list$Sigma <- Sigma_k

    if (identical(pf$pfamily, "dNormal")) {
      d_k <- pf$prior_list$dispersion
      if (isTRUE(pf$prior_list$ddef)) {
        warning(
          fn_name, ": pfamily_list[[\"", k, "\"]] uses the default ",
          "dispersion = 1 (none was supplied to dNormal()); the Block 1 ",
          "random-effect variance tau^2 for \"", k, "\" is therefore 1.",
          call. = FALSE
        )
      }
    } else {
      ## ING: both truncation bounds are required.  disp_lower doubles as
      ## the conservative tau^2 plug-in for the eigenvalue / TV calibration;
      ## together the bounds fix the tau^2_k truncation window across all
      ## inner Gibbs sweeps (one-sided specifications would fall back to a
      ## per-sweep surrogate-posterior window inside the envelope code).
      d_k <- pf$prior_list$disp_lower
      if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
          !is.finite(d_k) || d_k <= 0) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a positive scalar ",
          "'disp_lower' (lower dispersion truncation). It is used as the ",
          "conservative tau^2 plug-in for the convergence calibration.",
          call. = FALSE
        )
      }
      u_k <- pf$prior_list$disp_upper
      if (is.null(u_k) || !is.numeric(u_k) || length(u_k) != 1L ||
          !is.finite(u_k) || u_k <= as.numeric(d_k)) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a finite scalar ",
          "'disp_upper' > 'disp_lower' (upper dispersion truncation), so ",
          "the tau^2 truncation window is fixed across Gibbs sweeps. ",
          "pfamily_list(Prior_Setup_lmebayes(...)) sets both bounds to the ",
          "0.01/0.99 prior dispersion quantiles by default.",
          call. = FALSE
        )
      }
    }

    tau2[[k]] <- as.numeric(d_k)
    prior_list[[k]] <- list(
      mu_fixef         = mu_k,
      Sigma_fixef      = Sigma_k,
      dispersion_fixef = as.numeric(d_k)
    )
  }

  Sigma_ranef <- diag(unname(tau2), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(re_names, re_names)

  list(
    pfamily_list          = pfamily_list,
    dispersion_ranef      = dispersion_ranef,
    dispersion_mode       = disp_res$mode,
    dispersion_pfamily    = disp_res$dispersion_pfamily,
    dispersion_prior_list = disp_res$dispersion_prior_list,
    Sigma_ranef           = Sigma_ranef,
    prior_list            = prior_list,
    ptypes         = ptypes,
    any_non_normal = any(ptypes != "dNormal")
  )
}

#' Restructure \code{lme4} sparse \code{Z} to per-observation RE loadings
#'
#' Collapses the stacked \code{lme4} random-effects design (one column per
#' \code{group|coef|level}) into an \code{n_obs x p_re} matrix whose columns
#' are the random coefficient names (\code{(Intercept)}, slopes, etc.). Row
#' \code{i} contains the loadings on the within-group random vector for
#' observation \code{i}'s group.
#'
#' @param matrices_list List from \code{\link{get_lme4_components}}.
#' @param group_name Name of the grouping factor.
#' @param re_coef_names Character vector of random coefficient names
#'   (as in \code{reTrms$cnms}).
#' @return Numeric matrix \code{n_obs x length(re_coef_names)} with
#'   \code{colnames = re_coef_names}.
#' @keywords internal
extract_re_Z_obs <- function(matrices_list, group_name, re_coef_names) {
  Z_sparse <- matrices_list$Z_random_sparse
  col_map <- matrices_list$Z_random_column_map
  group_factor <- matrices_list$groups[[group_name]]

  if (is.null(Z_sparse)) {
    stop("matrices_list has no Z_random_sparse component.", call. = FALSE)
  }
  if (is.null(group_factor)) {
    stop(sprintf("Grouping factor '%s' not found.", group_name), call. = FALSE)
  }
  if (length(re_coef_names) < 1L) {
    stop("re_coef_names must be non-empty.", call. = FALSE)
  }

  n <- nrow(Z_sparse)
  if (length(group_factor) != n) {
    stop(
      "Length of groups (", length(group_factor),
      ") does not match nrow(Z) (", n, ").",
      call. = FALSE
    )
  }

  p <- length(re_coef_names)
  Z <- matrix(0, nrow = n, ncol = p)
  colnames(Z) <- re_coef_names

  g_chr <- as.character(group_factor)
  Zm <- as.matrix(Z_sparse)

  for (j in seq_along(re_coef_names)) {
    coef_nm <- re_coef_names[j]
    map_j <- col_map[col_map$coef == coef_nm, , drop = FALSE]
    if (nrow(map_j) < 1L) {
      stop(
        "No Z columns found for random coefficient '", coef_nm, "'.",
        call. = FALSE
      )
    }
    for (lev in unique(map_j$level)) {
      rows <- which(g_chr == lev)
      if (length(rows) == 0L) {
        next
      }
      col_idx <- map_j$column_index[map_j$level == lev]
      if (length(col_idx) != 1L) {
        stop(
          "Expected one Z column for coef '", coef_nm, "' at level '", lev,
          "', found ", length(col_idx), ".",
          call. = FALSE
        )
      }
      Z[rows, j] <- Zm[rows, col_idx]
    }
  }

  Z
}

#' Extract per-group submatrices from lme4 design components
#'
#' @param matrices_list List from \code{\link{get_lme4_components}}.
#' @param group_name Name of the grouping factor (e.g. \code{"school_id"}).
#' @param target_level Level of \code{group_name} to slice.
#' @return List with per-group \code{X_fixed_level1}, \code{Z_random_level1},
#'   and school-constant \code{X_fixed_level2}.
#' @keywords internal
extract_lme4_submatrices <- function(matrices_list, group_name, target_level) {
  X_fixed <- matrices_list$X_fixed
  Z_random_sparse <- matrices_list$Z_random_sparse
  group_factor <- matrices_list$groups[[group_name]]
  
  # 1. Random-effects columns are shared across levels; slice rows by group below.
  Z_level_sparse <- Z_random_sparse

  fixed_cls <- classify_lme4_fixed_columns(X_fixed, group_factor)
  level1_fixed_cols <- fixed_cls$level1_cols
  level2_fixed_cols <- fixed_cls$level2_cols
  
  # 3. Filter rows for target level
  group_rows <- which(group_factor == target_level)
  
  # Output structures reflect standard lme4 matrix groupings
  return(list(
    group_variable     = group_name,
    group_level        = target_level,
    n_observations     = length(group_rows),
    X_fixed_level1     = X_fixed[group_rows, level1_fixed_cols, drop = FALSE],
    Z_random_level1    = {
      Zb <- Z_level_sparse[group_rows, , drop = FALSE]
      Zm <- as.matrix(Zb)
      rownames(Zm) <- rownames(Zb)
      colnames(Zm) <- colnames(Zb)
      Zm
    },
    X_fixed_level2     = X_fixed[group_rows[1], level2_fixed_cols, drop = FALSE] # Single-row vector
  ))
}


#' Extract the group-constant fixed-effects matrix (level-2 X)
#'
#' @param matrices_list List from \code{\link{get_lme4_components}}.
#' @param group_name Name of the grouping factor.
#' @return Matrix with one row per group level and school-constant fixed columns.
#' @keywords internal
extract_lme4_fixed_group_matrix <- function(matrices_list, group_name) {
  X_fixed <- matrices_list$X_fixed
  group_factor <- matrices_list$groups[[group_name]]
  
  if (is.null(group_factor)) {
    stop(sprintf("Grouping factor '%s' not found.", group_name))
  }
  
  all_levels <- levels(group_factor)
  fixed_cls <- classify_lme4_fixed_columns(X_fixed, group_factor)
  level2_fixed_cols <- fixed_cls$level2_cols
  first_row_indices <- match(all_levels, group_factor)
  
  # Extract the row-per-group matrix from lme4's X
  X_group_level_matrix <- X_fixed[first_row_indices, level2_fixed_cols, drop = FALSE]
  rownames(X_group_level_matrix) <- all_levels
  
  return(X_group_level_matrix)
}
