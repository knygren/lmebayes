#' Small 5-school lmerb fixture for dGamma Block~1 smoke tests.
#'
#' Fixed effects: intercept + distracted_ppvt (required by Prior_Setup for RE slope).
#' Random effects: uncorrelated (Intercept) + distracted_ppvt per school (2 RE columns).
#' Independent Block~1 face grid: 5 schools x gs faces (gs = 9 with Gridtype 3, l1 = 2).
#'
#' @param n_schools Number of schools to keep (default 5).
#' @return List with \code{dat}, \code{form}, \code{design}.
#' @keywords internal
.prepare_small5_lmerb_manual <- function(n_schools = 5L) {
  if (!requireNamespace("bayesrules", quietly = TRUE)) {
    stop("Install bayesrules.", call. = FALSE)
  }
  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(
    dat,
    !is.na(score_ppvt) &
      !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, c(
        "score_ppvt", "distracted_ppvt", "school_id"
      )])
  )

  form <- score_ppvt ~ 1 + distracted_ppvt + (1 + distracted_ppvt || school_id)

  design_scr <- model_setup(form, data = dat)
  full_rank_schools <- names(design_scr$re_rank)[design_scr$re_rank]
  if (length(full_rank_schools) < n_schools) {
    stop(
      "Need at least ", n_schools, " full-rank schools; have ", length(full_rank_schools),
      call. = FALSE
    )
  }

  # Prefer schools with most pupils (stable rank, faster envelope build).
  tab <- sort(table(dat$school_id[dat$school_id %in% full_rank_schools]), decreasing = TRUE)
  keep <- names(tab)[seq_len(n_schools)]
  dat <- subset(dat, school_id %in% keep)
  dat$school_id <- droplevels(dat$school_id)

  design <- model_setup(form, data = dat)
  stopifnot(all(design$re_rank))
  stopifnot(identical(design$re_coef_names, c("(Intercept)", "distracted_ppvt")))
  stopifnot(nlevels(design$groups) == n_schools)

  message(sprintf(
    "Small fixture: %d schools, %d obs, RE = [%s]",
    nlevels(design$groups),
    nrow(dat),
    paste(design$re_coef_names, collapse = ", ")
  ))

  list(dat = dat, form = form, design = design)
}

#' All full-rank schools for the same 2-RE model as \code{.prepare_small5_lmerb_manual}.
#'
#' Same \code{form} and variable structure as the 5-school fixture; keeps every
#' \code{school_id} level that is algebraically full-rank on \code{Z_j} for
#' \code{p_re = 2} (Intercept + distracted_ppvt).
#'
#' @return List with \code{dat}, \code{form}, \code{design}.
#' @keywords internal
.prepare_small5_all_full_rank_manual <- function() {
  if (!requireNamespace("bayesrules", quietly = TRUE)) {
    stop("Install bayesrules.", call. = FALSE)
  }
  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(
    dat,
    !is.na(score_ppvt) &
      !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, c(
        "score_ppvt", "distracted_ppvt", "school_id"
      )])
  )

  form <- score_ppvt ~ 1 + distracted_ppvt + (1 + distracted_ppvt || school_id)

  design_scr <- model_setup(form, data = dat)
  full_rank_schools <- names(design_scr$re_rank)[design_scr$re_rank]
  if (length(full_rank_schools) < 1L) {
    stop("No full-rank schools remain for the 2-RE small model.", call. = FALSE)
  }
  if (length(full_rank_schools) < length(design_scr$re_rank)) {
    dropped <- names(design_scr$re_rank)[!design_scr$re_rank]
    message(
      "All-rank small fixture: dropping ", length(dropped), " rank-deficient ",
      "school_id level(s): ", paste(dropped, collapse = ", ")
    )
    dat <- subset(dat, school_id %in% full_rank_schools)
    dat$school_id <- droplevels(dat$school_id)
  }

  design <- model_setup(form, data = dat)
  stopifnot(all(design$re_rank))
  stopifnot(identical(design$re_coef_names, c("(Intercept)", "distracted_ppvt")))

  n_g <- as.numeric(table(dat$school_id))
  message(sprintf(
    paste0(
      "All-rank small fixture: %d schools, %d obs ",
      "(mean %.1f obs/school, RE = [%s])"
    ),
    nlevels(design$groups),
    nrow(dat),
    mean(n_g),
    paste(design$re_coef_names, collapse = ", ")
  ))

  list(dat = dat, form = form, design = design)
}
