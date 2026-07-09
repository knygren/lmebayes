#' Shared Big Word Club model for manual lmerb validation scripts.
#'
#' Same \code{dat}, \code{form}, and school levels as
#' \code{test_lmerb_mer_re_validation.R} and
#' \code{test_lmerb_dgamma_mer_re_validation.R} (only the dispersion /
#' Block~2 prior route differs across those scripts).
#'
#' After the usual NA / complete-case filter, drops \code{school_id} levels that
#' are not algebraically full-rank on \code{Z_j} (\code{model_setup()$re_rank})
#' before \code{\link{Prior_Setup_lmebayes}} — same rule as
#' \code{demo("Ex_24_lmerb_dGamma_BigWordClub")}.
#'
#' @param pwt Prior weight for Block~2 setup.
#' @param pwt_dispersion Optional prior weight for ING Block~2 calibration.
#' @return List with \code{dat}, \code{form}, \code{design}, \code{ps}.
#' @keywords internal
.prepare_bwc_lmerb_manual <- function(
    pwt = 0.01,
    pwt_dispersion = NULL
) {
  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(
    dat,
    !is.na(score_ppvt) &
      !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, c(
        "score_ppvt", "distracted_a1", "distracted_ppvt",
        "private_school", "title1", "free_reduced_lunch", "school_id"
      )])
  )

  form <- score_ppvt ~
    private_school + title1 + free_reduced_lunch +
    distracted_ppvt + distracted_a1 +
    free_reduced_lunch:distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)

  design_scr <- model_setup(form, data = dat)
  full_rank_schools <- names(design_scr$re_rank)[design_scr$re_rank]
  if (length(full_rank_schools) < length(design_scr$re_rank)) {
    dropped <- names(design_scr$re_rank)[!design_scr$re_rank]
    message(
      "BWC manual fixture: dropping ", length(dropped), " rank-deficient ",
      "school_id level(s): ", paste(dropped, collapse = ", ")
    )
    if (length(full_rank_schools) < 1L) {
      stop("No full-rank schools remain after filtering.", call. = FALSE)
    }
    dat <- subset(dat, school_id %in% full_rank_schools)
    dat$school_id <- droplevels(dat$school_id)
    design <- model_setup(form, data = dat)
  } else {
    design <- design_scr
  }
  stopifnot(all(design$re_rank))

  ps_args <- list(form = form, data = dat, pwt = pwt)
  if (!is.null(pwt_dispersion)) {
    ps_args$pwt_dispersion <- pwt_dispersion
  }
  ps <- do.call(Prior_Setup_lmebayes, ps_args)

  list(dat = dat, form = form, design = design, ps = ps)
}
