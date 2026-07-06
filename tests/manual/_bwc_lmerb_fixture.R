#' Shared Big Word Club model for manual lmerb validation scripts.
#'
#' Same \code{dat}, \code{form}, and school levels as
#' \code{test_lmerb_mer_re_validation.R} and
#' \code{test_lmerb_dgamma_mer_re_validation.R} (only the dispersion /
#' Block~2 prior route differs across those scripts).
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

  design <- model_setup(form, data = dat)

  ps_args <- list(form = form, data = dat, pwt = pwt)
  if (!is.null(pwt_dispersion)) {
    ps_args$pwt_dispersion <- pwt_dispersion
  }
  ps <- do.call(Prior_Setup_lmebayes, ps_args)

  list(dat = dat, form = form, design = design, ps = ps)
}
