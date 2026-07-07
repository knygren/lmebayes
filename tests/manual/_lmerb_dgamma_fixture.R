#' dGamma \code{dispersion_ranef} builders for manual lmerb validation.
#'
#' Model \code{dat} / \code{form} come from \code{\link{.prepare_bwc_lmerb_manual}}
#' in \code{_bwc_lmerb_fixture.R} (same as the main lmerb manual tests).
#'
#' @keywords internal

source("tests/manual/_bwc_lmerb_fixture.R")

#' Build dGamma dispersion_ranef from Prior_Setup ing_prior_measurement.
#'
#' @param ps Output of \code{Prior_Setup_lmebayes()}.
#' @param narrow_window If TRUE, shrink truncation toward the interior (Ex_24).
#' @keywords internal
.dgamma_dispersion_ranef <- function(ps, narrow_window = FALSE) {
  m_disp <- ps$ing_prior_measurement
  disp_lower <- m_disp$disp_lower
  disp_upper <- m_disp$disp_upper
  if (isTRUE(narrow_window)) {
    disp_lower <- disp_lower + 0.25 * (disp_upper - disp_lower)
    disp_upper <- disp_upper - 0.25 * (disp_upper - disp_lower)
  }
  dGamma(
    shape          = m_disp$shape,
    rate           = m_disp$rate,
    beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
    Inv_Dispersion = TRUE,
    disp_lower     = disp_lower,
    disp_upper     = disp_upper
  )
}
