#' Bind \code{lmebayes:::.validate_lmerb_re} after \code{load_all()}.
#'
#' Implementation lives in \code{R/lmebayes_mer_re_compare.R} so
#' \code{devtools::load_all()} picks up changes.  Do not duplicate the body here.
#'
#' @keywords internal
.bind_lmerb_re_validate <- function() {
  if (!"lmebayes" %in% loadedNamespaces()) {
    stop("Call .manual_test_load() or devtools::load_all() first.", call. = FALSE)
  }
  assign(
    ".validate_lmerb_re",
    get(".validate_lmerb_re", envir = asNamespace("lmebayes"), inherits = FALSE),
    envir = parent.frame()
  )
}
