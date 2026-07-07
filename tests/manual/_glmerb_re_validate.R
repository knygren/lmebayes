#' Bind \code{lmebayes:::.validate_glmerb_re} after \code{load_all()}.
#'
#' Implementation lives in \code{R/lmebayes_mer_re_compare.R}.
#'
#' @keywords internal
.bind_glmerb_re_validate <- function() {
  if (!"lmebayes" %in% loadedNamespaces()) {
    stop("Call .manual_test_load() or devtools::load_all() first.", call. = FALSE)
  }
  assign(
    ".validate_glmerb_re",
    get(".validate_glmerb_re", envir = asNamespace("lmebayes"), inherits = FALSE),
    envir = parent.frame()
  )
}
