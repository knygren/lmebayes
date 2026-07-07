#' Bind \code{lmebayes:::.validate_manual_block2_fixef} after \code{load_all()}.
#'
#' @keywords internal
.bind_manual_block2_fixef <- function() {
  if (!"lmebayes" %in% loadedNamespaces()) {
    stop("Call .manual_test_load() or devtools::load_all() first.", call. = FALSE)
  }
  ns <- asNamespace("lmebayes")
  assign(
    ".validate_manual_block2_fixef",
    get(".validate_manual_block2_fixef", envir = ns, inherits = FALSE),
    envir = parent.frame()
  )
}
