#' Bootstrap for tests/manual scripts (not sourced by R CMD check).
#'
#' @param require_bayesrules Require bayesrules to be installed.
#' @param load_glmbayes_core If TRUE, load sibling glmbayesCore when present.
#' @return Package root path, invisibly.
#' @keywords internal
.manual_test_load <- function(
    require_bayesrules = TRUE,
    load_glmbayes_core = FALSE
) {
  root <- Sys.getenv("LMEBAYES_ROOT", unset = NA_character_)
  if (is.na(root) || !file.exists(file.path(root, "DESCRIPTION"))) {
    root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
  desc <- file.path(root, "DESCRIPTION")
  if (!file.exists(desc) || read.dcf(desc)[1L, "Package"] != "lmebayes") {
    stop(
      "Run from the lmebayes package root, e.g.\n",
      "  cd .../lmebayes\n",
      "  Rscript tests/manual/<script>.R\n",
      "Or set LMEBAYES_ROOT to the package directory.",
      call. = FALSE
    )
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install pkgload.", call. = FALSE)
  }
  if (isTRUE(require_bayesrules) &&
      !requireNamespace("bayesrules", quietly = TRUE)) {
    stop("Install bayesrules.", call. = FALSE)
  }
  if (isTRUE(load_glmbayes_core)) {
    core <- Sys.getenv(
      "GLMBAYESCORE_ROOT",
      unset = file.path(dirname(root), "glmbayesCore")
    )
    if (file.exists(file.path(core, "DESCRIPTION"))) {
      pkgload::load_all(core, export_all = FALSE, quiet = TRUE)
    }
  }
  pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
  invisible(root)
}
