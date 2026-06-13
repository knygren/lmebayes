# Internal R → C++ bridges.  User-facing R code must call these dotted
# wrappers, not symbols from R/RcppExports.R (compileAttributes output).

#' @noRd
#' @keywords internal
.has_opencl_cpp <- function() {
  .Call(`_lmebayes_has_opencl_cpp_export`)
}
