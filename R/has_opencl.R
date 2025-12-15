#' Check if OpenCL support was compiled in
#'
#' This function reports whether the package was compiled with OpenCL support.
#' It calls a small C++ routine that checks the compile-time flag \code{USE_OPENCL}.
#'
#' @usage has_opencl()
#' @return A logical scalar: \code{TRUE} if OpenCL support is available, \code{FALSE} otherwise.
#' @examples
#' if (has_opencl()) {
#'   message("OpenCL is available")
#' } else {
#'   message("OpenCL not available")
#' }
#' @export
has_opencl <- function() {
  .Call(`_glmbayes_has_opencl`)  # call the registered C++ routine directly
}