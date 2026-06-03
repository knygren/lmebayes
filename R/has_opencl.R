#' OpenCL compile-time status for \pkg{lmebayes}
#'
#' @description
#' Returns whether **this** \pkg{lmebayes} installation was built with OpenCL
#' support (link-time / `USE_OPENCL`). This is independent of
#' \code{\link[opencltools:gpu_diagnostics]{has_opencl}()} in \pkg{opencltools},
#' which probes the host runtime (ICD, drivers, headers).
#'
#' For workstation diagnostics (drivers, PATH, ICD), use
#' \code{\link[opencltools:gpu_diagnostics]{diagnose_glmbayes}()} from
#' \pkg{opencltools}. For kernel loading from \code{inst/cl/}, use
#' \code{\link[opencltools:load_kernel_source]{load_kernel_source}()} with the
#' appropriate \code{package} argument once \pkg{opencltools} is installed.
#'
#' @return Logical scalar.
#' @export
#' @seealso \pkg{opencltools}, \pkg{glmbayes}
has_opencl <- function() {
  .has_opencl_cpp()
}
