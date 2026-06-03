#' GPU and OpenCL diagnostics
#'
#' @description
#' \pkg{lmebayes} reuses OpenCL host utilities from \pkg{opencltools} and
#' compile-time probing via \code{\link{has_opencl}()}.
#'
#' @section Package-specific (\pkg{lmebayes}):
#' \itemize{
#'   \item \code{\link{has_opencl}()} --- \code{TRUE} if this \pkg{lmebayes}
#'     build was compiled with OpenCL.
#' }
#'
#' @section Host / runtime (\pkg{opencltools}):
#' \itemize{
#'   \item \code{\link[opencltools:gpu_diagnostics]{diagnose_glmbayes}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{detect_environment_and_gpus}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{detect_compute_runtimes}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{verify_opencl_runtime}()}
#'   \item \code{\link[opencltools:gpu_diagnostics]{check_runtime_env}()}
#'   \item \code{\link[opencltools:load_kernel_source]{load_kernel_source}()},
#'     \code{\link[opencltools:load_kernel_source]{load_kernel_library}()}
#' }
#'
#' @section Modelling with OpenCL (\pkg{glmbayes}):
#' Envelope and GLM sampling with \code{use_opencl = TRUE} use the
#' \pkg{glmbayes} DLL; see \code{\link[glmbayes:gpu_diagnostics]{has_opencl}()}
#' there for that build flag.
#'
#' @keywords internal
#' @name gpu_diagnostics
NULL
