opencl_kernels_available <- function() {
  if (!exists("has_opencl", mode = "function") || !has_opencl()) {
    return(FALSE)
  }
  path <- system.file("cl", "OPENCL.cl", package = "glmbayes")
  nzchar(path) && file.exists(path)
}

skip_if_no_opencl_kernels <- function() {
  skip_if(!opencl_kernels_available(), "OpenCL or kernel sources not available")
}
