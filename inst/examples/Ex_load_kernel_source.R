############################### Start of load_kernel_source example ####################

\dontrun{
if (has_opencl()) {
  src <- load_kernel_source("nmath/bd0.cl")
  lib <- load_kernel_library("nmath")
}
}

###############################################################################
## End of load_kernel_source example
###############################################################################
