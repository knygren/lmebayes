## Scratch repro wrapper (not part of any test suite): runs
## group_dgamma_smoke_test_allrank.R with a non-interactive error handler that
## dumps the full call stack (traceback) instead of dropping into recover(),
## so we can see exactly which R frame invoked the C++ dispersion sampler
## with a non-positive disp_lower/disp_upper.
options(error = function() {
  cat("\n\n=== FULL TRACEBACK ===\n\n")
  calls <- sys.calls()
  for (i in seq_along(calls)) {
    cat(sprintf("%3d: ", i))
    print(calls[[i]])
  }
  cat("\n=== END TRACEBACK ===\n\n")
  quit(save = "no", status = 1)
})

source("data-raw/group_dgamma_smoke_test_allrank.R")
