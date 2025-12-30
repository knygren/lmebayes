# glmbayes (development version – 0.1.0)

## Major Enhancements

### OpenCL and GPU Acceleration
- Completed the OpenCL-based grid construction framework for large models.
- Added GPU-aware envelope sizing and improved OpenCL failure handling.
- Introduced diagnostic utilities to assess OpenCL availability and performance.
- Improved configure scripts to detect OpenCL and provide informative messages.
- Expanded OpenCL documentation and added a dedicated vignette chapter.

### Parallel CPU Sampling (RcppParallel)
- Enabled parallel envelope construction and parallel iid sampling.
- Added pilot functions for large-dimension grid estimation.
- Implemented thread-safe parallel sampling for independent normal-gamma models.

### Core Statistical Improvements
- Migrated to an improved independent normal-gamma simulation algorithm.
- Added theoretical derivations for independent normal-gamma regression.
- Improved UB2 and RSS minimization routines, including scaling corrections.
- Enhanced Prior_Setup() to support family-specific prior construction.
- Added dedicated envelope evaluation and sizing functions.

### Package Infrastructure
- Significant cleanup to remove NOTES and improve CRAN readiness.
- Improved configure and Makevars files for portability.
- Added testthat tests, including OpenCL-specific tests.
- Consolidated envelope-building functions into a cleaner structure.

### Documentation
- Major updates to README and package-level documentation.
- Added multiple new vignettes and expanded existing ones.
- Improved examples for `lmb()`, `rlmb()`, and OpenCL models.

## Bug Fixes
- Corrected scaling in UB2 minimization.
- Improved error handling for missing OpenCL functionality.
- Fixed various small issues uncovered during parallelization work.