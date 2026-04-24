# CRAN submission comments — glmbayes 0.9.0

## Package summary

glmbayes provides iid sampling for Bayesian Generalized Linear Models
(Gaussian, Poisson, Binomial, Gamma) via accept-reject methods based on
likelihood subgradients (Nygren & Nygren, 2006). It mirrors the interface
of base R's glm() and lm(), and optionally accelerates envelope
construction via OpenCL for high-dimensional models. OpenCL is an optional
capability; the package detects its absence at build time and disables that
code path gracefully — all checks pass on platforms without OpenCL.

## Test environments

### Local (developer machine)
- Windows 11, ASUS TUF F16, GeForce RTX GPU, OpenCL installed
- R 4.6.0 RC, glmbayes built with OpenCL enabled
- Rcpp 1.1.1.1
- Command: `devtools::check(vignettes = TRUE, args = "--as-cran", remote = TRUE, manual = TRUE)`

- 0 errors, 0 warnings, 3 notes

  1. New submission (see Notes)
  2. Rcpp workaround (see Notes)
  3. Long-running examples on OpenCL-enabled machine (see Notes)
   
### Win-builder

- R release 
    -R version 4.6.0 RC (2026-04-22 r89945 ucrt)
    -Rcpp 1.1.1.1    
    -0 errors, 0 warnings, 2 notes

- R-devel   
    -4.6.0 RC(2026-04-20 r89921 ucrt)
    -Rcpp 1.1.1.1    
    -0 errors, 0 warnings, 2 notes

- R-oldrelease 
    -R version 4.5.3 (2026-03-11 ucrt)
    -Rcpp 1.1.1    
    -0 errors, 0 warnings, 3 notes


  1. New submission (see Notes)
  2. Rcpp workaround (see Notes)
  3. Long-running non-OpenCL (see Notes - oldrelease only)

### Mac-builder
- macOS release (mac.R-project.org): 0 errors, 0 warnings, N notes
- macOS devel  (mac.R-project.org): 0 errors, 0 warnings, N notes

### R-universe
- All platforms pass except wasm (WebAssembly), which is expected:
  the package includes compiled C/C++ code that is not compatible
  with the wasm toolchain.

### rhub (via rhub::rhub_check())
- linux, macos-arm64, windows, m1-san, atlas, c23,
  clang16–clang22, gcc13–gcc16, intel, lto, mkl,
  nold, noremap, ubuntu-clang, ubuntu-gcc12,
  ubuntu-release, donttest:

-  0 errors, 0 warnings, N notes
  
  [Note: Rcpp was special handled see Rcpp note below]

- valgrind, clang-asan, clang-ubsan, gcc-asan:

  0 errors, 0 warnings, N notes

- rchk: [describe outcome and explain here]

### GPU / OpenCL on Linux (Vast.ai virtual machine)
- Ubuntu [version], OpenCL enabled, R [version]
- Confirms OpenCL code path builds and runs correctly outside Windows
- Result: 0 errors, 0 warnings, N notes


## Comments Related to Notes appearing on various systems

All checks produced 0 errors and 0 warnings. The following 3 notes were
observed on the local Windows machine (R 4.5.3, OpenCL enabled):

### Note: **New submission** 

       Maintainer: 'Kjell Nygren <kjell.a.nygren@gmail.com>'
       New submission

   Expected for an initial CRAN submission. No action required.

### Note: Rcpp listed in both Imports and Suggests

Rcpp 1.1.1-1 introduced a fix for `R_UnboundValue` required under R 4.6.0
but acknowledged by the Rcpp team (RcppCore/Rcpp#1466) to have backward
compatibility concerns on older R versions. Rcpp 1.1.1-1 is therefore not
available on older platforms. `Imports: Rcpp (>= 1.1.1)` ensures the
package installs on older R platforms; `Suggests: Rcpp (>= 1.1.1-1)` signals
the preference for the newer version where available and allows it to install on 
newer R platforms. This is a temporary workaround pending a stable Rcpp release 
that resolves the version boundary.


### Note: **OpenCL Examples with long CPU or elapsed time**

       Examples with CPU (user + system) or elapsed time > 5s
                        user  system elapsed
       Boston_centered 150.89  16.16  105.20
       Cleveland        42.25   3.00   29.34

   Boston_centered and Cleveland are GPU/OpenCL examples where part of code is guarded by
   `has_opencl()` that does not execute on machines without OpenCL installed.
   They will not appear on CRAN check servers. These examples 
   are used on OpenCL machines to demonstrate bigger models.

### Note: **Non-OpenCL Examples with long CPU or elapsed time**

       Examples with CPU (user + system) or elapsed time > 5s
                user  system elapsed
       rlmb    12.60    0.45   10.61

This appears only on select platforms/machines. On many, this note is never
triggered as elapsed time falls below the 5-second threshold.


### Note on rchk
[rchk checks for PROTECT issues in C code. Describe what rchk flagged,
whether it is a false positive, and what you did to investigate or
mitigate it. If the flag is in Rcpp-generated code rather than your
own C, say so explicitly.]

---
_This file is listed in `.Rbuildignore` and is not included in the built
source tarball. When submitting, paste the content above into the
"Optional comments" field on the CRAN submission form at
https://cran.r-project.org/submit.html._