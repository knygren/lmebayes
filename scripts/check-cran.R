# CRAN-oriented checks for glmbayes (maintainer script; excluded from package build).
#
# Prerequisites:
#   - Working directory = package root (directory containing DESCRIPTION).
#   - devtools installed.
#   - Win-builder calls need network access; ensure Maintainer email in DESCRIPTION is valid.
#
# Usage:
#   setwd("/path/to/glmbayes")   # or RStudio: project root
#   source("scripts/check-cran.R")

if (!file.exists("DESCRIPTION")) {
  stop("Set working directory to the glmbayes package root (where DESCRIPTION is).")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools: install.packages(\"devtools\")")
}

# Local check with CRAN incoming-style environment
devtools::check(cran = TRUE)

# Submit source package to win-builder (results arrive by email; runs can be long).
# Comment out any of the following if you only want a subset.
devtools::check_win_release()
devtools::check_win_devel()
devtools::check_win_oldrelease()

# rhub::rhub_check() needs a GitHub token — set GITHUB_PAT in .Renviron or the shell, never in source.
if (!nzchar(Sys.getenv("GITHUB_PAT", ""))) {
  message("Skipping rhub::rhub_check(): set environment variable GITHUB_PAT (do not commit tokens).")
} else {
  rhub::rhub_check(platform = "windows")
}
