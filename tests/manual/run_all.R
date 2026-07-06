# Run all manual MER / sampler validation scripts.
#
#   Rscript tests/manual/run_all.R
#   Rscript tests/manual/run_all.R small   # glmerb Poisson on airbnb_small n=1000

args <- commandArgs(trailingOnly = TRUE)
extra <- if (length(args)) paste(args, collapse = " ") else ""

root <- Sys.getenv("LMEBAYES_ROOT", unset = normalizePath(getwd(), winslash = "/"))
scripts <- c(
  file.path(root, "tests/manual/test_lmerb_mer_re_validation.R"),
  file.path(root, "tests/manual/test_lmerb_dgamma_mer_re_validation.R"),
  file.path(root, "tests/manual/test_glmerb_mer_re_validation.R")
)

for (scr in scripts) {
  cmd <- paste(c("Rscript", shQuote(scr), extra), collapse = " ")
  cat("\n=== ", cmd, " ===\n\n", sep = "")
  status <- system(cmd)
  if (!identical(status, 0L)) {
    stop("Manual validation failed: ", scr, call. = FALSE)
  }
}

cat("\nrun_all: OK\n")
