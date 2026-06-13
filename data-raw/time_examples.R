# Time every example file in inst/examples (one R session, like R CMD check):
# each file is sourced into a fresh environment with output suppressed, and
# wall-clock seconds are reported sorted from slowest to fastest.
#
#   Rscript data-raw/time_examples.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all("../glmbayesCore", export_all = FALSE, quiet = TRUE)
pkgload::load_all(export_all = FALSE, quiet = TRUE)

files <- list.files("inst/examples", pattern = "^Ex_.*\\.R$", full.names = TRUE)

results <- data.frame(
  file    = basename(files),
  seconds = NA_real_,
  status  = NA_character_,
  stringsAsFactors = FALSE
)

pdf(NULL)  # absorb any plotting without writing Rplots.pdf

for (i in seq_along(files)) {
  f <- files[i]
  cat(sprintf("[%2d/%d] %-40s ", i, length(files), basename(f)))
  env <- new.env(parent = globalenv())
  t0  <- proc.time()
  err <- NULL
  out <- utils::capture.output(
    tryCatch(
      source(f, local = env, echo = FALSE),
      error = function(e) err <<- conditionMessage(e)
    )
  )
  el <- (proc.time() - t0)[["elapsed"]]
  results$seconds[i] <- el
  results$status[i]  <- if (is.null(err)) "ok" else paste("ERROR:", err)
  cat(sprintf("%8.2f s  %s\n", el, if (is.null(err)) "" else "ERROR"))
  if (!is.null(err)) cat("        ", err, "\n")
}

invisible(dev.off())

cat("\n===== Examples by elapsed time (slowest first) =====\n")
ord <- order(-results$seconds)
print(
  data.frame(
    file    = results$file[ord],
    seconds = round(results$seconds[ord], 2),
    status  = results$status[ord]
  ),
  row.names = FALSE, right = FALSE
)
cat(sprintf("\nTotal: %.1f s across %d files\n",
            sum(results$seconds), nrow(results)))
