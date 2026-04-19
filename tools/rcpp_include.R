# Invoked by configure / configure.win. Prints a single -I"path" line to stdout (Rcpp include dir).
# Messages go to stderr.
#
# GLMBAYES_RCPP_LIB — optional: explicit library directory that contains Rcpp/ (e.g. CI user library).
# If unset and several libraries contain Rcpp, the newest packageVersion("Rcpp") wins.

ov <- Sys.getenv("GLMBAYES_RCPP_LIB", "")
lp <- .libPaths()
hits <- Filter(function(L) file.exists(file.path(L, "Rcpp", "DESCRIPTION")), lp)

if (!length(hits)) {
  writeLines("configure: no Rcpp under .libPaths(); install Rcpp before building", con = stderr())
  quit(status = 1L)
}

pick_lib <- function() {
  if (nzchar(ov)) {
    d <- ov
    if (!file.exists(file.path(d, "Rcpp", "DESCRIPTION"))) {
      writeLines(sprintf("configure: GLMBAYES_RCPP_LIB=%s does not contain Rcpp", ov), con = stderr())
      quit(status = 1L)
    }
    writeLines(sprintf("configure: Rcpp include from GLMBAYES_RCPP_LIB=%s", d), con = stderr())
    return(d)
  }
  if (length(hits) == 1L) {
    writeLines(sprintf("configure: single Rcpp installation: %s", hits[[1L]]), con = stderr())
    return(hits[[1L]])
  }
  best <- hits[[1L]]
  bv <- packageVersion("Rcpp", lib.loc = best)
  for (i in 2L:length(hits)) {
    v <- packageVersion("Rcpp", lib.loc = hits[[i]])
    if (v > bv) {
      best <- hits[[i]]
      bv <- v
    }
  }
  writeLines(sprintf("configure: multiple Rcpp — using newest (%s): %s", as.character(bv), best), con = stderr())
  for (L in hits) {
    if (!identical(L, best)) {
      writeLines(
        sprintf("configure:   (other) %s @ %s", L, as.character(packageVersion("Rcpp", lib.loc = L))),
        con = stderr()
      )
    }
  }
  best
}

lib <- pick_lib()
inc <- normalizePath(
  system.file("include", package = "Rcpp", lib.loc = lib),
  winslash = "/",
  mustWork = TRUE
)
p <- gsub("\\\\", "/", inc)
cat(sprintf('-I"%s"', p))
