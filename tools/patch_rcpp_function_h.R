# Post-install patch for Rcpp headers on R-devel/R 4.5.x when R no longer declares
# R_NamespaceRegistry but Rcpp still uses it in the R 4.5.* #elif branch.
# Replaces the R_getVarEx(..., R_NamespaceRegistry, ...) line with R_getRegisteredNamespace.
#
# Usage: Rscript tools/patch_rcpp_function_h.R
# Lib: GLMBAYES_RCPP_LIB or R_LIBS_USER or .libPaths()[1]

lib <- Sys.getenv("GLMBAYES_RCPP_LIB", "")
if (!nzchar(lib)) lib <- Sys.getenv("R_LIBS_USER", "")
if (!nzchar(lib)) lib <- .libPaths()[1L]

fh <- file.path(lib, "Rcpp", "include", "Rcpp", "Function.h")
if (!file.exists(fh)) {
  message("patch_rcpp_function_h: missing ", fh, " — skip")
  quit(status = 0L)
}

lines <- readLines(fh, warn = FALSE)
orig <- lines

# 1) Single-line form (CRAN / GitHub): entire call on one line
idx <- which(
  grepl("R_getVarEx", lines, fixed = TRUE) &
    grepl("R_NamespaceRegistry", lines, fixed = TRUE) &
    grepl("R_UnboundValue|R_NilValue", lines, perl = TRUE)
)
if (length(idx)) {
  i <- idx[[1L]]
  indent <- regmatches(lines[i], regexpr("^[[:space:]]*", lines[i], perl = TRUE))
  if (!length(indent) || !nzchar(indent[1L])) {
    indent <- ""
  } else {
    indent <- indent[1L]
  }
  lines[i] <- paste0(indent, "Shield env(R_getRegisteredNamespace(ns.c_str()));")
  writeLines(lines, fh)
  message("patch_rcpp_function_h: patched line ", i, " in ", fh)
  quit(status = 0L)
}

# 2) Fallback: whole-file regex (multiline / odd spacing)
txt <- paste(lines, collapse = "\n")
pat <- paste0(
  "(?s)Shield(?:<SEXP>)?[[:space:]]+env\\(",
  "R_getVarEx\\(",
  "Rf_install\\(ns\\.c_str\\(\\)\\),[[:space:]]*",
  "R_NamespaceRegistry,[[:space:]]*FALSE,[[:space:]]*",
  "R_(?:UnboundValue|NilValue)",
  "\\)\\)[[:space:]]*;"
)
repl <- "Shield env(R_getRegisteredNamespace(ns.c_str()));"

if (!grepl(pat, txt, perl = TRUE)) {
  message("patch_rcpp_function_h: no R_getVarEx/R_NamespaceRegistry pattern matched — skip")
  quit(status = 0L)
}

txt2 <- gsub(pat, repl, txt, perl = TRUE)
if (identical(txt2, txt)) {
  quit(status = 0L)
}
writeLines(strsplit(txt2, "\n", fixed = TRUE)[[1L]], fh)
message("patch_rcpp_function_h: patched (regex) ", fh)
