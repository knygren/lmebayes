# Scan non-exported R functions for reachability from package entry points.
# Run: Rscript data-raw/_dead_code_scan.R

files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
defs <- character()
def_file <- character()
for (f in files) {
  for (ln in readLines(f, warn = FALSE)) {
    if (grepl("^[.a-zA-Z][.a-zA-Z0-9_]*\\s*<-\\s*function", ln)) {
      nm <- sub("\\s*<-.*", "", ln)
      defs <- c(defs, nm)
      def_file <- c(def_file, f)
    }
  }
}
names(def_file) <- defs
defs <- unique(defs)

ns <- readLines("NAMESPACE")
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
s3 <- grep("^S3method\\(", ns, value = TRUE)
s3methods <- sub("^S3method\\(([^,]+),([^)]+)\\)$", "\\1.\\2", s3)
roots <- unique(c(exports, s3methods))

all_r <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
ref_dirs <- c("R", "src", "tests", "inst", "data-raw", "demo")
ref_code <- paste(
  unlist(lapply(ref_dirs, function(d) {
    if (!dir.exists(d)) return(character())
    fs <- list.files(d, pattern = "\\.[RrCc]$", recursive = TRUE, full.names = TRUE)
    unlist(lapply(fs, readLines, warn = FALSE))
  })),
  collapse = "\n"
)

is_called <- function(name, text) {
  grepl(paste0(name, "\\s*\\("), text, perl = TRUE) ||
    grepl(paste0(":::", name, "\\s*\\("), text, perl = TRUE) ||
    grepl(paste0("::", name, "\\s*\\("), text, perl = TRUE)
}

callees_of <- setNames(vector("list", length(defs)), defs)
for (caller in defs) {
  body <- readLines(def_file[[caller]], warn = FALSE)
  body <- paste(body, collapse = "\n")
  for (callee in defs) {
    if (identical(caller, callee)) next
    if (is_called(callee, body)) {
      callees_of[[caller]] <- c(callees_of[[caller]], callee)
    }
  }
}

reachable <- character()
queue <- intersect(roots, defs)
while (length(queue)) {
  fn <- queue[[1L]]
  queue <- queue[-1L]
  if (fn %in% reachable) next
  reachable <- c(reachable, fn)
  queue <- unique(c(queue, callees_of[[fn]]))
}

non_export <- setdiff(defs, c(exports, s3methods))
dead_from_exports <- setdiff(non_export, reachable)

# Any use anywhere in package tree (incl. tests, C++ via symbol names)?
used_anywhere <- vapply(non_export, function(nm) {
  is_called(nm, ref_code) && sum(vapply(
    files,
    function(f) any(grepl(paste0("^", nm, "\\s*<-\\s*function"), readLines(f, warn = FALSE))),
    logical(1L)
  )) # definition exists; check refs outside def line
}, logical(1L))

cat("Defined:", length(defs),
    "| Export/S3 roots in R/:", length(intersect(roots, defs)),
    "| Reachable from roots:", length(reachable), "\n\n")

cat("=== Non-exported, NOT reachable from exports (review for removal) ===\n")
for (nm in sort(dead_from_exports)) {
  refs <- is_called(nm, ref_code)
  tag <- if (!refs) "NO REFS" else if (is_called(nm, all_r)) "R-only refs (unreachable chain)" else "refs in tests/demo/data-raw only"
  cat(sprintf("  %-45s %s  [%s]\n", nm, def_file[[nm]], tag))
}

cat("\n=== Non-exported, reachable from exports (keep) ===\n")
for (nm in sort(setdiff(non_export, dead_from_exports))) {
  cat(sprintf("  %s\n", nm))
}
