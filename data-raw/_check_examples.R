pkg <- "lmebayes"
src <- normalizePath("..", winslash = "/")
tar <- file.path(tempdir(), paste0(pkg, "_0.1.0.tar.gz"))
cmd <- function(...) system2("R", c("CMD", ...), stdout = TRUE, stderr = TRUE)
stopifnot(identical(cmd("build", src, "-o", tempdir()), character(0)))
stopifnot(identical(cmd("INSTALL", "--no-multiarch", "--with-keep.source", tar), character(0)))
tools::testInstalledPackage(pkg, outDir = tempdir(), types = "examples")
cat("examples: OK\n")
