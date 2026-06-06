# Smoke tests for model_setup() - grouping factors and fixed-effect level checks.
#
#   0 factors:                      error
#   1 factor, level-2 fixed only:   OK (returns model_setup with X_hyper)
#   1 factor, level-1 fixed:        error
#   2 factors:                      error

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload: install.packages('pkgload')")
}

pkgload::load_all(export_all = FALSE)

data("iris", package = "datasets")

iris$region <- factor(rep(c("north", "south", "east"), each = 50))

form_ok <- Sepal.Length ~ region + (0 + Sepal.Width || Species)

expect_err <- function(expr, label) {
  cat("---", label, "---\n")
  tryCatch(
    expr,
    error = function(e) cat("OK:", conditionMessage(e), "\n")
  )
}

expect_err(
  model_setup(Sepal.Length ~ Sepal.Width, data = iris),
  "0 factors (expect error)"
)

cat("--- 1 factor, level-2 fixed + random slope (expect model_setup) ---\n")
design <- model_setup(form_ok, data = iris)
X_sw <- design$X_hyper[["Sepal.Width"]]
stopifnot(
  inherits(design, "model_setup"),
  identical(design$re_coef_names, "Sepal.Width"),
  nrow(design$Z) == nrow(iris),
  ncol(design$Z) == 1L,
  colnames(design$Z) == "Sepal.Width",
  isTRUE(all.equal(design$Z[, 1L], iris$Sepal.Width)),
  length(design$groups) == nrow(iris),
  length(design$y) == nrow(iris),
  isTRUE(all.equal(
    design$y,
    stats::model.response(lme4::lFormula(form_ok, data = iris)$fr)
  )),
  nrow(X_sw) == 3L,
  ncol(X_sw) == 1L,
  colnames(X_sw) == "(Intercept)",
  all(rownames(X_sw) == levels(iris$Species))
)
cat("OK: X_hyper[[Sepal.Width]] dim =", paste(dim(X_sw), collapse = " x "), "\n")
print(X_sw)

expect_err(
  model_setup(Sepal.Length ~ Sepal.Width + (1 | Species), data = iris),
  "1 factor, level-1 fixed (expect error)"
)

iris2 <- iris
iris2$plate <- factor(rep(1:2, length.out = nrow(iris2)))
expect_err(
  model_setup(
    Sepal.Length ~ region + (1 | Species) + (1 | plate),
    data = iris2
  ),
  "2 factors (expect error)"
)

cat("\nmodel_setup smoke tests: OK\n")
