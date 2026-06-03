# Sequential multi-response lmb (iris) — compare to glmbayes::lmb with cbind LHS / mlmb
#
#   Rscript data-raw/test_multi_lmb_sequential.R

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
pkgload::load_all(export_all = FALSE)

data("iris", package = "datasets")

set.seed(42)
n_draw <- 150L

form_multi <- cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species

ps_multi <- multi_prior_setup(
  form_multi,
  data = iris,
  family = gaussian()
)
resp <- names(ps_multi)

pfamily_list <- lapply(ps_multi, function(ps) {
  dNormal_Gamma(
    mu = ps$mu,
    Sigma_0 = ps$Sigma_0,
    shape = ps$shape,
    rate = ps$rate
  )
})
stopifnot(length(pfamily_list) == length(resp))

fits <- setNames(vector("list", length(resp)), resp)

for (nm in resp) {
  f_j <- stats::reformulate("Species", response = nm)
  message("lmb: ", deparse(f_j))
  fits[[nm]] <- lmb(
    formula = f_j,
    pfamily = pfamily_list[[nm]],
    data = iris,
    n = n_draw,
    use_parallel = FALSE
  )
}

stopifnot(
  all(vapply(fits, inherits, logical(1), "lmb")),
  nrow(fits[[1]]$coefficients) == n_draw
)

fit_mlm <- stats::lm(form_multi, data = iris)
cat("\n--- lm(cbind(...) ~ Species) coefficients ---\n")
print(stats::coef(fit_mlm))

coef_post <- do.call(cbind, lapply(fits, `[[`, "coef.means"))
rownames(coef_post) <- names(fits[[1]]$coef.means)
colnames(coef_post) <- resp
cat("\n--- Posterior mean coefficients (lmb) ---\n")
print(coef_post)

for (nm in resp) {
  cat("\n========== ", nm, " ==========\n", sep = "")
  print(fits[[nm]])
}

dic_tab <- cbind(
  pD = vapply(fits, `[[`, numeric(1), "pD"),
  DIC = vapply(fits, `[[`, numeric(1), "DIC")
)
rownames(dic_tab) <- resp
cat("\n--- DIC by response ---\n")
print(dic_tab)
cat("Sum DIC:", sum(dic_tab[, "DIC"]), "\n")

out_mlmb <- fits
attr(out_mlmb, "call") <- match.call()
attr(out_mlmb, "formula") <- form_multi
attr(out_mlmb, "coef_names") <- resp
class(out_mlmb) <- "mlmb"

cat("\nSequential multi-lmb (iris): OK\n")
