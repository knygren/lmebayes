## rglmerb() on demo Ex_17 now samples through .rGLMM_sweep_cpp.
##
## Same model as compare_sweep_indexed.R:
##   removed_i ~ violent_i + income_c + income_c:violent_i +
##     (1 + violent_i || state)
##   binomial(), Prior_Setup_GLMM(pop.pwt = 0.01), n = 1000, mode_gap_max = 1.
## rwt = 1 (compare_sweep_indexed.R uses 0.5) shortens the Corollary 5.1 main
## stage. The pilot length (Theorem 3) does not depend on rwt.
## Checks: draw_engine, pilot and main sweep_history, and both convergence
## plots for both stages.
##
##   Rscript C:/Rpackages/lmebayes/data-raw/check_glmerb_cpp_route.R

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("C:/Rpackages/lmebayesCore", quiet = TRUE)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

fail <- function(...) stop(paste0(...), call. = FALSE)

data(book_banning, package = "bayesrules")
dat <- book_banning[, c("state", "removed", "violent", "median_income")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
dat$income_c <- as.numeric(scale(dat$median_income))

form <- removed_i ~ violent_i + income_c + income_c:violent_i +
  (1 + violent_i || state)
design <- lmebayesCore::model_setup(form, data = dat, family = binomial())
ps <- lmebayesCore::Prior_Setup_GLMM(
  form, data = dat, family = binomial(), pop.pwt = 0.01
)
pf <- lmebayesCore::pfamily_list(ps)

t0 <- proc.time()[["elapsed"]]
fit <- lmebayesCore::rglmerb(
  n            = 1000L,
  design       = design,
  pfamily_list = pf,
  family       = binomial(),
  mode_gap_max = 1.0,
  verbose      = FALSE,
  progbar      = TRUE,
  rwt          = 1
)
elapsed <- proc.time()[["elapsed"]] - t0

show <- function(x) if (is.null(x)) "NULL" else format(x)
engine <- fit$convergence_info$draw_engine
cat(sprintf(
  "class=%s  draw_engine=%s  pilot n=%s  pilot m=%s  main m=%s  elapsed=%.1fs\n",
  paste(class(fit), collapse = "/"), show(engine),
  show(fit$pilot$n), show(fit$pilot$m_convergence),
  show(fit$m_convergence), elapsed
))
if (!identical(engine, "rGLMM_sweep_cpp")) {
  fail("rglmerb() did not use the C++ sweep (draw_engine = ", show(engine), ").")
}

for (stage in c("pilot", "main")) {
  hist <- if (stage == "pilot") fit$pilot$draws$sweep_history else fit$sweep_history
  if (!inherits(hist, "two_block_sweep_history")) {
    fail("No ", stage, " sweep_history on the fit.")
  }
  if (hist$n_sweeps < 1L || nrow(hist$table) == 0L) {
    fail("The ", stage, " sweep_history is empty.")
  }
  cat(sprintf("%s sweep_history: stage label '%s', %d sweeps, %d table rows\n",
              stage, show(hist$stage), hist$n_sweeps, nrow(hist$table)))
}

png_file <- file.path(tempdir(), "check_glmerb_cpp_route_%d.png")
grDevices::png(png_file)
for (stage in c("pilot", "main")) {
  lmebayesCore::plot_mean_convergence(fit, stage = stage)
  lmebayesCore::plot_var_convergence(fit, stage = stage)
}
grDevices::dev.off()
cat("Convergence plots for pilot and main written to", dirname(png_file), "\n")
cat("rglmerb() C++ route checks passed.\n")
