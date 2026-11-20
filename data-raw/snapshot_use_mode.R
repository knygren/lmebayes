## Deterministic before/after check for moving the sweep backend to glmbayesCore.
##
## Runs the indexed two-block sweep with every draw replaced by its posterior
## mode (use_mode = TRUE), on the R and C++ engines, plus the single groups and
## population posterior modes. None of these depend on random numbers.
##
##   Rscript snapshot_use_mode.R save      # before the change
##   Rscript snapshot_use_mode.R compare   # after the change
##
## Builds lmebayesCore with optimisation (debug = FALSE) and loads it without
## recompiling.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "compare"
snap_file <- file.path("C:/Rpackages/lmebayes/data-raw", "snapshot_use_mode.rds")

suppressPackageStartupMessages(library(pkgload))
pkgbuild::compile_dll("C:/Rpackages/lmebayesCore", debug = FALSE, quiet = TRUE)
pkgload::load_all("C:/Rpackages/lmebayesCore", compile = FALSE, quiet = TRUE)
ns <- asNamespace("lmebayesCore")
get_ns <- function(name) get(name, envir = ns)

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
ps <- lmebayesCore::Prior_Setup_GLMM(form, data = dat, family = binomial(),
                                     pop.pwt = 0.01)
pf <- lmebayesCore::pfamily_list(ps)
re <- design$groupef.names
gl <- levels(design$group)
block1_prior <- list(ddef = TRUE, P = get_ns(".rLMM_P_from_pfamily_list")(pf, re))
icm <- get_ns(".two_block_icm_at_start")(
  design = design, prior_list = block1_prior, pfamily_list = pf,
  re_names = re, family = binomial(), tol = 1e-10, maxit = 200L
)
map <- get_ns(".two_block_index_map")(design, re, gl, icm$start)
spec <- get_ns(".two_block_spec")(design, block1_prior, pf, binomial(), re, map)

n_chains <- 4L
state <- get_ns(".two_block_init")(
  n_chains, icm$start, icm$b_start, spec, design, re, gl
)
shift <- c(-0.30, -0.10, 0.15, 0.40)
for (k in seq_along(state$fixef)) {
  for (i in seq_len(n_chains)) {
    state$fixef[[k]][i, ] <- state$fixef[[k]][i, ] +
      shift[i] * seq_len(ncol(state$fixef[[k]])) / k
  }
}

mu <- get_ns(".two_block_groups_prior_mean")(spec, state$fixef, 3L)
y2 <- get_ns(".two_block_population_response")(
  spec, get_ns(".two_block_groups_place")(
    spec, get_ns(".two_block_groups_draw")(spec, mu)$coef_mode), 2L)

sw_r <- get_ns(".two_block_sweep_r")(spec, state, 3L, use_mode = TRUE)
sw_c <- get_ns(".two_block_sweep_cpp")(
  spec, state$b, state$fixef, state$iters, state$iters_ranef,
  3L, TRUE, FALSE, ""
)

out <- list(
  groups_mode_r     = get_ns(".two_block_groups_draw")(spec, mu)$coef_mode,
  groups_mode_cpp   = get_ns(".two_block_groups_draw_cpp")(spec, mu)$coef_mode,
  population_mode_r = get_ns(".two_block_population_draw")(spec, y2, 2L)$coef_mode,
  population_mode_cpp = get_ns(".two_block_population_draw_cpp")(spec, y2, 2L)$coef_mode,
  sweep_r_b         = sw_r$b,
  sweep_cpp_b       = sw_c$b,
  sweep_r_fixef     = sw_r$fixef,
  sweep_cpp_fixef   = sw_c$fixef
)

if (mode == "save") {
  saveRDS(out, snap_file)
  cat("Saved", length(out), "items to", snap_file, "\n")
} else {
  ref <- readRDS(snap_file)
  ok <- TRUE
  for (nm in names(ref)) {
    gap <- max(abs(unlist(ref[[nm]]) - unlist(out[[nm]])))
    pass <- is.finite(gap) && gap <= 1e-8
    ok <- ok && pass
    cat(sprintf("  %-22s max |before - after| = %-10s %s\n",
                nm, format(gap, digits = 3), if (pass) "ok" else "DIFFERENT"))
  }
  cat(if (ok) "All items match (tol 1e-8).\n" else "Some items differ.\n")
}
