## Time one compiled groups draw (all states) at likelihood weights 1 - rwt,
## and one C++ sweep step, to find where rwt = 0.75 slows down.
##
##   Rscript C:/Rpackages/lmebayes/data-raw/time_groups_draw_rwt.R

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("C:/Rpackages/lmebayesCore", quiet = TRUE, compile = FALSE)
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
state <- get_ns(".two_block_init")(1L, icm$start, icm$b_start, spec, design, re, gl)
mu_D <- get_ns(".two_block_groups_prior_mean")(spec, state$fixef, 1L)

for (w in c(1, 0.75, 0.5, 0.25, 0.1875)) {
  spec_w <- spec
  spec_w$wt1 <- spec$wt1 * w
  t0 <- proc.time()[["elapsed"]]
  d <- get_ns(".two_block_groups_draw_cpp")(spec_w, mu_D)
  cat(sprintf("weight %.4f: one groups draw %.2f s, mean candidates/state %.1f\n",
              w, proc.time()[["elapsed"]] - t0, d$iters))
  flush.console()
}
