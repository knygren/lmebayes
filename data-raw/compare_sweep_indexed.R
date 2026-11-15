## Current, indexed, and C++-heavy two-block Gibbs on demo Ex_17.
##
## demo/Ex_17_glmerb_book_banning_slopes.R, unchanged:
##   removed_i ~ violent_i + income_c + income_c:violent_i +
##     (1 + violent_i || state)
##   binomial(), Prior_Setup_GLMM(pop.pwt = 0.01), n = 1000, mode_gap_max = 1,
##   rwt = 0.5.
## The pilot length is the Theorem 3 count. The main-stage length is the
## Corollary 5.1 count at this rwt.
## The three pipelines are separate Monte Carlo samples. For every group and
## population coefficient the report gives the difference of means relative
## to the Monte Carlo SE of that difference.
## glmerb() / rglmerb() send this model through the non-Gaussian GLMM
## engine. The calls use that engine's matrix arguments. One uses
## rGLMM_sweep, one uses .rGLMM_sweep_indexed, and one uses
## .rGLMM_sweep_cpp.
##
##   Rscript C:/Rpackages/lmebayes/data-raw/compare_sweep_indexed.R

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

form_glmerb <- removed_i ~ violent_i + income_c + income_c:violent_i +
  (1 + violent_i || state)

## Larger than the glmerb() default of 0.01 so the Corollary 5.1 main
## stage is short enough to compare the two sweeps.
rwt <- 0.5

design <- lmebayesCore::model_setup(form_glmerb, data = dat, family = binomial())
ps <- lmebayesCore::Prior_Setup_GLMM(
  form_glmerb, data = dat, family = binomial(), pop.pwt = 0.01
)
pf <- lmebayesCore::pfamily_list(ps)

widths <- vapply(design$W, ncol, integer(1))
cat(
  "family=binomial  states=", nlevels(design$group),
  " obs=", nrow(design$D),
  "\nD: ", paste(colnames(design$D), collapse = " | "),
  "\nW: ", paste(names(widths), widths, sep = "=", collapse = ", "),
  "\n", sep = ""
)
flush.console()
if (!any(widths >= 2L)) {
  fail("Ex_17 group designs have no level-2 predictor.")
}

prep <- lmebayesCore:::.rlmerb_prepare_prior(
  design         = design,
  pfamily_list   = pf,
  dispprior_list = NULL,
  family         = binomial(),
  fn_name        = "rglmerb"
)
ma <- lmebayesCore:::.lmebayes_matrix_args_glmm(
  n              = 1000L,
  design         = design,
  prior          = prep$prior,
  family         = binomial(),
  mode_gap_max   = 1.0,
  verbose        = FALSE,
  progbar        = TRUE,
  collect_block1 = TRUE,
  rwt            = rwt
)

group_name <- lmebayesCore:::.lmebayes_resolve_group_name(
  ma$group, NULL, fn_name = "rGLMM_reg_known_vcov"
)
inp <- lmebayesCore:::.rGLMM_validate_matrix_inputs(
  n              = ma$n,
  y              = ma$y,
  D              = ma$D,
  group          = ma$group,
  W              = ma$W,
  tv_tol         = ma$tv_tol,
  group_name     = group_name,
  family         = ma$family,
  mode_gap_max   = ma$mode_gap_max,
  gap_tol        = ma$gap_tol,
  pfamily_list   = ma$pfamily_list,
  dispprior_list = ma$dispprior_list,
  rwt            = ma$rwt
)
if (!identical(inp$family$family, "binomial")) {
  fail("The engine family is not binomial.")
}
inp$offset <- ma$offset
inp$weights <- ma$weights

args <- list(
  inp              = inp,
  icm_tol          = 1e-10,
  icm_maxit        = 200L,
  mode_gap_max     = ma$mode_gap_max,
  verbose          = FALSE,
  progbar          = TRUE,
  stage_verbose    = FALSE,
  rate_calibration = NULL,
  collect_block1   = TRUE,
  engine_label     = "Ex_17_glmerb",
  result_class     = "rGLMM_reg_known_vcov",
  cl               = NULL
)

cat(sprintf(
  "rwt = %g (pilot = Theorem 3; main = Corollary 5.1 at this rwt)\n",
  inp$rwt
))
flush.console()
cat("starting current rGLMM_sweep\n")
flush.console()
cur <- do.call(lmebayesCore:::.rGLMM_reg_run_pipeline, args)
cat("starting indexed rGLMM_sweep\n")
flush.console()
idx <- do.call(lmebayesCore:::.rGLMM_reg_run_pipeline_indexed, args)
cat("starting C++ rGLMM_sweep\n")
flush.console()
cpp <- do.call(lmebayesCore:::.rGLMM_reg_run_pipeline_cpp, args)

pilot_n <- if (is.null(cur$pilot$n)) 0L else as.integer(cur$pilot$n)
cat(sprintf(
  "family=%s  pilot n=%d  main m=%d  engines=%s / %s / %s\n",
  cur$family$family, pilot_n, cur$m_convergence,
  cur$convergence_info$draw_engine,
  idx$convergence_info$draw_engine,
  cpp$convergence_info$draw_engine
))
if (pilot_n < 1L) fail("The non-Gaussian pilot did not run.")
if (!identical(idx$convergence_info$draw_engine, "rGLMM_sweep_indexed")) {
  fail("The indexed pipeline did not record its sweep.")
}
if (!identical(cpp$convergence_info$draw_engine, "rGLMM_sweep_cpp")) {
  fail("The C++ pipeline did not record its sweep.")
}

map <- lmebayesCore:::.two_block_index_map(
  design       = cur$design,
  re_names     = cur$design$groupef.names,
  group_levels = levels(cur$design$group),
  fixef_start  = cur$popef.mode
)
cat(
  "d_col:", paste(map$d_col, collapse = ","),
  " w_col identity:",
  all(vapply(map$w_col, function(z) identical(as.integer(z), seq_along(z)), logical(1))),
  " b_for_wrow identity:",
  all(vapply(map$b_for_wrow, function(z) identical(as.integer(z), seq_along(z)), logical(1))),
  "\n"
)

pair_stats <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  if (!length(y) || length(y) != n) {
    fail("The two samples do not have the same number of draws.")
  }
  sx <- stats::sd(x)
  sy <- stats::sd(y)
  se <- sqrt(sx^2 / n + sy^2 / n)
  diff <- mean(x) - mean(y)
  z <- if (is.finite(se) && se > 0) diff / se else NA_real_
  data.frame(
    n = n,
    diff = diff,
    sd = max(sx, sy, na.rm = TRUE),
    se = se,
    z = z,
    abs_draw = max(abs(x - y)),
    stringsAsFactors = FALSE
  )
}

align_groupef <- function(a, b, gname) {
  key_a <- paste(a$draw, as.character(a[[gname]]), sep = "\r")
  key_b <- paste(b$draw, as.character(b[[gname]]), sep = "\r")
  ord <- match(key_a, key_b)
  if (anyNA(ord)) fail("Group/draw keys do not match.")
  b[ord, , drop = FALSE]
}

groupef_rows <- function(a, b, gname, re) {
  b <- align_groupef(a, b, gname)
  pieces <- list()
  for (lev in unique(as.character(a[[gname]]))) {
    ia <- as.character(a[[gname]]) == lev
    ib <- as.character(b[[gname]]) == lev
    for (k in re) {
      st <- pair_stats(a[ia, k], b[ib, k])
      st$name <- paste(lev, k, sep = " / ")
      pieces[[length(pieces) + 1L]] <- st
    }
  }
  do.call(rbind, pieces)
}

popef_rows <- function(a, b, re) {
  pieces <- list()
  for (k in re) {
    A <- a[[k]]
    B <- b[[k]]
    cn <- colnames(A)
    if (is.null(cn) || is.null(colnames(B)) || !setequal(cn, colnames(B))) {
      fail("Population coefficient names differ for ", k, ".")
    }
    for (col in cn) {
      st <- pair_stats(A[, col], B[, col])
      st$name <- paste(k, col, sep = " / ")
      pieces[[length(pieces) + 1L]] <- st
    }
  }
  do.call(rbind, pieces)
}

report_means <- function(title, rows) {
  z <- rows$z
  finite <- is.finite(z)
  worst <- rows[order(-abs(z), na.last = TRUE), , drop = FALSE]
  worst <- head(worst, 5L)
  cat(sprintf(
    "%s: %d coefficients, n = %d\n  elementwise max |draw difference| = %.3e\n  max |z| = %s   |z| > 2: %d   |z| > 3: %d\n",
    title,
    nrow(rows),
    rows$n[1L],
    max(rows$abs_draw),
    if (any(finite)) sprintf("%.2f", max(abs(z[finite]))) else "NA",
    sum(finite & abs(z) > 2),
    sum(finite & abs(z) > 3)
  ))
  cat("  largest |z|:\n")
  for (i in seq_len(nrow(worst))) {
    cat(sprintf(
      "    %-40s  diff = % .4g  sd = %.4g  se = %.4g  z = %s\n",
      worst$name[i],
      worst$diff[i],
      worst$sd[i],
      worst$se[i],
      if (is.finite(worst$z[i])) sprintf("% .2f", worst$z[i]) else "NA"
    ))
  }
}

gname <- cur$design$group_name
re <- cur$design$groupef.names
n_main <- nrow(cur$popef[[re[1L]]])
cat(sprintf(
  paste0(
    "SE of one mean is sd/sqrt(n); ",
    "SE of a difference of two independent means is sd*sqrt(2/n).\n",
    "At n = %d those factors are %.5f and %.5f.\n"
  ),
  n_main, 1 / sqrt(n_main), sqrt(2 / n_main)
))
report_means(
  "main group coefficients, indexed vs current",
  groupef_rows(cur$groupef, idx$groupef, gname, re)
)
report_means(
  "main population coefficients, indexed vs current",
  popef_rows(cur$popef, idx$popef, re)
)
report_means(
  "pilot group coefficients, indexed vs current",
  groupef_rows(cur$pilot$draws$groupef, idx$pilot$draws$groupef, gname, re)
)
report_means(
  "pilot population coefficients, indexed vs current",
  popef_rows(cur$pilot$draws$popef, idx$pilot$draws$popef, re)
)
report_means(
  "main group coefficients, cpp vs current",
  groupef_rows(cur$groupef, cpp$groupef, gname, re)
)
report_means(
  "main population coefficients, cpp vs current",
  popef_rows(cur$popef, cpp$popef, re)
)
report_means(
  "pilot group coefficients, cpp vs current",
  groupef_rows(cur$pilot$draws$groupef, cpp$pilot$draws$groupef, gname, re)
)
report_means(
  "pilot population coefficients, cpp vs current",
  popef_rows(cur$pilot$draws$popef, cpp$pilot$draws$popef, re)
)

## Same components glmerb() stores after rglmerb().
pipeline_as_glmerb <- function(sampler, formula, design, prior) {
  sampler <- lmebayesCore:::.lmebayes_add_popef_summaries(sampler)
  run_pilot <- !is.null(sampler$pilot$n) && sampler$pilot$n > 0L
  out <- list(
    call            = quote(glmerb()),
    formula         = formula,
    family          = binomial(),
    glmer           = design$glmer,
    dispformula     = ~1,
    dispersion_fit  = NULL,
    prior           = prior,
    model_setup     = design,
    popef.mode      = sampler$popef.mode,
    popef.init      = if (run_pilot) sampler$popef.init else NULL,
    groupef.mode    = sampler$groupef.mode,
    popef.means     = sampler$popef.means,
    popef           = sampler$popef,
    groupef         = sampler$groupef,
    popef.dispersion = sampler$popef.dispersion,
    popef.dispersion.mean = sampler$popef.dispersion.mean,
    popef.iters     = sampler$popef.iters,
    popef.iters.mean = sampler$popef.iters.mean,
    groupef.iters   = sampler$groupef.iters,
    groupef.iters.mean = sampler$groupef.iters.mean,
    group.dispersion = sampler$group.dispersion,
    group.dispersion.mean = sampler$group.dispersion.mean,
    group.dispersion.iters = sampler$group.dispersion.iters,
    group.dispersion.iters.mean = sampler$group.dispersion.iters.mean,
    draw_engine     = sampler$convergence_info$draw_engine,
    sim_method_used = "TWO_BLOCK_GIBBS",
    m_convergence   = sampler$m_convergence,
    pilot           = sampler$pilot,
    gap_tol         = 0.0196,
    mode_gap_max    = 1,
    rwt             = sampler$rwt,
    convergence     = sampler$convergence_info,
    sweep_history   = list(
      pilot = if (run_pilot) sampler$pilot$draws$sweep_history else NULL,
      main  = sampler$sweep_history
    )
  )
  structure(out, class = c("glmerb", "list"))
}

glmerb_cur <- pipeline_as_glmerb(cur, form_glmerb, design, prep$prior)
glmerb_idx <- pipeline_as_glmerb(idx, form_glmerb, design, prep$prior)
glmerb_cpp <- pipeline_as_glmerb(cpp, form_glmerb, design, prep$prior)
cat("glmerb_cur, glmerb_idx, and glmerb_cpp are class glmerb.\n")
