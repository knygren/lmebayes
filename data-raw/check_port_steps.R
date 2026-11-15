## Step-by-step check of the indexed two-block sweep: R steps vs C++ steps.
##
## Model: demo Ex_17 (binomial logit, (1 + violent_i || state), income_c as
## the level-2 predictor), Prior_Setup_GLMM(pop.pwt = 0.01).
##
## Each R step (lmebayesCore R/two_block_groups_indexed.R,
## R/two_block_population_indexed.R, R/two_block_sweep_indexed.R) and its C++
## twin (glmbayes::sim::twoblock, src/two_block_*.cpp) get identical inputs.
## Outputs that do not depend on random numbers must agree exactly:
##   groups prior mean, groups posterior mode, column placement,
##   population response, population posterior mode, write-back,
##   and a full sweep with every draw replaced by its posterior mode.
## Random draws are not compared. Instead, one run's draws are written to a
## text file and replayed through both production loops (via the dev replay
## observer), so every sampler input (groups prior mean, population response)
## can be compared exactly for every sweep and chain, and against references
## built from names only.
##
##   Rscript C:/Rpackages/lmebayes/data-raw/check_port_steps.R

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("C:/Rpackages/lmebayesCore", quiet = TRUE)
ns <- asNamespace("lmebayesCore")
get_ns <- function(name) get(name, envir = ns)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

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

spec <- get_ns(".two_block_spec")(
  design, block1_prior, pf, binomial(), re, map
)

## Four chains with different population coefficients, so chain-specific
## indexing errors cannot hide behind a common start.
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

results <- list()
check <- function(step, r_val, c_val) {
  same_shape <- identical(dim(r_val), dim(c_val)) &&
    identical(length(r_val), length(c_val))
  gap <- if (same_shape) max(abs(as.numeric(r_val) - as.numeric(c_val))) else NA
  exact <- same_shape && identical(unname(as.numeric(r_val)),
                                   unname(as.numeric(c_val)))
  shape <- if (is.null(dim(r_val))) {
    sprintf("vector %d", length(r_val))
  } else {
    paste(dim(r_val), collapse = " x ")
  }
  cat(sprintf("  %-34s %-12s max |R - C++| = %-10s bitwise identical: %s\n",
              step, shape, format(gap, digits = 3), exact))
  results[[step]] <<- isTRUE(same_shape) && isTRUE(gap <= 1e-12)
  invisible(NULL)
}

i <- 3L
k <- 2L
cat("Steps on chain", i, "and random effect", re[k], "\n")

mu_r <- get_ns(".two_block_groups_prior_mean")(spec, state$fixef, i)
mu_c <- get_ns(".two_block_groups_prior_mean_cpp")(spec, state$fixef, i)
check("groups_prior_mean", mu_r, mu_c)

d1_r <- get_ns(".two_block_groups_draw")(spec, mu_r)
d1_c <- get_ns(".two_block_groups_draw_cpp")(spec, mu_r)
check("groups_draw: posterior mode", d1_r$coef_mode, d1_c$coef_mode)
cat(sprintf("  %-34s R %s, C++ %s\n", "groups_draw: draw shape",
            paste(dim(d1_r$coef), collapse = " x "),
            paste(dim(d1_c$coef), collapse = " x ")))

b_i_r <- get_ns(".two_block_groups_place")(spec, d1_r$coef_mode)
b_i_c <- get_ns(".two_block_groups_place_cpp")(spec, d1_r$coef_mode)
check("groups_place", b_i_r, b_i_c)

y_r <- get_ns(".two_block_population_response")(spec, b_i_r, k)
y_c <- get_ns(".two_block_population_response_cpp")(spec, b_i_r, k)
check("population_response", y_r, y_c)

d2_r <- get_ns(".two_block_population_draw")(spec, y_r, k)
d2_c <- get_ns(".two_block_population_draw_cpp")(spec, y_r, k)
check("population_draw: posterior mode", d2_r$coef_mode, d2_c$coef_mode)

f_r <- get_ns(".two_block_population_write")(spec, state$fixef[[k]], i, k,
                                             d2_r$coef_mode)
f_c <- get_ns(".two_block_population_write_cpp")(spec, state$fixef[[k]], i, k,
                                                 d2_r$coef_mode)
check("population_write", f_r, f_c)

cat("\nR steps against name-based references\n")
fixef_i <- lapply(seq_along(re), function(kk) {
  stats::setNames(state$fixef[[kk]][i, ], names(icm$start[[re[kk]]]))
})
names(fixef_i) <- re
mu_named <- get_ns(".two_block_index_mu_all")(map, fixef_i, design)
check("prior mean vs .two_block_index_mu_all",
      mu_r[spec$d_col, , drop = FALSE], unname(mu_named))
W_k <- spec$W[[k]]
V <- solve(crossprod(W_k) / spec$disp2[k] + spec$P2[[k]])
post_mean <- as.numeric(V %*% (crossprod(W_k, y_r) / spec$disp2[k] +
                                 spec$P2[[k]] %*% spec$mu2[[k]]))
gap_post <- max(abs(d2_r$coef_mode - post_mean))
cat(sprintf("  %-34s max |mode - closed form| = %.3g\n",
            "population mode vs posterior mean", gap_post))
results[["population mode vs posterior mean"]] <- gap_post < 1e-8

n_sweeps <- 3L
cat("\nFull sweep,", n_chains, "chains,", n_sweeps,
    "sweeps, every draw replaced by its posterior mode\n")
sw_r <- get_ns(".two_block_sweep_r")(spec, state, n_sweeps, use_mode = TRUE)
sw_c <- get_ns(".two_block_sweep_cpp")(
  spec, state$b, state$fixef, state$iters, state$iters_ranef,
  n_sweeps, TRUE, FALSE, ""
)
check("sweep: b", sw_r$b, sw_c$b)
for (kk in seq_along(re)) {
  check(paste0("sweep: fixef ", re[kk]), sw_r$fixef[[kk]], sw_c$fixef[[kk]])
}
check("sweep: iters", sw_r$iters, sw_c$iters)
check("sweep: sweep_mean", sw_r$sweep_mean, sw_c$sweep_mean)
check("sweep: sweep_cov", sw_r$sweep_cov, sw_c$sweep_cov)
## The dev observer (no tape, no record) must leave both loops unchanged.
ob_r <- get_ns(".dev_two_block_sweep_replay")(spec, state, n_sweeps,
                                              use_mode = TRUE, engine = "r")
ob_c <- get_ns(".dev_two_block_sweep_replay")(spec, state, n_sweeps,
                                              use_mode = TRUE, engine = "cpp")
check("observer: R loop b unchanged", sw_r$b, ob_r$b)
check("observer: C++ loop b unchanged", sw_c$b, ob_c$b)

## Replay. One R run with random draws records every draw; the draws go to a
## text file; both loops then read the same draws in place of the samplers
## and log every sampler input. Given identical groups draws, the population
## inputs must be identical, and vice versa, for every sweep and chain.
write_tape <- function(log, path) {
  parts <- list(block1 = log$block1, iters1 = log$iters1)
  for (kk in seq_along(log$block2)) {
    parts[[paste0("block2_", kk)]] <- log$block2[[kk]]
  }
  lines <- character(0)
  for (nm in names(parts)) {
    a <- parts[[nm]]
    lines <- c(lines, paste("#", nm, paste(dim(a), collapse = " ")),
               sprintf("%.17g", as.numeric(a)))
  }
  writeLines(lines, path)
}
read_tape <- function(path) {
  lines <- readLines(path)
  hdr <- grep("^#", lines)
  ends <- c(hdr[-1L] - 1L, length(lines))
  parts <- list()
  for (h in seq_along(hdr)) {
    tok <- strsplit(lines[hdr[h]], " ", fixed = TRUE)[[1L]]
    vals <- as.numeric(lines[(hdr[h] + 1L):ends[h]])
    parts[[tok[2L]]] <- array(vals, as.integer(tok[-(1:2)]))
  }
  n_b2 <- length(grep("^block2_", names(parts)))
  list(block1 = parts$block1, iters1 = parts$iters1,
       block2 = unname(parts[paste0("block2_", seq_len(n_b2))]))
}

n_replay <- 4L
cat("\nReplay,", n_chains, "chains,", n_replay,
    "sweeps, draws read from a text file\n")
replay <- get_ns(".dev_two_block_sweep_replay")
rec <- replay(spec, state, n_replay, record = TRUE, engine = "r")
tape_file <- file.path(tempdir(), "two_block_tape.txt")
write_tape(rec$log, tape_file)
tape <- read_tape(tape_file)
cat("  tape file:", tape_file, "\n")
check("tape round trip: block1", rec$log$block1, tape$block1)
check("tape round trip: iters1", rec$log$iters1, tape$iters1)
for (kk in seq_along(re)) {
  check(paste0("tape round trip: block2 ", re[kk]),
        rec$log$block2[[kk]], tape$block2[[kk]])
}

rp_r <- replay(spec, state, n_replay, tape = tape, record = TRUE, engine = "r")
rp_c <- replay(spec, state, n_replay, tape = tape, record = TRUE, engine = "cpp")
cat("  R replay vs C++ replay, all sweeps and chains\n")
check("groups input mu_D", rp_r$log$mu_D, rp_c$log$mu_D)
check("groups draw used", rp_r$log$block1, rp_c$log$block1)
for (kk in seq_along(re)) {
  check(paste0("population input y ", re[kk]), rp_r$log$y[[kk]], rp_c$log$y[[kk]])
  check(paste0("population draw used ", re[kk]),
        rp_r$log$block2[[kk]], rp_c$log$block2[[kk]])
}
check("final b", rp_r$b, rp_c$b)
for (kk in seq_along(re)) {
  check(paste0("final fixef ", re[kk]), rp_r$fixef[[kk]], rp_c$fixef[[kk]])
}
check("final iters", rp_r$iters, rp_c$iters)
check("final iters_ranef", rp_r$iters_ranef, rp_c$iters_ranef)
check("prior mean sum over chains", rp_r$prior_sum, rp_c$prior_sum)
check("prior mean sum of squares", rp_r$prior_sumsq, rp_c$prior_sumsq)

cat("  R replay vs the recorded random run (same inputs reproduced)\n")
check("groups input mu_D (recorded)", rec$log$mu_D, rp_r$log$mu_D)
for (kk in seq_along(re)) {
  check(paste0("population input y (recorded) ", re[kk]),
        rec$log$y[[kk]], rp_r$log$y[[kk]])
}

## Name-based references, built only from design names and the tape. Chain i
## at sweep m must use chain i's population draws from sweep m - 1 (or chain
## i's start) for the groups block, and chain i's groups draw from sweep m for
## the population block.
D_names <- colnames(design$D)
fixef_named <- function(m, i) {
  out <- lapply(seq_along(re), function(kk) {
    if (m == 1L) {
      stats::setNames(state$fixef[[kk]][i, ], names(icm$start[[re[kk]]]))
    } else {
      stats::setNames(tape$block2[[kk]][, i, m - 1L],
                      colnames(design$W[[re[kk]]]))
    }
  })
  names(out) <- re
  out
}
mu_by_name <- function(m, i) {
  g <- fixef_named(m, i)
  mu <- matrix(0, length(re), length(gl), dimnames = list(re, gl))
  for (kk in re) {
    Wk <- as.matrix(design$W[[kk]])
    for (grp in gl) {
      mu[kk, grp] <- sum(Wk[grp, ] * g[[kk]][colnames(Wk)])
    }
  }
  mu
}
y_by_name <- function(m, i, kk) {
  b_named <- matrix(tape$block1[, , i, m], length(gl), length(D_names),
                    dimnames = list(gl, D_names))
  b_named[rownames(design$W[[re[kk]]]), re[kk]]
}
ref_gap <- function(chain_of) {
  gap_mu <- 0
  gap_y <- 0
  for (m in seq_len(n_replay)) {
    for (i in seq_len(n_chains)) {
      mu_log <- rp_r$log$mu_D[match(re, D_names), , i, m, drop = FALSE]
      gap_mu <- max(gap_mu, abs(as.numeric(mu_log) -
                                  as.numeric(mu_by_name(m, chain_of(i)))))
      for (kk in seq_along(re)) {
        gap_y <- max(gap_y, abs(rp_r$log$y[[kk]][, i, m] -
                                  y_by_name(m, chain_of(i), kk)))
      }
    }
  }
  c(mu = gap_mu, y = gap_y)
}
cat("  Inputs vs name-based references, all sweeps and chains\n")
own <- ref_gap(function(i) i)
cat(sprintf("  %-34s max gap %.3g\n", "groups input mu_D, own chain", own[["mu"]]))
cat(sprintf("  %-34s max gap %.3g\n", "population input y, own chain", own[["y"]]))
results[["name reference mu_D"]] <- own[["mu"]] <= 1e-12
results[["name reference y"]] <- own[["y"]] <= 1e-12
other <- ref_gap(function(i) i %% n_chains + 1L)
cat(sprintf("  %-34s max gap %.3g (must be > 0)\n",
            "same, with chains shifted by one", max(other)))
results[["chain shift detected"]] <- min(other) > 1e-6

failed <- names(results)[!unlist(results)]
if (length(failed) > 0L) {
  stop("Steps that do not match: ", paste(failed, collapse = ", "), call. = FALSE)
}
cat("\nAll steps match.\n")
