# Smoke test: glmbBlock (Poisson) BY neighborhood on bayesrules::airbnb (full default)
#
# Default: full airbnb, n = 1000 draws per block. Poisson has no dispersion parameter;
# priors come from per-block Prior_Setup(family = poisson()).
#
#   Rscript data-raw/test_block_airbnb.R
#   Rscript data-raw/test_block_airbnb.R small   # airbnb_small
#   Rscript data-raw/test_block_airbnb.R quick     # n = 50

USE_FULL_DATA <- TRUE

args <- commandArgs(trailingOnly = TRUE)
use_full <- USE_FULL_DATA
if (any(tolower(args) %in% c("full", "--full"))) use_full <- TRUE
if (any(tolower(args) %in% c("small", "--small"))) use_full <- FALSE
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

if (use_full) {
  data("airbnb", package = "bayesrules")
  airbnb_dat <- airbnb
  message("Using bayesrules::airbnb (n = ", nrow(airbnb_dat), ")")
} else {
  data("airbnb_small", package = "bayesrules")
  airbnb_dat <- airbnb_small
  message("Using bayesrules::airbnb_small (n = ", nrow(airbnb_dat), ")")
}

airbnb_dat$rating_c <- airbnb_dat$rating - mean(airbnb_dat$rating)
airbnb_dat$room_type <- factor(airbnb_dat$room_type)
airbnb_dat <- airbnb_dat[complete.cases(airbnb_dat[, c(
  "reviews", "rating", "rating_c", "room_type", "neighborhood"
)]), ]

form <- reviews ~ rating_c

rank_info <- glmbayes:::.blmb_blocks_full_rank(
  formula = form,
  block = "neighborhood",
  data = airbnb_dat
)
n_nbhd_all <- nrow(rank_info$table)
drop_nbhd <- rank_info$drop
if (length(drop_nbhd)) {
  drop_tab <- rank_info$table[!rank_info$table$full_rank, , drop = FALSE]
  message(
    "Excluded ", nrow(drop_tab), " neighborhood(s) (rank-deficient design):",
    paste0(
      " ", drop_tab$id, " (n=", drop_tab$n, ", rank=", drop_tab$rank,
      ", p=", drop_tab$p, ")",
      collapse = ""
    )
  )
}
airbnb_dat <- airbnb_dat[airbnb_dat$neighborhood %in% rank_info$keep, , drop = FALSE]
airbnb_dat$neighborhood <- droplevels(factor(airbnb_dat$neighborhood))
k_expected <- length(rank_info$keep)
message(
  "Neighborhoods: ", k_expected, " kept (full-rank design), ",
  length(drop_nbhd), " excluded (of ", n_nbhd_all, " total)"
)

set.seed(42)
n_draw <- if (run_quick) 50L else 1000L
message("Posterior draws per block: n = ", n_draw)

ps_block <- Prior_SetupBlock(
  form,
  block = "neighborhood",
  data = airbnb_dat,
  family = poisson()
)
stopifnot(inherits(ps_block, "block_PriorSetup"), length(ps_block) == k_expected)

pfamily_list <- lapply(ps_block, function(ps) {
  dNormal(mu = ps$mu, Sigma = ps$Sigma)
})

out_glmb <- glmbBlock(
  form,
  block = "neighborhood",
  family = poisson(),
  pfamily_list = pfamily_list,
  data = airbnb_dat,
  n = n_draw,
  use_parallel = FALSE
)

stopifnot(inherits(out_glmb, "bglmb"), length(out_glmb) == k_expected)
stopifnot(inherits(out_glmb[[1L]], "glmb"))
stopifnot(nrow(out_glmb[[1L]]$coefficients) == n_draw)
stopifnot(identical(names(out_glmb), names(ps_block)))

print(out_glmb)
s_glmb <- summary(out_glmb)
stopifnot(inherits(s_glmb, "summary.bglmb"), length(s_glmb) == k_expected)
stopifnot(!is.null(attr(s_glmb, "coef_means")))

cat("\nblock_airbnb (Poisson, full-rank BY): OK (k = ", k_expected, ")\n", sep = "")
