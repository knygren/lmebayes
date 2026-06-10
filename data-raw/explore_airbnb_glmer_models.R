# Explore glmer (Poisson) models on bayesrules::airbnb — find RE with signal.
# Plain lme4 only; no lmebayes.
#
#   Rscript data-raw/explore_airbnb_glmer_models.R
#
# Selected for inst/examples/Ex_glmerb.R (m22):
#   reviews ~ walk_c + transit_c + rating_c + log_price_c +
#     walk_c:rating_c + transit_c:log_price_c +
#     (1 + rating_c + log_price_c || neighborhood)
# Level-2 covariates on intercept; cross-level moderation of random slopes.

library(lme4)
library(bayesrules)

data(airbnb)
dat <- airbnb
dat$rating_c <- dat$rating - mean(dat$rating)
dat$room_type <- factor(dat$room_type)
dat$neighborhood <- factor(dat$neighborhood)

vars <- c("reviews", "rating", "rating_c", "room_type", "neighborhood",
          "walk_score", "transit_score", "bike_score", "price", "bedrooms",
          "bathrooms", "beds", "number_of_reviews")
vars <- intersect(vars, names(dat))
dat <- dat[complete.cases(dat[, vars]), ]

cat("n =", nrow(dat), "  neighborhoods =", nlevels(dat$neighborhood), "\n\n")

# Center / scale listing-level predictors where useful
if ("price" %in% names(dat)) dat$log_price_c <- scale(log(dat$price + 1))[, 1]
if ("bedrooms" %in% names(dat)) dat$bedrooms_c <- dat$bedrooms - mean(dat$bedrooms)
if ("beds" %in% names(dat)) dat$beds_c <- dat$beds - mean(dat$beds)

summarize_fit <- function(fit, label) {
  vc <- VarCorr(fit)
  re_sd <- as.data.frame(vc, row.names = NULL)
  re_sd <- re_sd[re_sd$grp == "neighborhood", , drop = FALSE]
  fe <- fixef(fit)
  sing <- isSingular(fit)
  conv <- fit@optinfo$conv$opt == 0
  cat("=== ", label, " ===\n", sep = "")
  cat("  converged:", conv, "  singular:", sing, "\n")
  cat("  fixef (", length(fe), "):\n", sep = "")
  print(round(fe, 4))
  cat("  RE Std.Dev (all terms):\n")
  if (nrow(re_sd)) print(re_sd[, c("var1", "var2", "sdcor"), drop = FALSE])
  cat("\n")
  re_vars <- re_sd$sdcor^2
  names(re_vars) <- re_sd$var1
  invisible(list(vc = vc, singular = sing, converged = conv, re_sd = re_sd,
                 re_vars = re_vars, total_re_var = sum(re_vars, na.rm = TRUE)))
}

models <- list(
  m01 = reviews ~ rating_c + (1 | neighborhood),
  m02 = reviews ~ rating_c + (1 + rating_c || neighborhood),
  m03 = reviews ~ rating_c + room_type + (1 + rating_c || neighborhood),
  m04 = reviews ~ rating_c + room_type + (1 | neighborhood),
  m05 = reviews ~ rating_c + log_price_c + (1 + rating_c || neighborhood),
  m06 = reviews ~ rating_c + log_price_c + (1 + log_price_c || neighborhood),
  m07 = reviews ~ rating_c + log_price_c + (1 + rating_c + log_price_c || neighborhood),
  m08 = reviews ~ rating_c + bedrooms_c + (1 + rating_c || neighborhood),
  m09 = reviews ~ rating_c + beds_c + (1 + rating_c || neighborhood),
  m10 = reviews ~ rating_c + room_type + log_price_c + (1 + rating_c || neighborhood),
  m11 = reviews ~ rating_c + room_type + (1 + rating_c + log_price_c || neighborhood),
  m12 = reviews ~ rating_c + walk_score + (1 + rating_c || neighborhood)
)

results <- list()
for (nm in names(models)) {
  fit <- tryCatch(
    glmer(models[[nm]], data = dat, family = poisson(),
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl = list(maxfun = 2e5))),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    cat("=== ", nm, " FAILED ===\n", conditionMessage(fit), "\n\n", sep = "")
    results[[nm]] <- list(error = conditionMessage(fit))
  } else {
    results[[nm]] <- summarize_fit(fit, nm)
  }
}

cat("--- Rank by total RE variance (non-singular, converged) ---\n")
scores <- lapply(names(results), function(nm) {
  r <- results[[nm]]
  if (!is.null(r$error) || isTRUE(r$singular) || !isTRUE(r$converged)) {
    return(data.frame(model = nm, score = NA, max_sd = NA, n_re = NA,
                      slope_sd = NA, note = "skip", stringsAsFactors = FALSE))
  }
  sd <- r$re_sd$sdcor
  slope_sd <- sd[r$re_sd$var1 != "(Intercept)"]
  data.frame(
    model = nm,
    score = sum(sd^2, na.rm = TRUE),
    max_sd = max(sd, na.rm = TRUE),
    n_re = length(sd),
    slope_sd = if (length(slope_sd)) max(slope_sd, na.rm = TRUE) else 0,
    note = "ok",
    stringsAsFactors = FALSE
  )
})
score_df <- do.call(rbind, scores)
score_df <- score_df[order(-score_df$score), ]
print(score_df)

cat("\n--- Best candidate formulas (deparse) ---\n")
for (nm in score_df$model[score_df$note == "ok"][1:min(5, sum(score_df$note == "ok"))]) {
  cat(nm, ":", deparse(models[[nm]]), "\n")
}
