## glmerb: neighborhood random effects on bayesrules::airbnb
##
## Poisson model for review counts. Listing-level predictors (rating, log price)
## have fixed main effects and neighborhood random slopes. Neighborhood-level
## walkability and transit scores are level-2 covariates on the random intercept;
## they also moderate the random slopes via cross-level interactions
## (walk_c:rating_c, transit_c:log_price_c), analogous to
## free_reduced_lunch:distracted_a1 in inst/examples/Ex_lmerb.R.
##
## Workflow: model_setup(), then Prior_Setup_lmebayes(), then glmerb().

data(airbnb, package = "bayesrules")

dat <- airbnb
dat$rating_c    <- dat$rating - mean(dat$rating)
dat$log_price_c <- scale(log(dat$price + 1))[, 1]
dat$walk_c      <- dat$walk_score - mean(dat$walk_score)
dat$transit_c   <- dat$transit_score - mean(dat$transit_score)
dat <- dat[complete.cases(dat[, c(
  "reviews", "rating", "rating_c", "price", "log_price_c",
  "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
)]), ]

form_glmer <- reviews ~
  walk_c + transit_c +
  rating_c + log_price_c +
  walk_c:rating_c + transit_c:log_price_c +
  (1 + rating_c + log_price_c || neighborhood)

design <- model_setup(form_glmer, data = dat, family = poisson())
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(form_glmer, data = dat, family = poisson(), pwt = 0.01)
cat("\n=== Prior_Setup_lmebayes ===\n\n")
print(ps)

fit <- glmerb(
  form_glmer,
  data = dat,
  family = poisson(),
  measurement_prior_list = ps,
  n = 200L,
  seed = 42L
)

cat("\n=== summary(fit) ===\n\n")
print(summary(fit))

cat("\n=== Block 2 posterior mode vs glmer fixed effects ===\n\n")
fe_glmer <- lme4::fixef(fit$glmer)
for (k in names(fit$coef.mode)) {
  cat(sprintf("\n--- RE component: %s ---\n", k))
  cmp <- data.frame(
    parameter = names(fit$coef.mode[[k]]),
    glmer     = unname(fe_glmer[names(fit$coef.mode[[k]])]),
    glmerb    = unname(fit$coef.mode[[k]]),
    row.names = NULL
  )
  cmp$diff <- cmp$glmerb - cmp$glmer
  print(cmp)
}

grp_col  <- fit$model_setup$group_name
re_names <- fit$model_setup$re_coef_names
cat("\n=== mu_all varies by neighborhood (rating_c slope prior mean) ===\n\n")
mu_rating <- fit$mu_all["rating_c", , drop = TRUE]
walk_by_nbhd <- tapply(dat$walk_c, dat$neighborhood, function(x) x[1L])
walk_by_nbhd <- walk_by_nbhd[names(mu_rating)]
cor_mu_walk <- cor(mu_rating, walk_by_nbhd, use = "complete.obs")
cat(sprintf(
  "  Cor(mu_all[rating_c], neighborhood walk_c): %.3f\n",
  cor_mu_walk
))
cat("  (Positive correlation: higher walkability => higher prior mean rating slope.)\n")

cat("\nNeighborhood random effects (first 6 levels):\n")
ri_glmer <- coef(fit$glmer)[[grp_col]][, re_names, drop = FALSE]
ri_mode  <- fit$ranef.mode[, re_names, drop = FALSE]
print(head(cbind(glmer = ri_glmer, glmerb = ri_mode)))
