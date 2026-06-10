## glmerb: neighborhood random intercepts on bayesrules::airbnb
##
## Listing-level review counts (Poisson) with neighborhood-level walkability
## covariates (constant within neighborhood). Random intercepts capture
## unmodeled between-neighborhood heterogeneity.
##
## Workflow: model_setup(), then Prior_Setup_lmebayes(), then glmerb().

data(airbnb, package = "bayesrules")

dat <- airbnb
dat <- dat[complete.cases(dat[, c(
  "reviews", "neighborhood", "walk_score", "transit_score", "bike_score"
)]), ]
dat$walk_c    <- dat$walk_score    - mean(dat$walk_score)
dat$transit_c <- dat$transit_score - mean(dat$transit_score)
dat$bike_c    <- dat$bike_score    - mean(dat$bike_score)

form_glmer <- reviews ~ walk_c + transit_c + bike_c + (1 | neighborhood)

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

cat("\n=== glmer vs glmerb Block 2 posterior mode (log scale) ===\n\n")
fe_glmer <- lme4::fixef(fit$glmer)
fe_mode  <- unlist(fit$coef.mode)
cmp <- data.frame(
  parameter = names(fe_mode),
  glmer     = unname(fe_glmer[names(fe_mode)]),
  glmerb    = unname(fe_mode),
  row.names = NULL
)
cmp$diff <- cmp$glmerb - cmp$glmer
print(cmp)

grp_col  <- fit$model_setup$group_name
re_names <- fit$model_setup$re_coef_names
cat("\nNeighborhood random intercepts on log scale (first 6 levels):\n")
ri_glmer <- coef(fit$glmer)[[grp_col]][, re_names, drop = FALSE]
ri_mode  <- fit$ranef.mode[, re_names, drop = FALSE]
print(head(cbind(glmer = ri_glmer, glmerb = ri_mode)))
