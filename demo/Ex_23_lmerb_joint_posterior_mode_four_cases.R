## Demo: joint posterior mode -- four Gaussian LMM subcases (big_word_club)
##
## Same data and model as demo/Ex_12_lmerb_BigWordClub.R.
## Workflow matches Ex_12 / Ex_14:
##   model_setup() -> Prior_Setup_GLMM() -> pfamily_list(ps) -> lmerb()
##
## Here each case uses lmerb(simulate = FALSE) (joint posterior mode when
## Block~2 is ING and/or dispersion_ranef is dGamma(), else ICM at plug-ins).
##
## Four subcases (same 2x2 as rlmerb):
##   1  fixed sigma^2, fixed tau^2     (dNormal Block~2, ps$dispersion_ranef)
##   2  fixed sigma^2, random tau^2    (ING Block~2)
##   3  random sigma^2, fixed tau^2    (dGamma dispersion_ranef)
##   4  random sigma^2, random tau^2   (ING + dGamma)
##
##   demo("Ex_23_lmerb_joint_posterior_mode_four_cases", package = "lmebayes")

for (pkg in c("bayesrules", "lme4")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("This demo requires the '%s' package.", pkg), call. = FALSE)
  }
}

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_GLMM(form_lmer, data = dat, pwt = 0.01)
cat("\n=== Prior_Setup_GLMM ===\n\n")
print(ps)

pf <- pfamily_list(ps)
pf_ing <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

## Cases 3--4: dGamma() on observation dispersion -- all hyperparameters from ps.
m_disp <- ps$ing_prior_measurement
disp_pf <- dGamma(
  shape          = m_disp$shape,
  rate           = m_disp$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m_disp$disp_lower,
  disp_upper     = m_disp$disp_upper
)

cat("\n=== Case 1: fixed sigma^2, fixed tau^2 ===\n\n")
fit1 <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = ps$dispersion_ranef,
  simulate         = FALSE
)

cat("\n=== Case 2: fixed sigma^2, random tau^2 (ING Block~2) ===\n\n")
fit2 <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf_ing,
  dispersion_ranef = ps$dispersion_ranef,
  simulate         = FALSE
)

cat("\n=== Case 3: random sigma^2, fixed tau^2 (dGamma dispersion_ranef) ===\n\n")
fit3 <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = disp_pf,
  simulate         = FALSE
)

cat("\n=== Case 4: random sigma^2, random tau^2 (ING + dGamma) ===\n\n")
fit4 <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf_ing,
  dispersion_ranef = disp_pf,
  simulate         = FALSE
)

cat("\n=== lmer reference (classical REML) ===\n\n")
print(summary(lme4::lmer(form_lmer, data = dat, REML = TRUE)))

cat("\n=== Case 1 vs lmer (ICM at prior plug-ins; should be close) ===\n\n")
lmebayes:::print_coef_means(fit1)

invisible(list(
  prior_setup = ps,
  case1       = fit1,
  case2       = fit2,
  case3       = fit3,
  case4       = fit4
))