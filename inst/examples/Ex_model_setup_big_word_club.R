## model_setup() on bayesrules::big_word_club (private-school moderation)
##
## Likelihood (per student): random slopes age_c, female by school
## Hyper (per school):
##   (Intercept) ~ private_school + title1 + free_reduced_lunch
##   age_c         ~ 1
##   female        ~ private_school  (cross-level moderation via female:private_school)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Example requires the 'bayesrules' package.", call. = FALSE)
}

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$age_c <- dat$age_months - mean(dat$age_months, na.rm = TRUE)
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "age_c", "female",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  female:private_school +
  (1 + age_c + female || school_id)

design <- model_setup(form_lmer, data = dat)
print(design)

X_int <- design$X_hyper[["(Intercept)"]]
X_age <- design$X_hyper[["age_c"]]
X_fem <- design$X_hyper[["female"]]

cat("\nHyper design matrices (head):\n")
print(utils::head(X_int))
print(utils::head(X_age))
print(utils::head(X_fem))

sid <- levels(design$groups)[2L]
idx <- design$groups == sid
cat("\nSchool", sid, ": y and per-observation Z\n")
print(cbind(y = design$y[idx], design$Z[idx, , drop = FALSE]))

cat("\nvcov_formula:\n")
print(design$vcov_formula)
cat("\nVariance components:\n")
print(design$varcorr)

fe <- lme4::fixef(design$lmer_fit)
gamma_int <- fe[colnames(X_int)]
gamma_fem <- c("(Intercept)" = 0, private_school = unname(fe["private_school:female"]))

cat("\nSchool-level prior means (school 1) from X_hyper %*% gamma:\n")
cat("  intercept:", as.numeric(X_int[1L, ] %*% gamma_int), "\n")
cat("  age_c:     0\n")
cat("  female:   ", as.numeric(X_fem[1L, ] %*% gamma_fem), "\n")
