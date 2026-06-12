## Quick check: Prior_Setup_lmebayes calibrates from ALL groups (full-rank
## status is a design check only).  Run: Rscript data-raw/check_prior_setup_all_groups.R
pkgload::load_all("C:/Rpackages/lmebayes", export_all = FALSE, quiet = TRUE)

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c("score_ppvt", "distracted_a1", "distracted_ppvt",
                           "private_school", "title1", "free_reduced_lunch",
                           "school_id")])
)
form <- score_ppvt ~ private_school + title1 + free_reduced_lunch +
  distracted_a1 + distracted_ppvt + free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
print(ps)

## fit_ref must be the all-groups reference fit from model_setup
stopifnot(identical(ps$fit_ref, ps$design$lmer_fit))
stopifnot(nlevels(ps$design$groups) ==
            nlevels(lme4::getME(ps$fit_ref, "flist")[[1L]]))
cat("\nfit_ref uses all", nlevels(ps$design$groups), "schools: OK\n")
