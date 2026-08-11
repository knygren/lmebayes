## Pre-generate the expensive vignette outputs shipped in inst/extdata.
##
## The curriculum chapters that need a suggested package (bayesrules,
## glmmTMB) or a slow GLMM fit do not run their sampler at build time.  They
## show the code in an `eval = FALSE` chunk and load the stored result here,
## following the glmbayes convention (see glmbayes/data-raw/make_Chapter13_*).
##
## Run from the package root, with lmebayes installed:
##
##   Rscript data-raw/make_vignette_outputs.R
##
## Every artifact stores *captured console text* plus a handful of small
## numeric summaries, never a whole fit object: a 1000-draw lmerb fit is
## several megabytes, while its printed summary is a few kilobytes.

suppressMessages(library(lmebayes))

out_dir <- file.path("inst", "extdata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Capture everything a chunk would have printed, sampler chatter included,
## as a character vector ready for cat(x, sep = "\n").
capture_print <- function(expr) {
  txt <- utils::capture.output(suppressWarnings(print(expr)))
  txt[!grepl("^\\s*$", txt) | seq_along(txt) %in% which(nzchar(txt))]
}

quiet_fit <- function(expr) {
  utils::capture.output(fit <- suppressWarnings(expr))
  fit
}

message("Writing artifacts to ", normalizePath(out_dir, mustWork = FALSE))


# ---------------------------------------------------------------------------
# Chapter 04 -- linear mixed model on bayesrules::big_word_club
# ---------------------------------------------------------------------------

make_chapter04 <- function() {
  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  keep <- c("score_ppvt", "distracted_a1", "distracted_ppvt",
            "private_school", "title1", "free_reduced_lunch", "school_id")
  dat <- subset(
    dat,
    !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, keep])
  )
  dat$school_id <- droplevels(dat$school_id)

  form_simple <- score_ppvt ~ private_school + title1 +
    distracted_ppvt + distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)

  form_full <- score_ppvt ~ private_school + title1 + free_reduced_lunch +
    distracted_ppvt + distracted_a1 +
    free_reduced_lunch:distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)

  design_simple <- model_setup(form_simple, data = dat)
  design_full   <- model_setup(form_full,   data = dat)

  ps <- Prior_Setup_GLMM(form_full, data = dat, pop.pwt = 0.01)

  set.seed(2024)
  fit <- quiet_fit(lmerb(
    form_full,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$group.dispersion,
    n                = 1000L,
    progbar          = FALSE
  ))

  list(
    n_obs           = nrow(dat),
    n_schools       = nlevels(dat$school_id),
    design_simple   = capture_print(design_simple),
    design_full     = capture_print(design_full),
    prior_setup     = capture_print(ps),
    pfamily_list    = capture_print(pfamily_list(ps)),
    fit_print       = capture_print(fit),
    fit_summary     = capture_print(summary(fit)),
    ranef_head      = capture_print(utils::head(as.data.frame(lme4::ranef(fit)), 12L)),
    popef_mode      = fit$popef.mode,
    popef_means     = fit$popef.means,
    m_convergence   = fit$m_convergence,
    sim_method_used = fit$sim_method_used
  )
}


# ---------------------------------------------------------------------------
# Chapter 09 -- Poisson GLMM on bayesrules::airbnb
# ---------------------------------------------------------------------------

make_chapter09 <- function() {
  data(airbnb, package = "bayesrules")
  ab <- airbnb[, c("neighborhood", "reviews", "rating", "room_type")]
  ab <- ab[complete.cases(ab), ]
  ab$neighborhood <- droplevels(factor(ab$neighborhood))
  ab$rating_c <- ab$rating - mean(ab$rating)

  form <- reviews ~ rating_c + (1 + rating_c || neighborhood)

  design <- model_setup(form, data = ab, family = poisson())
  ps <- Prior_Setup_GLMM(form, data = ab, family = poisson(), pop.pwt = 0.01)

  set.seed(2024)
  fit <- quiet_fit(glmerb(
    form,
    data         = ab,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    n            = 200L,
    progbar      = FALSE
  ))

  list(
    n_obs           = nrow(ab),
    n_groups        = nlevels(ab$neighborhood),
    design          = capture_print(design),
    prior_setup     = capture_print(ps),
    fit_print       = capture_print(fit),
    fit_summary     = capture_print(summary(fit)),
    popef_mode      = fit$popef.mode,
    popef_means     = fit$popef.means,
    m_convergence   = fit$m_convergence,
    pilot           = fit$pilot,
    sim_method_used = fit$sim_method_used
  )
}


# ---------------------------------------------------------------------------
# Chapter 10 -- binomial GLMM on bayesrules::book_banning
# ---------------------------------------------------------------------------

make_chapter10 <- function() {
  data(book_banning, package = "bayesrules")
  bb <- book_banning[, c("state", "removed", "violent")]
  bb <- bb[complete.cases(bb), ]
  bb$removed_i <- as.integer(bb$removed == 1L | bb$removed == "1")
  bb$violent_i <- as.integer(
    bb$violent == TRUE | bb$violent == 1L | bb$violent == "TRUE"
  )
  bb$state <- factor(bb$state)

  ## Keep the twelve largest states so the stored fit stays small.
  keep <- names(sort(table(bb$state), decreasing = TRUE))[seq_len(12L)]
  bb <- droplevels(subset(bb, state %in% keep))

  form <- removed_i ~ violent_i + (1 + violent_i || state)

  design <- model_setup(form, data = bb, family = binomial())
  ps <- Prior_Setup_GLMM(form, data = bb, family = binomial(), pop.pwt = 0.01)

  set.seed(2024)
  fit <- quiet_fit(glmerb(
    form,
    data         = bb,
    family       = binomial(),
    pfamily_list = pfamily_list(ps),
    n            = 200L,
    progbar      = FALSE
  ))

  list(
    n_obs           = nrow(bb),
    n_groups        = nlevels(bb$state),
    states          = levels(bb$state),
    design          = capture_print(design),
    prior_setup     = capture_print(ps),
    fit_print       = capture_print(fit),
    fit_summary     = capture_print(summary(fit)),
    ranef_head      = capture_print(utils::head(as.data.frame(lme4::ranef(fit)), 12L)),
    popef_mode      = fit$popef.mode,
    popef_means     = fit$popef.means,
    m_convergence   = fit$m_convergence,
    sim_method_used = fit$sim_method_used
  )
}


# ---------------------------------------------------------------------------
# Chapter 08 -- per-group observation dispersion via dGamma_list()
# ---------------------------------------------------------------------------
# Sampling a separate sigma^2_j for each of the 18 subjects runs the ING
# dispersion envelope once per group per sweep, which takes about a minute.

make_chapter08 <- function() {
  data(sleepstudy, package = "lme4")
  dat <- sleepstudy
  dat$Days_c <- dat$Days - mean(dat$Days)
  form <- Reaction ~ Days_c + (1 + Days_c || Subject)

  ps_scalar <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
  ps_group  <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01,
                                dispformula = ~Subject)
  dl <- dGamma_list(ps_group)

  set.seed(2024)
  fit <- quiet_fit(lmerb(
    form, data = dat,
    pfamily_list     = pfamily_list(ps_group),
    dispersion_ranef = dl,
    dispformula      = ~Subject,
    n                = 200L,
    progbar          = FALSE
  ))

  s2 <- summary_sigma2(fit)

  list(
    scalar_dispersion = ps_scalar$group.dispersion,
    dGamma_list_first = capture_print(dl[[1L]]),
    dGamma_list_names = names(dl),
    fit_print         = capture_print(fit),
    sigma2_summary    = capture_print(s2),
    m_convergence     = fit$m_convergence,
    sim_method_used   = fit$sim_method_used
  )
}


# ---------------------------------------------------------------------------
# Chapter 13 -- lme4 / glmmTMB comparison on lme4::sleepstudy
# ---------------------------------------------------------------------------

make_chapter13 <- function() {
  data(sleepstudy, package = "lme4")
  dat <- sleepstudy
  dat$Days_c <- dat$Days - mean(dat$Days)
  form <- Reaction ~ Days_c + (1 + Days_c || Subject)

  ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)

  set.seed(2024)
  fit <- quiet_fit(lmerb(
    form, data = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$group.dispersion,
    n                = 1000L,
    progbar          = FALSE
  ))

  m_lmer <- lme4::lmer(form, data = dat, REML = TRUE)
  m_tmb  <- glmmTMB::glmmTMB(form, data = dat, family = gaussian())

  cmp <- data.frame(
    term      = names(lme4::fixef(m_lmer)),
    lmer      = unname(lme4::fixef(m_lmer)),
    glmmTMB   = unname(glmmTMB::fixef(m_tmb)$cond),
    lmerb     = unname(unlist(fit$popef.means)),
    row.names = NULL,
    check.names = FALSE
  )

  vc <- data.frame(
    component = c("(Intercept)", "Days_c", "Residual"),
    lmer_sd   = c(
      sqrt(unlist(lapply(lme4::VarCorr(m_lmer), function(z) z[1L, 1L]))),
      stats::sigma(m_lmer)
    ),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    lmer_summary    = capture_print(summary(m_lmer)),
    tmb_summary     = capture_print(summary(m_tmb)),
    lmerb_print     = capture_print(fit),
    fixef_compare   = cmp,
    varcomp         = vc,
    m_convergence   = fit$m_convergence,
    sim_method_used = fit$sim_method_used
  )
}


# ---------------------------------------------------------------------------

artifacts <- list(
  Chapter04_big_word_club        = make_chapter04,
  Chapter08_group_dispersion     = make_chapter08,
  Chapter09_airbnb_poisson       = make_chapter09,
  Chapter10_book_banning_binom   = make_chapter10,
  Chapter13_mer_comparison       = make_chapter13
)

for (nm in names(artifacts)) {
  message("  building ", nm, " ...")
  obj <- artifacts[[nm]]()
  path <- file.path(out_dir, paste0(nm, ".rds"))
  saveRDS(obj, path, version = 2L)
  message("    wrote ", path, " (",
          format(file.size(path) / 1024, digits = 4), " KB)")
}

message("Done.")
