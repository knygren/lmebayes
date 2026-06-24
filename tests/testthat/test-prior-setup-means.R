test_that("Prior_Setup_lmebayes: null intercept only; other hyperprior means are 0", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(book_banning, package = "bayesrules", envir = environment())
  dat <- book_banning[, c("state", "removed", "explicit", "language", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == TRUE | dat$removed == 1L)
  for (v in c("explicit", "language", "violent")) {
    dat[[paste0(v, "_i")]] <- as.integer(
      dat[[v]] == TRUE | dat[[v]] == 1L | dat[[v]] == "1"
    )
  }

  form <- removed_i ~
    explicit_i + language_i + violent_i +
    (1 + explicit_i + language_i + violent_i || state)

  ps <- Prior_Setup_lmebayes(form, dat, family = binomial(), pwt = 0.01)

  expect_equal(
    unname(ps$prior_list[["(Intercept)"]]$mu_fixef["(Intercept)"]),
    unname(lme4::fixef(lme4::glmer(
      removed_i ~ 1 + (1 | state),
      data = dat,
      family = binomial()
    ))["(Intercept)"])
  )
  expect_equal(unname(ps$prior_list[["explicit_i"]]$mu_fixef["(Intercept)"]), 0)
  expect_equal(unname(ps$prior_list[["language_i"]]$mu_fixef["(Intercept)"]), 0)
  expect_equal(unname(ps$prior_list[["violent_i"]]$mu_fixef["(Intercept)"]), 0)
})

test_that("Prior_Setup_lmebayes: effects_source = full_model uses glmer slopes", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(book_banning, package = "bayesrules", envir = environment())
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == TRUE | dat$removed == 1L)
  dat$violent_i <- as.integer(dat$violent == TRUE | dat$violent == 1L)

  form <- removed_i ~ violent_i + (1 + violent_i || state)
  ps <- Prior_Setup_lmebayes(
    form, dat, family = binomial(), pwt = 0.01,
    effects_source = "full_model"
  )

  fe <- lme4::fixef(ps$fit_ref)
  expect_equal(
    unname(ps$prior_list[["violent_i"]]$mu_fixef["(Intercept)"]),
    unname(fe["violent_i"])
  )
})
