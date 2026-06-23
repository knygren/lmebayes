test_that("Prior_Setup_lmebayes stops when glmer checkConv fails (Ex_18)", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(book_banning, package = "bayesrules", envir = environment())
  reasons <- c(
    "explicit", "antifamily", "occult",
    "language", "lgbtq", "violent"
  )
  dat <- book_banning[, c("state", "removed", reasons)]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == TRUE | dat$removed == 1L)
  for (v in reasons) {
    dat[[paste0(v, "_i")]] <- as.integer(dat[[v]] == TRUE | dat[[v]] == 1L)
  }

  fix <- paste(paste0(reasons, "_i"), collapse = " + ")
  re <- paste0(
    "(1 + ", paste(paste0(reasons, "_i"), collapse = " + "), " || state)"
  )
  form <- stats::as.formula(paste("removed_i ~", fix, "+", re))

  expect_error(
    suppressWarnings(
      Prior_Setup_lmebayes(form, dat, family = binomial(), pwt = 0.01)
    ),
    "requires converged glmer reference fits",
    fixed = FALSE
  )
})

test_that("Prior_Setup_lmebayes accepts converged glmer (violent only, Ex_16)", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(book_banning, package = "bayesrules", envir = environment())
  dat <- book_banning[, c("state", "removed", "violent")]
  dat <- dat[stats::complete.cases(dat), ]
  dat$removed_i <- as.integer(dat$removed == TRUE | dat$removed == 1L)
  dat$violent_i <- as.integer(dat$violent == TRUE | dat$violent == 1L)

  form <- removed_i ~ violent_i + (1 + violent_i || state)
  ps <- Prior_Setup_lmebayes(form, dat, family = binomial(), pwt = 0.01)

  expect_s3_class(ps, "lmebayes_prior_setup")
  expect_true(inherits(ps$fit_ref, "glmerMod"))
})
