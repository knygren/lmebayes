test_that("Bayesian binomial-logit regression-OpenCL", {
  skip_if(!has_opencl(), "OpenCL not available")
    # Skip on CRAN to avoid long runtime / external data fetch
  skip_on_cran()
  
  # Load Cleveland dataset (you may want to ship a small version in data/)
  url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data"
  df  <- read.csv(url, header = FALSE, na.strings = "?")
  names(df) <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
                 "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num")
  df <- na.omit(df)
  df$hd <- factor(ifelse(df$num > 0, "Yes", "No"), levels = c("No", "Yes"))
  
  # Prior setup
  ps <- Prior_Setup(hd ~ age + sex + cp + trestbps + chol +
                      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
                    family = binomial(logit),
                    data = df)
  
  mu <- ps$mu
  sigma <- ps$Sigma
  
  # Fit Bayesian GLM with parallel + OpenCL
  glmb_hd <- glmb(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family  = binomial(link = "logit"),
    pfamily = dNormal(mu = mu, Sigma = sigma),
    data    = df,
    n       = 10000,
    Gridtype = 2,
    use_parallel = TRUE,
    use_opencl   = TRUE
    ,verbose      = FALSE
  )
  
  # Test condition: average number of candidates per acceptance < 6
  avg_candidates <- mean(glmb_hd$iters)
  expect_true(avg_candidates < 6)
  

  
})

