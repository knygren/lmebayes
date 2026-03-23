############################### Start of EnvelopeDispersionBuild example ####################

# This example mirrors the current C++ algorithm path for Gaussian regression
# with an independent Normal-Gamma prior:
#   rIndepNormalGammaReg_std:
#     - compute initial dispersion (weighted lm.wfit residual variance)
#     - iterate: sample beta using rNormal_reg at current dispersion, compute
#       RSS_Post2, and update dispersion2 via Gamma posterior update
#     - optimize posterior mode for coefficients (optim + f2/f3)
#     - standardize the model (glmb_Standardize_Model)
#     - build coefficient envelope (EnvelopeBuild)
#     - build dispersion-aware envelope (EnvelopeDispersionBuild)
#     - sort envelope components (EnvelopeSort)
# It stops after envelope construction (no standardized-envelope sampling).

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

ps <- Prior_Setup(weight ~ group, gaussian())

x <- as.matrix(ps$x)
y <- as.vector(ps$y)
mu <- ps$mu
Sigma <- ps$Sigma
shape <- ps$shape
rate <- ps$rate

n_obs <- length(y)
wt <- rep(1, n_obs)
offset2 <- rep(0, n_obs)

# Reconstruct coefficient precision P (matches rindepNormalGamma_reg)
Rchol <- chol(Sigma)
Pinv <- chol2inv(Rchol)
P <- 0.5 * (Pinv + t(Pinv))

famfunc <- glmbfamfunc(gaussian())
f2 <- famfunc$f2
f3 <- famfunc$f3

Gridtype_core <- as.integer(2)  # C++ uses this in the RSS anchoring loop

###############################################################################
# Step A: Initial dispersion2 via weighted lm.wfit residual variance
###############################################################################
y_star <- y - offset2
fit0 <- lm.wfit(x = x, y = y_star, w = wt)
res0 <- fit0$residuals
RSS0 <- sum(res0^2)
p_rank <- as.integer(fit0$rank)

dispersion2 <- RSS0 / (n_obs - p_rank)

###############################################################################
# Step B: Iterate dispersion anchoring: compute RSS_Post2 and update dispersion
###############################################################################
n_beta_draws <- as.integer(10000)
n_rss_iter <- as.integer(10)

n_w <- sum(wt)

RSS_Post2 <- NA_real_
cpp_out <- NULL

for (j in seq_len(n_rss_iter)) {
  # Match C++: call rNormalReg at current dispersion2
  prior_list_loop <- list(
    mu = mu,
    P = P,
    dispersion = dispersion2
  )

  cpp_out <- rNormal_reg(
    n = n_beta_draws,
    y = y,
    x = x,
    prior_list = prior_list_loop,
    offset = offset2,
    weights = wt,
    family = gaussian(),
    Gridtype = Gridtype_core,
    use_parallel = FALSE,
    use_opencl = FALSE,
    verbose = FALSE,
    progbar = FALSE
  )

  beta_draws <- cpp_out$coefficients # n_beta_draws x p

  # Match C++ algebra (gaussian identity):
  #   lp_mat  = beta_draws %*% t(x)
  #   eta_mat = lp_mat + offset
  #   diff    = eta_mat - y
  #   RSS_temp = rowSums( (diff^2) * wt )
  #   RSS_Post2 = mean(RSS_temp)
  lp_mat <- beta_draws %*% t(x)
  eta_mat <- lp_mat + matrix(offset2, nrow = n_beta_draws, ncol = n_obs, byrow = TRUE)
  diff_mat <- eta_mat - matrix(y, nrow = n_beta_draws, ncol = n_obs, byrow = TRUE)
  res_sq <- diff_mat * diff_mat

  res_sq_weighted <- res_sq * matrix(wt, nrow = n_beta_draws, ncol = n_obs, byrow = TRUE)
  RSS_temp <- rowSums(res_sq_weighted)
  RSS_Post2 <- mean(RSS_temp)

  # Update shape2, rate2 and dispersion2 exactly as in C++
  shape2 <- shape + n_w / 2.0
  rate2 <- rate + RSS_Post2 / 2.0
  dispersion2 <- rate2 / (shape2 - 1.0)
}

###############################################################################
# Step C: Coefficient posterior mode optimization (optim + f2/f3)
###############################################################################
bstar_mode <- as.vector(cpp_out$coef.mode)
dispstar <- dispersion2

wt2_opt <- wt / dispstar
alpha <- as.vector(x %*% as.vector(mu) + offset2)

mu2 <- rep(0, length(as.vector(mu)))   # mu2 = 0 * mu (as in C++)
parin <- rep(0, length(as.vector(mu))) # parin = 0 vector (mu - mu)

opt_out <- optim(
  par = parin,
  fn = f2,
  gr = f3,
  y = as.vector(y),
  x = as.matrix(x),
  mu = as.vector(mu2),
  P = as.matrix(P),
  alpha = as.vector(alpha),
  wt = as.vector(wt2_opt),
  method = "BFGS",
  hessian = TRUE
)

bstar <- opt_out$par
A1 <- opt_out$hessian

###############################################################################
# Step D: Standardize model (glmb_Standardize_Model)
###############################################################################
Standard_Mod <- glmb_Standardize_Model(
  y = as.vector(y),
  x = as.matrix(x),
  P = as.matrix(P),
  bstar = as.matrix(bstar, ncol = 1),
  A1 = as.matrix(A1)
)

bstar2 <- Standard_Mod$bstar2
A <- Standard_Mod$A
x2_std <- Standard_Mod$x2
mu2_std <- Standard_Mod$mu2
P2_std <- Standard_Mod$P2

###############################################################################
# Step E: EnvelopeBuild (coefficient envelope at Gridtype = 3)
###############################################################################
max_disp_perc <- 0.99
n_env <- as.integer(200) # used by EnvelopeBuild for diagnostics/overhead
Gridtype_env <- as.integer(3) # EnvelopeOrchestrator overrides to 3

shape2_env <- shape + n_w / 2.0
rate3_env <- rate + RSS_Post2 / 2.0
d1_star <- rate3_env / (shape2_env - 1.0)

wt2_env <- wt / d1_star

Env2 <- EnvelopeBuild(
  bStar = as.vector(bstar2),
  A = as.matrix(A),
  y = as.vector(y),
  x = as.matrix(x2_std),
  mu = as.matrix(mu2_std, ncol = 1),
  P = as.matrix(P2_std),
  alpha = as.vector(alpha),
  wt = as.vector(wt2_env),
  family = "gaussian",
  link = "identity",
  Gridtype = Gridtype_env,
  n = n_env,
  n_envopt = as.integer(1),
  sortgrid = FALSE,
  use_opencl = FALSE,
  verbose = FALSE
)

###############################################################################
# Step F: EnvelopeDispersionBuild (dispersion-aware envelope)
###############################################################################
disp_env_out <- EnvelopeDispersionBuild(
  Env = Env2,
  Shape = shape,
  Rate = rate,
  P = as.matrix(P2_std),
  y = as.vector(y),
  x = as.matrix(x2_std),
  alpha = as.vector(alpha),
  n_obs = as.integer(n_obs),
  RSS_post = RSS_Post2,
  RSS_ML = NA_real_,
  mu = as.matrix(mu2_std, ncol = 1),
  wt = as.vector(wt),
  max_disp_perc = max_disp_perc,
  disp_lower = NULL,
  disp_upper = NULL,
  verbose = FALSE,
  use_parallel = TRUE
)

###############################################################################
# Step G: EnvelopeSort (mirror EnvelopeOrchestrator: disp_grid_type = 2)
###############################################################################
Env3_raw <- disp_env_out$Env_out
UB_list_new <- disp_env_out$UB_list
gamma_list_new <- disp_env_out$gamma_list

cbars <- Env3_raw$cbars
l1 <- ncol(cbars)
l2 <- nrow(cbars)

logP_vec <- Env3_raw$logP
logP_mat <- matrix(logP_vec, nrow = length(logP_vec), ncol = 1)

Env3 <- EnvelopeSort(
  l1 = l1,
  l2 = l2,
  GIndex = Env3_raw$GridIndex,
  G3 = Env3_raw$thetabars,
  cbars = cbars,
  logU = Env3_raw$logU,
  logrt = Env3_raw$logrt,
  loglt = Env3_raw$loglt,
  logP = logP_mat,
  LLconst = Env3_raw$LLconst,
  PLSD = Env3_raw$PLSD,
  a1 = Env3_raw$a1,
  E_draws = Env3_raw$E_draws,
  lg_prob_factor = UB_list_new$lg_prob_factor,
  UB2min = UB_list_new$UB2min
)

UB_list_final <- UB_list_new
UB_list_final$lg_prob_factor <- Env3$lg_prob_factor
UB_list_final$UB2min <- Env3$UB2min

env_final <- list(
  Env = Env3,
  gamma_list = gamma_list_new,
  UB_list = UB_list_final,
  diagnostics = disp_env_out$diagnostics,
  low = gamma_list_new$disp_lower,
  upp = gamma_list_new$disp_upper
)

print(env_final$low)
print(env_final$upp)
print(env_final$gamma_list[c("shape3", "rate2")])

env_final

###############################################################################
# End: envelope construction only
###############################################################################

