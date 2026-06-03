## Dobson (1990) Page 93: Randomized Controlled Trial :
counts    <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome   <- gl(3, 1, 9)
treatment <- gl(3, 3)
print(d.AD <- data.frame(treatment, outcome, counts))

ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
ps

rglmb.D93 <- rglmb(
  n = 200,
  y = ps$y,
  x = as.matrix(ps$x),
  pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
  family = poisson(),
  weights = rep(1, nrow(ps$x))
)
pfamily(rglmb.D93)

## Plant weight data (Gaussian)
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group  <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

ps2 <- Prior_Setup(weight ~ group, family = gaussian())
ps2

pf_norm <- dNormal(mu = ps2$mu, ps2$Sigma, dispersion = ps2$dispersion)
pf_ng <- dNormal_Gamma(
  ps2$mu,
  Sigma_0 = ps2$Sigma_0,
  shape = ps2$shape,
  rate = ps2$rate
)
pf_ing <- dIndependent_Normal_Gamma(
  ps2$mu,
  ps2$Sigma,
  shape = ps2$shape_ING,
  rate = ps2$rate
)

rlmb.D9 <- rlmb(
  n = 200,
  y = ps2$y,
  x = as.matrix(ps2$x),
  pfamily = pf_ing,
  weights = rep(1, length(ps2$y))
)
pfamily(rlmb.D9)
