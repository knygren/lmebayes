set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
print(d.AD <- data.frame(treatment, outcome, counts))

glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())

## Plain iid GLM fitting engine (rglmb) is not yet re-exported from lmebayes:
## it will eventually live in the stripped-down glmbayesCore, but for now the
## up-to-date copy is in lmebayesCore, so we source it from there directly.
ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())

rglmb.D93 <- lmebayesCore::rglmb(
  n = 200,
  y = ps$y,
  x = as.matrix(ps$x),
  pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
  family = poisson(),
  weights = rep(1, nrow(ps$x))
)
print(rglmb.D93)
summary(rglmb.D93)
