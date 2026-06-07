set.seed(333)
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
print(d.AD <- data.frame(treatment, outcome, counts))

glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())

ps <- glmbayesCore::Prior_Setup(counts ~ outcome + treatment, family = poisson())

rglmb.D93 <- glmbayesCore::rglmb(
  n = 200,
  y = ps$y,
  x = as.matrix(ps$x),
  pfamily = glmbayesCore::dNormal(mu = ps$mu, Sigma = ps$Sigma),
  family = poisson(),
  weights = rep(1, nrow(ps$x))
)
print(rglmb.D93)
summary(rglmb.D93)
