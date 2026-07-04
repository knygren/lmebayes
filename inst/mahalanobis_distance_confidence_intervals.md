# Confidence Intervals for Mahalanobis Distance from the Mean in Multivariate Normals

## Setup

For $X \sim N_p(\mu, \Sigma)$, the squared Mahalanobis distance from the mean is

$$D^2 = (X - \mu)'\Sigma^{-1}(X - \mu)$$

The exact distribution of $D^2$ (and the resulting confidence/tolerance thresholds) depends on whether $\mu$ and $\Sigma$ are known or estimated, and — if estimated — whether the point being evaluated is part of the estimation sample or a new independent observation.

---

## Case 1: μ and Σ known

This is the clean case: $D^2 \sim \chi^2_p$ exactly, regardless of $p$. A $(1-\alpha)$ region is

$$D^2 \le \chi^2_{p,\,1-\alpha}$$

If you want an interval on $D$ itself (not $D^2$), take square roots of the endpoints — it's a monotonic transform, so no extra work is needed.

---

## Case 2: μ, Σ estimated from a sample; distance of a new independent observation

This is the practically relevant case. Say $\bar{X}$ and $S$ are estimated from $n$ iid draws, and $X_0$ is a **new**, independent observation from the same population. Then

$$D_0^2 = (X_0 - \bar{X})'S^{-1}(X_0 - \bar{X})$$

is related to Hotelling's $T^2$:

$$T^2 = \frac{n}{n+1}D_0^2$$

and the standard result is

$$\frac{n-p}{p(n-1)}\,T^2 \sim F_{p,\,n-p}$$

Substituting through:

$$\frac{n(n-p)}{p(n-1)(n+1)}\,D_0^2 \sim F_{p,\,n-p}$$

which gives an exact $(1-\alpha)$ threshold:

$$D_0^2 \le \frac{p(n-1)(n+1)}{n(n-p)}\,F_{p,\,n-p,\,1-\alpha}$$

This is what's used for prediction/tolerance regions for a new observation relative to a reference sample (e.g., "is this new patient/sample an outlier relative to the training set?").

---

## Case 3: distance of a point that is part of the estimation sample

This is subtly different — if $X_i$ was one of the $n$ points used to compute $\bar{X}$ and $S$, its distance is no longer independent of the estimates, and the distribution changes. The classical result (Gnanadesikan & Kettenring) is

$$\frac{n}{(n-1)^2}D_i^2 \sim \text{Beta}\!\left(\frac{p}{2},\,\frac{n-p-1}{2}\right)$$

This beta-based threshold is the standard basis for internal multivariate outlier/leverage tests. It matters practically because using the chi-square or F approximation on in-sample points **understates** the tail probability — a common bug in ad hoc outlier-flagging code.

---

## Practical notes

- As $n \to \infty$, both the F- and Beta-based results collapse back to $\chi^2_p$, so for large $n$ relative to $p$ the naive chi-square cutoff is a fine approximation. The corrections mainly matter when $n$ isn't much bigger than $p$.
- If $\Sigma$ is estimated with a robust covariance estimator (e.g., MCD) rather than the sample $S$, these exact finite-sample results no longer hold. People typically fall back on bootstrap resampling of $D^2$ or an asymptotic chi-square with a finite-sample correction factor.
- For non-normal data, none of the exact distributions apply, and bootstrap confidence intervals on $D^2$ (or $D$) are the usual fallback.

### Bayesian framing

For a fully Bayesian version, rather than plugging in $\hat\mu, \hat\Sigma$, one can obtain the posterior predictive distribution of $D^2$ directly by simulation — drawing $\mu, \Sigma$ from their posterior and $X_0$ from the resulting predictive distribution, then computing $D^2$ for each draw to build an empirical interval.

---

## References

- Gnanadesikan, R. & Kettenring, J. R. (1972). *Robust Estimates, Residuals, and Outlier Detection with Multiresponse Data.* Biometrics.
- Johnson, R. A. & Wichern, D. W. *Applied Multivariate Statistical Analysis.* (Standard reference for the Hotelling's $T^2$ / F-distribution result.)
