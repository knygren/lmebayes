# An empirical-Bayes hyperprior for group-level dispersion (σ²ⱼ)

Exploratory derivation of a shared population-level hyperprior for
group-level (school-level) measurement dispersion, fit to the classical
per-group variance estimates. As shown in §8, the *prior* stays shared
across groups; group-specific behavior comes from the envelope's
truncation window and from each group's own likelihood at draw time, not
from a group-specific prior mean. Not yet wired into
`Prior_Setup_GLMM()`; see `data-raw/scale_school_smoke_test.R` for the
runnable version this note describes.

## 1. Motivation

`Prior_Setup_GLMM()` currently calibrates **one pooled** measurement
dispersion, `dispersion_ranef = sigma(lmer_fit)^2`, shared by every group.
Classical per-group estimates on the `big_word_club` all-full-rank fixture
(39 schools, `p_re = 2`) show substantial heterogeneity around that pooled
value:

```
Pooled dispersion_ranef:  365.63
Per-school sigma2_hat_g:  min = 114.95, median = 353.29, mean = 388.47,
                          max = 996.91, sd = 204.84
ratio range vs. pooled:   [0.31, 2.73]
```

The single shared σ² forces the Block~1 accept/reject envelope to be
evaluated at a candidate that is a compromise for *every* group at once;
when the true σ²ⱼ vary this widely, no single candidate is a good fit for
all groups simultaneously, and per-group envelope slack accumulates
multiplicatively over blocks — a likely driver of the poor Block~1
acceptance rate at larger group counts. See `data-raw/scale_school_smoke_test.R`
for the full per-school table and the derivation below for how to build a
hyperprior that captures this heterogeneity, as a first step toward
group-specific σ²ⱼ priors.

## 2. Notation

| Symbol | Meaning | Purely data-based? |
|--------|---------|:---:|
| `g = 1..k` | Group (school) index | — |
| `n_g` | Observations in group g | — |
| `df_g = n_g - p_re` | Residual degrees of freedom for group g | — |
| `RSS_g` | Residual sum of squares for group g (from the reference `lmer` fit's `residuals()`, grouped by school) | Yes |
| `s_g² = RSS_g / df_g` | **Group g's own classical dispersion estimate** | **Yes — no pooling** |
| `s₀² = Σ RSS_g / Σ df_g` | **Pooled (grand) classical dispersion estimate** (= `dispersion_ranef`, already computed by `model_setup()`) | **Yes — pooled across all groups, no prior** |
| `σ²ⱼ` | True (unobserved) dispersion for group g | — |
| `τ²` | Between-group heterogeneity: `Var(σ²ⱼ)` across groups | estimated below |
| `a₀, b₀` | Hyperprior `1/σ²ⱼ ~ Gamma(a₀, rate = b₀)` | fit from `(s₀², τ²)` |
| `n₀` | Hyperprior expressed as an equivalent prior sample size (same units as `n_prior_measurement`) | derived from `a₀` |

`s₀²` and `s_g²` are both **pure data summaries** — no prior information
enters either one. `s₀²` is used below purely because it happens to equal
the exact precision-weighted mean of the hierarchical model in §4, which
makes it the natural choice for the hyperprior's mean.

## 3. Model

Treat each group's true dispersion as an unobserved draw from a shared
population distribution, using the same Inverse-Gamma family `dGamma()`
already uses for the pooled prior:

```
1/σ²ⱼ            ~  Gamma(shape = a0, rate = b0)      iid across g = 1..k
RSS_g | sigma2_g ~  sigma2_g * ChiSq(df_g)
```

Equivalently, writing `W_g = RSS_g/2`:

```
W_g | sigma2_g  ~  Gamma(shape = df_g/2, rate = 1/sigma2_g)
```

This is the standard **Gamma–Gamma conjugate pair**: a Gamma prior on the
*rate* of a Gamma likelihood (same mechanism as Poisson–Gamma). Conjugate
updates combine shape with shape and rate with rate:

```
1/sigma2_g | data  ~  Gamma(a0 + df_g/2,  rate = b0 + RSS_g/2)
sigma2_g   | data  ~  InverseGamma(a0 + df_g/2,  rate = b0 + RSS_g/2)
```

## 4. Estimating the hyperprior mean (`s₀²`) — already done

The natural weights for this model are `w_g = df_g / (2 * s0hat^2)`
(inverse chi-square sampling variance). Using them for a weighted mean:

```
sum(w_g * s_g^2) / sum(w_g)
  = sum( [df_g/(2*s0hat^2)] * (RSS_g/df_g) ) / sum( df_g/(2*s0hat^2) )
  = sum(RSS_g) / sum(df_g)
```

which is exactly `s0hat^2` — the pooled estimate. **The precision-weighted
mean of this hierarchical model collapses algebraically to the pooled
`dispersion_ranef` already computed by `model_setup()`.** No separate mean
estimation step is required; only the *second* moment (the between-group
heterogeneity `tau2`) is new.

## 5. Estimating the between-group heterogeneity (`τ²`)

Use a **DerSimonian–Laird**-style moment estimator (the standard
random-effects meta-analysis heterogeneity estimator; here applied to
variances rather than means). It subtracts the *known* within-group
sampling variance of a chi-square estimator, `2 * sigma2_g^2 / df_g`
(exact), from the observed spread of `s_g²`, attributing only the
remainder to real between-group heterogeneity:

```
mu_hat <- s0hat2                                 # pooled estimate, from step 4
w_g    <- df_g / (2 * mu_hat^2)
Q      <- sum(w_g * (s_g2 - mu_hat)^2)
c      <- sum(w_g) - sum(w_g^2) / sum(w_g)
tau2_hat <- max(0, (Q - (k - 1)) / c)
```

`Q` vs. `df = k - 1` is a heterogeneity test: `Q >> df` means the spread in
`s_g²` is larger than sampling noise alone can explain.

**Worked numbers** (39-school fixture, `p_re = 2`):

```
Q = 64.30   (df = 38 under homogeneity)
tau2_hat = 17963.92
```

## 6. Mapping `(mu_hat, tau2_hat)` to Gamma hyperparameters

Same mean-matching convention already used by `ing_prior` /
`ing_prior_measurement` in `Prior_Setup_GLMM()`:

```
a0 <- 2 + mu_hat^2 / tau2_hat
b0 <- mu_hat * (a0 - 1)
n0 <- 2 * a0 - 1 - p_re      # equivalent prior sample size, same units as n_prior_measurement
```

**Worked numbers:**

```
1/sigma2_g ~ Gamma(shape = 9.442, rate = 3087)
Inverse-Gamma mean = 365.63 [= mu_hat], sd = 134.03
n0 (equivalent prior sample size) = 15.9
```

## 7. The per-group shrinkage estimator — and why it is exactly linear

From §3, the posterior of `sigma2_g` is `InverseGamma(A, rate = B)` with
`A = a0 + df_g/2`, `B = b0 + RSS_g/2`. The **Inverse-Gamma mean** is
`rate / (shape - 1)`:

```
E[sigma2_g | data] = (b0 + RSS_g/2) / (a0 - 1 + df_g/2)
```

Substituting `RSS_g = df_g * s_g^2` and `b0 = (a0-1) * mu_hat`:

```
E[sigma2_g | data] = [ (a0-1)*mu_hat + (df_g/2)*s_g^2 ] / [ (a0-1) + df_g/2 ]
```

This is a **weighted average by construction** — a ratio of a sum of two
numerator terms over the sum of their two weights (the "mediant" identity:
for any positive `a, b, c, d`, `(a+c)/(b+d)` equals a weighted average of
`a/b` and `c/d`, weighted by `b` and `d`). Writing it out:

```
sigma2_g_shrunk = w_g * s_g^2 + (1 - w_g) * mu_hat

w_g = (df_g/2) / ((a0-1) + df_g/2) = df_g / (df_g + n0 + p_re - 1)
```

**This is not special to Gaussian models.** It holds for any conjugate
exponential-family pair where the posterior mean is a ratio of two
additively-combined sufficient statistics — true for Normal–Normal
(precision-weighted mean), Beta–Binomial (weighted proportion), and
Gamma–Gamma / Inverse-Gamma (this case). The identical structure underlies
`limma`'s empirical-Bayes variance shrinkage in genomics
(Smyth, 2004): `s_g,post^2 = (d0*s0^2 + d_g*s_g^2) / (d0 + d_g)`.

### Worked numeric check (school 44: n_g = 13, df_g = 11, s_g² = 114.95)

```
A = a0 + df_g/2 = 14.942
B = b0 + RSS_g/2 = 3087 + 632.25 = 3719.25
Posterior mean = B / (A - 1) = 3719.25 / 13.942 = 266.7

w_g = (df_g/2) / (a0 - 1 + df_g/2) = 5.5 / 13.942 = 0.3945
0.3945 * 114.95 + 0.6055 * 365.63 = 45.35 + 221.35 = 266.7   [matches exactly]
```

### A note on Jensen's inequality

It would be **wrong** to claim `E[sigma2_g | data] = 1 / E[1/sigma2_g | data]`
— inversion is nonlinear, so `E[1/X] != 1/E[X]` in general. The derivation
above avoids this by using the Inverse-Gamma's own mean formula,
`rate/(shape-1)`, directly; the `-1` in `shape-1` is precisely the
correction for the nonlinearity of inversion (consistent with
`E[1/X] > 1/E[X]` for `X` a positive random variable). The weighted-average
identity in §7 is exact because it is applied to the *correct* Inverse-Gamma
mean, not to a naive reciprocal of the Gamma-scale mean.

## 8. Next step (not implemented) — one shared prior, per-group windows

A naive design would build **one `dGamma()` per group using
`sigma2_g_shrunk` (§7) as that group's own prior mean**. This is wrong: it
double-counts the group's own data. `sigma2_g_shrunk` already incorporates
group g's own `RSS_g` (via §7's weighted average); feeding it back in as
group g's *prior* mean, and then having Block~1's own envelope sampler
evaluate group g's actual `y_g, Z_g` likelihood on top of that, uses the
same data twice — once baked into the prior, once again in the live
likelihood. That inflates the effective weight of the group's own data well
beyond what `pwt_measurement` is meant to control, and biases the
implied prior information.

**The consistent design keeps the prior shared and only lets the
truncation window vary:**

| Quantity | Shared across groups? |
|---|:---:|
| Prior fed to the sampler: `dGamma(shape = a0, rate = b0)` | **Yes** — identical for every group |
| Posterior mean the sampler *produces* | No — equals `sigma2_g_shrunk` from §7, but obtained by combining the shared prior with group g's own likelihood inside the sampler, not pre-computed and re-injected |
| Envelope truncation window `[disp_lower_g, disp_upper_g]` | No — see below |

Every group is given the **same** `dGamma(a0, rate = b0)` prior object.
Nothing about the prior's information content is group-specific; the
shrinkage in §7 falls out for free once this shared prior meets group g's
own data during sampling — exactly mirroring how one shared ING prior on
`tau2_k` combines with each component's own `b_k` in
`TAU2_ING_FORMULAS.md`.

What *can* legitimately vary per group is the **envelope's truncation
window**, using the identical mean-matched convention
`ING_TRUNCATION_WINDOW.md` already applies to `tau2_k` — mean-matched to
the group's own classical `s_g^2`, with width set by the group's own `df_g`:

```
a_inf,g = (df_g + 1) / 2
b_inf,g = s_g^2 * (df_g - 1) / 2        # mean-matched to s_g^2, not to a0/b0
disp_lower_g = 1 / qgamma(0.99, a_inf,g, rate = b_inf,g)
disp_upper_g = 1 / qgamma(0.01, a_inf,g, rate = b_inf,g)
```

This is legitimate because the window only bounds where the envelope
sampler searches; it does not change what the prior asserts about
`sigma2_g`. It directly targets the mismatch mechanism in §1: instead of
one pooled window that has to cover every group's true σ²ⱼ at once, each
group's window only needs to cover *that* group's own plausible range —
exactly the fix the acceptance-rate problem in §1 calls for.

One caveat: `(a0, b0)` were fit from all `k` groups' `RSS_g/df_g`,
including group g itself, so there is a mild empirical-Bayes "double
dipping" effect in the *hyperprior* fit — standard and generally accepted
practice (Efron, *Large-Scale Inference*, 2010) since any single group's
contribution to a `k`-group hyperprior fit is `O(1/k)`, unlike the full
double-count of the naive per-group-prior design above.

## References

- DerSimonian, R. and Laird, N. (1986). *Meta-analysis in clinical trials.*
  Controlled Clinical Trials, 7(3), 177–188. — heterogeneity-variance moment
  estimator (§5), applied here to variances instead of means.
- Smyth, G.K. (2004). *Linear models and empirical Bayes methods for
  assessing differential expression in microarray experiments.* Statistical
  Applications in Genetics and Molecular Biology, 3(1). — the same
  Gamma–Gamma / Inverse-Gamma variance-shrinkage formula (§7), widely used
  and validated in the `limma` Bioconductor package.
- `inst/TAU2_ING_FORMULAS.md`, `inst/ING_TRUNCATION_WINDOW.md` — the same
  mean-matched Gamma/Inverse-Gamma calibration convention already used
  package-wide for `tau2_k` and the pooled measurement `sigma2`.

## Runnable check

```r
Rscript data-raw/scale_school_smoke_test.R
```

The exploratory block near the top (before the `K_SCHOOLS` loop) computes
the per-school table, the DerSimonian–Laird heterogeneity test, and the
resulting `(a0, b0, n0)` hyperprior from the live `big_word_club` all-rank
fixture.
