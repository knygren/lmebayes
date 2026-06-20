# Block Gibbs Ergodicity and the Full-Rank Condition

Maintainer note on why `block_prior_setup()`, `block_lmb()`, and
`block_rNormalGLM()` require each block's design matrix to be **full column
rank**, and what happens to geometric ergodicity of a coupled block Gibbs chain
when that condition fails.

---

## 1. The scalar two-block model

The cleanest illustration uses **univariate** β and μ with **fixed** variance
components.

| Symbol | Role |
|--------|------|
| μ | Population mean (hyper / block 1) |
| β_b | Neighborhood coefficient (data / block 2) |
| P₁₁ | Prior precision on μ |
| A = 1/τ² | Coupling precision (β \| μ prior) |
| C_b = x_b²/σ² | Likelihood precision contributed by block b |

Two-block Gibbs alternates:

1. **Draw** μ \| β, y
2. **Draw** β_b \| y_b, μ  for each block b

### Joint precision matrix (posterior)

Expanding the joint log-density of (μ, β) as a quadratic form:

```
        ┌ P₁₁ + A     −A    ┐
P  =    │                   │
        └   −A       A + C_b┘
```

- Off-diagonal: P₁₂ = −A (coupling only, no data)
- β diagonal:   P₂₂ = A + C_b
- μ diagonal:   P₁₁_total = P₁₁ + A

---

## 2. What happens when there is no data for β (C_b = 0)?

A rank-deficient X_b contributes **zero** likelihood precision in the deficient
direction: C_b = 0.  The precision matrix becomes:

```
        ┌ P₁₁ + A     −A ┐
P  =    │                │      det(P) = A · P₁₁
        └   −A          A┘
```

**Key structural fact: P₂₂ = |P₁₂| = A.**

The data block's precision equals (in magnitude) the cross-precision.  This is
the signature that β receives *all* its precision from the coupling to μ, none
from data.

### Covariance matrix Σ = P⁻¹

```
         1       ┌  A      A    ┐
Σ = ─────────  × │              │
     A · P₁₁    └  A    P₁₁+A  ┘
```

Reading off:

```
Var(μ)      = 1/P₁₁
Var(β)      = (P₁₁ + A) / (A · P₁₁)
Cov(μ, β)   = 1/P₁₁  =  Var(μ)
```

**Correlation:**

```
           ┌     A    ┐^(1/2)
Cor(μ,β) = │ ─────── │         ∈ [0, 1)
           └  P₁₁+A  ┘
```

---

## 3. Spectral radius of the Gibbs chain

For a 2D Gaussian Gibbs sampler the contraction factor (spectral radius of the
Markov operator on the centred subspace) is:

```
       P₁₂²              A²                A
r = ─────────────  =  ───────────────  =  ─────────    (C_b = 0 case)
    P₁₁_tot · P₂₂    (P₁₁+A)(A+C_b)      P₁₁+A
```

| Regime | r | Chain behaviour |
|--------|---|-----------------|
| P₁₁ → ∞ (tight μ prior) | r → 0 | Instant mixing |
| P₁₁ > 0, C_b = 0 | 0 < r < 1 | Geometrically ergodic |
| P₁₁ → 0, C_b = 0 | r → 1 | Mixing time → ∞ |
| P₁₁ = 0, C_b = 0 | det(P) = 0 | Joint improper, **null recurrent** |

### Data rescues the chain

With C_b > 0 (even one observation):

```
       A²
r = ──────────────  <  ──── A ────   for any P₁₁ ≥ 0
    (P₁₁+A)(A+C_b)     P₁₁+A
```

Even a flat prior on μ (P₁₁ ≈ 0) gives `r ≈ A/(A+C_b) < 1` — the chain is
still geometrically ergodic because **data pins down β**, breaking the
deterministic dependence on μ.

---

## 4. Two neighborhoods: one with data, one without

Extending to two neighborhoods (β₁ with data C₁ > 0, β₂ with no data C₂ = 0)
using the simple intercept-only hyper prior β_b ~ N(μ, τ²).

### Joint 3×3 precision in order (μ, β₁, β₂)

```
        ┌ P₁₁ + 2A     −A        −A  ┐
P  =    │   −A        A + C₁     0   │
        └   −A          0         A  ┘
```

β₁ and β₂ are **conditionally independent given μ** (zero off-diagonal between
them).

### Integrate out β₂

β₂ has no data so its full conditional is β₂ | μ ~ N(μ, τ²) — just the
prior.  The Schur complement for β₂ contributes correction A to P_μμ:

```
P_marg(μ, β₁)  =  [P₁₁ + 2A   −A  ]  −  [A  0]  =  [P₁₁ + A   −A  ]
                   [  −A      A+C₁  ]     [0  0]     [  −A     A+C₁ ]
```

The 2A coupling collapses back to A.  **Integrating out a data-free β₂ is
informationally equivalent to its never having existed** — it carries no
evidence about μ beyond the prior coupling already assigned to it.

### Marginal precision of μ and ergodicity

```
det P_marg  =  (P₁₁ + A)(A + C₁) − A²  =  A·P₁₁ + C₁(P₁₁ + A)
```

With P₁₁ = 0 (flat μ prior):

```
det  =  A · C₁  > 0    because A > 0 and C₁ > 0
```

The joint is **proper** — and the chain geometrically ergodic — purely because
β₁ has data, regardless of β₂.

Spectral radius (P₁₁ = 0, C₁ > 0):

```
       A²             A
r  =  ────────  =  ──────  < 1
      A(A+C₁)      A+C₁
```

Compare with the all-no-data case (Section 2): det = A·P₁₁ → 0 as P₁₁ → 0.
A single data-bearing neighborhood rescues the chain.

---

## 5. Hyper regression: β_b ~ N(X_nbhd[b,] μ, τ²)

The Airbnb neighborhood benchmark uses a richer hyper model: the mean of each
β_b is not a common scalar μ but a **linear function of neighborhood-level
covariates** (walk score, transit score, bike score):

```
β_b | μ  ~  N(x_b^(nbhd) μ,  τ²)
```

where x_b^(nbhd) is a q-vector of centered neighborhood scores (intercept +
walk_c + transit_c + bike_c, so q = 4) and μ is the q-vector of hyper
regression coefficients.

### Joint precision for (μ, β₁, β₂) with this prior

```
        ┌ P₁₁ + (x₁x₁'+x₂x₂')A    −x₁A      −x₂A ┐
P  =    │       −x₁A              A + C₁       0   │
        └       −x₂A                0           A  ┘
```

### Integrate out β₂ (no data, C₂ = 0)

Schur complement correction for β₂: x₂ (1/A) x₂' A = x₂x₂'A.

```
P_marg(μ, β₁)  =  [P₁₁ + x₁x₁'A    −x₁A  ]
                   [    −x₁A'        A+C₁  ]
```

x₂x₂'A drops out exactly — **the non-identified neighborhood's hyper
covariate x₂ is irrelevant**.  Its direction in the μ-space remains identified
(or not) solely through the data-bearing blocks.

### Conditional precision of μ with P₁₁ = 0

```
P_{μ | β₁}  =  x₁x₁'A  −  x₁A (A+C₁)⁻¹ Ax₁'  =  x₁x₁' · A·C₁/(A+C₁)
```

This is positive definite iff **x₁ ≠ 0**, i.e. neighborhood 1's score vector
is non-zero.  More generally, for k neighborhoods the marginal precision of μ
after integrating out all non-identified blocks is:

```
P_{μ|marg}  =  P₁₁  +  ∑_{b: C_b > 0}  x_b^(nbhd) (x_b^(nbhd))' · A·C_b/(A+C_b)
```

For P₁₁ = 0 this is **positive definite** iff the neighborhood-score vectors
of the **data-bearing** blocks collectively span R^q:

```
rank( X_nbhd[ data-bearing rows, ] )  ==  q
```

### The two-level full-rank criterion

| Level | Condition | What it ensures |
|-------|-----------|-----------------|
| **Data (per block)** | `qr(X_b)$rank == l1` | C_b > 0: block β is identified by its data |
| **Hyper (across blocks)** | `qr(X_nbhd[kept rows,])$rank == q` | μ is identified by the data-bearing β's through the neighborhood covariate regression |

**Both levels must hold for geometric ergodicity with an improper or
near-flat prior on μ.**  If either fails, there exists a direction in (β, μ)
space along which the chain is null recurrent.

### The Airbnb benchmarks

`benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block_covariates.R` uses:

```r
X_nbhd <- cbind(
  `(Intercept)` = 1,
  walk_c   = nbhd_unique$walk_score   - mean(nbhd_unique$walk_score),
  transit_c = nbhd_unique$transit_score - mean(nbhd_unique$transit_score),
  bike_c   = nbhd_unique$bike_score   - mean(nbhd_unique$bike_score)
)
```

The simpler `benchmark_airbnb_neighborhood_rNormalGLM_reg_block.R` uses:

```r
X_hyper <- kronecker(diag(l1), matrix(1, k, 1))
```

which is the intercept-only special case (x_b = 1 for all b).  There step 2
reduces to: at least one full-rank data block exists.

Both benchmarks satisfy step 2 because the retained neighborhoods have genuine
variation in walk/transit/bike scores (and the score vectors span R⁴), but the
condition is non-trivial: adding a neighborhood-level covariate that is
**constant across all retained neighborhoods** would introduce a degenerate
direction in μ exactly as rank deficiency in X_b introduces one in β_b.

---

## 6. From scalar to the multivariate / Poisson case

### Full rank in each block

In the **p-dimensional Gaussian** case with fixed Σ the precision structure
is the same block-by-block, but with matrices instead of scalars.  Full column
rank of X_b means C_b = X_b' X_b / σ² is **positive definite**, so each
eigenvalue of the data block contributes positively to P₂₂ − P₁₂ P₁₁⁻¹ P₂₁
(the Schur complement in the β block).  Rank deficiency → at least one zero
eigenvalue in C_b → at least one direction with C_b = 0, landing exactly in
the degenerate scalar picture above.

**Full column rank of every X_b is therefore necessary for the joint posterior
to be proper in all directions, and necessary for geometric ergodicity of the
two-block Gibbs chain when P₁₁ is small or improper.**

### Poisson (and other log-concave GLMs)

For non-Gaussian likelihoods, C_b is replaced by the **expected Fisher
information** at the current draw.  The intuition carries over:

- Rank-deficient X_b → zero information in some β direction → the β
  conditional is improper along that direction (with a proper but non-informative
  prior on μ) → Gibbs chain is not ergodic in that direction.
- Full-rank X_b is **necessary** for the block conditional to be proper.
- **Sufficient** for geometric ergodicity of the *full* coupled chain requires
  additional conditions (bounded log-concave conditionals, coupling through the
  hyper block, etc.) — see `DESIGN_RGLM_BLOCKS.md §4b`.

---

## 7. Package implementation

### Level 1 (per-block): enforced automatically

`block_prior_setup()`, `block_lmb()`, and the helper `.blmb_blocks_full_rank()`
enforce full column rank of X_b for every retained block.  The filter in
`data-raw/test_block_airbnb.R` illustrates this:

```r
rank_info <- glmbayes:::.blmb_blocks_full_rank(
  formula = reviews ~ rating_c,
  block   = "neighborhood",
  data    = airbnb_dat
)
# keep only neighborhoods where qr(X_b)$rank == ncol(X_b)
airbnb_dat <- airbnb_dat[airbnb_dat$neighborhood %in% rank_info$keep, ]
```

Saturated blocks (n = p, full rank) are retained: β is identified by the data;
only the residual variance is not (for Gaussian); the Poisson path needs no
dispersion, so n = p blocks are usable.

Rank-deficient blocks (rank < p) are excluded — they correspond exactly to the
C_b = 0 failure mode above, where the Gibbs chain degenerates.

### Level 2 (hyper design): caller's responsibility

When the hyper prior is β_b ~ N(X_nbhd[b,] μ, Σ) with neighborhood-level
covariates (e.g. walk/transit/bike scores in the Airbnb benchmarks), the caller
must verify that the retained blocks collectively identify μ:

```r
# After applying the level-1 rank filter:
stopifnot(qr(X_nbhd[rank_info$keep, ])$rank == ncol(X_nbhd))
```

For the intercept-only hyper design (X_hyper = kronecker(diag(l1), ones(k,1))),
level 2 is automatic once any full-rank data block exists.  For richer hyper
designs it is a genuine condition: a neighborhood-level covariate that is
constant across all retained blocks introduces a degenerate direction in μ
exactly as rank deficiency in X_b introduces one in β_b.

---

## 8. Recommended identifiability workflow

The two-level analysis above suggests the following preflight procedure before
running a coupled block Gibbs sampler.  It is implemented by
`block_check_identifiability()`.

### Algorithm

```
Step 1.  Temporarily subset to blocks with full-rank data design:
         identified ← { b : qr(X_b)$rank == l1 }

Step 2.  Check that μ is identified by those blocks through the hyper design:
         is_level2 ← qr( X_nbhd[identified, ] )$rank == q

Step 3.  If is_level2 is FALSE:
           stop (or warn) — the chain will be null recurrent in some
           direction of μ even with a proper prior; the full model is
           not identifiable.

Step 4.  If is_level2 is TRUE:
           return to the full model (all blocks, including non-identified).
           Non-identified β_b's draw from their prior N(X_nbhd[b,] μ, τ²)
           each iteration; they are harmless passengers — they cannot
           disrupt ergodicity as long as Step 2 holds.
```

### Why Step 4 includes non-identified blocks

A block with rank-deficient X_b contributes C_b = 0.  From Section 4, its
β_b integrates out cleanly and leaves the marginal precision of μ unchanged.
In the running Gibbs chain that means:

- β_b | μ draws from the prior (no likelihood information) — proper because
  the prior is proper.
- μ | β receives β_b's contribution, but since β_b ~ prior(μ), it carries
  no additional signal about μ beyond what the data-bearing blocks already
  provide.

The chain remains geometrically ergodic (Step 2 ensures μ is identified through
the data-bearing β's) and the non-identified β_b marginals correctly reflect
prior-level uncertainty — which is the honest posterior answer when there are no
listings in that neighborhood.

### Contrast with BY models (block_lmb / block_glmb)

For independent BY fits (no coupling across blocks) non-identified blocks cannot
be "carried along" — each neighborhood is fitted in isolation, so a rank-deficient
block has an improper posterior.  Those are dropped (Step 1 only), and there is
no Step 2 / Step 4.

### R usage

```r
# Simple case: intercept-only hyper (X_nbhd = NULL → all-ones assumed)
chk <- lmebayes:::block_check_identifiability(
  formula = reviews ~ rating_c,
  block   = "neighborhood",
  data    = airbnb_dat
)

# Richer hyper design: neighborhood-level covariates
chk <- lmebayes:::block_check_identifiability(
  formula = reviews ~ rating_c + room_type,
  block   = "neighborhood",
  data    = airbnb_dat,
  X_nbhd  = X_nbhd          # k × q matrix (intercept + walk_c + ...)
)

chk$level1_table    # data.frame: id, n, rank, p, full_rank
chk$level1_keep     # character: block ids passing step 1
chk$level1_drop     # character: rank-deficient block ids
chk$level2_rank     # integer: rank of X_nbhd[keep, ]
chk$level2_ok       # logical: step 2 satisfied?
chk$action          # "proceed" | "warn" | "stop"
```

The function emits a **warning** when level-2 fails (non-identified μ
direction) and an **error** when `on_failure = "stop"` is passed.  With
`on_failure = "warn"` (default) the caller receives the result and decides.

---

## References

- Liu, J.S., Wong, W.H., Kong, A. (1994). Covariance structure of the Gibbs
  sampler with applications to the comparisons of estimators and augmentation
  schemes. *Biometrika* **81**(1), 27–40.
- Roberts, G.O., Rosenthal, J.S. (1997). Geometric ergodicity and hybrid Markov
  chains. *Electronic Communications in Probability* **2**, 13–25.
- See also `inst/DESIGN_RGLM_BLOCKS.md` for the full block-sampler API contract.
