#' Bayesian linear mixed-effects model fit
#'
#' Fits a Bayesian linear mixed model from an \code{lmer}-style formula,
#' returning posterior draws. It is the mixed-model counterpart of
#' \code{\link{lmb}} for fixed-effects models.
#'
#' \code{lmerb} builds the design with \code{\link{model_setup}}, embeds the
#' companion \code{\link[lme4]{lmer}} fit for reference, and passes both to
#' the \pkg{lmebayesCore} sampler. Priors are supplied as
#' \code{pfamily_list} (population level) and \code{dispersion_ranef}
#' (residual variance), normally built from one
#' \code{\link{Prior_Setup_GLMM}} object.
#'
#' @details
#' ## The model being estimated
#'
#' For group \eqn{j = 1, \dots, J}, \code{lmerb} estimates
#' \deqn{y_j \mid \beta_j, \sigma^2 \sim N(D_j \beta_j,\; \sigma^2 I_{n_j})
#'   \qquad \textrm{(within group)}}
#' \deqn{\beta_j \mid \gamma, \Psi \sim N(\mathcal{W}_j \gamma,\; \Psi)
#'   \qquad \textrm{(across groups)}}
#'
#' The second line is the model for how groups differ: each group's own
#' coefficient vector \eqn{\beta_j} varies around a population expectation
#' \eqn{\mathcal{W}_j\gamma}, and \eqn{\Psi} measures how widely. This is a
#' \emph{centered} parameterization --- the full coefficient \eqn{\beta_j}
#' enters the likelihood, not a mean-zero deviation --- so \code{coef()} on
#' the fit reports \eqn{\beta_j} directly and \code{ranef()} reports the
#' deviation \eqn{\beta_j - \mathcal{W}_j\gamma} that \pkg{lme4} calls
#' \eqn{b_j}. Marginalizing recovers the usual \pkg{lme4} model, with
#' \eqn{\gamma} playing the role of the fixed effects.
#'
#' \eqn{\Psi} is diagonal because \code{formula} may only contain
#' uncorrelated (\code{||}) random-effect terms and a single grouping
#' factor; there is no random-effect correlation to estimate. See
#' \code{vignette("Chapter-02", package = "lmebayes")} for the formula rules.
#'
#' ## What the priors decide
#'
#' The two prior arguments are where you declare which variances are known
#' and which are estimated, and that choice determines what comes back.
#'
#' \describe{
#'   \item{\code{pfamily_list}}{One component per random-effect coefficient.
#'     A \code{dNormal} component fixes the between-group variance
#'     \eqn{\tau^2_k} at the supplied value; a
#'     \code{dIndependent_Normal_Gamma} component estimates it, under a
#'     Gamma prior on \eqn{1/\tau^2_k} truncated to
#'     \code{[disp_lower, disp_upper]}.}
#'   \item{\code{dispersion_ranef}}{The residual variance \eqn{\sigma^2}:
#'     a fixed scalar, a pooled \code{dGamma()} prior to estimate one
#'     common value, a \code{dGamma_list()} to estimate a separate
#'     \eqn{\sigma^2_j} per group, or a named vector of fixed per-group
#'     values. \code{dispformula} must agree with the shape supplied.}
#' }
#'
#' When every variance is fixed, the posterior is exactly multivariate
#' normal and the draws are exact. Estimating any variance makes the
#' posterior non-Gaussian and the draws approximate in the sense described
#' next.
#'
#' ## What the returned draws are
#'
#' \code{lmerb} does not return one long chain. It runs \code{n}
#' \strong{independent replicate chains} and keeps \strong{one draw from
#' each}, taken at that chain's final sweep. Two things follow:
#'
#' \enumerate{
#'   \item \strong{The draws are mutually independent.} They form an iid
#'     sample, so ordinary Monte Carlo standard errors apply. Nothing needs
#'     thinning, there is no effective sample size to compute, and trace
#'     plots or \eqn{\hat{R}}{R-hat} on the returned draws are not
#'     meaningful.
#'   \item \strong{Each draw comes from a distribution deliberately made
#'     close to the posterior.} The sweep count is not a burn-in guess: it
#'     is calibrated before sampling so that the distribution reached at the
#'     final sweep is within \code{tv_tol} of the exact posterior in total
#'     variation. The resulting \code{m_convergence} is reported on the fit.
#' }
#'
#' \code{n} and \code{tv_tol} therefore control different errors --- Monte
#' Carlo noise and per-draw bias --- and can be set independently. To
#' certify all \code{n} draws jointly at level \eqn{\alpha}, pass
#' \code{tv_tol = alpha / n}; cost grows only logarithmically in
#' \eqn{1/}\code{tv_tol}, so this usually adds a few sweeps.
#'
#' The calibration itself, and what to do if a warning reports slow
#' convergence, are covered in
#' \code{vignette("Chapter-06", package = "lmebayes")}.
#'
#' @references
#' Nygren, K. (2020). \emph{On the total variation distance between multivariate
#' normal densities with applications to two-block Gibbs samplers.}
#' Unpublished manuscript.
#'
#' Jones, G. L. and Hobert, J. P. (2001). Honest exploration of intractable
#' probability distributions via Markov chain Monte Carlo.
#' \emph{Statistical Science} \bold{16}, 312--334.
#'
#' @param formula Mixed-model formula (single grouping factor; same constraints
#'   as \code{\link{model_setup}}).
#' @param data Data frame containing all variables in \code{formula}.
#' @param pfamily_list Required named list of
#'   \code{\link[glmbayesCore]{pfamily}} objects, one per random-effect
#'   coefficient (names must match the random-effect coefficient names, any
#'   order).  Supplies the Block~2 hyperpriors (\code{mu}, \code{Sigma}) and
#'   the Block~1 random-effect variances \eqn{\tau^2_k}.  \code{dNormal}
#'   components treat \eqn{\tau^2_k} (the pfamily \code{dispersion}) as
#'   known and make conjugate \eqn{\gamma_k} draws.
#'   \code{dIndependent_Normal_Gamma} components place a Gamma prior on the
#'   Block~2 precision \eqn{1/\tau^2_k}: Block~2 then makes a joint
#'   \eqn{(\gamma_k, \tau^2_k)} draw via the likelihood-subgradient envelope
#'   sampler (\code{\link[glmbayesCore]{rindepNormalGamma_reg}}), and the
#'   sampled \eqn{\tau^2_k} feeds back into the Block~1 prior precision.
#'   ING components must supply both truncation bounds: each
#'   \eqn{\tau^2_k} draw is hard-truncated to
#'   \code{[disp_lower, disp_upper]}, fixed across all inner Gibbs sweeps,
#'   with \code{disp_lower} doubling as the conservative \eqn{\tau^2_k}
#'   plug-in for the eigenvalue / TV calibration (smaller \eqn{\tau^2}
#'   increases the contraction rate \eqn{\lambda^*}, so the bound holds
#'   for every dispersion in the truncated support).  They must also
#'   satisfy the prior-vs-data guard \eqn{n_{\mathrm{prior}} \le J}
#'   (\code{pwt_dispersion} \eqn{\le 0.5}).  Typically built with
#'   \code{\link[lmebayesCore:pfamily_list.Prior_Setup_GLMM]{pfamily_list}} from a
#'   \code{\link{Prior_Setup_GLMM}} object.
#' @param dispersion_ranef Observation-level measurement dispersion
#'   \eqn{\sigma^2} for Block~1.  One of: a positive scalar (treated as
#'   known; typically \code{Prior_Setup_GLMM(...)$dispersion_ranef}), a
#'   single \code{\link[glmbayesCore]{dGamma}()} \code{pfamily} (pooled
#'   \eqn{\sigma^2} across groups), a named list of \code{dGamma()}
#'   objects (one per group level) from
#'   \code{\link[lmebayesCore:dGamma_list.Prior_Setup_GLMM]{dGamma_list}(Prior_Setup_GLMM(...))},
#'   or a named numeric vector of positive, fixed per-group values (names
#'   must match the random-effects grouping factor's levels exactly; each
#'   group's \eqn{\sigma^2_j} is then treated as known, like the pooled
#'   scalar case but allowed to vary by group). Which of these four shapes
#'   is accepted depends on \code{dispformula} (see below).
#' @param dispformula One-sided formula selecting the measurement-dispersion
#'   structure: \code{~1} (default, pooled) requires \code{dispersion_ranef}
#'   to be a fixed scalar or a single (pooled) \code{dGamma()};
#'   \code{~<group_name>}, matching the random-effects grouping factor
#'   exactly, requires \code{dispersion_ranef} to be a \code{dGamma_list(...)}
#'   (one \code{dGamma()} per group level) or a named numeric vector (one
#'   fixed value per group level). Any other formula is an error. \code{~1}
#'   never fits an extra reference model. \code{~<group_name>} with a
#'   \code{dGamma_list(...)} additionally requires a \code{glmmTMB} reference
#'   fit (\pkg{glmmTMB} must be installed), stored as \code{dispersion_fit};
#'   \code{~<group_name>} with a fixed numeric vector never fits one, since
#'   the per-group dispersion is directly user-supplied, not a prior to
#'   calibrate. When \code{dispersion_ranef = dGamma_list(Prior_Setup_GLMM(...,
#'   dispformula = dispformula))}, that call already fit the glmmTMB
#'   reference model to calibrate the priors, and \code{lmerb()} reuses it
#'   here rather than fitting \code{glmmTMB} a second time; keep
#'   \code{dispformula} identical between the two calls, since it is not
#'   re-validated against the reused fit. \code{lmer} is always the plain
#'   \code{\link[lme4]{lmer}} fit regardless of \code{dispformula}; the
#'   sampler route (pooled vs. per-group) already follows from
#'   \code{dispersion_ranef}'s shape alone.
#' @param n Number of iid draws per group (default \code{1000L}, as in \code{\link{lmb}}).
#' @param tv_tol Total variation tolerance per stored draw, in (0, 1)
#'   (default \code{0.01}, the conventional threshold of the honest-burn-in
#'   literature; Jones and Hobert 2001).  The number of inner Gibbs sweeps
#'   per stored draw is derived so that each draw is within \code{tv_tol} of
#'   the exact joint posterior in total variation (Nygren 2020, Theorem 3;
#'   see Details).  To certify the whole \code{n}-draw sample at level
#'   \eqn{\alpha} pass \code{tv_tol = alpha / n}; the cost grows only
#'   logarithmically in \code{1/tv_tol}.
#' @param gap_tol Legacy mode--mean gap tolerance for the pilot stage when
#'   any Block~2 component uses \code{dIndependent_Normal_Gamma} and
#'   \code{tv_tol} is \code{NULL}.  When \code{tv_tol} is set (default),
#'   the pilot chain count is chosen by cost optimization instead.  Set
#'   \code{NULL} to skip the pilot unless \code{tv_tol} is set.  Ignored for
#'   all-\code{dNormal} models.
#' @param mode_gap_max Maximum per-coordinate mode--mean gap (in posterior
#'   standard deviation units) used to calibrate pilot inner sweeps when ING
#'   Block~2 components run a pilot stage (default \code{1.0}).  Ignored for
#'   all-\code{dNormal} models.
#' @param REML Logical; passed to \code{\link{model_setup}}.
#' @param control \code{\link[lme4]{lmerControl}} settings; passed to \code{model_setup}.
#' @param start Optional starting values; passed to \code{model_setup}.
#' @param verbose Verbosity flag; passed to \code{model_setup}.
#' @param subset Optional subset; passed to \code{model_setup}.
#' @param weights Optional weights; passed to \code{model_setup}.
#' @param na.action Missing-data handler; passed to \code{model_setup}.
#' @param offset Optional offset; passed to \code{model_setup}.
#' @param contrasts Optional contrasts; passed to \code{model_setup}.
#' @param devFunOnly If \code{TRUE}, return deviance function only; passed to \code{model_setup}.
#' @param simulate Logical (default \code{TRUE}). Passed through to
#'   \code{\link[lmebayesCore]{rlmerb}}. When \code{TRUE} the two-block Gibbs
#'   sampler is run for \code{n} iterations and posterior draws are stored.
#'   When \code{FALSE}, \code{rlmerb} returns fixed-VC plug-in point estimates
#'   only (exact joint Gaussian posterior mean). Simulation-only fields
#'   (\code{groupef}, \code{popef.means}, \code{popef}) are \code{NULL} when
#'   \code{simulate = FALSE}.
#' @param fixef Optional named list of hyper-parameter vectors (Block 2 state).
#'   When \code{NULL} (default), iter-0 means are taken from the
#'   \code{pfamily_list} prior means.
#' @param diag_sweeps Temporary diagnostic flag (ING models with pilot).
#'   Non-\code{dNormal} sampling already runs via
#'   \code{lmebayesCore::run_sweep_outer_chains_v6()} (R sweep-outer;
#'   pilot then main).  When \code{TRUE}, each stage auto-prints one combined
#'   Block~2 chain-mean table when that stage finishes (same layout as
#'   \code{print()} on \code{$sweep_history}).  \code{sweep_history} is
#'   collected regardless.  Default \code{FALSE}.
#' @param progbar Logical. Show text progress bars during sampling (passed to
#'   \code{\link{rlmerb}}). Default \code{NULL}: \code{FALSE} when
#'   \code{diag_sweeps = TRUE}, otherwise \code{TRUE}.
#' @param sim_method Sampling engine: \code{"DEFAULT"} or
#'   \code{"TWO_BLOCK_GIBBS"}. Only changes behavior when
#'   \code{dispersion_ranef} is fixed (a scalar or a named per-group vector)
#'   \strong{and} every \code{pfamily_list} component is \code{dNormal()}
#'   (known variance components) -- the joint posterior is then exactly
#'   multivariate normal, and \code{"DEFAULT"} draws directly, iid, from that
#'   closed form (no Gibbs sweeps, no burn-in, no autocorrelation between
#'   draws); \code{"TWO_BLOCK_GIBBS"} forces the two-block Gibbs sampler
#'   described above instead. Every other model (any
#'   \code{dIndependent_Normal_Gamma} component, or a sampled/estimated
#'   variance component) only has the two-block Gibbs engine, so both values
#'   behave identically there. See
#'   \code{\link[lmebayesCore]{rLMMNormal_reg_known_vcov}}.
#' @param ... Reserved for future use.
#' @return Object of class \code{"lmerb"}: a list with the following
#'   components (parallel to \code{\link{glmb}} and \code{\link{lmb}}).
#'
#'   Every draw slot holds one row per replicate chain, taken at that
#'   chain's final sweep, so the rows are mutually independent and can be
#'   averaged directly. Accessors \code{fixef()}, \code{ranef()},
#'   \code{coef()} and \code{summary()} are usually more convenient than
#'   reaching into these fields; see
#'   \code{vignette("Chapter-05", package = "lmebayes")}.
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{formula}}{The formula supplied.}
#'     \item{\code{lmer}}{\code{\link[lme4]{lmer}} fit from
#'       \code{model_setup} (full \code{formula}), embedded as a sub-object —
#'       analogous to \code{glmb$glm} and \code{lmb$lm}.  Use
#'       \code{coef(fit$lmer)} for per-group classical coefficients.  Always
#'       the plain pooled-dispersion \code{lmer} fit, regardless of
#'       \code{dispformula}.}
#'     \item{\code{dispformula}}{The \code{dispformula} supplied.}
#'     \item{\code{dispersion_fit}}{\code{NULL} when \code{dispformula = ~1};
#'       otherwise the diagnostic-only \code{\link[glmmTMB]{glmmTMB}} fit with
#'       per-group residual dispersion (\code{dispformula}). Identical to
#'       \code{glmmTMB} below; kept for backward compatibility.}
#'     \item{\code{glmmTMB}}{Same value as \code{dispersion_fit} (an alias):
#'       \code{NULL} when \code{dispformula = ~1}, otherwise the
#'       \code{\link[glmmTMB]{glmmTMB}} per-group-dispersion reference fit.
#'       Prefer this field name in new code; downstream consumers that need
#'       a per-group-dispersion-aware reference fit (e.g.\
#'       \code{summary.lmerb()}) use \code{glmmTMB} when non-\code{NULL} and
#'       fall back to \code{lmer} otherwise.}
#'     \item{\code{prior}}{Normalized prior container: \code{pfamily_list}
#'       (as supplied, reordered to the RE coefficient names),
#'       \code{dispersion_ranef}, the reconstructed \code{Sigma_ranef}, and
#'       the per-component \code{prior_list} (\code{mu_fixef},
#'       \code{Sigma_fixef}, \code{dispersion_fixef}) — analogous to
#'       \code{glmb$Prior}.}
#'     \item{\code{model_setup}}{The \code{\link{model_setup}} object built
#'       inside \code{lmerb} from \code{formula} and \code{data}.}
#'     \item{\code{popef.mode}}{Named list of exact posterior mode (= mean,
#'       since the joint posterior is Gaussian) vectors for the level-2 fixed
#'       effects \eqn{\gamma_k}, computed by
#'       \code{\link[lmebayesCore]{lmerb_posterior_mean}} (ICM).}
#'     \item{\code{groupef.mode}}{\eqn{J \times p_{\mathrm{re}}} numeric matrix
#'       of exact posterior mode random effects from ICM.  Rows are group
#'       levels (\code{levels(design$group)}); columns are
#'       \code{design$re_coef_names}.}
#'     \item{\code{popef}}{Named list of \eqn{n \times q_k} matrices of
#'       population draws \eqn{\gamma_k}, one per random-effect coefficient.
#'       \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{popef.means}}{Named list of posterior mean vectors, the
#'       column means of \code{popef}.  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{groupef}}{\code{data.frame} with \code{n * J} rows:
#'       \code{draw}, the grouping-factor column, and one column per RE
#'       variable.  These are full group coefficients \eqn{\beta_j},
#'       comparable to \code{lme4::coef()}, not mean-zero deviations.
#'       Average over \code{draw} within each group for posterior
#'       means (see Examples).  \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{popef.dispersion}}{\eqn{n \times p_{\mathrm{re}}} matrix of
#'       the between-group variances \eqn{\tau^2_k} at each stored draw:
#'       sampled values for \code{dIndependent_Normal_Gamma} components,
#'       constant columns for \code{dNormal} components, where the value was
#'       treated as known.  \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{popef.dispersion.mean}}{Named vector of posterior means of
#'       \eqn{\tau^2_k} (\code{colMeans(fixef.dispersion)}).  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{group.dispersion}}{Observation-level residual variance \eqn{\sigma^2}:
#'       a scalar when \code{dispersion_ranef} is fixed, a length-\code{n} vector
#'       of final-sweep draws when \code{dispersion_ranef} is \code{dGamma()};
#'       \code{NULL} when not applicable.  When \code{simulate = FALSE}, the
#'       prior plug-in scalar (\code{prior$dispersion_ranef}) is returned for
#'       Gaussian models.}
#'     \item{\code{group.dispersion.mean}}{Posterior mean of \eqn{\sigma^2}; equals
#'       \code{group.dispersion} when fixed.  \code{NULL} when \code{group.dispersion} is
#'       \code{NULL}.}
#'     \item{\code{draw_engine}}{Name of the Block~1 sampling engine used
#'       (e.g. \code{"rGLMM_sweep_ing_block1_ind"}). \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{fixef.iters}}{\eqn{n \times p_{\mathrm{re}}} matrix of the
#'       total number of Block~2 candidates generated per stored draw,
#'       summed over the \code{m_convergence} inner sweeps.
#'       \code{dIndependent_Normal_Gamma} components count envelope
#'       accept-reject candidates until acceptance (as \code{iters} in
#'       \pkg{glmbayes} samplers); \code{dNormal} components count exactly
#'       one conjugate draw per sweep.  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{fixef.iters.mean}}{Named vector of average candidates per
#'       accepted Block~2 draw
#'       (\code{colMeans(fixef.iters)/m_convergence}; equals 1 for
#'       \code{dNormal} components, and approximately the reciprocal
#'       acceptance rate of the ING envelope sampler otherwise).
#'       \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{fixef.mu}}{Numeric matrix \code{p_re x J} of Block 1 prior
#'       means at the final Gibbs state (from
#'       \code{\link[lmebayesCore]{build_mu_all}}).}
#'     \item{\code{convergence}}{List describing the sweep-count calibration:
#'       \code{method} (\code{"exact"}, \code{"exact_iid"} when
#'       \code{sim_method_used = "DEFAULT"}, or \code{"local_gaussian_mode"}
#'       for non-Gaussian \code{\link{glmerb}}), \code{tv_tol},
#'       \code{lambda_star}, \code{eigenvalues}, \code{m_min} (derived
#'       minimum sweeps), and \code{m_convergence} (sweeps actually used).
#'       \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{sim_method_used}}{\code{"DEFAULT"} (exact iid draws) or
#'       \code{"TWO_BLOCK_GIBBS"} (two-block Gibbs sweeps), whichever engine
#'       actually ran. \code{NULL} when \code{simulate = FALSE}.}
#'   }
#' @example inst/examples/Ex_lmerb.R
#' @seealso \code{\link{Prior_Setup_GLMM}}, \code{\link{model_setup}},
#'   \code{\link[lmebayesCore]{build_mu_all}},
#'   \code{\link[lmebayesCore]{two_block_rNormal_reg}},
#'   \code{\link[lmebayesCore]{lmerb_posterior_mean}},
#'   \code{\link[lmebayesCore]{rNormal_reg_group}},
#'   \code{\link{lmb}}, \code{\link{glmb}}
#' @param digits Number of significant digits to use when printing.
#' @title Fit a Bayesian linear mixed-effects model (LMM) to data, via two-Block Gibbs sampling
#' @aliases lmerb print.lmerb
#' @export
lmerb <- function(
    formula,
    data = NULL,
    pfamily_list,
    dispersion_ranef,
    dispformula = ~1,
    n = 1000L,
    tv_tol = 0.01,
    gap_tol = 0.0196,
    mode_gap_max = 1.0,
    diag_sweeps = FALSE,
    progbar = NULL,
    sim_method = "DEFAULT",
    simulate = TRUE,
    REML = TRUE,
    control = lme4::lmerControl(),
    start = NULL,
    verbose = 0L,
    subset,
    weights,
    na.action,
    offset,
    contrasts = NULL,
    devFunOnly = FALSE,
    fixef = NULL,
    ...
) {
  cl <- match.call()
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (missing(pfamily_list) || is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required. Build it with ",
      "pfamily_list(Prior_Setup_GLMM(...)) and pass the result to lmerb().",
      call. = FALSE
    )
  }
  if (missing(dispersion_ranef)) {
    stop(
      "'dispersion_ranef' is required for lmerb(). Typically ",
      "Prior_Setup_GLMM(...)$group.dispersion.",
      call. = FALSE
    )
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }
  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }
  sim_method <- lmebayesCore:::.rLMM_validate_sim_method(sim_method, fn_name = "lmerb")

  setup_args <- list(
    formula = formula,
    data = data,
    REML = REML,
    control = control,
    verbose = verbose,
    devFunOnly = devFunOnly,
    dispformula = dispformula
  )
  if (!missing(start) && !is.null(start)) {
    setup_args$start <- start
  }
  if (!missing(subset)) {
    setup_args$subset <- subset
  }
  if (!missing(weights)) {
    setup_args$weights <- weights
  }
  if (!missing(na.action)) {
    setup_args$na.action <- na.action
  }
  if (!missing(offset)) {
    setup_args$offset <- offset
  }
  if (!missing(contrasts)) {
    setup_args$contrasts <- contrasts
  }

  design <- do.call(model_setup, c(setup_args, list(...)))
  if (!inherits(design, "model_setup")) {
    stop("model_setup() must return a model_setup object.", call. = FALSE)
  }

  ## Dispersion mode only (for dispformula / glmmTMB). Full prior unpacking
  ## lives inside rlmerb().
  gd <- lmebayesCore:::.lmebayes_dispprior_list_as_group_dispersion(
    dispersion_ranef
  )
  disp_mode <- lmebayesCore:::.lmebayes_resolve_group_dispersion(
    group.dispersion = gd,
    family           = gaussian(),
    design           = design,
    fn_name          = "lmerb"
  )$mode

  dispformula_kind <- .lmebayes_validate_dispformula(
    dispformula = dispformula,
    group_name  = design$group_name,
    family      = gaussian(),
    disp_mode   = disp_mode
  )
  dispersion_fit <- NULL
  if (identical(disp_mode, "gamma_list")) {
    ## dGamma_list(Prior_Setup_GLMM(..., dispformula = dispformula))
    ## already carries its glmmTMB reference fit forward as an attribute;
    ## reuse it instead of re-fitting glmmTMB here. Failing that, the
    ## model_setup() call above already fit and stored the same reference
    ## (design$glmmTMB_fit) whenever dispformula requests per-group
    ## dispersion, so only fit a third copy if both are unavailable. A
    ## "fixed_vector" dispersion_ranef is a directly user-supplied constant,
    ## not a prior to calibrate, so it never needs a glmmTMB reference fit.
    dispersion_fit <- attr(dispersion_ranef, "group.dispersion.fit")
    if (is.null(dispersion_fit)) {
      dispersion_fit <- design$glmmTMB_fit
    }
    if (is.null(dispersion_fit)) {
      dispersion_fit <- .lmebayes_fit_glmmtmb_dispersion(
        formula           = formula,
        data              = data,
        family            = gaussian(),
        dispformula       = dispformula,
        REML              = REML,
        mer_optional_args = lmebayesCore:::.lmebayes_mer_optional_args(
          start     = start,
          subset    = subset,
          weights   = weights,
          na.action = na.action,
          offset    = offset,
          contrasts = contrasts
        )
      )
    }
  }

  lmer_fit <- design$lmer

  sampler <- lmebayesCore::rlmerb(
    n              = n,
    design         = design,
    pfamily_list   = pfamily_list,
    dispprior_list = dispersion_ranef,
    tv_tol         = tv_tol,
    progbar        = progbar,
    verbose        = TRUE,
    gap_tol        = gap_tol,
    mode_gap_max   = mode_gap_max,
    diag_sweeps    = diag_sweeps,
    sim_method     = sim_method,
    simulate       = simulate
  )

  prior <- sampler$prior
  convergence_info <- sampler$convergence
  m_convergence    <- sampler$m_convergence
  run_pilot <- !is.null(sampler$pilot$n) && sampler$pilot$n > 0L

  out <- list(
    call                  = cl,
    formula               = formula,
    lmer                  = lmer_fit,
    glmmTMB               = dispersion_fit,
    dispformula           = dispformula,
    dispersion_fit        = dispersion_fit,
    prior                 = prior,
    model_setup           = design,
    popef.mode            = sampler$popef.mode,
    popef.init            = if (!isTRUE(simulate) || run_pilot) {
      sampler$popef.init
    } else {
      NULL
    },
    groupef.mode          = sampler$groupef.mode,
    popef.means           = sampler$popef.means,
    popef                 = sampler$popef,
    groupef               = sampler$groupef,
    popef.dispersion      = sampler$popef.dispersion,
    popef.dispersion.mean = sampler$popef.dispersion.mean,
    popef.iters           = sampler$popef.iters,
    popef.iters.mean      = sampler$popef.iters.mean,
    groupef.iters         = sampler$groupef.iters,
    groupef.iters.mean    = sampler$groupef.iters.mean,
    group.dispersion      = sampler$group.dispersion,
    group.dispersion.mean = sampler$group.dispersion.mean,
    group.dispersion.iters = sampler$group.dispersion.iters,
    group.dispersion.iters.mean = sampler$group.dispersion.iters.mean,
    draw_engine           = if (!is.null(sampler$convergence_info)) {
      sampler$convergence_info$draw_engine
    } else {
      NULL
    },
    sim_method_used       = if (!is.null(sampler$convergence_info)) {
      sampler$convergence_info$sim_method_used
    } else {
      NULL
    },
    m_convergence         = m_convergence,
    pilot                 = sampler$pilot,
    gap_tol               = if (isTRUE(prior$any_non_normal)) gap_tol else NULL,
    mode_gap_max          = if (isTRUE(prior$any_non_normal)) mode_gap_max else NULL,
    convergence           = convergence_info,
    sweep_history         = list(
      pilot = if (run_pilot && !is.null(sampler$pilot$draws$sweep_history)) {
        sampler$pilot$draws$sweep_history
      } else {
        NULL
      },
      main = sampler$sweep_history
    )
  )
  if (!isTRUE(simulate)) {
    out$joint_mode <- sampler$joint_mode
    out$tau2.mode <- sampler$tau2.mode
    out$group.dispersion.mode <- sampler$group.dispersion.mode
  }
  structure(out, class = c("lmerb", "list"))
}

#' @rdname lmerb
#' @method print lmerb
#' @param x Object of class \code{"lmerb"}.
#' @param sweep_history If \code{TRUE}, print stored Block~2 sweep history
#'   (see \code{$sweep_history}) after the usual summary: one combined table
#'   per stage with mode plus all inner sweeps.
#' @param sweep_history_stage Which stage to print when \code{sweep_history =
#'   TRUE}: \code{"main"}, \code{"pilot"}, or \code{"both"}.
#' @param max_sweeps Passed to \code{print()} on sweep-history objects; limits
#'   how many inner sweeps are shown per stage.
#' @param components Optional RE components passed to sweep-history \code{print()}.
#' @param covariate Optional covariate names passed to sweep-history \code{print()}.
#' @export
print.lmerb <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    sweep_history = FALSE,
    sweep_history_stage = c("main", "pilot", "both"),
    max_sweeps = Inf,
    components = NULL,
    covariate = NULL,
    ...
) {

  re_names <- x$model_setup$groupef.names
  grp      <- x$model_setup$group_name
  n_obs    <- length(x$model_setup$y)
  n_grp    <- nlevels(x$model_setup$group)
  simulated <- !is.null(x$groupef)

  # --- Call ---
  cat("Call:\n  ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  # --- Header line ---
  if (simulated) {
    n_draws <- nrow(x$popef[[re_names[1L]]])
    cat(sprintf(
      "Bayesian linear mixed model  [%d draws, %s]\n",
      n_draws, .lmerb_engine_label(x$sim_method_used)))
  } else {
    cat("Bayesian linear mixed model  [ICM only; use simulate = TRUE for draws]\n")
  }
  cat("Formula:", deparse1(x$formula), "\n\n")

  # --- Variance components ---
  any_non_normal <- isTRUE(x$any_non_normal) || isTRUE(x$prior$any_non_normal)
  if (any_non_normal) {
    cat("Random effects (lmer reference; tau^2 sampled for non-dNormal components):\n")
  } else {
    cat("Random effects (variance components fixed at lmer estimates):\n")
  }
  vc_fit <- if (!is.null(x$glmmTMB)) x$glmmTMB else x$lmer
  print(VarCorr(vc_fit), comp = "Std.Dev.", digits = digits)
  cat(sprintf("Number of obs: %d,  groups: %s, %d\n\n", n_obs, grp, n_grp))
  if (any_non_normal && !is.null(x$popef.dispersion.mean)) {
    cat("Posterior mean tau^2_k: ",
        paste(sprintf("%s = %.4g", names(x$popef.dispersion.mean),
                      x$popef.dispersion.mean),
              collapse = ", "),
        "\n\n", sep = "")
  }

  # --- Block 2 hyperparameter reference (ICM) vs MCMC means when simulated ---
  mode_col <- if (any_non_normal) {
    cat("--- Block 2 hyperparameters (gamma at lmer tau^2 plug-in; MCMC means when simulated) ---\n\n")
    "gamma @ lmer tau2"
  } else {
    cat("--- Posterior means (ICM exact, under fixed variance components) ---\n\n")
    "fixef.mode"
  }

  rows <- do.call(rbind, lapply(re_names, function(k) {
    nms <- names(x$popef.mode[[k]])
    data.frame(
      re  = k,
      par = nms,
      mode = unname(x$popef.mode[[k]]),
      stringsAsFactors = FALSE
    )
  }))

  w_re  <- max(nchar(rows$re),  nchar("RE component"))
  w_par <- max(nchar(rows$par), nchar("parameter"))

  if (!simulated) {
    cat(sprintf("  %-*s  %-*s  %12s\n",
                w_re, "RE component", w_par, "parameter", mode_col))
    cat(sprintf("  %s  %s  %s\n",
                strrep("-", w_re), strrep("-", w_par), strrep("-", 12L)))
    for (i in seq_len(nrow(rows))) {
      cat(sprintf("  %-*s  %-*s  %12.*f\n",
                  w_re, rows$re[i], w_par, rows$par[i],
                  digits, rows$mode[i]))
    }
    cat("\n")

  } else {
    rows$means <- unlist(lapply(re_names, function(k) unname(x$popef.means[[k]])))
    rows$sd    <- unlist(lapply(re_names, function(k) {
      apply(x$popef[[k]], 2L, sd)
    }))

    cat(sprintf("  %-*s  %-*s  %12s  %12s  %10s\n",
                w_re, "RE component", w_par, "parameter",
                mode_col, "fixef.means", "draws SD"))
    cat(sprintf("  %s  %s  %s  %s  %s\n",
                strrep("-", w_re), strrep("-", w_par),
                strrep("-", 12L), strrep("-", 12L), strrep("-", 10L)))
    for (i in seq_len(nrow(rows))) {
      cat(sprintf("  %-*s  %-*s  %12.*f  %12.*f  %10.*f\n",
                  w_re, rows$re[i], w_par, rows$par[i],
                  digits, rows$mode[i],
                  digits, rows$means[i],
                  digits, rows$sd[i]))
    }
    cat("\n")
  }

  if (isTRUE(sweep_history) && !is.null(x$sweep_history)) {
    stage <- match.arg(sweep_history_stage)
    stages <- switch(
      stage,
      main = "main",
      pilot = "pilot",
      both = c("pilot", "main")
    )
    for (st in stages) {
      hist <- x$sweep_history[[st]]
      if (!is.null(hist)) {
        print(
          hist,
          max_sweeps = max_sweeps,
          components = components,
          covariate  = covariate,
          digits     = digits,
          ...
        )
      }
    }
  }

  invisible(x)
}
