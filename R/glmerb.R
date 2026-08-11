#' Bayesian generalized linear mixed-effects model fit
#'
#' Draft entry point for \pkg{lmebayes} GLMM models with a \code{glmer}-like
#' interface, analogous to \code{\link{lmerb}} for Gaussian responses and to
#' \code{\link{glmb}} for fixed-effects GLMs.
#'
#' It takes the same arguments as \code{\link{lmerb}} plus \code{family}.
#' With \code{family = gaussian()} it fits the same model as
#' \code{\link{lmerb}}, differing only in that the embedded reference fit
#' comes from \code{\link[lme4]{glmer}} rather than \code{\link[lme4]{lmer}}.
#'
#' @details
#' \strong{The model being estimated.}
#' For group \eqn{j = 1, \dots, J} with link \eqn{g},
#' \deqn{y_j \mid \beta_j \sim \mathrm{family}\big(g^{-1}(D_j \beta_j)\big),
#'       \qquad \beta_j \mid \gamma, \Psi \sim N(\mathcal{W}_j \gamma,\; \Psi),}
#' where \eqn{\gamma} collects the level-2 fixed effects, \eqn{\Psi} is
#' diagonal with entries \eqn{\tau^2_k} (one per random-effect coefficient),
#' and \eqn{\mathcal{W}_j} is built by \code{\link{model_setup}} from the
#' level-2 predictors.  As in \code{\link{lmerb}}, only a single grouping
#' factor and uncorrelated (\code{||}) random-effect terms are supported.
#' Marginalizing recovers the usual \code{\link[lme4]{glmer}} form with
#' \eqn{X_j = D_j \mathcal{W}_j}, \eqn{Z_j = D_j}, and mean-zero random
#' effects \eqn{b_j = \beta_j - \mathcal{W}_j \gamma}.
#'
#' Supported families are \code{gaussian()}, \code{poisson()},
#' \code{binomial()}, and \code{Gamma()} with the standard links.  What they
#' have in common is a \emph{log-concave} likelihood, which is the condition
#' the envelope accept-reject sampler requires (Nygren and Nygren 2006).
#'
#' Because the link is applied within groups, \eqn{\beta_j} and
#' \eqn{\gamma} live on the link scale --- log rates for \code{poisson()},
#' log-odds for \code{binomial()} --- so any summary on the response scale
#' must apply the inverse link to the draws rather than to their averages.
#'
#' \strong{What the priors decide.}
#' \code{pfamily_list} carries one component per random-effect coefficient
#' and declares whether the between-group variance \eqn{\tau^2_k} is known
#' (\code{dNormal}) or estimated (\code{dIndependent_Normal_Gamma}, under a
#' Gamma prior on \eqn{1/\tau^2_k} truncated to
#' \code{[disp_lower, disp_upper]}).  \code{dispersion_ranef} is the
#' residual variance and applies to \code{gaussian()} only; it must be
#' \code{NULL} for \code{poisson()} and \code{binomial()}, which have no
#' free dispersion parameter.  Both are normally built from one
#' \code{\link{Prior_Setup_GLMM}} object.
#'
#' \strong{What the returned draws are.}
#' As in \code{\link{lmerb}}, \code{glmerb} runs \code{n} independent
#' replicate chains and stores one draw from each, taken at that chain's
#' final sweep.  The draws are therefore mutually independent --- an iid
#' sample, not an autocorrelated trace --- so ordinary Monte Carlo standard
#' errors apply and there is nothing to thin or diagnose for
#' autocorrelation.
#'
#' The accompanying guarantee is weaker here than for a Gaussian response,
#' and the difference matters.  With \code{family = gaussian()} and fixed
#' variances the posterior is exactly normal and the sweep count is
#' certified.  A Poisson or binomial posterior is not normal, so the same
#' calibration is applied at local curvature near the posterior mode.  The
#' resulting \code{m_convergence} is a well-motivated sweep \emph{budget}
#' rather than a proof, and \code{tv_tol} should be read that way.
#' Tightening it is cheap: cost grows only logarithmically in
#' \eqn{1/}\code{tv_tol}.
#'
#' The same fact --- that the mode of a skewed posterior is not its mean ---
#' is why a short \strong{pilot stage} runs first.  Its cross-chain average
#' gives the main chains a better starting point than the mode does, which
#' shortens the main stage.  Pilot draws are diagnostic only and are never
#' part of the returned sample.  \code{gap_tol} and \code{mode_gap_max} are
#' legacy controls used when \code{tv_tol} is \code{NULL}.
#'
#' \strong{Point estimates (simulate = FALSE).}
#' \code{simulate = FALSE} returns the joint posterior \emph{mode} by
#' iterated conditional modes (the exact Gaussian mean when
#' \code{family = gaussian()}), at fixed variance-component plug-ins.  No
#' sweeps run and no RNG is consumed, so the result is deterministic ---
#' useful while iterating on model specification.
#'
#' @inheritParams lmerb
#' @param family A \code{\link[stats]{family}} object describing the response
#'   distribution and link. Defaults to \code{gaussian()}.
#' @param dispersion_ranef Observation-level measurement dispersion, treated
#'   as known during sampling.  For families with a dispersion parameter
#'   (e.g. \code{gaussian()}): a positive scalar, a single pooled
#'   \code{dGamma()}, or (with \code{dispformula = ~<group_name>}) a
#'   \code{dGamma_list(...)} or a named numeric vector of positive, fixed
#'   per-group values (names must match the random-effects grouping factor's
#'   levels exactly).  Must be \code{NULL} (default) for \code{poisson()} and
#'   \code{binomial()}.  Typically \code{Prior_Setup_GLMM(...)$dispersion_ranef}.
#'   Which shapes are accepted depends on \code{dispformula} (see below). The
#'   fixed numeric vector shape currently only takes effect for
#'   \code{family = gaussian()} (shared code path with \code{\link{lmerb}});
#'   observation-level dispersion is not yet wired into the sampler for the
#'   other dispersion families (\code{Gamma()}, \code{quasipoisson()},
#'   \code{quasibinomial()}).
#' @param dispformula One-sided formula selecting the measurement-dispersion
#'   structure: \code{~1} (default, pooled) requires \code{dispersion_ranef}
#'   to be a fixed scalar, \code{NULL}, or a single (pooled) \code{dGamma()};
#'   \code{~<group_name>}, matching the random-effects grouping factor
#'   exactly, requires \code{dispersion_ranef} to be a \code{dGamma_list(...)}
#'   or a named numeric vector, and errors for families without a dispersion
#'   parameter. \code{~1} never fits an extra reference model.
#'   \code{~<group_name>} with a \code{dGamma_list(...)} additionally
#'   requires a \code{glmmTMB} reference fit (\pkg{glmmTMB} must be
#'   installed), stored as \code{dispersion_fit} (reused from
#'   \code{dispersion_ranef}'s \code{"dispersion_fit"} attribute when
#'   \code{dispersion_ranef} was built via
#'   \code{dGamma_list(Prior_Setup_GLMM(..., dispformula = dispformula))},
#'   rather than re-fitting \code{glmmTMB}); \code{~<group_name>} with a fixed
#'   numeric vector never fits one, since the per-group dispersion is
#'   directly user-supplied, not a prior to calibrate. \code{glmer} is always
#'   the plain \code{\link[lme4]{glmer}} fit regardless of \code{dispformula}.
#' @param gap_tol Legacy mode--mean gap tolerance. When \code{tv_tol} is
#'   \code{NULL}, the number of pilot chains is derived as
#'   \code{ceiling((qnorm(0.975) / gap_tol)^2)} (default \code{gap_tol = 0.0196}
#'   gives \code{n_pilot = 10000}). When \code{tv_tol} is set (default),
#'   \code{n_pilot} is instead chosen by
#'   \code{\link[lmebayesCore]{two_block_optimize_pilot_cost}} to minimize total
#'   inner-sweep cost. Set \code{NULL} to skip the pilot unless \code{tv_tol}
#'   is set. Ignored for \code{gaussian()} without ING Block~2 components.
#' @param tv_tol Total variation tolerance per stored draw, in (0, 1)
#'   (default \code{0.01}).  For \code{family = gaussian()} the joint
#'   posterior is exactly multivariate normal and the number of inner Gibbs
#'   sweeps per stored draw is calibrated exactly as in \code{\link{lmerb}}
#'   (Nygren 2020, Theorem 3).  For non-Gaussian families the same
#'   calibration is applied to the \emph{local-Gaussian approximation of the
#'   posterior at its mode}: per-observation likelihood precisions are
#'   evaluated at the ICM posterior mode
#'   (\code{two_block_mode_weights()} in \pkg{lmebayesCore}) and fed to
#'   \code{\link[lmebayesCore]{two_block_rate_from_pfamily_list}}.  The derived sweep count is
#'   then the \emph{minimum} number of iterations required to converge to
#'   that hypothetical multivariate normal approximation -- a lower bound
#'   for the true (non-normal) posterior, not a guarantee.
#' @param mode_gap_max Maximum per-coordinate mode--mean gap (in posterior
#'   standard deviation units) used to calibrate pilot inner sweeps
#'   (default \code{1.0}).  Applies only to non-Gaussian families
#'   when \code{gap_tol} is not \code{NULL}.  The pilot chains start at the ICM mode, which is at
#'   Mahalanobis distance \eqn{D_{\max} = \sqrt{p}\,\times\,\texttt{mode\_gap\_max}}
#'   from the posterior mean (assuming \code{mode_gap_max} SDs per coordinate
#'   across \eqn{p} fixed-effect dimensions).  The number of pilot sweeps is
#'   the smallest \eqn{l} satisfying
#'   \eqn{\mathrm{erf}_1(0.5\,\lambda^{*l}\,D_{\max}/\sqrt{2}) \le
#'   \texttt{tv\_tol}} (Nygren 2020, Theorem 3 mean-shift term), floored at
#'   \code{m_min}.  Set \code{mode_gap_max = NULL} to fall back to the
#'   Theorem~3 minimum sweep count.  Ignored for
#'   \code{family = gaussian()} without ING Block~2 components.
#' @param control Optional \code{\link[lme4]{glmerControl}} settings passed to
#'   the reference \code{\link[lme4]{glmer}} fit. Defaults to \code{NULL}
#'   (lme4 defaults). When \code{family = gaussian()}, lme4's \code{glmer}
#'   shortcut to \code{lmer} does not accept an explicit \code{glmerControl};
#'   leave \code{control = NULL} or pass \code{\link[lme4]{lmerControl}}.
#' @param progbar Logical. When \code{TRUE}, show text progress bars during
#'   pilot and main replicate sampling inside \code{\link{rglmerb}}. Default
#'   \code{FALSE}.
#' @param sim_method Sampling engine: \code{"DEFAULT"} or
#'   \code{"TWO_BLOCK_GIBBS"}. Only changes behavior for
#'   \code{family = gaussian()} with fixed \code{dispersion_ranef} (a scalar
#'   or a named per-group vector) \strong{and} all-\code{dNormal()}
#'   \code{pfamily_list} components (known variance components) -- see
#'   \code{\link{lmerb}}'s \code{sim_method} for details. Every other
#'   \code{glmerb} model (non-Gaussian families, or any
#'   \code{dIndependent_Normal_Gamma} component, or a sampled variance
#'   component) only has the two-block Gibbs engine, so both values behave
#'   identically there.
#' @return Object of class \code{"glmerb"}: same \code{fixef.*} structure as
#'   \code{"lmerb"}, including \code{sigma2} and \code{sigma2.mean} for
#'   \code{family = gaussian()} (see \code{\link{lmerb}}), with additional
#'   \code{family}, \code{glmer} (reference
#'   \code{\link[lme4]{glmer}} fit, always the plain pooled-dispersion fit
#'   regardless of \code{dispformula}), \code{dispformula} (as supplied),
#'   \code{dispersion_fit} (\code{NULL} unless \code{dispformula} requests
#'   per-group dispersion, in which case the diagnostic-only
#'   \code{\link[glmmTMB]{glmmTMB}} fit), \code{fixef.init} (main-chain start from
#'   pilot colMeans when a pilot runs; \code{NULL} when no pilot runs),
#'   \code{pilot_chisq} (Hotelling chi-squared test of
#'   pilot mean vs ICM mode), \code{gap_tol}, and \code{mode_gap_max}.
#' @references
#' Nygren, K. (2020). \emph{On the total variation distance between multivariate
#' normal densities with applications to two-block Gibbs samplers.}
#' Unpublished manuscript.
#'
#' Nygren, K. N. and Nygren, L. M. (2006). Likelihood Subgradient Densities.
#' \emph{Journal of the American Statistical Association} \bold{101}(475),
#' 1144--1156.
#'
#' Jones, G. L. and Hobert, J. P. (2001). Honest exploration of intractable
#' probability distributions via Markov chain Monte Carlo.
#' \emph{Statistical Science} \bold{16}, 312--334.
#' @seealso \code{\link{lmerb}}, \code{\link[lmebayesCore]{glmerb_posterior_mode}},
#'   \code{\link{glmb}}; \code{\link[utils]{demo}} for the full sampling workflow
#'   (\code{demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")}).
#' @param digits Number of significant digits to use when printing.
#' @example inst/examples/Ex_glmerb.R
#' @title Fit a Bayesian generalized linear mixed-effects model (GLMM) to data, via two-Block Gibbs sampling
#' @aliases glmerb print.glmerb
#' @export
glmerb <- function(
    formula,
    data = NULL,
    family = gaussian(),
    pfamily_list,
    dispersion_ranef = NULL,
    dispformula = ~1,
    n = 1000L,
    gap_tol = 0.0196,
    mode_gap_max = 1.0,
    tv_tol = 0.01,
    simulate = TRUE,
    REML = TRUE,
    control = NULL,
    start = NULL,
    verbose = 0L,
    subset,
    weights,
    na.action,
    offset,
    contrasts = NULL,
    devFunOnly = FALSE,
    fixef = NULL,
    progbar = FALSE,
    sim_method = "DEFAULT",
    ...
) {
  cl <- match.call()
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (missing(family) || is.null(family)) {
    family <- gaussian()
  }
  if (!inherits(family, "family")) {
    stop("'family' must be a family object.", call. = FALSE)
  }
  if (missing(pfamily_list) || is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required. Build it with ",
      "pfamily_list(Prior_Setup_GLMM(...)) and pass the result to glmerb().",
      call. = FALSE
    )
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }
  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }
  sim_method <- lmebayesCore:::.rLMM_validate_sim_method(sim_method, fn_name = "glmerb")
  setup_args <- list(
    formula = formula,
    data = data,
    family = family,
    fit_mer = FALSE
  )

  design <- do.call(model_setup, c(setup_args, list(...)))
  if (!inherits(design, "model_setup")) {
    stop("model_setup() must return a model_setup object.", call. = FALSE)
  }

  ## Dispersion mode only (for dispformula / glmmTMB). Full prior unpacking
  ## lives inside rglmerb().
  gd <- lmebayesCore:::.lmebayes_dispprior_list_as_group_dispersion(
    dispersion_ranef
  )
  disp_mode <- lmebayesCore:::.lmebayes_resolve_group_dispersion(
    group.dispersion = gd,
    family           = family,
    design           = design,
    fn_name          = "glmerb"
  )$mode

  dispformula_kind <- .lmebayes_validate_dispformula(
    dispformula = dispformula,
    group_name  = design$group_name,
    family      = family,
    disp_mode   = disp_mode
  )

  mer_optional_args <- lmebayesCore:::.lmebayes_mer_optional_args(
    start = start,
    subset = subset,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts
  )

  dispersion_fit <- NULL
  if (identical(disp_mode, "gamma_list")) {
    ## dGamma_list(Prior_Setup_GLMM(..., dispformula = dispformula))
    ## already carries its glmmTMB reference fit forward as an attribute;
    ## reuse it instead of re-fitting glmmTMB here. A "fixed_vector"
    ## dispersion_ranef is a directly user-supplied constant, not a prior to
    ## calibrate, so it never needs a glmmTMB reference fit.
    dispersion_fit <- attr(dispersion_ranef, "group.dispersion.fit")
    if (is.null(dispersion_fit)) {
      dispersion_fit <- .lmebayes_fit_glmmtmb_dispersion(
        formula           = formula,
        data              = data,
        family            = family,
        dispformula       = dispformula,
        REML              = REML,
        mer_optional_args = mer_optional_args
      )
    }
  }

  glmer_args <- c(
    list(
      formula = formula,
      data = data,
      family = family,
      verbose = verbose
    ),
    if (!is.null(control)) list(control = control),
    mer_optional_args,
    list(...)
  )
  glmer_fit <- do.call(lme4::glmer, glmer_args)

  sampler <- lmebayesCore::rglmerb(
    n              = n,
    design         = design,
    pfamily_list   = pfamily_list,
    family         = family,
    dispprior_list = dispersion_ranef,
    gap_tol        = gap_tol,
    tv_tol         = tv_tol,
    mode_gap_max   = mode_gap_max,
    collect_block1 = TRUE,
    verbose        = TRUE,
    progbar        = progbar,
    sim_method     = sim_method,
    simulate       = simulate
  )

  prior <- sampler$prior
  run_pilot <- !is.null(sampler$pilot$n) && sampler$pilot$n > 0L

  out <- list(
    call                  = cl,
    formula               = formula,
    family                = family,
    glmer                 = glmer_fit,
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
    m_convergence         = sampler$m_convergence,
    pilot                 = sampler$pilot,
    gap_tol               = gap_tol,
    mode_gap_max          = mode_gap_max,
    convergence           = sampler$convergence,
    sweep_history         = list(
      pilot = if (run_pilot) sampler$pilot$draws$sweep_history else NULL,
      main  = sampler$sweep_history
    )
  )
  if (!isTRUE(simulate)) {
    out$joint_mode <- sampler$joint_mode
    out$tau2.mode <- sampler$tau2.mode
    out$group.dispersion.mode <- sampler$group.dispersion.mode
  }
  structure(out, class = c("glmerb", "list"))
}

#' Print posterior estimates by RE component for a glmerb / lmerb fit
#'
#' Displays a side-by-side table of the lmer/glmer MLE reference, posterior
#' mode or ICM mean (\code{fixef.mode}), and posterior mean (\code{fixef.means})
#' for every (RE-component, parameter) pair.  When \code{x} is a bare
#' \code{fixef.means} list rather than a full fit object, only the posterior
#' mean column is shown.
#'
#' For \code{lmerb} objects the reference column is labelled \code{"lmer"} and
#' the \code{fixef.mode} column is labelled \code{"ICM.mean"} (the Gaussian
#' posterior mean and mode coincide exactly).  For \code{glmerb} objects the
#' reference column is labelled \code{"glmer"} and the \code{fixef.mode} column
#' is labelled \code{"post.mode"}.
#'
#' @param x A \code{glmerb} or \code{lmerb} object, or a bare \code{fixef.means}
#'   list.
#' @param digits Number of decimal places for numeric columns.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @keywords internal
print_coef_means <- function(x, digits = 4L, ...) {
  is_fit    <- inherits(x, c("glmerb", "lmerb"))
  is_lmerb  <- inherits(x, "lmerb")
  cm        <- if (is_fit) .lmerb_popef_means(x) else x
  if (is.null(cm)) {
    cat("popef.means: NULL (simulation not yet run)\n")
    return(invisible(x))
  }

  rows <- do.call(rbind, lapply(names(cm), function(k) {
    v <- cm[[k]]
    data.frame(component = k, parameter = names(v),
               post_mean = unname(v), stringsAsFactors = FALSE)
  }))

  # Reference MLE column: glmer for glmerb, lmer for lmerb.
  mer_fit   <- if (is_fit) (if (is_lmerb) x$lmer else x$glmer) else NULL
  has_mer   <- !is.null(mer_fit)
  mer_label <- if (is_lmerb) "lmer" else "glmer"

  if (has_mer) {
    mer_v <- lme4::fixef(mer_fit)
    # Map (component, parameter) -> fixef name using the same convention as
    # fe_name_for() in Prior_Setup_GLMM:
    #   (Intercept) component, col X  -> fixef["X"]
    #   component K, (Intercept) col  -> fixef["K"]
    #   component K, col X            -> fixef["X:K"] or fixef["K:X"]
    rows[[mer_label]] <- mapply(function(k, col) {
      nm <- if (k == "(Intercept)") {
        col
      } else if (col == "(Intercept)") {
        k
      } else {
        cand <- c(paste0(col, ":", k), paste0(k, ":", col))
        hit  <- cand[cand %in% names(mer_v)]
        if (length(hit)) hit[1L] else NA_character_
      }
      if (!is.na(nm) && nm %in% names(mer_v)) unname(mer_v[nm]) else NA_real_
    }, rows$component, rows$parameter)
  }

  # fixef.mode column: "ICM.mean" for lmerb (exact posterior mean), "post.mode"
  # for glmerb (posterior mode from ICM optimisation).
  popef_mode <- if (is_fit) .lmerb_popef_mode(x) else NULL
  has_mode   <- !is.null(popef_mode)
  mode_label <- if (is_lmerb) "ICM.mean" else "post.mode"
  if (has_mode) {
    rows[[mode_label]] <- unlist(lapply(names(popef_mode), function(k) {
      unname(popef_mode[[k]])
    }))
  }

  w_c <- max(nchar(rows$component), nchar("RE component"))
  w_p <- max(nchar(rows$parameter),  nchar("parameter"))
  w_v <- digits + 4L

  cols <- character(0)
  if (has_mer)  cols <- c(cols, mer_label)
  if (has_mode) cols <- c(cols, mode_label)
  cols <- c(cols, "post.mean")
  n_val <- length(cols)

  # Pre-format numeric values so the outer sprintf only ever sees %s.
  num_fmt <- sprintf("%%%d.%df", w_v, digits)   # e.g. "%8.4f"
  val_hdr <- paste(formatC(cols, width = w_v, flag = " "), collapse = "  ")
  val_sep <- paste(rep(strrep("-", w_v), n_val), collapse = "  ")

  cat(sprintf("  %-*s  %-*s  %s\n", w_c, "RE component", w_p, "parameter", val_hdr))
  cat(sprintf("  %-*s  %-*s  %s\n", w_c, strrep("-", w_c), w_p, strrep("-", w_p), val_sep))
  for (i in seq_len(nrow(rows))) {
    vals <- character(0L)
    if (has_mer)  vals <- c(vals, sprintf(num_fmt, rows[[mer_label]][i]))
    if (has_mode) vals <- c(vals, sprintf(num_fmt, rows[[mode_label]][i]))
    vals <- c(vals, sprintf(num_fmt, rows$post_mean[i]))
    cat(sprintf("  %-*s  %-*s  %s\n",
                w_c, rows$component[i],
                w_p, rows$parameter[i],
                paste(vals, collapse = "  ")))
  }
  invisible(x)
}

#' @rdname glmerb
#' @method print glmerb
#' @param x Object of class \code{"glmerb"}.
#' @param sweep_history If \code{TRUE}, print stored Block~2 sweep history
#'   (see \code{$sweep_history}) after the usual summary.
#' @param sweep_history_stage Which stage to print when \code{sweep_history =
#'   TRUE}: \code{"main"}, \code{"pilot"}, or \code{"both"}.
#' @param max_sweeps Passed to \code{print()} on sweep-history objects; limits
#'   how many inner sweeps are shown.
#' @export
print.glmerb <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    sweep_history = FALSE,
    sweep_history_stage = c("main", "pilot", "both"),
    max_sweeps = Inf,
    ...
) {

  re_names <- x$model_setup$groupef.names
  grp      <- x$model_setup$group_name
  n_obs    <- length(x$model_setup$y)
  n_grp    <- nlevels(x$model_setup$group)
  popef_draws <- if (!is.null(x$popef)) x$popef else x$fixef
  popef_mode  <- .lmerb_popef_mode(x)
  popef_means <- .lmerb_popef_means(x)
  simulated <- !is.null(.lmerb_groupef_draws(x))
  fam      <- if (!is.null(x$family)) x$family$family else "gaussian"

  cat("Call:\n  ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  if (simulated) {
    n_draws <- nrow(popef_draws[[re_names[1L]]])
    cat(sprintf(
      "Bayesian generalized linear mixed model  [%s; %d draws, %s]\n",
      fam, n_draws, .lmerb_engine_label(x$sim_method_used)))
  } else {
    cat(sprintf(
      "Bayesian generalized linear mixed model  [%s; ICM only]\n", fam))
  }
  cat("Formula:", deparse1(x$formula), "\n\n")

  any_non_normal <- isTRUE(x$any_non_normal) || isTRUE(x$prior$any_non_normal)
  if (any_non_normal) {
    cat("Random effects (glmer reference; tau^2 sampled for non-dNormal components):\n")
  } else {
    cat("Random effects (variance components fixed at glmer estimates):\n")
  }
  print(lme4::VarCorr(x$glmer), comp = "Std.Dev.", digits = digits)
  cat(sprintf("Number of obs: %d,  groups: %s, %d\n\n", n_obs, grp, n_grp))
  disp_mean <- if (!is.null(x$popef.dispersion.mean)) {
    x$popef.dispersion.mean
  } else {
    x$fixef.dispersion.mean
  }
  if (any_non_normal && !is.null(disp_mean)) {
    cat("Posterior mean tau^2_k: ",
        paste(sprintf("%s = %.4g", names(disp_mean), disp_mean),
              collapse = ", "),
        "\n\n", sep = "")
  }

  mode_col <- if (any_non_normal) {
    cat("--- Block 2 hyperparameters (gamma at lmer tau^2 plug-in; MCMC means when simulated) ---\n\n")
    "gamma @ lmer tau2"
  } else if (identical(fam, "gaussian")) {
    cat("--- Posterior means (ICM exact, under fixed variance components) ---\n\n")
    "popef.mode"
  } else {
    cat("--- Block 2 hyperparameters (ICM posterior mode; MCMC means when simulated) ---\n\n")
    "ICM mode"
  }

  rows <- do.call(rbind, lapply(re_names, function(k) {
    nms <- names(popef_mode[[k]])
    data.frame(
      re  = k,
      par = nms,
      mode = unname(popef_mode[[k]]),
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
    rows$means <- unlist(lapply(re_names, function(k) unname(popef_means[[k]])))
    rows$sd    <- unlist(lapply(re_names, function(k) {
      apply(popef_draws[[k]], 2L, sd)
    }))

    cat(sprintf("  %-*s  %-*s  %12s  %12s  %10s\n",
                w_re, "RE component", w_par, "parameter",
                mode_col, "popef.means", "draws SD"))
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
        print(hist, max_sweeps = max_sweeps, digits = digits, ...)
      }
    }
  }

  invisible(x)
}
