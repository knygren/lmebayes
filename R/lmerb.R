#' Bayesian linear mixed model fit (draft)
#'
#' Entry point for \pkg{lmebayes} models with an \code{lmer}-like interface,
#' analogous to \code{\link{lmb}} and \code{\link{glmb}} for fixed-effects models.
#'
#' Calls \code{\link[lme4]{lmer}} on \code{formula} and \code{data} inside
#' \code{lmerb} so the returned \code{lmer} fit always matches the model
#' arguments passed here. Measurement priors (\code{dispersion_ranef},
#' \code{Sigma_ranef}, \code{design}, iter-0 \code{fixef}) must be supplied
#' via \code{measurement_prior_list} from a prior call to
#' \code{\link{Prior_Setup_lmebayes}}; \code{lmerb} does not run prior setup
#' internally.
#'
#' Runs a two-block Gibbs sampler for \code{n} iterations. Block 1 draws
#' group-level random effects \eqn{b_j} given the current hyper means; Block 2
#' updates the hyper means (level-2 fixed effects \eqn{\boldsymbol{\gamma}_k})
#' given the current \eqn{b_j} draw, using \code{\link{multi_rNormal_reg_v2}}
#' with the hyper design matrices from \code{design$X_hyper}.
#'
#' @details
#' \strong{Exact posterior and convergence characterisation.}
#' When variance components are treated as fixed at their \code{lmer} estimates
#' (as done here), the joint posterior over the random-effect coefficients and
#' the level-2 fixed effects is \emph{exactly} multivariate normal.  In this
#' setting the convergence of the two-block Gibbs sampler is fully
#' characterised: Corollary 1 of Nygren (2020) shows that the total variation
#' (TV) distance between the \eqn{l}-step kernel and the target density is
#' bounded above by a geometrically decreasing function of \eqn{l}, with
#' contraction rate \eqn{\lambda^* \in [0, 1)}, the maximal eigenvalue of the
#' matrix
#' \deqn{A \;=\; P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2}}
#' where \eqn{P_{11}}, \eqn{P_{22}}, and \eqn{P_{12}} are the corresponding
#' blocks of the joint precision matrix.  Because the bound is explicit and
#' computable, the required number of iterations to reach a pre-specified TV
#' tolerance \eqn{\varepsilon} can be derived analytically once \eqn{\lambda^*}
#' is known.
#'
#' \strong{\code{m_convergence}.}
#' \code{lmerb} maintains an internal constant \code{m_convergence} as a proxy
#' for this required number of iterations.  It is initialised to \code{100L}
#' and will be replaced by the formula-derived value (a function of
#' \eqn{\lambda^*} and \eqn{\varepsilon}) once Block 2 Gibbs sampling is
#' implemented and \eqn{\lambda^*} can be computed from the fitted model
#' parameters.  In the current Block-1-only implementation every draw is
#' already an exact iid draw from the conditional posterior
#' \eqn{p(b_j \mid y, \mathrm{fixef}, \sigma^2, \Sigma_b)}, so
#' \code{m_convergence} is reserved for future use.
#'
#' @references
#' Nygren, K. (2020). \emph{On the total variation distance between multivariate
#' normal densities with applications to two-block Gibbs samplers.}
#' Unpublished manuscript.
#'
#' @param formula Mixed-model formula (single grouping factor; same constraints
#'   as \code{\link{model_setup}}). Must match
#'   \code{measurement_prior_list$formula}.
#' @param data Data frame containing all variables in \code{formula}.
#' @param measurement_prior_list Required object from
#'   \code{\link{Prior_Setup_lmebayes}}. Supplies \code{design},
#'   \code{dispersion_ranef}, \code{Sigma_ranef}, and iter-0 \code{fixef}
#'   means via \code{prior_list}. Call \code{\link{Prior_Setup_lmebayes}}
#'   explicitly before \code{lmerb}.
#' @param n Number of iid draws per group (default \code{1000L}, as in \code{\link{lmb}}).
#' @param REML Logical; passed to \code{\link[lme4]{lmer}}.
#' @param control \code{\link[lme4]{lmerControl}} settings; passed to \code{lmer}.
#' @param start Optional starting values; passed to \code{lmer}.
#' @param verbose Verbosity flag; passed to \code{lmer}.
#' @param subset Optional subset; passed to \code{lmer}.
#' @param weights Optional weights; passed to \code{lmer}.
#' @param na.action Missing-data handler; passed to \code{lmer}.
#' @param offset Optional offset; passed to \code{lmer}.
#' @param contrasts Optional contrasts; passed to \code{lmer}.
#' @param devFunOnly If \code{TRUE}, return deviance function only; passed to \code{lmer}.
#' @param fixef Optional named list of hyper-parameter vectors (Block 2 state).
#'   When \code{NULL} (default), iter-0 means are taken from
#'   \code{measurement_prior_list$prior_list$mu_fixef}.
#' @param seed Optional; sets the RNG seed before sampling.
#' @param ... Reserved for future use.
#' @return Object of class \code{"lmerb"}: a list with components:
#'   \describe{
#'     \item{\code{model_setup}}{The \code{\link{model_setup}} object (from
#'       \code{measurement_prior_list$design}).}
#'     \item{\code{lmer}}{\code{\link[lme4]{lmer}} fit on \code{formula} and
#'       \code{data} as passed to \code{lmerb}. Use \code{coef(fit$lmer)} for
#'       per-group coefficients on the same scale as \code{fit$coefficients}.}
#'     \item{\code{mu_all}}{Numeric matrix \code{p_re x J} of Block 1 prior
#'       means from the final Gibbs iteration (from \code{\link{build_mu_all}}
#'       at the final \code{fixef}). Rows are \code{design$re_coef_names};
#'       columns are grouping levels.}
#'     \item{\code{fixef}}{Named list of hyper-parameter vectors: the Block 2
#'       draw from the final Gibbs iteration.}
#'     \item{\code{coefficients}}{\code{data.frame} with \code{n * J} rows and
#'       \code{2 + p_re} columns: \code{draw}, the grouping-factor column
#'       (\code{design$group_name}), and one column per random-effect variable
#'       (\code{design$re_coef_names}). Average over \code{draw} within each
#'       grouping level for factor-level posterior means (see Examples).}
#'   }
#' @examples
#' \donttest{
#'   source(system.file("examples", "Ex_lmerb.R", package = "lmebayes"))
#' }
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{build_mu_all}},
#'   \code{\link[glmbayesCore]{block_rNormalReg}},
#'   \code{\link{lmb}}, \code{\link{glmb}}
#' @export
lmerb <- function(
    formula,
    data = NULL,
    measurement_prior_list,
    n = 1000L,
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
    seed = NULL,
    ...
) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (missing(measurement_prior_list) || is.null(measurement_prior_list)) {
    stop(
      "'measurement_prior_list' is required. ",
      "Build it with Prior_Setup_lmebayes() and pass the result to lmerb().",
      call. = FALSE
    )
  }
  if (!inherits(measurement_prior_list, "lmebayes_prior_setup")) {
    stop(
      "'measurement_prior_list' must be an object from Prior_Setup_lmebayes().",
      call. = FALSE
    )
  }
  if (!identical(deparse(formula), deparse(measurement_prior_list$formula))) {
    stop(
      "'formula' does not match measurement_prior_list$formula.",
      call. = FALSE
    )
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }

  design <- measurement_prior_list$design
  if (!inherits(design, "model_setup")) {
    stop("measurement_prior_list$design must be a model_setup object.", call. = FALSE)
  }

  lmer_args <- list(
    formula = formula,
    data = data,
    REML = REML,
    control = control,
    verbose = verbose,
    devFunOnly = devFunOnly
  )
  if (!missing(start) && !is.null(start)) {
    lmer_args$start <- start
  }
  if (!missing(subset)) {
    lmer_args$subset <- subset
  }
  if (!missing(weights)) {
    lmer_args$weights <- weights
  }
  if (!missing(na.action)) {
    lmer_args$na.action <- na.action
  }
  if (!missing(offset)) {
    lmer_args$offset <- offset
  }
  if (!missing(contrasts)) {
    lmer_args$contrasts <- contrasts
  }

  lmer_fit <- do.call(lme4::lmer, c(lmer_args, list(...)))

  if (is.null(fixef)) {
    fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
    names(fixef) <- design$re_coef_names
  }

  Sigma_ranef <- measurement_prior_list$Sigma_ranef
  if (is.null(Sigma_ranef)) {
    stop("measurement_prior_list must contain 'Sigma_ranef'.", call. = FALSE)
  }
  P <- solve(Sigma_ranef)

  dispersion <- measurement_prior_list$dispersion_ranef
  if (is.null(dispersion)) {
    stop("measurement_prior_list must contain 'dispersion_ranef'.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Convergence proxy: target number of two-block Gibbs iterations required to
  # bring the sampler within a pre-specified TV distance of the exact joint
  # posterior.  When variance components are fixed the joint posterior is
  # exactly multivariate normal and the TV bound decays geometrically at rate
  # (lambda*)^l (Nygren 2020, Corollary 1), where lambda* is the maximal
  # eigenvalue of A = P11^{-1/2} P12 P22^{-1} P21 P11^{-1/2}.  Initialized to
  # 100L as a placeholder; will be derived from the model parameters once
  # lambda* is computed from fitted model parameters.
  m_convergence <- 100L

  grp_col  <- design$group_name
  re_names <- design$re_coef_names
  coef_cols <- c("draw", grp_col, re_names)

  # Block 1 argument template; prior_list$mu is updated at each Gibbs step.
  mu_all <- as.matrix(build_mu_all(design, fixef)$mu_all)
  block1_args <- list(
    n          = 1L,
    y          = design$y,
    x          = design$Z,
    block      = design$groups,
    prior_list = list(
      mu         = mu_all,
      P          = P,
      dispersion = dispersion,
      ddef       = FALSE
    )
  )

  # Block 2 prior: fixed hyperpriors from Prior_Setup_lmebayes.
  # The prior on fixef_k does not change across iterations; only the response
  # y (the Block 1 draws b_j) changes.
  block2_prior_list <- stats::setNames(
    lapply(re_names, function(k) {
      pl_k <- measurement_prior_list$prior_list[[k]]
      list(
        mu         = pl_k$mu_fixef,
        Sigma      = pl_k$Sigma_fixef,
        dispersion = pl_k$dispersion_fixef
      )
    }),
    re_names
  )

  draw_rows <- vector("list", n)

  for (i in seq_len(n)) {

    # -- Block 1: b_j | fixef, sigma^2, Sigma_b -------------------------
    mu_all <- as.matrix(build_mu_all(design, fixef)$mu_all)
    block1_args$prior_list$mu <- mu_all
    block_i <- do.call(block_rNormalReg, block1_args)
    b_i <- block_i$coefficients
    if (is.null(rownames(b_i))) {
      rownames(b_i) <- block_i$block_info$ids
    }
    colnames(b_i) <- re_names

    # -- Block 2: fixef_k | b_j, tau^2_k --------------------------------
    fixef_draw <- multi_rNormal_reg_v2(
      n          = 1L,
      y          = b_i,
      x          = design$X_hyper,
      prior_list = block2_prior_list,
      progbar    = FALSE
    )
    fixef <- stats::setNames(
      lapply(re_names, function(k) fixef_draw[[k]]$coefficients[1L, ]),
      re_names
    )

    # -- Store Block 1 draw ----------------------------------------------
    J_i <- nrow(b_i)
    draw_df <- data.frame(
      draw = rep(i, J_i),
      stringsAsFactors = FALSE
    )
    draw_df[[grp_col]] <- rownames(b_i)
    for (nm in re_names) {
      draw_df[[nm]] <- b_i[, nm]
    }
    draw_rows[[i]] <- draw_df
  }

  coefficients <- do.call(rbind, draw_rows)
  rownames(coefficients) <- NULL
  coefficients <- coefficients[, coef_cols, drop = FALSE]

  structure(
    list(
      model_setup  = design,
      lmer         = lmer_fit,
      mu_all       = mu_all,
      fixef        = fixef,
      coefficients = coefficients
    ),
    class = c("lmerb", "list")
  )
}
