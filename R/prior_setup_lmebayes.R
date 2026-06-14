#' Prior setup for the two-block Gibbs lmebayes sampler
#'
#' Calibrates priors for the level-2 fixed effects (\code{fixef}) of a
#' hierarchical mixed model using the reference \code{lmer}/\code{glmer} fits
#' on \strong{all} groups (from \code{\link{model_setup}}).  Per-group design
#' rank (\code{re_rank}) is a diagnostic check only and does not subset the
#' data.  Random-effect variances are treated as fixed at their mixed-model
#' estimates.  The returned object provides all inputs needed for the
#' two-block Gibbs sampler:
#'
#' \strong{Block 1} (per-group, independent):
#' \deqn{p(\mathbf{b}_j \mid \mathbf{y}, \mathrm{fixef}, \sigma^2, \Sigma_b)
#'       = \mathcal{N}(\boldsymbol{\mu}_{b,j}^*, \boldsymbol{\Sigma}_{b,j}^*)}
#' \deqn{\boldsymbol{\Sigma}_{b,j}^{*-1}
#'       = \mathbf{Z}_j'\mathbf{Z}_j / \sigma^2
#'         + \mathrm{diag}(1/\tau^2_k)}
#' when \code{family = gaussian()}.  For non-Gaussian families there is no
#' observation-level dispersion; Block~1 uses \code{dNormal} with
#' \code{ddef = TRUE} (see \code{\link[glmbayesCore]{dNormal}}).
#'
#' \strong{Block 2} (per-RE coefficient \eqn{k}, independent):
#' \deqn{p(\mathrm{fixef}_k \mid \mathbf{b}_k, \tau^2_k)
#'       = \mathcal{N}(\boldsymbol{\mu}_{\mathrm{fixef},k}^*,
#'                     \boldsymbol{\Sigma}_{\mathrm{fixef},k}^*)}
#' \deqn{\boldsymbol{\Sigma}_{\mathrm{fixef},k}^{*-1}
#'       = \mathbf{X}_k'\mathbf{X}_k / \tau^2_k
#'         + \boldsymbol{\Sigma}_{\mathrm{fixef},k}^{-1}}
#'
#' @param formula Mixed-model formula passed to \code{\link{model_setup}},
#'   whose reference \code{lmer}/\code{glmer} fits (all groups) supply the
#'   calibration quantities.
#' @param data Data frame containing all variables in \code{formula}.
#' @param family Model \code{\link[stats]{family}}.  Default \code{gaussian()}.
#'   Non-Gaussian families use \code{\link[lme4]{glmer}} for calibration;
#'   \code{dispersion_ranef} is omitted (analogous to
#'   \code{\link[glmbayesCore]{Prior_Setup}} for flat GLMs).
#' @param pwt Prior weight(s) in \eqn{(0, 1)}.  Either a \strong{scalar}
#'   (applied to every random-effect component and every Block~2 predictor),
#'   or a \strong{list with one element per random-effect component} (named
#'   with the RE coefficient names in any order, or unnamed positional).
#'   Each list element is a scalar (recycled over that component's Block~2
#'   predictors) or a vector of length \eqn{p_k} (optionally named with the
#'   predictor column names of \code{X_hyper[[k]]}, reordered to match).
#'   The prior covariance for each \code{fixef_k} block is scaled relative to
#'   \code{vcov(fit_ref)} following the \code{\link[glmbayesCore]{Prior_Setup}}
#'   convention: \eqn{(1-\mathrm{pwt})/\mathrm{pwt}} for a scalar, and
#'   elementwise
#'   \eqn{\sqrt{(1-\mathrm{pwt}_i)/\mathrm{pwt}_i}\,
#'        \sqrt{(1-\mathrm{pwt}_j)/\mathrm{pwt}_j}} for vectors.
#' @param pwt_dispersion Optional \emph{relative} prior weight(s) in
#'   \eqn{(0, 1)} for the Block~2 dispersion (precision) prior, decoupled from
#'   \code{pwt}.  A scalar, or a list / numeric vector with one value per
#'   random-effect component (named or positional).  Converted internally to
#'   an effective prior sample size \eqn{n_k = J\,w_k/(1-w_k)} where \eqn{J}
#'   is the number of groups.  At most one of \code{pwt_dispersion} and
#'   \code{n_prior_dispersion} may be supplied; when neither is, the value
#'   is derived from \code{pwt} (the mean across a component's predictors),
#'   keeping the dispersion prior consistent with the coefficient prior
#'   strength.  Weak values carry no computational penalty for
#'   \code{dIndependent_Normal_Gamma} sampling: the \eqn{\tau^2} truncation
#'   window comes from limiting-posterior quantiles independent of the
#'   prior strength (see \code{ing_prior} below).
#' @param n_prior_dispersion Optional \emph{absolute} effective prior sample
#'   size(s) (in group units) for the Block~2 dispersion prior.  A positive
#'   scalar, or a list / numeric vector with one value per random-effect
#'   component (named or positional).  See \code{pwt_dispersion}.
#'
#' @return Object of class \code{"lmebayes_prior_setup"} with fields:
#'   \describe{
#'     \item{\code{formula}}{Model formula.}
#'     \item{\code{family}}{Family object.}
#'     \item{\code{pwt}}{Prior weight(s) used: the scalar as supplied, or the
#'       canonical named list of per-predictor weight vectors when a list was
#'       supplied.}
#'     \item{\code{pwt_dispersion}}{Named per-component vector of relative
#'       dispersion prior weights (always present; consistent with
#'       \code{n_prior_dispersion} via \eqn{w_k = n_k/(n_k + J)}).}
#'     \item{\code{n_prior_dispersion}}{Named per-component vector of
#'       effective prior sample sizes for the Block~2 dispersion prior
#'       (always present; used by
#'       \code{\link[=pfamily_list.lmebayes_prior_setup]{pfamily_list}} to
#'       calibrate \code{dIndependent_Normal_Gamma} components).}
#'     \item{\code{design}}{Full \code{\link{model_setup}} object (all groups).}
#'     \item{\code{fit_ref}}{Reference \code{lmer}/\code{glmer} fit on all
#'       groups (the full-formula fit from \code{\link{model_setup}}).}
#'     \item{\code{dispersion_ranef}}{Scalar \eqn{\sigma^2} for Gaussian models
#'       only; \code{NULL} otherwise.}
#'     \item{\code{Sigma_ranef}}{Diagonal RE covariance matrix (Block~1).}
#'     \item{\code{prior_list}}{Named Block~2 prior list per RE coefficient.}
#'     \item{\code{ing_prior}}{Named per-component list of the prospective
#'       \code{dIndependent_Normal_Gamma} calibration: Gamma precision-prior
#'       \code{shape} \eqn{= (n_0 + 1 + p_k)/2} and \code{rate}
#'       \eqn{= \hat\tau^2_k (n_0 + p_k - 1)/2} (the glmbayesCore default
#'       calibration with \eqn{n_0 =} \code{n_prior_dispersion}; since
#'       \code{rate} \eqn{= \hat\tau^2_k (\code{shape} - 1)}, the implied
#'       inverse-Gamma prior on \eqn{\tau^2_k} has mean exactly
#'       \eqn{\hat\tau^2_k}), and the default \eqn{\tau^2_k} truncation
#'       window \code{disp_lower} / \code{disp_upper}: the 0.01 / 0.99
#'       quantiles of the \emph{limiting posterior}
#'       \eqn{\Gamma((J+1)/2,\; \hat\tau^2_k (J-1)/2)} -- the weak-prior
#'       (\eqn{n_0 \to 0}) limit of the Block~2 posterior Gamma for the
#'       precision (glmbayesCore Chapter A12, Theorem 2; inverted to a
#'       \eqn{\tau^2} interval).  This window is identical for all
#'       \eqn{n_0}, covers \eqn{\ge} ~98\% of the exact posterior for every
#'       prior strength, and keeps the envelope sampler's cost stable as
#'       priors weaken; see \code{inst/ING_TRUNCATION_WINDOW.md}.  Used by
#'       \code{\link{pfamily_list}} when
#'       \code{ptypes = "dIndependent_Normal_Gamma"}; ignored for
#'       \code{dNormal} priors.}
#'   }
#' @details
#' \strong{Why default calibration depends on classical estimates.}
#' \code{Prior_Setup_lmebayes} anchors each Block~2 mean at the classical
#' mixed-model estimate and scales covariances from \code{vcov(fit_ref)} by
#' \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}.  This requires:
#' \enumerate{
#'   \item Every \code{X_hyper[[k]]} column maps to a \code{fixef(fit_ref)} term.
#'   \item Each RE variance \eqn{\tau^2_k} from \code{fit_ref} is strictly positive.
#' }
#' @seealso \code{\link{model_setup}}, \code{\link[glmbayesCore]{Prior_Setup}},
#'   \code{\link[glmbayesCore]{build_mu_all}},
#'   \code{\link{print.lmebayes_prior_setup}}
#' @export
Prior_Setup_lmebayes <- function(formula,
                                 data,
                                 family = gaussian(),
                                 pwt    = 0.01,
                                 pwt_dispersion = NULL,
                                 n_prior_dispersion = NULL,
                                 intercept_source = c("null_model", "full_model"),
                                 effects_source   = c("null_effects", "full_model")) {

  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (is.numeric(pwt)) {
    if (length(pwt) != 1L || is.na(pwt) || pwt <= 0 || pwt >= 1) {
      stop(
        "'pwt' must be a scalar in (0, 1) or a list with one element per ",
        "random-effect component.",
        call. = FALSE
      )
    }
  } else if (!is.list(pwt)) {
    stop(
      "'pwt' must be a scalar in (0, 1) or a list with one element per ",
      "random-effect component.",
      call. = FALSE
    )
  }
  if (!is.null(pwt_dispersion) && !is.null(n_prior_dispersion)) {
    stop(
      "Supply at most one of 'pwt_dispersion' and 'n_prior_dispersion'.",
      call. = FALSE
    )
  }
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "family") || is.null(family$family)) {
    stop("'family' must be a family object.", call. = FALSE)
  }

  is_gaussian <- identical(family$family, "gaussian")
  mer_label   <- if (is_gaussian) "lmer" else "glmer"

  if (is_gaussian) {
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  } else {
    # Explicit glmerControl() can crash some lme4 builds; let glmer defaults apply.
    ctrl <- NULL
  }

  design <- model_setup(
    formula = formula,
    data = data,
    family = family,
    control = ctrl
  )

  ## Full-rank status is a per-group DESIGN CHECK only (reported by print();
  ## groups with rank-deficient Z_j are still fully used below).  All
  ## calibration quantities come from the reference fits on ALL groups,
  ## already computed by model_setup(): lmer_fit/glmer_fit (full formula) for
  ## fixed effects and their covariance; the matching vcov fit for the
  ## per-component RE variances and the residual variance.
  fit_ref <- if (is_gaussian) design$lmer_fit else design$glmer_fit
  if (is.null(fit_ref)) {
    stop(
      "model_setup() did not return a reference ", mer_label, " fit.",
      call. = FALSE
    )
  }

  dispersion_ranef <- if (is_gaussian) design$residual_var else NULL
  tau2_vec         <- design$vcov_re

  p_re        <- length(design$re_coef_names)
  Sigma_ranef <- diag(unname(tau2_vec), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(design$re_coef_names, design$re_coef_names)

  fe   <- lme4::fixef(fit_ref)
  V_fe <- as.matrix(stats::vcov(fit_ref))

  fe_name_for <- function(k, col) {
    if (k == "(Intercept)") {
      if (col %in% names(fe)) col else NA_character_
    } else if (col == "(Intercept)") {
      if (k %in% names(fe)) k else NA_character_
    } else {
      cand <- c(paste0(col, ":", k), paste0(k, ":", col))
      hit  <- cand[cand %in% names(fe)]
      if (length(hit)) hit[1L] else NA_character_
    }
  }

  re_names  <- design$re_coef_names
  tau_tol   <- sqrt(.Machine$double.eps)
  re_issues <- character(0)

  for (k in re_names) {
    X_k    <- design$X_hyper[[k]]
    cols_k <- colnames(X_k)
    fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
    miss_idx <- is.na(fe_nms) | !fe_nms %in% names(fe)

    if (any(miss_idx)) {
      if (k != "(Intercept)" &&
          length(cols_k) == 1L &&
          identical(cols_k, "(Intercept)")) {
        re_issues <- c(
          re_issues,
          sprintf(
            paste0(
              "%s: random slope has no fixed main effect in ", mer_label,
              " (add '%s' to the fixed part of the formula)"
            ),
            k, k
          )
        )
      } else {
        expected_fe <- vapply(seq_along(cols_k), function(i) {
          col <- cols_k[i]
          if (k == "(Intercept)") {
            col
          } else if (col == "(Intercept)") {
            k
          } else {
            paste0(col, ":", k)
          }
        }, character(1L))
        re_issues <- c(
          re_issues,
          sprintf(
            "%s: no %s fixed effect for %s",
            k, mer_label,
            paste(expected_fe[miss_idx], collapse = ", ")
          )
        )
      }
    }

    tau2_k <- unname(tau2_vec[[k]])
    if (is.na(tau2_k) || tau2_k <= tau_tol) {
      re_issues <- c(
        re_issues,
        sprintf(
          paste0(
            "%s: random-effect variance is zero or on the boundary ",
            "(singular fit); group-level variation is not identified"
          ),
          k
        )
      )
    }
  }

  if (length(re_issues) > 0L) {
    stop(
      "Prior_Setup_lmebayes() cannot calibrate default hyperpriors:\n  - ",
      paste(re_issues, collapse = "\n  - "),
      "\n\nRevise the formula (e.g. add a fixed main effect for each random ",
      "slope and avoid RE terms with zero estimated variance), or supply ",
      "hyperpriors manually without Prior_Setup_lmebayes().",
      call. = FALSE
    )
  }

  pwt_list <- .lmebayes_resolve_pwt(pwt, design)

  J_groups <- nlevels(design$groups)
  disp     <- .lmebayes_resolve_disp_prior(
    pwt_dispersion     = pwt_dispersion,
    n_prior_dispersion = n_prior_dispersion,
    J                  = J_groups,
    re_names           = re_names,
    pwt_list           = pwt_list
  )

  # ---- null model for prior means ------------------------------------------
  # Level-2 predictors: non-intercept columns of X_hyper[["(Intercept)"]] that
  # are NOT individual-level random slopes (random slopes appear in re_names).
  random_slopes <- setdiff(re_names, "(Intercept)")
  lvl2_preds <- if ("(Intercept)" %in% re_names) {
    setdiff(
      colnames(design$X_hyper[["(Intercept)"]]),
      c("(Intercept)", random_slopes)
    )
  } else {
    character(0L)
  }

  fe_null <- fe   # default: fall back to full model
  if (intercept_source == "null_model" && length(lvl2_preds) > 0L) {
    null_formula <- formula
    for (v in lvl2_preds) {
      null_formula <- stats::update(
        null_formula,
        stats::as.formula(paste(". ~ . -", v))
      )
    }
    null_fit <- if (is_gaussian) {
      lme4::lmer(null_formula, data = data, REML = TRUE, control = ctrl)
    } else {
      lme4::glmer(null_formula, family = family, data = data)
    }
    fe_null <- lme4::fixef(null_fit)
  }

  prior_list <- stats::setNames(
    lapply(re_names, function(k) {
      X_k    <- design$X_hyper[[k]]
      cols_k <- colnames(X_k)
      p_k    <- length(cols_k)
      tau2_k <- tau2_vec[[k]]

      fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
      fe_idx <- fe_nms

      mu_fixef <- vapply(seq_len(p_k), function(i) {
        col <- cols_k[i]
        if (col == "(Intercept)") {
          nm_i <- fe_name_for(k, col)
          if (intercept_source == "null_model" && !is.na(nm_i) &&
              nm_i %in% names(fe_null)) {
            unname(fe_null[nm_i])
          } else {
            unname(fe[fe_nms[i]])
          }
        } else {
          if (effects_source == "null_effects") 0
          else unname(fe[fe_nms[i]])
        }
      }, numeric(1L))
      names(mu_fixef) <- cols_k

      ## Elementwise scaling sqrt(s_i) * sqrt(s_j) with s_i = (1-w_i)/w_i;
      ## reduces to the scalar (1-pwt)/pwt factor when all weights are equal.
      sc_k <- sqrt((1 - pwt_list[[k]]) / pwt_list[[k]])
      Sigma_fixef <- V_fe[fe_idx, fe_idx, drop = FALSE] * outer(sc_k, sc_k)
      dimnames(Sigma_fixef) <- list(cols_k, cols_k)

      list(
        mu_fixef         = mu_fixef,
        Sigma_fixef      = Sigma_fixef,
        dispersion_fixef = tau2_k
      )
    }),
    re_names
  )

  pwt_out <- if (is.numeric(pwt)) pwt else pwt_list

  ## Prospective dIndependent_Normal_Gamma calibration per component (used
  ## only when pfamily_list(ptypes = "dIndependent_Normal_Gamma") is chosen):
  ## Gamma precision prior shape/rate from the glmbayesCore default
  ## calibration (compute_gaussian_prior() with k = 1):
  ##   shape_ING = (n0 + 1 + p_k)/2,  b_0 = tau2_k * (n0 + p_k - 1)/2.
  ## Since b_0 = tau2_k * (shape_ING - 1), the implied inverse-Gamma prior on
  ## tau^2_k has mean exactly tau2_k for every n0 and p_k.
  ##
  ## The tau^2 truncation window (disp_lower / disp_upper) uses the
  ## *limiting posterior* of glmbayesCore Chapter A12, Theorem 2 -- the
  ## weak-prior (n0 -> 0) limit of the Block 2 posterior Gamma:
  ##   a_inf = (J + 1)/2,  b_inf = tau2_k * (J - 1)/2
  ## (so b_inf/(a_inf - 1) = tau2_k: mean-matched, like the prior).  The
  ## window is its central 98% mass (0.01/0.99 quantiles).  Quantiles of the
  ## *prior* would stretch without bound as n0 -> 0 (posterior coverage ->
  ## 100%, envelope acceptance -> 0); the limiting-posterior window instead
  ## has coverage >= ~98% of the exact posterior for every n0 (the finite-n0
  ## posterior is strictly more concentrated than the limit), is identical
  ## for all n0, and keeps the envelope sampler's candidates-per-draw
  ## roughly constant as priors weaken.  See inst/ING_TRUNCATION_WINDOW.md.
  ## Stored here so print() can display the window and pfamily_list()
  ## consumes one source of truth.
  ing_prior <- stats::setNames(
    lapply(re_names, function(k) {
      n0_k    <- unname(disp$n_prior_dispersion[[k]])
      p_k     <- length(prior_list[[k]]$mu_fixef)
      tau2_k  <- unname(prior_list[[k]]$dispersion_fixef)
      shape_k <- (n0_k + 1) / 2 + p_k / 2
      rate_k  <- tau2_k * (n0_k + p_k - 1) / 2
      a_inf   <- (J_groups + 1) / 2
      b_inf   <- tau2_k * (J_groups - 1) / 2
      list(
        shape      = shape_k,
        rate       = rate_k,
        disp_lower = 1 / stats::qgamma(0.99, shape = a_inf, rate = b_inf),
        disp_upper = 1 / stats::qgamma(0.01, shape = a_inf, rate = b_inf)
      )
    }),
    re_names
  )

  structure(
    list(
      formula            = formula,
      family             = family,
      pwt                = pwt_out,
      pwt_dispersion     = disp$pwt_dispersion,
      n_prior_dispersion = disp$n_prior_dispersion,
      intercept_source   = intercept_source,
      effects_source     = effects_source,
      design             = design,
      fit_ref            = fit_ref,
      dispersion_ranef   = dispersion_ranef,
      Sigma_ranef        = Sigma_ranef,
      prior_list         = prior_list,
      ing_prior          = ing_prior
    ),
    class = "lmebayes_prior_setup"
  )
}

## Resolve 'pwt' (scalar or list) into a canonical named list with one named
## numeric vector per random-effect component, ordered like X_hyper[[k]].
#' @keywords internal
#' @noRd
.lmebayes_resolve_pwt <- function(pwt, design) {

  re_names <- design$re_coef_names
  p_re     <- length(re_names)

  check_range <- function(v, what) {
    if (!is.numeric(v) || anyNA(v) || any(v <= 0) || any(v >= 1)) {
      stop(sprintf("%s must be numeric with all values in (0, 1).", what),
           call. = FALSE)
    }
  }

  if (is.numeric(pwt)) {
    check_range(pwt, "'pwt'")
    return(stats::setNames(
      lapply(re_names, function(k) {
        cols_k <- colnames(design$X_hyper[[k]])
        stats::setNames(rep(as.numeric(pwt), length(cols_k)), cols_k)
      }),
      re_names
    ))
  }

  if (length(pwt) != p_re) {
    stop(sprintf(
      paste0("'pwt' list has length %d but there are %d random-effect ",
             "components (%s)."),
      length(pwt), p_re, paste(re_names, collapse = ", ")
    ), call. = FALSE)
  }

  nms <- names(pwt)
  if (!is.null(nms) && any(nzchar(nms))) {
    if (!setequal(nms, re_names)) {
      stop(
        "Names of 'pwt' must match the random-effect coefficient names: ",
        paste(re_names, collapse = ", "),
        call. = FALSE
      )
    }
    pwt <- pwt[re_names]
  } else {
    names(pwt) <- re_names
  }

  out <- stats::setNames(vector("list", p_re), re_names)
  for (k in re_names) {
    cols_k <- colnames(design$X_hyper[[k]])
    p_k    <- length(cols_k)
    v      <- pwt[[k]]
    what   <- sprintf("'pwt[[\"%s\"]]'", k)
    check_range(v, what)

    if (length(v) == 1L) {
      v <- rep(as.numeric(v), p_k)
    } else if (length(v) == p_k) {
      vn <- names(v)
      if (!is.null(vn) && any(nzchar(vn))) {
        if (!setequal(vn, cols_k)) {
          stop(sprintf(
            "Names of %s must match the Block 2 predictors: %s.",
            what, paste(cols_k, collapse = ", ")
          ), call. = FALSE)
        }
        v <- v[cols_k]
      }
      v <- as.numeric(v)
    } else {
      stop(sprintf(
        "%s must have length 1 or %d (one value per Block 2 predictor).",
        what, p_k
      ), call. = FALSE)
    }
    out[[k]] <- stats::setNames(v, cols_k)
  }
  out
}

## Resolve the Block 2 dispersion-prior weights into mutually consistent
## per-component vectors: n_k = J * w_k / (1 - w_k)  <=>  w_k = n_k / (n_k + J).
#' @keywords internal
#' @noRd
.lmebayes_resolve_disp_prior <- function(pwt_dispersion,
                                         n_prior_dispersion,
                                         J,
                                         re_names,
                                         pwt_list) {

  p_re <- length(re_names)

  expand <- function(x, what) {
    if (is.list(x)) {
      ok <- vapply(
        x, function(e) is.numeric(e) && length(e) == 1L && !is.na(e),
        logical(1L)
      )
      if (!all(ok)) {
        stop(sprintf(
          "'%s' list elements must each be a single numeric value.", what
        ), call. = FALSE)
      }
      nms <- names(x)
      x <- vapply(x, as.numeric, numeric(1L))
      names(x) <- nms
    }
    if (!is.numeric(x) || anyNA(x)) {
      stop(sprintf("'%s' must be numeric without missing values.", what),
           call. = FALSE)
    }
    if (length(x) == 1L) {
      x <- rep(unname(x), p_re)
    } else if (length(x) == p_re) {
      nms <- names(x)
      if (!is.null(nms) && any(nzchar(nms))) {
        if (!setequal(nms, re_names)) {
          stop(sprintf(
            "Names of '%s' must match the random-effect coefficient names: %s.",
            what, paste(re_names, collapse = ", ")
          ), call. = FALSE)
        }
        x <- x[re_names]
      }
    } else {
      stop(sprintf(
        paste0("'%s' must have length 1 or %d (one value per random-effect ",
               "component)."),
        what, p_re
      ), call. = FALSE)
    }
    stats::setNames(as.numeric(x), re_names)
  }

  if (!is.null(pwt_dispersion)) {
    w <- expand(pwt_dispersion, "pwt_dispersion")
    if (any(w <= 0) || any(w >= 1)) {
      stop("'pwt_dispersion' values must be in (0, 1).", call. = FALSE)
    }
    n   <- J * w / (1 - w)
    src <- "user-supplied (pwt_dispersion)"
  } else if (!is.null(n_prior_dispersion)) {
    n <- expand(n_prior_dispersion, "n_prior_dispersion")
    if (any(n <= 0) || any(!is.finite(n))) {
      stop("'n_prior_dispersion' values must be positive and finite.",
           call. = FALSE)
    }
    w   <- n / (n + J)
    src <- "user-supplied (n_prior_dispersion)"
  } else {
    ## Default: derive from the coefficient pwt (mean across predictors per
    ## component), keeping the dispersion prior consistent with the
    ## coefficient prior strength.  This was briefly replaced by a fixed
    ## 0.2 when the ING tau^2 truncation window came from *prior*
    ## quantiles (weak priors widened the window and collapsed envelope
    ## acceptance); the window now uses limiting-posterior quantiles
    ## independent of n0 (inst/ING_TRUNCATION_WINDOW.md), so weak
    ## dispersion priors no longer carry a computational penalty.
    w <- vapply(re_names, function(k) mean(pwt_list[[k]]), numeric(1L))
    n   <- J * w / (1 - w)
    src <- "derived from pwt"
  }

  w <- stats::setNames(w, re_names)
  n <- stats::setNames(n, re_names)
  attr(w, "source") <- src
  attr(n, "source") <- src
  list(pwt_dispersion = w, n_prior_dispersion = n, source = src)
}

#' Print method for \code{lmebayes_prior_setup} objects
#'
#' @param x Object of class \code{"lmebayes_prior_setup"}.
#' @param digits Number of decimal places for numeric output.  Default 4.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.lmebayes_prior_setup <- function(x, digits = 4L, ...) {

  re_names <- x$design$re_coef_names
  n_fr     <- sum(x$design$re_rank)
  n_all    <- nlevels(x$design$groups)

  disp_src <- attr(x$pwt_dispersion, "source")

  cat("Call: Prior_Setup_lmebayes()\n\n")
  cat(sprintf("  family           : %s (%s link)\n",
              x$family$family, x$family$link))
  cat(sprintf("  intercept_source : %s\n",
              if (!is.null(x$intercept_source)) x$intercept_source else "full_model"))
  cat(sprintf("  effects_source   : %s\n",
              if (!is.null(x$effects_source)) x$effects_source else "full_model"))
  if (!is.null(x$dispersion_ranef)) {
    cat(sprintf(
      "  dispersion_ranef : %.4f  (sigma2, fixed from all %d %s)\n",
      x$dispersion_ranef, n_all, x$design$group_name
    ))
  } else {
    cat("  dispersion_ranef : NULL  (no observation-level dispersion)\n")
  }
  cat(sprintf("  Full-rank groups : %d of %d %s  (design check only)\n\n",
              n_fr, n_all, x$design$group_name))

  cat("--- Sigma_ranef (diagonal RE covariance) ---\n")
  print(round(x$Sigma_ranef, digits))

  cat("\n--- prior_list: mu_fixef / Sigma_fixef / dispersion_fixef (Block 2) ---\n")
  for (nm in re_names) {
    pl <- x$prior_list[[nm]]
    cat(sprintf("\n  [%s]\n", nm))
    pwt_k <- if (is.numeric(x$pwt)) x$pwt else x$pwt[[nm]]
    pwt_str <- if (length(unique(pwt_k)) == 1L) {
      sprintf("%.4g", pwt_k[1L])
    } else {
      paste(sprintf("%s=%.4g", names(pwt_k), pwt_k), collapse = ", ")
    }
    cat(sprintf("  pwt             : %s\n", pwt_str))
    if (!is.null(x$pwt_dispersion)) {
      cat(sprintf(
        "  pwt_disp        : %.4g  [%s]\n",
        x$pwt_dispersion[[nm]],
        if (is.null(disp_src)) "unknown" else disp_src
      ))
    }
    if (!is.null(x$n_prior_dispersion)) {
      cat(sprintf(
        "  n_prior_disp    : %.4g  (= J * pwt_disp / (1 - pwt_disp))\n",
        x$n_prior_dispersion[[nm]]
      ))
    }
    cat("  mu_fixef:\n")
    print(round(pl$mu_fixef, digits))
    cat("  Sigma_fixef:\n")
    print(round(pl$Sigma_fixef, digits))
    cat(sprintf(
      "  dispersion_fixef: %.4f  (RE variance tau^2_k; Block 2 scale)\n",
      pl$dispersion_fixef))
    ing_k <- x$ing_prior[[nm]]
    if (!is.null(ing_k)) {
      cat(sprintf(
        "  ING tau^2 window: [%.4g, %.4g]  (0.01/0.99 limiting-posterior quantiles; upper/tau2 = %.3g)\n",
        ing_k$disp_lower, ing_k$disp_upper,
        ing_k$disp_upper / unname(pl$dispersion_fixef)
      ))
      cat(sprintf(
        "  ING shape, rate : %.4g, %.4g  (Gamma prior on 1/tau^2_k; used only with ptypes = \"dIndependent_Normal_Gamma\")\n",
        ing_k$shape, ing_k$rate
      ))
    }
  }

  invisible(x)
}
