#' Prior setup for the two-block Gibbs lmebayes sampler
#'
#' Calibrates priors for the level-2 fixed effects (\code{fixef}) of a
#' hierarchical mixed model using an \code{lmer}/\code{glmer} fit on the
#' full-rank subset of groups.  Random-effect variances are treated as fixed
#' at their mixed-model estimates.  The returned object provides all inputs
#' needed for the two-block Gibbs sampler:
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
#' @param formula Mixed-model formula passed to \code{\link{model_setup}} and
#'   the full-rank reference \code{lmer}/\code{glmer} fit.
#' @param data Data frame containing all variables in \code{formula}.
#' @param family Model \code{\link[stats]{family}}.  Default \code{gaussian()}.
#'   Non-Gaussian families use \code{\link[lme4]{glmer}} for calibration;
#'   \code{dispersion_ranef} is omitted (analogous to
#'   \code{\link[glmbayesCore]{Prior_Setup}} for flat GLMs).
#' @param pwt Scalar prior weight in \eqn{(0, 1)}.  The prior covariance for
#'   each \code{fixef_k} block is scaled by \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}
#'   relative to \code{vcov(fit_fr)}, matching the
#'   \code{glmbayesCore::compute_gaussian_prior} convention.
#'
#' @return Object of class \code{"lmebayes_prior_setup"} with fields:
#'   \describe{
#'     \item{\code{formula}}{Model formula.}
#'     \item{\code{family}}{Family object.}
#'     \item{\code{pwt}}{Prior weight used.}
#'     \item{\code{design}}{Full \code{\link{model_setup}} object (all groups).}
#'     \item{\code{fit_fr}}{Reference \code{lmer}/\code{glmer} fit on full-rank
#'       groups only.}
#'     \item{\code{dispersion_ranef}}{Scalar \eqn{\sigma^2} for Gaussian models
#'       only; \code{NULL} otherwise.}
#'     \item{\code{Sigma_ranef}}{Diagonal RE covariance matrix (Block~1).}
#'     \item{\code{prior_list}}{Named Block~2 prior list per RE coefficient.}
#'   }
#' @details
#' \strong{Why default calibration depends on classical estimates.}
#' \code{Prior_Setup_lmebayes} anchors each Block~2 mean at the classical
#' mixed-model estimate and scales covariances from \code{vcov(fit_fr)} by
#' \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}.  This requires:
#' \enumerate{
#'   \item Every \code{X_hyper[[k]]} column maps to a \code{fixef(fit_fr)} term.
#'   \item Each RE variance \eqn{\tau^2_k} from \code{fit_fr} is strictly positive.
#' }
#' @seealso \code{\link{model_setup}}, \code{\link[glmbayesCore]{Prior_Setup}},
#'   \code{\link[glmbayesCore]{build_mu_all}},
#'   \code{\link{print.lmebayes_prior_setup}}
#' @export
Prior_Setup_lmebayes <- function(formula,
                                 data,
                                 family = gaussian(),
                                 pwt    = 0.01) {

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(pwt) || length(pwt) != 1L || pwt <= 0 || pwt >= 1) {
    stop("'pwt' must be a scalar in (0, 1).", call. = FALSE)
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
    control = if (is_gaussian) ctrl else lme4::lmerControl()
  )

  fr_levs  <- names(design$re_rank)[design$re_rank]
  grp_col  <- design$group_name
  grp_vals <- as.character(data[[grp_col]])
  data_fr  <- data[grp_vals %in% fr_levs, , drop = FALSE]
  data_fr[[grp_col]] <- droplevels(factor(data_fr[[grp_col]]))

  fit_fr <- if (is_gaussian) {
    lme4::lmer(formula, data = data_fr, control = ctrl)
  } else {
    lme4::glmer(formula, data = data_fr, family = family)
  }

  vc_fr <- extract_mer_variance_components(fit_fr, design$re_coef_names)
  dispersion_ranef <- if (is_gaussian) vc_fr$residual_var else NULL
  tau2_vec         <- vc_fr$vcov_re

  p_re        <- length(design$re_coef_names)
  Sigma_ranef <- diag(unname(tau2_vec), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(design$re_coef_names, design$re_coef_names)

  fe   <- lme4::fixef(fit_fr)
  V_fe <- as.matrix(stats::vcov(fit_fr))

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

  scale <- (1 - pwt) / pwt

  prior_list <- stats::setNames(
    lapply(re_names, function(k) {
      X_k    <- design$X_hyper[[k]]
      cols_k <- colnames(X_k)
      p_k    <- length(cols_k)
      tau2_k <- tau2_vec[[k]]

      fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
      fe_idx <- fe_nms

      mu_fixef <- vapply(seq_len(p_k), function(i) {
        unname(fe[fe_nms[i]])
      }, numeric(1L))
      names(mu_fixef) <- cols_k

      Sigma_fixef <- V_fe[fe_idx, fe_idx, drop = FALSE] * scale
      dimnames(Sigma_fixef) <- list(cols_k, cols_k)

      list(
        mu_fixef         = mu_fixef,
        Sigma_fixef      = Sigma_fixef,
        dispersion_fixef = tau2_k
      )
    }),
    re_names
  )

  structure(
    list(
      formula          = formula,
      family           = family,
      pwt              = pwt,
      design           = design,
      fit_fr           = fit_fr,
      dispersion_ranef = dispersion_ranef,
      Sigma_ranef      = Sigma_ranef,
      prior_list       = prior_list
    ),
    class = "lmebayes_prior_setup"
  )
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

  cat("Call: Prior_Setup_lmebayes()\n\n")
  cat(sprintf("  family           : %s (%s link)\n",
              x$family$family, x$family$link))
  cat(sprintf("  pwt              : %.4g\n", x$pwt))
  if (!is.null(x$dispersion_ranef)) {
    cat(sprintf(
      "  dispersion_ranef : %.4f  (sigma2, fixed from %d full-rank %s)\n",
      x$dispersion_ranef, n_fr, x$design$group_name
    ))
  } else {
    cat("  dispersion_ranef : NULL  (no observation-level dispersion)\n")
  }
  cat(sprintf("  Full-rank groups : %d of %d %s\n\n",
              n_fr, n_all, x$design$group_name))

  cat("--- Sigma_ranef (diagonal RE covariance) ---\n")
  print(round(x$Sigma_ranef, digits))

  cat("\n--- prior_list: mu_fixef / Sigma_fixef / dispersion_fixef (Block 2) ---\n")
  for (nm in re_names) {
    pl <- x$prior_list[[nm]]
    cat(sprintf("\n  [%s]\n", nm))
    cat("  mu_fixef:\n")
    print(round(pl$mu_fixef, digits))
    cat("  Sigma_fixef:\n")
    print(round(pl$Sigma_fixef, digits))
    cat(sprintf(
      "  dispersion_fixef: %.4f  (RE variance tau^2_k; Block 2 scale)\n",
      pl$dispersion_fixef))
  }

  invisible(x)
}
