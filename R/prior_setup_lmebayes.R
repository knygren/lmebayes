#' Prior setup for the two-block Gibbs lmebayes sampler
#'
#' Calibrates priors for the level-2 fixed effects (\code{fixef}) of a
#' hierarchical linear mixed model using a \code{lmer} fit on the full-rank
#' subset of groups.  Variance components are treated as fixed at their lmer
#' estimates (Gaussian / \code{dNormal} analog).  The returned object provides
#' all inputs needed for the two-block Gibbs sampler:
#'
#' \strong{Block 1} (per-group, independent):
#' \deqn{p(\mathbf{b}_j \mid \mathbf{y}, \mathrm{fixef}, \sigma^2, \Sigma_b)
#'       = \mathcal{N}(\boldsymbol{\mu}_{b,j}^*, \boldsymbol{\Sigma}_{b,j}^*)}
#' \deqn{\boldsymbol{\Sigma}_{b,j}^{*-1}
#'       = \mathbf{Z}_j'\mathbf{Z}_j / \sigma^2
#'         + \mathrm{diag}(1/\tau^2_k)}
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
#'   \code{\link[lme4]{lmer}}.
#' @param data Data frame containing all variables in \code{formula}.
#' @param family Model family.  Default \code{gaussian()}.  Reserved for future
#'   extension to \code{dNormal_Gamma} and \code{dIndependent_Normal_Gamma};
#'   only Gaussian is implemented in this version.
#' @param pwt Scalar prior weight in \eqn{(0, 1)}.  The prior covariance for
#'   each \code{fixef_k} block is scaled by \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}
#'   relative to \code{vcov(lmer_fit)}, matching the
#'   \code{glmbayes::compute_gaussian_prior} convention.  A small \code{pwt}
#'   (default 0.01) gives a diffuse prior; a large \code{pwt} gives a
#'   more informative prior centred on the lmer estimates.
#'
#' @return Object of class \code{"lmebayes_prior_setup"} with fields:
#'   \describe{
#'     \item{\code{formula}}{Model formula.}
#'     \item{\code{family}}{Family object.}
#'     \item{\code{pwt}}{Prior weight used.}
#'     \item{\code{design}}{Full \code{\link{model_setup}} object (all groups).
#'       Contains \code{re_rank}, \code{hyper_rank}, \code{rank_ok}, \code{Z},
#'       \code{X_hyper}, etc.}
#'     \item{\code{fit_fr}}{lmer fit on full-rank groups only.}
#'     \item{\code{dispersion_ranef}}{Scalar \eqn{\sigma^2}: residual variance
#'       fixed at the full-rank lmer estimate.  Used in Block 1.}
#'     \item{\code{Sigma_ranef}}{Diagonal \eqn{p_\mathrm{re} \times
#'       p_\mathrm{re}} matrix with \eqn{\tau^2_k} on the diagonal (one per RE
#'       coefficient).  Used in Block 1.  Off-diagonals are zero under the
#'       \code{||} zero-correlation structure; will become a full matrix when
#'       correlated RE (\code{|}) are supported.}
#'     \item{\code{prior_list}}{Named list, one entry per RE coefficient \eqn{k}
#'       (names match \code{design$re_coef_names}), each containing:
#'       \describe{
#'         \item{\code{mu_fixef}}{Prior mean vector for \code{fixef_k} (length
#'           \code{ncol(X_hyper[[k]])}).  Set from \code{fixef(fit_fr)} for every
#'           hyper parameter; see Details.}
#'         \item{\code{Sigma_fixef}}{Prior covariance matrix for \code{fixef_k}
#'           (\code{ncol x ncol}).  Scaled from \code{vcov(fit_fr)} submatrix
#'           by \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}.}
#'         \item{\code{dispersion_fixef}}{\eqn{\tau^2_k} scalar (equals
#'           \code{Sigma_ranef[k,k]}).  Used in Block 2 for the dNormal prior;
#'           will be replaced by \code{shape}/\code{rate} for
#'           \code{dNormal_Gamma}.}
#'       }}
#'   }
#' @details
#' \strong{Why default calibration depends on classical estimates.}
#' \code{Prior_Setup_lmebayes} anchors the prior mean for each Block 2
#' parameter vector \eqn{\mathrm{fixef}_k} at the classical \code{lmer}
#' estimate and scales the prior covariance from \code{vcov(fit_fr)} by
#' \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}.  This data-driven approach requires
#' the classical estimates to be well defined, which imposes two conditions on
#' each random-effect coefficient \eqn{k}:
#'
#' \enumerate{
#'   \item \strong{Fixed-effect requirement.}  Every column of
#'     \code{X_hyper[[k]]} must correspond to an entry in
#'     \code{fixef(fit_fr)}.
#'     \itemize{
#'       \item For the random intercept (\code{"(Intercept)"}), each level-2
#'         covariate in the hyper model (e.g., \code{private_school},
#'         \code{title1}) must appear as a fixed main effect in the formula.
#'       \item For each random slope \eqn{k}, the formula must include a fixed
#'         main effect for \eqn{k} (its population mean \eqn{\gamma_{10}}).
#'         If cross-level moderation is specified via
#'         \code{moderator:k} (e.g., \code{free_reduced_lunch:distracted_a1}),
#'         that interaction term must also appear in the fixed part so that
#'         \code{lmer} produces an estimate for it.
#'     }
#'     A slope that appears only in the random part---without a matching fixed
#'     main effect---violates this requirement.  \code{lmer} can fit such
#'     models (the slope has no fixed population mean), but the classical
#'     \code{fixef} vector contains no value to anchor the prior, so
#'     automatic calibration is not possible.
#'
#'   \item \strong{Non-singularity requirement.}  The estimated random-effect
#'     variance \eqn{\tau^2_k} from \code{fit_fr} must be strictly positive.
#'     \eqn{\tau^2_k} plays two roles: as the \eqn{k}-th diagonal entry of
#'     \code{Sigma_ranef} (Block 1 prior precision on \eqn{b_j[k]}) and as
#'     \code{dispersion_fixef} (the scale linking Block 1 and Block 2 in the
#'     full Gibbs).  A boundary estimate (\eqn{\tau^2_k = 0}) means the
#'     \code{lmer} fit found no evidence of group-level variation in
#'     coefficient \eqn{k}: Block 1 would collapse to a point-mass prior, and
#'     Block 2 would receive a degenerate likelihood.
#' }
#'
#' \strong{Models where default calibration fails.}
#' Both conditions encode the requirement that the classical mixed-model
#' estimates are well defined for the specified formula.  When either fails,
#' \code{Prior_Setup_lmebayes()} stops with a diagnostic message explaining
#' which coefficients are problematic.  The usual remedies are to revise the
#' formula (add the missing fixed main effect; remove or reparameterise an RE
#' term with zero estimated variance).
#'
#' Such models can still be fitted with \code{\link{lmerb}}, but the user
#' must supply the measurement prior directly, without calling
#' \code{Prior_Setup_lmebayes}.  The user-constructed object must contain
#' \code{dispersion_ranef}, \code{Sigma_ranef}, and a \code{prior_list} (named
#' by \code{re_coef_names}, each element providing \code{mu_fixef},
#' \code{Sigma_fixef}, and \code{dispersion_fixef}) based on substantive
#' knowledge or an alternative calibration strategy, and must be passed as
#' \code{measurement_prior_list} to \code{\link{lmerb}}.
#' @seealso \code{\link{model_setup}}, \code{\link{build_mu_all}},
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
    family <- get(family, mode = "function", envir = parent.frame())()
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!identical(family$family, "gaussian")) {
    stop("Only Gaussian family is implemented in this version.", call. = FALSE)
  }

  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

  # ------------------------------------------------------------------
  # Step 1: model_setup on full data -- populates re_rank, X_hyper, Z
  # ------------------------------------------------------------------
  design <- model_setup(formula, data = data, control = ctrl)

  # ------------------------------------------------------------------
  # Step 2: refit lmer on full-rank groups only
  # ------------------------------------------------------------------
  fr_levs  <- names(design$re_rank)[design$re_rank]
  grp_col  <- design$group_name
  grp_vals <- as.character(data[[grp_col]])
  data_fr  <- data[grp_vals %in% fr_levs, , drop = FALSE]
  data_fr[[grp_col]] <- droplevels(factor(data_fr[[grp_col]]))

  fit_fr <- lme4::lmer(formula, data = data_fr, control = ctrl)

  # ------------------------------------------------------------------
  # Step 3: variance components from the full-rank fit (treated as fixed)
  # ------------------------------------------------------------------
  vc_fr            <- extract_lmer_variance_components(fit_fr, design$re_coef_names)
  dispersion_ranef <- vc_fr$residual_var    # sigma2  (scalar)
  tau2_vec         <- vc_fr$vcov_re         # named numeric vector, one per RE

  # Sigma_ranef: diagonal p_re x p_re matrix with tau2_k on the diagonal.
  # Off-diagonal elements are zero (|| zero-correlation structure).
  # When extended to full | correlated RE, this becomes a full covariance matrix.
  p_re        <- length(design$re_coef_names)
  Sigma_ranef <- diag(unname(tau2_vec), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(design$re_coef_names, design$re_coef_names)

  # ------------------------------------------------------------------
  # Step 4: fixed effects and their vcov from the full-rank fit
  # ------------------------------------------------------------------
  fe   <- lme4::fixef(fit_fr)
  V_fe <- as.matrix(stats::vcov(fit_fr))

  # Helper: map one X_hyper column name to its lmer fixef name given RE k
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

  # ------------------------------------------------------------------
  # Step 4b: calibration checks (each RE needs fixef + positive tau^2)
  # ------------------------------------------------------------------
  re_names  <- design$re_coef_names
  tau_tol   <- sqrt(.Machine$double.eps)
  re_issues <- character(0)

  for (k in re_names) {
    X_k    <- design$X_hyper[[k]]
    cols_k <- colnames(X_k)
    fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
    miss_idx <- is.na(fe_nms) | !fe_nms %in% names(fe)

    if (any(miss_idx)) {
      miss_idx <- is.na(fe_nms) | !fe_nms %in% names(fe)
      if (k != "(Intercept)" &&
          length(cols_k) == 1L &&
          identical(cols_k, "(Intercept)")) {
        # Hyper ~ 1 for a random slope: the missing piece is fixef[k], not
        # a literal (Intercept) fixed effect.
        re_issues <- c(
          re_issues,
          sprintf(
            paste0(
              "%s: random slope has no fixed main effect in lmer ",
              "(add '%s' to the fixed part of the formula)"
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
            "%s: no lmer fixed effect for %s",
            k,
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
            "(singular fit); school-level variation is not identified"
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

  # ------------------------------------------------------------------
  # Step 5: build prior_list for each RE k
  # ------------------------------------------------------------------

  # Scaling: (1-pwt)/pwt matches glmbayes::compute_gaussian_prior which uses
  # Sigma = (n_eff/n_prior) * dispersion * (X'X)^{-1} and n_eff/n_prior = (1-pwt)/pwt
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
  w        <- max(nchar(re_names))

  cat("Call: Prior_Setup_lmebayes()\n\n")
  cat(sprintf("  pwt              : %.4g\n", x$pwt))
  cat(sprintf("  dispersion_ranef : %.4f  (sigma2, fixed from %d full-rank %s)\n",
              x$dispersion_ranef, n_fr, x$design$group_name))
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
      "  dispersion_fixef: %.4f  (dNormal; replaced by shape/rate for dNormal_Gamma)\n",
      pl$dispersion_fixef))
  }

  invisible(x)
}
