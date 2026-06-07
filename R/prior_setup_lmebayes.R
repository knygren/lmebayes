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
#'           \code{ncol(X_hyper[[k]])}).  Set from \code{fixef(fit_fr)} where
#'           available; for slope-only REs with no explicit fixed effect the
#'           mean of per-group \code{coef()} values is used.}
#'         \item{\code{Sigma_fixef}}{Prior covariance matrix for \code{fixef_k}
#'           (\code{ncol x ncol}).  Scaled from \code{vcov(fit_fr)} submatrix
#'           by \eqn{(1-\mathrm{pwt})/\mathrm{pwt}}; dimensions not present in
#'           \code{vcov} use \eqn{\tau^2_k \times (1-\mathrm{pwt})/\mathrm{pwt}}.}
#'         \item{\code{dispersion_fixef}}{\eqn{\tau^2_k} scalar (equals
#'           \code{Sigma_ranef[k,k]}).  Used in Block 2 for the dNormal prior;
#'           will be replaced by \code{shape}/\code{rate} for
#'           \code{dNormal_Gamma}.}
#'       }}
#'   }
#' @seealso \code{\link{model_setup}}, \code{\link{print.lmebayes_prior_setup}}
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
  fe         <- lme4::fixef(fit_fr)
  V_fe       <- as.matrix(stats::vcov(fit_fr))
  coef_df    <- coef(fit_fr)[[grp_col]]
  coef_means <- colMeans(coef_df)

  # ------------------------------------------------------------------
  # Step 5: build prior_list for each RE k
  # ------------------------------------------------------------------

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

  # Scaling: (1-pwt)/pwt matches glmbayes::compute_gaussian_prior which uses
  # Sigma = (n_eff/n_prior) * dispersion * (X'X)^{-1} and n_eff/n_prior = (1-pwt)/pwt
  scale    <- (1 - pwt) / pwt
  re_names <- design$re_coef_names

  prior_list <- stats::setNames(
    lapply(re_names, function(k) {

      X_k    <- design$X_hyper[[k]]
      cols_k <- colnames(X_k)
      p_k    <- length(cols_k)
      tau2_k <- tau2_vec[[k]]

      # Map X_hyper column names -> lmer fixef names (NA if not in fixef)
      fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)

      # mu_fixef: use lmer fixef where available; fall back to mean of
      # per-group coef() for slope-only REs with no explicit fixed effect
      mu_fixef <- vapply(seq_len(p_k), function(i) {
        nm <- fe_nms[i]
        if (!is.na(nm) && nm %in% names(fe)) {
          unname(fe[nm])
        } else if (cols_k[i] == "(Intercept)" && k %in% names(coef_means)) {
          unname(coef_means[k])
        } else {
          0
        }
      }, numeric(1L))
      names(mu_fixef) <- cols_k

      # Sigma_fixef: vcov(fit_fr) submatrix * scale for known fixef dimensions;
      # tau2_k * scale for slope-only dimensions not present in vcov
      known       <- !is.na(fe_nms) & fe_nms %in% rownames(V_fe)
      Sigma_fixef <- matrix(0, p_k, p_k, dimnames = list(cols_k, cols_k))

      if (any(known)) {
        fe_idx <- fe_nms[known]
        Sigma_fixef[known, known] <- V_fe[fe_idx, fe_idx, drop = FALSE] * scale
      }
      if (any(!known)) {
        diag(Sigma_fixef)[!known] <- tau2_k * scale
      }

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
