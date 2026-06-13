#' Build pfamily objects from a Prior_Setup_lmebayes object
#'
#' Converts the per-component Block-2 hyperprior parameters stored in a
#' \code{\link{Prior_Setup_lmebayes}} object into a named list of
#' \code{\link[glmbayesCore]{pfamily}} objects, one per random-effect
#' coefficient (e.g. \code{"(Intercept)"}, slope names).
#'
#' For each random-effect coefficient \eqn{k}, the prior parameters come
#' from \code{object$prior_list[[k]]}:
#' \itemize{
#'   \item \code{"dNormal"}: \code{dNormal(mu = mu_fixef, Sigma =
#'     Sigma_fixef, dispersion = dispersion_fixef)}.  The Block-2
#'     dispersion (the random-effect variance \eqn{\tau^2_k}) is treated
#'     as known.
#'   \item \code{"dIndependent_Normal_Gamma"}: the same \code{mu} and
#'     \code{Sigma}, plus a Gamma prior on the Block-2 precision
#'     \eqn{1/\tau^2_k} calibrated with the same convention as
#'     \code{\link[glmbayesCore]{Prior_Setup}}.  The per-component
#'     effective prior sample size \eqn{n_0} is taken from
#'     \code{object$n_prior_dispersion[[k]]} (set by
#'     \code{\link{Prior_Setup_lmebayes}} via \code{pwt_dispersion} /
#'     \code{n_prior_dispersion}, derived from \code{pwt} by default).
#'     Then
#'     \deqn{shape = (n_0 + 1 + p_k)/2, \qquad
#'           rate = \tau^2_k \, (n_0 + p_k - 1)/2,}
#'     where \eqn{p_k} is the number of Block-2 coefficients for
#'     component \eqn{k} (the \code{shape_ING} convention with the
#'     glmbayesCore default rate \eqn{b_0}).  Because
#'     \eqn{rate = \tau^2_k (shape - 1)}, the implied inverse-Gamma prior
#'     on the dispersion has mean exactly \eqn{\tau^2_k} for every
#'     \eqn{n_0} and \eqn{p_k}, while small \code{pwt_dispersion} keeps it
#'     deliberately diffuse.
#'
#'     The dispersion prior must not outweigh the data: \eqn{n_0 \le J}
#'     (equivalently \code{pwt_dispersion} \eqn{\le 0.5}) is required,
#'     mirroring the sampler-side guard in
#'     \code{\link[glmbayesCore]{two_block_rNormal_reg_v2}} (the ING
#'     dispersion envelope caps its log-tilt at the data contribution
#'     \eqn{J/2}; a prior-dominated calibration would invalidate it).
#'
#'     \code{disp_lower} and \code{disp_upper} default to the 0.01 and
#'     0.99 quantiles of the \emph{limiting posterior} for \eqn{\tau^2_k}
#'     -- the weak-prior (\eqn{n_0 \to 0}) limit of the Block-2 posterior
#'     Gamma for the precision (glmbayesCore Chapter A12, Theorem 2),
#'     \eqn{\Gamma(a_\infty, b_\infty)} with
#'     \deqn{a_\infty = (J+1)/2, \qquad b_\infty = \tau^2_k\,(J-1)/2,}
#'     inverted to a \eqn{\tau^2} interval:
#'     \deqn{disp\_lower = 1 / q_{\Gamma}(0.99;\; a_\infty, b_\infty),
#'           \qquad
#'           disp\_upper = 1 / q_{\Gamma}(0.01;\; a_\infty, b_\infty).}
#'     Quantiles of the limiting posterior -- rather than of the prior --
#'     make the window independent of \eqn{n_0}: prior quantiles stretch
#'     without bound as the dispersion prior weakens (collapsing the
#'     envelope acceptance rate), whereas this window covers at least
#'     ~98\% of the exact posterior for every \eqn{n_0} (the finite-
#'     \eqn{n_0} posterior is strictly more concentrated than the limit)
#'     and keeps sampling cost stable.  See
#'     \code{inst/ING_TRUNCATION_WINDOW.md} for the derivation.  The
#'     values are computed once by \code{\link{Prior_Setup_lmebayes}}
#'     (stored in its \code{ing_prior} field and shown by its print
#'     method); this function reads them from the object.  During
#'     sampling each \eqn{\tau^2_k} draw is hard-truncated to this window
#'     (renormalized inverse-CDF), and because both bounds are supplied
#'     the window is \emph{fixed} across all inner Gibbs sweeps, rather
#'     than re-derived per sweep from a surrogate posterior.
#'     \code{\link{lmerb}} and \code{\link{glmerb}} use \code{disp_lower}
#'     as the conservative \eqn{\tau^2_k} plug-in for their eigenvalue /
#'     TV convergence calibration, so the resulting bound holds over the
#'     chain's entire (truncated) dispersion support.
#' }
#'
#' @param object An object of class \code{"lmebayes_prior_setup"} as
#'   returned by \code{\link{Prior_Setup_lmebayes}}.
#' @param ptypes Character: either a single string applied to every
#'   random-effect component, or a character vector / list with one
#'   string per component.  Allowed values are \code{"dNormal"} and
#'   \code{"dIndependent_Normal_Gamma"}.  A vector may be named with the
#'   random-effect coefficient names (any order); unnamed vectors are
#'   matched positionally against \code{names(object$prior_list)}.
#' @param ... Currently ignored.
#'
#' @return A named list of \code{"pfamily"} objects, with names equal to
#'   \code{names(object$prior_list)} (the random-effect coefficient
#'   names).
#'
#' @seealso \code{\link{Prior_Setup_lmebayes}},
#'   \code{\link[glmbayesCore]{dNormal}},
#'   \code{\link[glmbayesCore]{dIndependent_Normal_Gamma}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("bayesrules", quietly = TRUE)) {
#'   data(big_word_club, package = "bayesrules")
#'   dat <- big_word_club
#'   dat$school_id <- factor(dat$school_id)
#'   dat <- subset(dat, !is.na(score_ppvt))
#'
#'   ps <- Prior_Setup_lmebayes(
#'     score_ppvt ~ private_school + (1 | school_id),
#'     data = dat
#'   )
#'
#'   ## All components as dNormal (known Block-2 dispersion)
#'   pf1 <- pfamily_list(ps)
#'   print(pf1[["(Intercept)"]])
#'
#'   ## All components with a Gamma prior on the Block-2 precision
#'   pf2 <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
#' }
#' }
#'
#' @export
#' @method pfamily_list lmebayes_prior_setup
pfamily_list.lmebayes_prior_setup <- function(object,
                                              ptypes = "dNormal",
                                              ...) {

  allowed <- c("dNormal", "dIndependent_Normal_Gamma")

  re_names <- names(object$prior_list)
  p_re     <- length(re_names)

  ## --- validate ptypes -----------------------------------------------------
  if (is.list(ptypes)) {
    ok <- vapply(
      ptypes,
      function(p) is.character(p) && length(p) == 1L && !is.na(p),
      logical(1L)
    )
    if (!all(ok)) {
      stop("'ptypes' list elements must each be a single string.",
           call. = FALSE)
    }
    nms    <- names(ptypes)
    ptypes <- vapply(ptypes, identity, character(1L))
    names(ptypes) <- nms
  }
  if (!is.character(ptypes) || length(ptypes) < 1L || anyNA(ptypes)) {
    stop("'ptypes' must be a character vector or list of strings.",
         call. = FALSE)
  }
  bad <- setdiff(unique(ptypes), allowed)
  if (length(bad) > 0L) {
    stop(
      "Invalid 'ptypes' value(s): ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (length(ptypes) == 1L) {
    ptypes <- stats::setNames(rep(unname(ptypes), p_re), re_names)
  } else {
    if (length(ptypes) != p_re) {
      stop(
        sprintf(
          "'ptypes' has length %d but the prior setup has %d random-effect component(s): %s.",
          length(ptypes), p_re, paste(re_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(names(ptypes)) && any(nzchar(names(ptypes)))) {
      if (!setequal(names(ptypes), re_names)) {
        stop(
          "Names of 'ptypes' must match the random-effect coefficient names: ",
          paste(re_names, collapse = ", "), ".",
          call. = FALSE
        )
      }
      ptypes <- ptypes[re_names]
    } else {
      names(ptypes) <- re_names
    }
  }

  ## --- Gamma hyperparameters (shape_ING convention from glmbayesCore) ------
  ## Effective prior sample size on the Block-2 scale (per component):
  ## taken from object$n_prior_dispersion when present (new objects), else
  ## derived from pwt as n0 = J * pwt/(1-pwt), where J is the number of
  ## groups (Block-2 "observations").
  J   <- nlevels(object$design$groups)
  npd <- object$n_prior_dispersion

  n_prior_for <- function(k) {
    if (!is.null(npd)) {
      return(unname(npd[[k]]))
    }
    w <- if (is.list(object$pwt)) mean(object$pwt[[k]]) else object$pwt
    (w / (1 - w)) * J
  }

  out <- stats::setNames(vector("list", p_re), re_names)

  for (k in re_names) {
    pl    <- object$prior_list[[k]]
    mu_k  <- pl$mu_fixef
    Sig_k <- pl$Sigma_fixef
    d_k   <- unname(pl$dispersion_fixef)
    p_k   <- length(mu_k)

    out[[k]] <- switch(
      ptypes[[k]],
      dNormal = glmbayesCore::dNormal(
        mu         = mu_k,
        Sigma      = Sig_k,
        dispersion = d_k
      ),
      dIndependent_Normal_Gamma = {
        n_prior_k <- n_prior_for(k)
        ## Prior-vs-data balance guard (mirrors the two_block_rNormal_reg_v2
        ## sampler check): the Block-2 hyper-regression has J group-level
        ## observations and the ING dispersion envelope caps its log-tilt at
        ## J/2, so the dispersion prior may not carry more effective
        ## observations than the groups supply.
        if (n_prior_k > J) {
          stop(
            "Component \"", k, "\": the dispersion prior has effective ",
            "prior sample size n_prior_dispersion = ", signif(n_prior_k, 4),
            ", but there are only J = ", J, " groups. ",
            "dIndependent_Normal_Gamma sampling requires ",
            "n_prior_dispersion <= J (pwt_dispersion <= 0.5); lower ",
            "'pwt_dispersion'/'n_prior_dispersion' in Prior_Setup_lmebayes().",
            call. = FALSE
          )
        }
        ## Gamma shape/rate and the default tau^2 truncation window (central
        ## 98% mass of the limiting n0 -> 0 posterior; see
        ## inst/ING_TRUNCATION_WINDOW.md) come from object$ing_prior,
        ## computed once by Prior_Setup_lmebayes() and shown by its print()
        ## method; the formulas are re-derived here only for legacy objects
        ## that predate the ing_prior field.  Fixed bounds keep the tau^2_k
        ## truncation identical across all inner Gibbs sweeps, so the
        ## disp_lower-based lambda* upper bound holds over the chain's
        ## entire dispersion support.
        ing_k <- object$ing_prior[[k]]
        if (is.null(ing_k)) {
          shape_k <- (n_prior_k + 1) / 2 + p_k / 2
          rate_k  <- d_k * (n_prior_k + p_k - 1) / 2
          ## Truncation window from the limiting (n0 -> 0) posterior Gamma
          ## (glmbayesCore Ch. A12, Thm 2), matching Prior_Setup_lmebayes().
          a_inf <- (J + 1) / 2
          b_inf <- d_k * (J - 1) / 2
          ing_k <- list(
            shape      = shape_k,
            rate       = rate_k,
            disp_lower = 1 / stats::qgamma(0.99, shape = a_inf,
                                           rate = b_inf),
            disp_upper = 1 / stats::qgamma(0.01, shape = a_inf,
                                           rate = b_inf)
          )
        }
        glmbayesCore::dIndependent_Normal_Gamma(
          mu         = mu_k,
          Sigma      = Sig_k,
          shape      = ing_k$shape,
          rate       = ing_k$rate,
          disp_lower = ing_k$disp_lower,
          disp_upper = ing_k$disp_upper
        )
      }
    )
  }

  out
}
