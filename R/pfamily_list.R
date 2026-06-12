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
#'     \code{n_prior_dispersion}, or derived from \code{pwt} as
#'     \eqn{n_0 = J \cdot pwt/(1 - pwt)} with \eqn{J} groups).  Then
#'     \deqn{shape = (n_0 + 1)/2 + p_k/2, \qquad
#'           rate = \tau^2_k \, n_0/2,}
#'     where \eqn{p_k} is the number of Block-2 coefficients for
#'     component \eqn{k} (the \code{shape_ING} convention).  The prior
#'     mean of the precision is \eqn{(n_0 + 1 + p_k)/(n_0 \tau^2_k)},
#'     i.e. centered near \eqn{1/\tau^2_k} when \eqn{n_0} is moderate and
#'     deliberately diffuse for small \code{pwt}.
#'
#'     \code{disp_lower} defaults to the 0.01 quantile of the implied
#'     inverse-Gamma dispersion prior, i.e. the value below which the
#'     dispersion falls with prior probability 1\%:
#'     \deqn{disp\_lower = 1 / q_{\Gamma}(0.99;\; shape,\; rate),}
#'     equivalently the reciprocal of the 99th percentile of the Gamma
#'     precision prior.  \code{\link{lmerb}} and \code{\link{glmerb}} use
#'     this value as the conservative \eqn{\tau^2_k} plug-in for their
#'     eigenvalue / TV convergence calibration, so the resulting bound
#'     holds over 99\% of the prior dispersion mass.  Note that with the
#'     diffuse default calibration (small \code{pwt}) this quantile can be
#'     far below \eqn{\hat\tau^2_k}, giving conservative (large) sweep
#'     counts.
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
        shape_k <- (n_prior_k + 1) / 2 + p_k / 2
        rate_k  <- d_k * (n_prior_k / 2)
        glmbayesCore::dIndependent_Normal_Gamma(
          mu    = mu_k,
          Sigma = Sig_k,
          shape = shape_k,
          rate  = rate_k,
          ## Default lower truncation: 0.01 quantile of the inverse-Gamma
          ## dispersion prior (= 1 / 99th percentile of the Gamma precision).
          disp_lower = 1 / stats::qgamma(0.99, shape = shape_k, rate = rate_k)
        )
      }
    )
  }

  out
}
