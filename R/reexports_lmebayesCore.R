## Re-export symbols from glmbayesCore and lmebayesCore.

#' @inherit glmbayesCore::Prior_Setup title description details params return seealso references
#' @examples
#' \dontrun{
#' ## Full runnable examples are maintained in \pkg{glmbayesCore}:
#' example(Prior_Setup, package = "glmbayesCore", ask = FALSE, echo = TRUE)
#' }
#' @export
Prior_Setup <- glmbayesCore::Prior_Setup

#' @inherit glmbayesCore::pfamily title description details params return seealso references
#' @name pfamily
#' @aliases dNormal dNormal_Gamma dIndependent_Normal_Gamma dGamma
NULL

#' @rdname pfamily
#' @export
dNormal <- glmbayesCore::dNormal

#' @rdname pfamily
#' @export
dNormal_Gamma <- glmbayesCore::dNormal_Gamma

#' @rdname pfamily
#' @export
dIndependent_Normal_Gamma <- glmbayesCore::dIndependent_Normal_Gamma

#' @rdname pfamily
#' @export
dGamma <- glmbayesCore::dGamma

## Proper re-export (not a copy) so S3 methods register against the
## lmebayesCore generic.  The Prior_Setup_GLMM method lives in Core;
## see ?lmebayesCore::pfamily_list.Prior_Setup_GLMM.
#' @inherit lmebayesCore::pfamily_list title description details return seealso
#' @param object A prior-specification object (typically from
#'   [Prior_Setup_GLMM()]).
#' @param ptypes Character: either a single string applied to every
#'   group-effect coefficient, or a character vector / list with one
#'   string per coefficient.  `NULL` (the generic's default) lets the
#'   method choose; for [Prior_Setup_GLMM()] that resolves to
#'   `"dNormal"` (known \eqn{\tau^2_k}), and the other allowed value is
#'   `"dIndependent_Normal_Gamma"` (Gamma prior on precision
#'   \eqn{1/\tau^2_k}).  A vector may be named with the group-effect
#'   coefficient names (any order); unnamed vectors are matched
#'   positionally against `names(object$pop.prior_list)`.
#' @param ... Additional arguments passed to methods.
#' @usage pfamily_list(object, ptypes = NULL, ...)
#' @export
pfamily_list <- lmebayesCore::pfamily_list

#' @inherit lmebayesCore::dGamma_list title description details params return seealso
#' @name dGamma_list
#' @importFrom lmebayesCore dGamma_list
#' @export
lmebayesCore::dGamma_list

#' @inherit lmebayesCore::plot_mean_convergence title description details params return seealso
#' @param hist Object of class \code{"two_block_sweep_history"}, as carried
#'   by a fitted object's \code{sweep_history} (see
#'   \code{\link[lmebayesCore]{print.two_block_sweep_history}}).
#' @export
plot_mean_convergence <- lmebayesCore::plot_mean_convergence

#' @inherit lmebayesCore::plot_var_convergence title description details params return seealso
#' @export
plot_var_convergence <- lmebayesCore::plot_var_convergence

## Replicate-chain and block Gibbs engines below rlmerb()/rglmerb() (e.g.
## rGLMM, rLMMNormal_reg, rNormalGLM_reg_group) are lmebayesCore-only; lmebayes
## calls them with lmebayesCore:: internally.  C++ callbacks (EnvelopeOpt,
## EnvelopeSort, glmbfamfunc, rNormal_reg.wfit, rgamma_ct) resolve from the
## glmbayesCore namespace inside lmebayesCore's compiled code; they are not
## re-exported from lmebayes.

#' @inherit lmebayesCore::rlmerb title description details params return seealso
#' @export
rlmerb <- lmebayesCore::rlmerb

#' @inherit lmebayesCore::rglmerb title description details params return seealso
#' @export
rglmerb <- lmebayesCore::rglmerb

#' @inherit lmebayesCore::model_setup title description details params return seealso
#' @export
model_setup <- lmebayesCore::model_setup

#' @inherit lmebayesCore::Prior_Setup_GLMM title description details params return seealso
#' @export
Prior_Setup_GLMM <- lmebayesCore::Prior_Setup_GLMM

#' @inherit lmebayesCore::Prior_SetupGroup title description details params return seealso
#' @export
Prior_SetupGroup <- lmebayesCore::Prior_SetupGroup
