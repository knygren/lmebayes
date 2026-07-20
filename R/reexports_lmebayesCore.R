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
## lmebayesCore generic.  The lmebayes_prior_setup method lives in Core;
## see ?lmebayesCore::pfamily_list.lmebayes_prior_setup.
#' @inherit lmebayesCore::pfamily_list title description details params return seealso
#' @name pfamily_list
#' @importFrom lmebayesCore pfamily_list
#' @export
lmebayesCore::pfamily_list

#' @inherit lmebayesCore::dGamma_list title description details params return seealso
#' @name dGamma_list
#' @importFrom lmebayesCore dGamma_list
#' @export
lmebayesCore::dGamma_list

#' @inherit lmebayesCore::plot_sweep_history_diag title description details params return seealso examples
#' @export
plot_sweep_history_diag <- lmebayesCore::plot_sweep_history_diag

## Replicate-chain and block Gibbs engines below rlmerb()/rglmerb() (e.g.
## rGLMM, rLMMNormal_reg, block_rNormalGLM) are lmebayesCore-only; lmebayes
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

#' @inherit lmebayesCore::Prior_Setup_lmebayes title description details params return seealso
#' @export
Prior_Setup_lmebayes <- lmebayesCore::Prior_Setup_lmebayes
