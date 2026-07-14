## Re-export symbols from glmbayesCore.

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
## glmbayesCore generic.  The lmebayes_prior_setup method lives in Core;
## see ?glmbayesCore::pfamily_list.lmebayes_prior_setup.
#' @inherit glmbayesCore::pfamily_list title description details params return seealso
#' @name pfamily_list
#' @importFrom glmbayesCore pfamily_list
#' @export
glmbayesCore::pfamily_list

#' @inherit glmbayesCore::dGamma_list title description details params return seealso
#' @name dGamma_list
#' @importFrom glmbayesCore dGamma_list
#' @export
glmbayesCore::dGamma_list

#' @inherit glmbayesCore::plot_sweep_history_diag title description details params return seealso examples
#' @export
plot_sweep_history_diag <- glmbayesCore::plot_sweep_history_diag

## Replicate-chain and block Gibbs engines below rlmerb()/rglmerb() (e.g.
## rGLMM, rLMMNormal_reg, block_rNormalGLM) are glmbayesCore-only; lmebayes
## calls them with glmbayesCore:: internally.  C++ callbacks (EnvelopeOpt,
## EnvelopeSort, glmbfamfunc, rNormal_reg.wfit, rgamma_ct) resolve from the
## glmbayesCore namespace in Core; they are not re-exported from lmebayes.

#' @inherit glmbayesCore::rlmerb title description details params return seealso
#' @export
rlmerb <- glmbayesCore::rlmerb

#' @inherit glmbayesCore::rglmerb title description details params return seealso
#' @export
rglmerb <- glmbayesCore::rglmerb

#' @inherit glmbayesCore::model_setup title description details params return seealso
#' @export
model_setup <- glmbayesCore::model_setup

#' @inherit glmbayesCore::Prior_Setup_lmebayes title description details params return seealso
#' @export
Prior_Setup_lmebayes <- glmbayesCore::Prior_Setup_lmebayes
