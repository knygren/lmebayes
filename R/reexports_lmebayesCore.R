## Re-export symbols from lmebayesCore (mixed-model setup and matrix samplers).

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
