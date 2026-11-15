## Re-export symbols from lmebayesCore (mixed-model setup and matrix samplers).

#' @inherit glmbayesCore::dGamma_list title description details params return seealso examples
#' @export
dGamma_list <- glmbayesCore::dGamma_list

#' @inherit glmbayesCore::plot_mean_convergence title description details params return seealso
#' @param hist Object of class \code{"two_block_sweep_history"}, as carried
#'   by a fitted object's \code{sweep_history} (see
#'   \code{\link[glmbayesCore]{print.two_block_sweep_history}}).
#' @export
plot_mean_convergence <- glmbayesCore::plot_mean_convergence

#' @inherit glmbayesCore::plot_var_convergence title description details params return seealso
#' @export
plot_var_convergence <- glmbayesCore::plot_var_convergence

## Replicate-chain and block Gibbs engines below rlmerb()/rglmerb() (e.g.
## rGLMM, rLMMNormal_reg, rNormalGLM_reg_group) are lmebayesCore-only; lmebayes
## calls them with lmebayesCore:: internally.  C++ callbacks (EnvelopeOpt,
## EnvelopeSort, glmbfamfunc, rNormal_reg.wfit, rgamma_ct) resolve from the
## glmbayesCore namespace inside lmebayesCore's compiled code; they are not
## re-exported from lmebayes.

#' @inherit glmbayesCore::rlmerb title description details params return seealso
#' @export
rlmerb <- glmbayesCore::rlmerb

#' @inherit lmebayesCore::rglmerb title description details params return seealso
#' @export
rglmerb <- lmebayesCore::rglmerb

#' @inherit glmbayesCore::model_setup title description details params return seealso examples
#' @export
model_setup <- glmbayesCore::model_setup

#' @inherit glmbayesCore::Prior_Setup_GLMM title description details params return seealso examples
#' @export
Prior_Setup_GLMM <- glmbayesCore::Prior_Setup_GLMM
