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
## glmbayesCore generic.
#' @importFrom glmbayesCore pfamily_list
#' @export
glmbayesCore::pfamily_list

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
