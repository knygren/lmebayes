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
## rGLMM, rLMMNormal_reg, block_rNormalGLM, run_sweep_outer_chains_v6) are
## glmbayesCore-only; lmebayes calls them with glmbayesCore:: internally.

## C++ callback symbols re-exported for search-path lookup when lmebayes is
## attached. glmbayesCore C++ resolves these by unqualified name
## (Rcpp::Function("EnvelopeOpt"), etc.).

#' C++ callback symbols re-exported for search-path lookup
#'
#' Re-exported from \pkg{glmbayesCore} so compiled code can resolve them by
#' unqualified name when \pkg{lmebayes} is attached. Not part of the
#' user-facing \pkg{lmebayes} API.
#'
#' @name glmbayesCore-callbacks
#' @aliases EnvelopeOpt EnvelopeSort rNormal_reg.wfit glmbfamfunc rgamma_ct
#' @keywords internal
NULL

#' @rdname glmbayesCore-callbacks
#' @export
EnvelopeOpt <- glmbayesCore::EnvelopeOpt

#' @rdname glmbayesCore-callbacks
#' @export
EnvelopeSort <- glmbayesCore::EnvelopeSort

#' @rdname glmbayesCore-callbacks
#' @export
rNormal_reg.wfit <- glmbayesCore::rNormal_reg.wfit

#' @rdname glmbayesCore-callbacks
#' @export
glmbfamfunc <- glmbayesCore::glmbfamfunc

#' @rdname glmbayesCore-callbacks
#' @export
rgamma_ct <- glmbayesCore::rgamma_ct
