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

#' @inherit glmbayesCore::block_rNormalGLM title description params return seealso
#' @export
block_rNormalGLM <- glmbayesCore::block_rNormalGLM

#' @inherit glmbayesCore::block_rNormalGLM_update title description params return seealso
#' @export
block_rNormalGLM_update <- glmbayesCore::block_rNormalGLM_update

#' @inherit glmbayesCore::block_rNormalReg title description params return seealso
#' @export
block_rNormalReg <- glmbayesCore::block_rNormalReg

#' @inherit glmbayesCore::block_rNormalReg_update title description params return seealso
#' @export
block_rNormalReg_update <- glmbayesCore::block_rNormalReg_update

#' @inherit glmbayesCore::normalize_block
#' @export
normalize_block <- glmbayesCore::normalize_block

#' @inherit glmbayesCore::build_mu_all
#' @export
build_mu_all <- glmbayesCore::build_mu_all

#' @inherit glmbayesCore::lmerb_posterior_mean
#' @export
lmerb_posterior_mean <- glmbayesCore::lmerb_posterior_mean

#' @inherit glmbayesCore::glmerb_posterior_mode
#' @export
glmerb_posterior_mode <- glmbayesCore::glmerb_posterior_mode

#' @inherit glmbayesCore::two_block_rNormal_reg
#' @export
two_block_rNormal_reg <- glmbayesCore::two_block_rNormal_reg

#' @inherit glmbayesCore::multi_rNormal_reg
#' @export
multi_rNormal_reg <- glmbayesCore::multi_rNormal_reg

#' @inherit glmbayesCore::multi_rlmb
#' @examples
#' \dontrun{
#' ## Full runnable examples are maintained in \pkg{glmbayesCore}:
#' example(multi_rlmb, package = "glmbayesCore", ask = FALSE, echo = TRUE)
#' }
#' @export
multi_rlmb <- glmbayesCore::multi_rlmb

#' @inherit glmbayesCore::multi_rNormalGamma_reg
#' @export
multi_rNormalGamma_reg <- glmbayesCore::multi_rNormalGamma_reg

#' @inherit glmbayesCore::multi_rindepNormalGamma_reg
#' @export
multi_rindepNormalGamma_reg <- glmbayesCore::multi_rindepNormalGamma_reg

#' @inherit glmbayesCore::summary.mrglmb
#' @param ... Additional optional arguments passed to
#'   \code{\link[glmbayesCore]{summary.rglmb}}.
#' @export
#' @method summary mrglmb
summary.mrglmb <- getFromNamespace("summary.mrglmb", "glmbayesCore")

#' @inherit glmbayesCore::print.summary.mrglmb
#' @param ... Additional optional arguments passed to
#'   \code{\link[glmbayesCore]{print.summary.rglmb}}.
#' @export
#' @method print summary.mrglmb
print.summary.mrglmb <- getFromNamespace("print.summary.mrglmb", "glmbayesCore")

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
