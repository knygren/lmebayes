#' @aliases lmebayes
#'
#' @title lmebayes: Bayesian Linear Mixed-Effects Models via Two-Block Gibbs Sampling
#'
#' @description
#' Two-block Gibbs samplers for Bayesian linear and generalized linear mixed-effects
#' models, following \pkg{lme4} notation. Builds on \pkg{glmbayes} for iid GLM
#' sampling within blocks.
#'
#' @details
#' Row-block interfaces include \code{\link{block_lmb}} and \code{\link{block_glmb}};
#' mixed-model setup from \pkg{lme4} formulas via \code{\link{model_setup}}.
#' Lower-level simulation uses \code{\link[glmbayesCore]{simfunction}} and envelope
#' utilities from \pkg{glmbayesCore}.
#'
#' See the package README at \url{https://github.com/knygren/lmebayes} for examples.
#'
#' @section OpenCL startup checks:
#' In interactive sessions, attaching the package with \code{library(lmebayes)}
#' may emit a short \code{\link{packageStartupMessage}}
#' when \code{has_opencl()} is \code{FALSE} (typical for CRAN binaries) but a
#' GPU or OpenCL stack appears available on the host. OpenCL modelling paths
#' require a source install with OpenCL at compile time;
#' \code{has_opencl()} then reports whether that build succeeded.
#' Set \code{options(glmbayes.quiet_opencl_startup = TRUE)} to suppress attach
#' notes (recommended for CI and \command{R CMD check}).
#'
#' @example inst/examples/Ex_lmebayes-package.R
#'
#' @seealso
#' \code{\link{lmerb}}, \code{\link{model_setup}}, \code{\link{block_lmb}}, \code{\link{block_glmb}};
#' \code{\link[glmbayesCore]{simfunction}}, \code{\link[glmbayesCore]{EnvelopeBuild}};
#' \code{\link[glmbayes]{lmb}} and \code{\link[glmbayes]{glmb}} for fixed-effects-only Bayesian
#' linear and generalized linear models (from \pkg{glmbayes});
#' \code{\link[glmbayesCore]{rlmb}} and \code{\link[glmbayesCore]{rglmb}} for iid
#' posterior draws.
#'
#' Useful links:
#' \itemize{
#'   \item GitHub: \url{https://github.com/knygren/lmebayes}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @author
#' Kjell Nygren
#'
#' @import stats Rcpp glmbayesCore
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom glmbayes lmb glmb
#' @importFrom glmbayes glmb.covratio glmb.dffits glmb.influence.measures
#' @importFrom glmbayes extractDIC directional_tail
#' @importFrom glmbayes has_opencl get_opencl_core_count
#' @importFrom glmbayesCore Prior_Setup dNormal dNormal_Gamma multi_prior_setup rlmb rglmb
#' @useDynLib lmebayes, .registration = TRUE
"_PACKAGE"
