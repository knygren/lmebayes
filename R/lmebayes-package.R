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
#' Lower-level simulation uses \code{\link{simfuncs}} and envelope utilities from
#' the \pkg{glmbayes} lineage.
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
#' \code{\link{model_setup}}, \code{\link{block_lmb}}, \code{\link{block_glmb}},
#' \code{\link{simfuncs}}, \code{\link{EnvelopeBuild}};
#' \pkg{glmbayes} for \code{\link[glmbayes]{glmb}}, \code{\link[glmbayes]{lmb}},
#' \code{\link[glmbayes]{rglmb}}, and \code{\link[glmbayes]{rlmb}}.
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
#' @import stats Rcpp glmbayes
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib lmebayes, .registration = TRUE
"_PACKAGE"
