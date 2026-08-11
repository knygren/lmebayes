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
#' Row-block interfaces include \code{\link{lmbBlock}} and \code{\link{glmbBlock}};
#' mixed-model setup from \pkg{lme4} formulas via \code{\link{model_setup}}.
#' Lower-level simulation uses \code{\link[glmbayesCore]{simfunction}} and envelope
#' utilities from \pkg{lmebayesCore}.
#'
#' See the package README at \url{https://github.com/knygren/lmebayes} for examples.
#'
#' @section Vignette index:
#' The vignettes are a curriculum, meant to be read roughly in order. Each
#' assumes the previous ones.
#' \describe{
#'   \item{\code{vignette("Chapter-00")}}{Roadmap: what the package does and
#'     which chapter answers which question.}
#'   \item{\code{vignette("Chapter-01")}}{Getting started: one model end to
#'     end on \code{lme4::sleepstudy}.}
#'   \item{\code{vignette("Chapter-02")}}{From \code{lm}/\code{glm} to mixed
#'     models: pooling, the centered parameterization, and the model-class
#'     restrictions.}
#'   \item{\code{vignette("Chapter-03")}}{The prior workflow:
#'     \code{model_setup}, \code{Prior_Setup_GLMM}, \code{pfamily_list}.}
#'   \item{\code{vignette("Chapter-04")}}{Linear mixed models with
#'     group-level covariates and cross-level moderation.}
#'   \item{\code{vignette("Chapter-05")}}{Reading the output: every accessor
#'     and container.}
#'   \item{\code{vignette("Chapter-06")}}{Convergence and \code{tv_tol}.}
#'   \item{\code{vignette("Chapter-07")}}{Estimated variance components and
#'     the pilot stage.}
#'   \item{\code{vignette("Chapter-08")}}{Observation dispersion, including
#'     per-group \eqn{\sigma^2_j}.}
#'   \item{\code{vignette("Chapter-09")}}{Poisson GLMMs.}
#'   \item{\code{vignette("Chapter-10")}}{Binomial GLMMs.}
#'   \item{\code{vignette("Chapter-11")}}{Point estimates without
#'     simulation.}
#'   \item{\code{vignette("Chapter-12")}}{Row-block models.}
#'   \item{\code{vignette("Chapter-13")}}{Comparison with \pkg{lme4} and
#'     \pkg{glmmTMB}.}
#' }
#' Engine internals, derivations, and convergence theory are documented in
#' the \code{Chapter-B*} vignettes of \pkg{lmebayesCore}; the iid GLM stack
#' underneath is documented in \pkg{glmbayes}.
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
#' \code{\link{lmerb}}, \code{\link{model_setup}}, \code{\link{lmbBlock}}, \code{\link{glmbBlock}};
#' \code{\link[glmbayesCore]{simfunction}}, \code{\link[glmbayesCore]{EnvelopeBuild}};
#' \code{\link{lmb}} and \code{\link{glmb}} for fixed-effects-only Bayesian
#' linear and generalized linear models;
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
#' @import stats Rcpp lmebayesCore
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom utils flush.console
#' @importFrom glmbayes glmb.covratio glmb.dffits glmb.influence.measures
#' @importFrom glmbayes extractDIC
#' @importFrom glmbayes has_opencl get_opencl_core_count
#' @importFrom glmbayesCore Prior_Setup dNormal dNormal_Gamma multi_prior_setup rlmb rglmb
#' @useDynLib lmebayes, .registration = TRUE
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables("yint", package = "lmebayes", add = FALSE)
}
