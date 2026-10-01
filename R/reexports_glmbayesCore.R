## Re-export symbols from glmbayesCore (priors, calibration, row-block setup).

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

#' @inherit glmbayesCore::Prior_SetupGroup title description details params return seealso
#' @export
Prior_SetupGroup <- glmbayesCore::Prior_SetupGroup

## Row-block dispatch uses glmbayesCore's generic; GLMM uses lmebayesCore's
## (separate S3 generic objects — do not assign glmbayesCore::pfamily_list alone).
#' @inherit glmbayesCore::pfamily_list title description details return seealso
#' @section Mixed-model method:
#' Objects from \code{\link{Prior_Setup_GLMM}()} dispatch via \pkg{lmebayesCore}; see
#' \code{\link[lmebayesCore:pfamily_list.Prior_Setup_GLMM]{pfamily_list.Prior_Setup_GLMM}}.
#' @param object A prior-specification object (from \code{\link{Prior_SetupGroup}()}
#'   or \code{\link{Prior_Setup_GLMM}()}).
#' @param ptypes For \code{\link{Prior_Setup_GLMM}()}: one string recycled to every
#'   group-effect coefficient, or a named vector / list per coefficient.
#'   \code{NULL} resolves to \code{"dNormal"} for GLMM setup. For row blocks,
#'   see \code{\link{Prior_SetupGroup}()}.
#' @param ... Passed to methods.
#' @usage pfamily_list(object, ptypes = NULL, ...)
#' @export
pfamily_list <- function(object, ptypes = NULL, ...) {
  if (inherits(object, "Prior_Setup_GLMM")) {
    return(lmebayesCore::pfamily_list(object, ptypes = ptypes, ...))
  }
  glmbayesCore::pfamily_list(object, ptypes = ptypes, ...)
}
