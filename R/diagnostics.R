#' Extract Diagnostic Quantities of \pkg{OncoBayes2} Models
#'
#' Extract quantities that can be used to diagnose sampling behavior
#' of the algorithms applied by \pkg{Stan} at the back-end of
#' \pkg{OncoBayes2}.
#'
#' @name diagnostic-quantities
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param object A \code{blrmfit} or \code{blrmtrial} object.
#' @param pars An optional character vector of parameter names.
#'   For \code{nuts_params} these will be NUTS sampler parameter
#'   names rather than model parameters. If \code{pars} is omitted
#'   all parameters are included.
#' @param ... Arguments passed to individual methods.
#'
#' @return The exact form of the output depends on the method.
#'
#' @details For more details see
#'   \code{\link[bayesplot:bayesplot-extractors]{bayesplot-extractors}}.
#'
#' @template start-example
#' @examples
#' example_model("single_agent", silent=TRUE)
#'
#' head(log_posterior(blrmfit))
#'
#' np <- nuts_params(blrmfit)
#' str(np)
#' # extract the number of divergence transitions
#' sum(subset(np, Parameter == "divergent__")$Value)
#'
#' head(rhat(blrmfit))
#' head(neff_ratio(blrmfit))
#'
#' @template stop-example
#'
NULL

#' @rdname diagnostic-quantities
#' @method log_posterior blrmfit
#' @importFrom bayesplot log_posterior
#' @export log_posterior
#' @export
log_posterior.blrmfit <- function(object, ...) {
  .contains_draws(object)
  bayesplot::log_posterior(object$stanfit, ...)
}

#' @rdname diagnostic-quantities
#' @method nuts_params blrmfit
#' @importFrom bayesplot nuts_params
#' @export nuts_params
#' @export
nuts_params.blrmfit <- function(object, pars = NULL, ...) {
  .contains_draws(object)
  bayesplot::nuts_params(object$stanfit, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @method rhat blrmfit
#' @importFrom bayesplot rhat
#' @export rhat
#' @export
rhat.blrmfit <- function(object, pars = NULL, ...) {
  .contains_draws(object)
  bayesplot::rhat(object$stanfit, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @method neff_ratio blrmfit
#' @importFrom bayesplot neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.blrmfit <- function(object, pars = NULL, ...) {
  .contains_draws(object)
  bayesplot::neff_ratio(object$stanfit, pars = pars, ...)
}

## --- internal

.contains_draws <- function(object) {
    assert_that(nsamples(object) > 0, msg="The model does not contain posterior draws.")
}
