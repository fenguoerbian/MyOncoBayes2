#' Posterior predictive intervals
#'
#' Posterior predictive intervals of the model.
#'
#' @template args-methods
#' @template args-prob
#' @template args-dots-ignored
#'
#' @details
#'
#' Reports for each row of the input data set the predictive interval
#' according to the fitted model.
#'
#' @return Matrix with as many rows as the input data set and two
#'     columns which contain the lower and upper quantile
#'     corresponding to the central probability mass \code{prob} for
#'     the number of responses of the predictive distribution.
#'
#' @template example-start
#' @examples
#' example_model("single_agent")
#'
#' predictive_interval(blrmfit)
#'
#' @template example-stop
#'
#' @method predictive_interval blrmfit
#' @aliases predictive_interval
#' @export
predictive_interval.blrmfit <- function(object, prob=0.95, newdata, ...) {
    yrep <- posterior_predict(object, newdata=newdata)
    rstantools::predictive_interval(yrep, prob=prob)
}
