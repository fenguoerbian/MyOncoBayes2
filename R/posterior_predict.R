#' Posterior of predictive
#'
#' Simulation of the predictive distribution.
#'
#' @template args-methods
#' @template args-posterior
#' @template args-dots-ignored
#'
#' @details
#'
#' Simulates the posterior predictive of the model \code{object} for
#' the specified data set.
#'
#' @template note-groups
#' @template return-samples
#'
#' @template example-start
#' @examples
#'
#' example_model("single_agent")
#'
#' post_pred  <- posterior_predict(blrmfit)
#' ## turn DLT counts into DLT rates
#' post_pred_rate <- sweep(post_pred, 2, hist_SA$Npat, "/")
#'
#' library(bayesplot)
#' library(ggplot2)
#'
#' ## compare posterior predictive of the model for the response rates
#' ## with observed data
#' with(hist_SA, ppc_intervals(Ntox / Npat, post_pred_rate, x=DosesAdm1, prob_outer=0.95)) +
#'     xlab("Dose [mg]")
#'
#' @template example-stop
#'
#' @method posterior_predict blrmfit
#' @aliases posterior_predict
#' @export
posterior_predict.blrmfit <- function(object, newdata, draws, ...) {
    dat <- inv_logit(pp_data(object, newdata=newdata, draws=draws))
    num_sims <- nrow(dat)
    num_obs <- ncol(dat)
    num_trials <- pp_binomial_trials(object, newdata)

    ##pr <- t(apply(dat, 1, rbinom, n=num_obs, size=num_trials))
    pr <- matrix(rbinom(num_sims*num_obs, matrix(num_trials, nrow=num_sims, ncol=num_obs, byrow=TRUE), dat), nrow=num_sims, ncol=num_obs)
    colnames(pr) <- colnames(dat)
    pr
}
