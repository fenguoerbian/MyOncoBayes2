#' Posterior of linear predictor
#'
#' Calculates the posterior of the linear predictor.
#'
#' @template args-methods
#' @template args-posterior
#' @template args-transform
#' @template args-dots-ignored
#'
#' @details
#'
#' Simulates the posterior of the linear predictor of the model
#' \code{object} for the specified data set.
#'
#' @template note-groups
#' @template return-samples
#'
#' @template example-start
#' @examples
#'
#' ## run single-agent analysis which defines blrmfit model object
#' example_model("single_agent")
#'
#' ## obtain posterior of linear prediction on 0-1 scale, but first name
#' ## rows of input data to obtain nice labels with bayesplot
#' trial_design <- hist_SA
#' row.names(trial_design) <- hist_SA$DosesAdm1
#' post_prob_dlt <- posterior_linpred(blrmfit, TRUE, newdata=trial_design)
#'
#' library(bayesplot)
#' library(ggplot2)
#' mcmc_intervals(post_prob_dlt, prob=0.5, prob_outer=0.95) +
#'     coord_flip() +
#'     vline_at(c(0.16, 0.33), linetype=2) +
#'     ylab("Dose [mg]") +
#'     ggtitle("Posterior Probability of a DLT") +
#'     scale_x_continuous(breaks=c(0.1,0.16,0.33, 0.5, 0.75))
#'
#' @template example-stop
#'
#' @method posterior_linpred blrmfit
#' @aliases posterior_linpred
#' @export
posterior_linpred.blrmfit <- function(object, transform=FALSE, newdata, draws, ...) {
    dat <- pp_data(object, newdata=newdata, draws=draws)
    if (transform)
        dat <- inv_logit(dat)
    return(dat)
}
