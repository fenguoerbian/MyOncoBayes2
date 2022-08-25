#' Summarise model results
#'
#' Provides model summaries for \code{\link{blrm_exnex}} and
#' \code{\link{blrm_trial}} analyses.
#'
#' @template args-methods
#' @template args-prob
#' @param interval_prob optional vector of sorted quantiles for which
#'     the interval probabilities are calculated
#' @param predictive logical indicates if the posterior predictive is
#'     being summarized. Defaults to \code{FALSE}.
#' @template args-transform
#' @param ... not used in this function
#'
#' @details
#' The calculated posterior summaries are returned as a
#' \code{data.frame} and contain optional interval probabilites for
#' the specified vector of sorted quantiles. These summaries are
#' calculated on the response scale by default and can be obtained on
#' the link scale when setting \code{transform=FALSE}.
#'
#' When the results are requested for the predictive distribution with
#' \code{predictive=TRUE}, then the link scale refers to the total
#' counts while the transformed scale divides the (predictive) counts
#' by the number of trials such that results are on the 0-1 scale.
#'
#' @return Returns a \code{data.frame} of the key summaries of the
#'     posterior mean, standard deviation, central probability
#'     interval, median and optional interval probabilities. Each row
#'     of the \code{data.frame} corresponds to the respective input
#'     data which is by default the same data set as used for the
#'     \code{\link{blrm_exnex}} analysis or the data specified in the
#'     \code{newdata} argument.
#'
#' @template start-example
#' @examples
#' example_model("single_agent", silent=TRUE)
#'
#' ## obtain underdosing (0-0.16), target dosing (0.16-0.33) and
#' ## overdosing (0.33-1) probabilities
#' summary(blrmfit, interval_prob=c(0,0.16,0.33,1))
#'
#' ## obtain predictive distribution for respective cohorts and
#' ## calculate probability for no event, 1 event or >1 event
#' ## note that this does the calculation for the cohort sizes
#' ## as put into the data-set
#' summary(blrmfit, interval_prob=c(-1,0,1,10), predictive=TRUE)
#'
#' ## to obtain the predictive for a cohort-size of 6 for all patients
#' ## in the data-set one would need to use the newdata argument, e.g.
#' summary(blrmfit, newdata=transform(hist_SA, num_patients=6),
#'                  interval_prob=c(-1,0,1,10), predictive=TRUE)
#'
#' @template stop-example
#'
#' @method summary blrmfit
#' @export
summary.blrmfit <- function(object, newdata, transform=!predictive, prob=0.95, interval_prob, predictive=FALSE, ...) {
    assert_numeric(prob, lower=0, upper=1, finite=TRUE, any.missing=FALSE, min.len=1)

    data <- object$data
    if(!missing(newdata))
        data <- newdata
    rn <- rownames(data)
    nr <- nrow(data)

    probs <- sort(unique(c(0.5, 0.5- prob/2, 0.5+ prob/2)))

    ## use the SAS defintion for quantiles whenever the output are
    ## discrete as for the predictive which is type 2. Type 7 is the
    ## default in R which gives non-integer outputs for integer only
    ## data.
    quantile_type <- ifelse(predictive, 2, 7)

    if(predictive) {
        sizes <- pp_binomial_trials(object, newdata=newdata)

        post_prob <- posterior_linpred(object, transform=TRUE, newdata=newdata)

        ## integrate binomial sampling density per outcome and then
        ## summarize the obtained predictive p(y)

        draws <- nrow(post_prob)
        max_size <- max(sizes)
        pred_lpmf <- matrix(-Inf, nrow=max_size+1, ncol=nr)
        pred_pmf <- matrix(0, nrow=max_size+1, ncol=nr)
        pred_cpmf <- matrix(1.0, nrow=max_size+1, ncol=nr)
        for(r in seq_len(nr)) {
            ## create matrix of outcomes x draws of prob
            outcomes_row <- 0:sizes[r]
            nor <- length(outcomes_row)
            Pr <- matrix(post_prob[,r], nrow=nor, ncol=draws, byrow=TRUE)
            pred_lpmf[seq_len(nor),r] <- apply(dbinom(outcomes_row, sizes[r], prob=Pr, log=TRUE), 1, log_mean_exp)
            ## normalize to sum-to-one per column
            pred_lpmf[seq_len(nor),r] <- pred_lpmf[seq_len(nor),r] - matrixStats::logSumExp(pred_lpmf[seq_len(nor),r])
            pred_pmf[,r] <- exp(pred_lpmf[,r])
            pred_cpmf[seq_len(nor),r] <- cumsum(pred_pmf[seq_len(nor),r])
            ## ensure that largest number is truly 1.0
            pred_cpmf[,r]  <- pmin(pred_cpmf[,r], 1.0)
        }

        outcomes <- 0:max_size
        all_outcomes <- matrix(outcomes, nrow=max_size+1, ncol=nr)

        if(transform) {
            ## transform outcomes to 0-1 scaled response rates
            all_outcomes <- sweep(all_outcomes, 2, sizes, "/")
            if(any(sizes == 0))
                warning("Data rows ", paste(which(sizes == 0), collapse=", "), " contain trials of size 0.")
        }

        m <- colSums(all_outcomes * pred_pmf)
        s <- sqrt(colSums(sweep(all_outcomes, 2, m)^2 * pred_pmf))

        quants <- matrix(NA, nrow=nr, ncol=length(probs))
        for(r in seq_len(nr)) {
            quants[r,]  <- findInterval(probs, pred_cpmf[seq_len(sizes[r]+1),r])
        }
        if(transform)
            quants <- sweep(quants, 1, sizes, "/")

        out_sum <- cbind(data.frame(mean=m, sd=s), quants)

        if(!missing(interval_prob)) {
            if(transform) {
                assert_numeric(interval_prob, lower=-1, upper=1, any.missing=FALSE, sorted=TRUE, finite=TRUE)
            } else {
                assert_integerish(interval_prob, any.missing=FALSE, sorted=TRUE)
            }

            nip <- length(interval_prob)
            out_interval <- matrix(NA, nrow=nr, ncol=nip+1)
            for(r in seq_len(nr)) {
                all_outcomes_r <- all_outcomes[seq_len(sizes[r]+1),r]
                outcome_idx <- findInterval(c(-Inf, interval_prob, Inf), all_outcomes_r)
                cpmf <- c(0,pred_cpmf[,r])
                out_interval[r,] <- diff(cpmf[outcome_idx+1])
            }
            out_interval <- out_interval[,-c(1,nip+1),drop=FALSE]
            colnames(out_interval) <- paste0("(", interval_prob[seq_len(nip-1)], ",", interval_prob[1 + seq_len(nip-1)], "]")
            out_sum <- cbind(out_sum, as.data.frame(out_interval))
        }

    } else {
        post <- posterior_linpred(object, transform=transform, newdata=newdata)
        S <- nrow(post)

        if (ncol(post) > 0) {
          out_sum <- as.data.frame(t(apply(post, 2, function(l) {
            if(all(is.na(l)))
              return(rep(NA, times=2+length(probs)))
            c(mean=mean(l), sd=sd(l), quantile(l, probs=probs, type=quantile_type))
          })))
        } else {
          out_sum <- data.frame(0.0 * matrix(ncol = (2+length(probs)), nrow = 0))
        }

        if(!missing(interval_prob)) {
            if(transform) {
                ## negative values are admissable as we otherwise cannot
                ## include 0
                assert_numeric(interval_prob, lower=-1, upper=1, any.missing=FALSE, sorted=TRUE, finite=TRUE)
            } else {
                assert_numeric(interval_prob, any.missing=FALSE, sorted=TRUE, finite=TRUE)
            }

            if (ncol(post) > 0)
            {
              out_interval <- as.data.frame(t(apply(post, 2, function(l) table(cut(l, breaks=c(-Inf, interval_prob, Inf)))/S )))
            } else {
              out_interval <- data.frame(0.0 * matrix(ncol = (1+length(interval_prob)), nrow = 0))
              colnames(out_interval) <- levels(cut(numeric(0), breaks=c(-Inf, interval_prob, Inf)))
            }

            out_interval <- out_interval[,-c(1,ncol(out_interval)), drop=FALSE]
            out_sum <- cbind(out_sum, out_interval)
        }
    }

    names(out_sum)[seq_len(2+length(probs))] <- c("mean", "sd", paste0(100 * probs, "%"))

    rownames(out_sum) <- rn
    out_sum
}

#' Summarise trial
#'
#' Provides model summaries for \code{\link{blrm_trial}} analyses.
#' @param object \code{\link{blrm_trial}} object
#' @param summarize one of the following options:
#' \itemize{
#'   \item{}{\code{blrmfit}: summary of the underlying blrmfit object with further arguments ...}
#'   \item{}{\code{blrm_exnex_call}: blrm_exnex call used to create the \code{blrmfit} object}
#'   \item{}{\code{drug_info}: drug_info for the trial, contains drugs, reference doses and units}
#'   \item{}{\code{dose_info}: dose_info that were defined}
#'   \item{}{\code{dose_prediction} prediction for the defined \code{dose_info}}
#'   \item{}{\code{data}: data that were observed}
#'   \item{}{\code{data_prediction}: prediction for the observed data}
#'   \item{}{\code{newdata_prediction}: prediction for data provided with the \code{newdata} argument}
#'   \item{}{\code{dimensionality}: numeric vector with entries \code{"num_components"},
#'     \code{"num_interaction_terms"}, \code{"num_groups"}, \code{"num_strata"}}
#'   \item{}{\code{interval_prob}: interval probabilities reported in the standard outputs}
#'   \item{}{\code{interval_max_mass}: named vector defining for each interval of the
#'     \code{interval_prob} vector a maxmimal admissable
#'     probability mass for a given dose level}
#'   \item{}{\code{ewoc_check}: MCMC diagnostic and precision estimates of ewoc defining metrics for the defined doses in \code{dose_info} (default) or for the doses provided in the \code{newdata} argument. Please refer to the details for reported diagnostics.}
#' }
#' @param ... further arguments for \code{\link{summary.blrmfit}}
#'
#' @details The \code{ewoc_check} summary routine allows to assess the
#'     accuracy and reliability of the ewoc criterion with respect to
#'     MCMC sampling noise. The returned summary provides detailled
#'     MCMC convergence and precision estimates for all criteria
#'     defined by \code{interval_prob} and \code{interval_max_mass}
#'     which contribute to EWOC metric. That is, for each interval
#'     probability with a maximal mass of less than unity the routine
#'     will return these columns:
#' \describe{
#' \item{\code{est}}{the MCMC estimate defining the critical value. For intervals defined by a tail probability this corresponds to the respective critical quantile while for interval probabilites this is equal to the interval probability.}
#' \item{\code{stat}}{centered and standardized test quantity. The estimate is centered by the critical value and scaled by the Monte-Carlo standard error (MCSE) of the estimate. Hence, negative (positive) values correspond to the constraint being (not) fulfilled. The standardization with the MCSE allows to compare the values to standard normal quantiles accordingly.}
#' \item{\code{mcse}}{the Monte-Carlo standard error of the estimate determined with \code{\link[posterior]{mcse_quantile}} (tail probability) or \code{\link[posterior]{mcse_mean}} (interval probability) functions.}
#' \item{\code{ess}}{the Monte-Carlo effective sample size of the estimate determined with \code{\link[posterior]{ess_quantile}} (tail probability) or \code{\link[posterior]{ess_mean}} (interval probability) functions.}
#' \item{\code{rhat}}{the Monte-Carlo non-convergence diagnostic Rhat as determined with the \code{\link[posterior]{rhat} function}.}
#' }
#'
#' For the common case of requiring that 33\% DLT probability is not
#' exceeded by more than 25\% of the posterior probability mass, the
#' estimate column \code{est} contains the 75\% quantile
#' \eqn{q_{75\%}}{q75\%} and the standardized statistic \code{stat} is
#' defined as:
#'
#' \deqn{\mbox{stat} = \frac{q_{75\%} - 33\%}{\mbox{mcse}_{q_{75\%}}}}{stat = (q75\% - 33\%)/mcse_q75\%}
#'
#' The statistic is approximately distributed as a standard normal
#' variate. The \code{ewoc_check} summary can be used to ensure that
#' the MCMC estimation accuracy is sufficient.
#'
#' @template start-example
#' @examples
#' # construct initial blrm_trial object from built-in example datasets
#' combo2_trial_setup <- blrm_trial(
#'   data = hist_combo2,
#'   dose_info = dose_info_combo2,
#'   drug_info = drug_info_combo2,
#'   simplified_prior = TRUE
#' )
#'
#' # extract blrm_call to see setup of the prior as passed to blrm_exnex
#' summary(combo2_trial_setup, "blrm_exnex_call")
#'
#' # extract ewoc precision accuracy
#' ec <- summary(combo2_trial_setup, "ewoc_check")
#'
#' # find any ewoc metrics which are within 95% MCMC error of the threshold
#' # these are counted as "imprecise" when printing blrm_trial objects
#' subset(ec, abs(prob_overdose_stat) < qnorm(0.975))
#'
#' # ensure that the ewoc metric only flags "ok" whenever the MCMC error
#' # is with 95% below the threshold
#' ewoc_ok <- ec$prob_overdose_stat < qnorm(0.025)
#'
#' @template stop-example
#'
#' @method summary blrm_trial
#' @export
summary.blrm_trial <- function(
  object,
  summarize = c(
    "blrmfit",
    "blrm_exnex_call",
    "data",
    "drug_info",
    "dose_info",
    "dose_prediction",
    "data_prediction",
    "newdata_prediction",
    "dimensionality",
    "interval_prob",
    "interval_max_mass",
    "ewoc_check"
  ),
  ...
)
{
  args <- list(...)
  if (missing(summarize)) {
    if ("newdata" %in% names(args)) {
      summarize <- "newdata_prediction"
    } else {
      summarize <- "data_prediction"
    }
  } else {
    summarize <- match.arg(summarize)
  }



  if (summarize %in% c("blrmfit", "blrm_exnex_call", "dose_prediction", "data_prediction", "newdata_prediction", "ewoc_check")) {
    .assert_is_blrm_trial_and_prior_is_set(object)
  } else {
    .assert_is_blrm_trial(object)
  }

  newdata_predict <- function() {
    args[["newdata"]] <- .blrm_trial_sanitize_data(
      trial=object,
      data=args[["newdata"]],
      warning_if_dose_not_prespecified = FALSE,
      require_num_patients = FALSE,
      require_num_toxicities = FALSE
    )

    args_without_newdata <- args
    args_without_newdata[["newdata"]] <- NULL
    do.call(".blrm_trial_predict", c(list(object, args[["newdata"]]), args_without_newdata))
  }

  ewoc_check <- function() {
      if(!("newdata" %in% names(args))) {
          .blrm_trial_ewoc_check(object, object$ewoc_check)
          return(object$ewoc_check)
      }

      args[["newdata"]] <- .blrm_trial_sanitize_data(
          trial=object,
          data=args[["newdata"]],
          warning_if_dose_not_prespecified = FALSE,
          require_num_patients = FALSE,
          require_num_toxicities = FALSE
      )

      args_without_newdata <- args
      args_without_newdata[["newdata"]] <- NULL
      ewoc_check <- do.call(".blrm_trial_compute_ewoc_check", c(list(object, newdata=args[["newdata"]]), args_without_newdata))
      .blrm_trial_ewoc_check(object, ewoc_check)
      ewoc_check
  }

  switch(summarize,
    blrmfit = summary(object$blrmfit, ...),
    blrm_exnex_call = object$blrmfit$call,
    drug_info = object$drug_info,
    dose_info = object$dose_info,
    dose_prediction = .blrm_trial_predict(object, object$dose_info, ...),
    data = object$data,
    data_prediction = .blrm_trial_predict(object, object$data, ...),
    newdata_prediction = newdata_predict(),
    dimensionality = object[c("num_components", "num_interaction_terms", "num_groups", "num_strata")],
    interval_prob = object$interval_prob,
    interval_max_mass = object$interval_max_mass,
    ewoc_check = ewoc_check()
  )
}
