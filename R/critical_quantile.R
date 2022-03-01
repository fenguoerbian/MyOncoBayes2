#' Critical quantile
#'
#' @description Calculates the critical value of a model covariate for
#'     which a fixed quantile of the crude rate is equal to the
#'     specified (tail or interval) probability. The specified
#'     quantile can relate to the posterior or the predictive
#'     distribution of the response rate.
#'
#' @name critical_quantile
#'
#' @template args-methods
#' @param x character giving the covariate name which is being search
#'     for to meet the critical conditions.  This also supports 'tidy'
#'     parameter selection by specifying \code{x = vars(...)}, where
#'     \code{...}  is specified the same way as in
#'     \code{\link[dplyr:select]{dplyr::select()}} and similar
#'     functions. Examples of using \code{x} in this way can be found
#'     in the examples. For \code{blrm_trial} methods, it defaults to
#'     the first entry in \code{summary(blrm_trial,
#'     "drug_info")$drug_name}.
#' @param p cumulative probability corresponding to the critical
#'     quantile range.
#' @param qc if given as a sorted vector of length 2, then the two
#'     entries define the interval on the response rate scale which
#'     must correspond to the cumulative probability \code{p}. In this
#'     case \code{lower.tail} must be \code{NULL} or left missing. If
#'     \code{qc} is only a single number, then \code{lower.tail} must
#'     be set to complete the specification of the tail area.
#' @param lower.tail defines if probabilities are lower or upper tail
#'     areas. Must be defined if \code{qc} is a single number and must
#'     not be defined if \code{qc} denotes an interval.
#' @param interval.x interval the covariate \code{x} is
#'     searched. Defaults to the minimal and maximal value of \code{x}
#'     in the data set used. Whenever \code{lower.tail} is set such
#'     that it is known that a tail area is used, then the function
#'     will automatically attempt to enlarge the search interval since
#'     the direction of the inverted function is determined around the
#'     critical value.
#' @param extendInt.x controls if the interval given in
#'     \code{interval.x} should be extended during the root
#'     finding. The default \code{"auto"} attempts to guess an
#'     adequate setting in cases thats possible. Other possible values
#'     are the same as for the \code{extendInt} argument of the
#'     \code{\link[stats:uniroot]{stats::uniroot()}} function.
#' @param log.x determines if during the numerical search for the
#'     critical value the covariate space is logarithmized. By default
#'     \code{x} is log transformed whenever all values are positive.
#' @param predictive logical indicates if the critical quantile
#'     relates to the posterior (\code{FALSE}) or the predictive
#'     (\code{TRUE}) distribution. Defaults to \code{FALSE}.
#' @param maxiter maximal number of iterations for root
#'     finding. Defaults to \code{100}.
#' @template args-dots-ignored
#'
#'
#' @details
#'
#' The function searches for a critical value \eqn{x_c} of the covariate
#' \eqn{x} (\code{x}) at which the integral over the model response in the
#' specified interval \eqn{\pi_1} to \eqn{\pi_2} (\code{qc}) is equal to the
#' given probability \eqn{p} (\code{p}). The given interval is considered
#' a tail area when \code{lower.tail} is used such that \eqn{\pi_1 = 0}
#' for \code{TRUE} and \eqn{\pi_2=1} for \code{FALSE}. At the
#' critical value \eqn{x_c} the equality holds
#'
#' \deqn{p = \int_{\pi_1}^{\pi_2} p(\pi | x_c) \, d\pi .}
#'
#' Note that a solution not guranteed to exist under all
#' circumstances. However, for a single agent model and when
#' requesting a tail probability then a solution does exist. In case
#' an interval probability is specified, then the solution may not
#' necessarily exist and the root finding is confined to the range
#' specified in \code{interval.x}.
#'
#' In case the solution is requested for the predictive distribution
#' (\code{predictive=TRUE}), then the respective problem solved leads
#' to the equality of
#'
#' \deqn{p = \sum_{y= r_1 = \lceil n \, \pi_1 \rceil }^{r_2 = \lceil n \, \pi_2 \rceil } \int p(y|\pi, n) \, p(\pi | x_c) \, d\pi .}
#'
#' Furthermore, the covariate space is log transformed by default
#' whenever all values of the covariate \eqn{x} are positive in the data
#' set. Values of \eqn{0} are shifted into the positive domain by the
#' machine precision to avoid issues with an ill-defined \eqn{\log(0)}.
#'
#' For \code{blrm_trial} objects the default arguments for \code{p},
#' \code{qc} and \code{lower.tail} are set to correspond to the
#' highest interval of \code{interval_prob} to be constrained by the
#' respective \code{interval_max_mass} (which are defined as part of
#' the \code{blrm_trial} objects). This will in the default case
#' correspond to the EWOC metric. These defaults are only applied if
#' \code{p}, \code{qc} and \code{lower.tail} are missing.
#'
#' @return Vector of length equal to the number of rows of the input
#'     data set with the crticial value for the covariate specified as
#'     \code{x} which fullfills the requirements as detailled
#'     above. May return \code{NA} for cases where no solution is
#'     found on the specified interval \code{interval.x}.
#'
#'
#' @template start-example
#' @examples
#' # fit an example model. See documentation for "combo2" example
#' example_model("combo2", silent=TRUE)
#'
#' # Find dose of drug_A at which EWOC criterium is just fulfilled
#' data_trial_ab <- subset(codata_combo2, group_id=="trial_AB")
#' drug_A_crit <- critical_quantile(blrmfit, newdata=data_trial_ab,
#'                                  x="drug_A", p=0.25, qc=0.33,
#'                                  lower.tail=FALSE)
#' data_trial_ab$drug_A <- drug_A_crit
#' summary(blrmfit, newdata=data_trial_ab, interval_prob=c(0,0.16,0.33,1))
#'
#' @template stop-example
#'
NULL

#' @rdname critical_quantile
#' @export
critical_quantile <- function(object, ...) UseMethod("critical_quantile")

#' @method critical_quantile default
#' @noRd
critical_quantile.default <- function(object, ...) { stop("object must inherit blrmfit or blrm_trial class") }

#' @rdname critical_quantile
#' @method critical_quantile blrmfit
#' @export
critical_quantile.blrmfit <- function(object, newdata, x, p, qc, lower.tail, interval.x, extendInt.x=c("auto", "no", "yes", "downX", "upX"), log.x, predictive=FALSE, maxiter=100, ...) {
    if(missing(newdata)) {
        data <- object$data
    } else {
        data <- newdata
    }
    xvar <- check_plot_variables(x, newdata=data)$x
    assert_that(length(xvar) == 1, msg="Can only have a single x variable.")
    if(missing(interval.x)) {
        interval.x <- c(min(data[,xvar], na.rm=TRUE), max(data[,xvar], na.rm=TRUE))
    }
    checkmate::assert_numeric(interval.x, finite=TRUE, len=2, any.missing=FALSE, sorted=TRUE)
    assert_that(interval.x[1] < interval.x[2], msg=paste("Lower interval.x[1] =", interval.x[1], "is not smaller than the upper interval.x[2] =", interval.x[2]))
    checkmate::assert_number(p, lower=0, upper=1)
    if(missing(lower.tail) || is.null(lower.tail)) {
        checkmate::assert_numeric(qc, lower=0, upper=1, len=2, sorted=TRUE, any.missing=FALSE)
        extendInt_auto <- "no"
        lower.tail <- NULL
    } else {
        checkmate::assert_logical(lower.tail, len=1, any.missing=FALSE)
        checkmate::assert_number(qc, lower=0, upper=1)
        extendInt_auto <- if(lower.tail) "downX" else "upX"
    }
    if(missing(log.x)) {
        log.x <- all(data[,xvar] >= 0)
    }
    checkmate::assert_logical(log.x, len=1, any.missing=FALSE)
    checkmate::assert_logical(predictive, len=1, any.missing=FALSE)
    extendInt.x <- match.arg(extendInt.x)
    if(extendInt.x == "auto")
        extendInt <- extendInt_auto
    else
        extendInt <- extendInt.x
    num_results <- nrow(data)
    qx <- rep(NA, times=num_results)
    machine_low <- .Machine$double.eps
    if(log.x)
        interval.x <- log(interval.x + machine_low) ## avoid trouble with 0
    Lp <- logit(max(min(p, 1-machine_low), machine_low))
    for(i in seq_len(num_results)) {
        data_row <- data[i,,drop=FALSE]
        root <- function(xc) {
            data_row[1,xvar] <- if(log.x) exp(xc) else xc
            logit(min(max(.model_distribution(object, data_row, qc, lower.tail=lower.tail, predictive=predictive), machine_low), 1-machine_low)) - Lp
        }
        tryCatch({
            qx[i] <- uniroot(root, interval.x, extendInt=extendInt, check.conv=TRUE, maxiter=maxiter)$root
        },
        error=function(root_error) {
            warning("No critical value found for row ", i,". Error message from uniroot:\n", root_error)
        })
    }
    return(if(log.x) exp(qx) else qx)
}

#' @rdname critical_quantile
#' @method critical_quantile blrm_trial
#' @export
critical_quantile.blrm_trial <- function(object, newdata, x, p, qc, lower.tail, interval.x, extendInt.x=c("auto", "no", "yes", "downX", "upX"), log.x, predictive=FALSE, maxiter=100, ...) {
    .assert_is_blrm_trial_and_prior_is_set(object)

    drug_info <- summary(object, "drug_info")

    if(missing(newdata)) newdata <- summary(object, "dose_info")

    if(missing(x)) x <- drug_info$drug_name[1]

    if(missing(p) && missing(qc) && missing(lower.tail)) {
        interval_max_mass <- summary(object, "interval_max_mass")
        p <- interval_max_mass[length(interval_max_mass)]
        ## ensure that all other categories are not constraining the
        ## criteria
        assert_that(all(interval_max_mass[-length(interval_max_mass)] == 1), msg="Non-standard EWOC detected. Please set p, qc and lower.tail.")

        interval_prob <- summary(object, summarize="interval_prob")
        qc <- interval_prob[(length(interval_prob)-1):length(interval_prob)]
        if(length(qc) == 2 && qc[2] == 1) {
            lower.tail <- FALSE
            qc <- qc[1]
        }
        if(length(qc) == 2 && qc[1] == 0) {
            lower.tail <- TRUE
            qc <- qc[2]
        }
    }

    return(critical_quantile(object$blrmfit, newdata, x, p, qc, lower.tail, interval.x, extendInt.x, log.x, predictive, maxiter))
}

## given a model return Pr( pi < q ) if lower.tail=TRUE. This is
## evaluted within the data set used for fitting (newdata left empty)
## or within the data.frame of newdata. In case the predictive
## distribution is used, then the crude rate is used for the
## quantiles.
.model_distribution <- function(model, newdata, q, lower.tail, predictive=FALSE) {
    if(missing(lower.tail) || is.null(lower.tail)) {
        ##cat("Missing: lower.tail\n")
        checkmate::assert_numeric(q, lower=0, upper=1, len=2, sorted=TRUE, any.missing=FALSE)
        interval <- q
    } else {
        ##cat("Not missing: lower.tail =", lower.tail, "\n")
        checkmate::assert_logical(lower.tail, len=1, any.missing=FALSE)
        checkmate::assert_number(q, lower=0, upper=1)
        if(lower.tail) {
            if (predictive) {
                ## for the predictive case we have to start from a
                ## negative value in order to include 0 in the
                ## interval
                interval <- c(-0.01, q)
            } else {
                interval <- c(0, q)
            }
        } else {
            interval <- c(q, 1)
        }
    }

    prob_col <- levels(cut(numeric(0), breaks=c(-Inf, interval , Inf)))[2]
    unname(summary(model, newdata=newdata, interval_prob=interval, predictive=predictive, transform=TRUE)[,prob_col])
}

