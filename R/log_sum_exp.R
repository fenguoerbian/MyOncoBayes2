#' Numerically stable mean of logs
#' @keywords internal
log_mean_exp <- function(x) {
    return(matrixStats::logSumExp(x) - length(x))
}
