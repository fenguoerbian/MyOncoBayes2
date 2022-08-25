#' OncoBayes2
#'
#' Bayesian logistic regression model with optional
#' EXchangeability-NonEXchangeability parameter modelling for flexible
#' borrowing from historical or concurrent data-sources. The safety
#' model can guide dose-escalation decisions for adaptive Oncology
#' phase I dose-escalation trials which involve an arbitrary number of
#' drugs.
#'
#' @section Global Options:
#'
#' \tabular{lcl}{
#' Option \tab Default \tab Description \cr
#' \code{OncoBayes2.MC.warmup} \tab 1000 \tab MCMC warmup iterations \cr
#' \code{OncoBayes2.MC.iter} \tab 2000 \tab total MCMC iterations \cr
#' \code{OncoBayes2.MC.save_warmup} \tab TRUE \tab save warmup samples \cr
#' \code{OncoBayes2.MC.chains} \tab 4 \tab MCMC chains\cr
#' \code{OncoBayes2.MC.thin} \tab 1 \tab MCMC thinning \cr
#' \code{OncoBayes2.MC.control} \tab \code{list(adapt_delta=0.99,} \tab sets \code{control} argument for Stan call\cr
#'  \tab \code{stepsize=0.1}) \tab \cr
#' \code{OncoBayes2.MC.backend} \tab rstan \tab Backend used to run Stan (\code{rstan} or \code{cmdstanr}) \cr
#' \code{OncoBayes2.abbreviate.min} \tab 0 \tab Minimal length of variable names \cr
#'    \tab \tab when abbreviating variable names. \cr
#'    \tab \tab The default 0 disables abbreviation.\cr
#' }
#'
#' @template ref-mac
#' @template ref-exnex
#' @template ref-critical_aspects
#' @template ref-bayesindustry
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
#' @name OncoBayes2
#' @aliases OncoBayes2
#' @docType package
#' @useDynLib OncoBayes2, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom rstan sampling get_sampler_params get_adaptation_info summary stanc_builder get_stancode get_divergent_iterations
#' @importFrom posterior as_draws_array as_draws_rvars as_draws_matrix summarise_draws subset_draws resample_draws default_convergence_measures rvar as_rvar ndraws variables draws_of variables<- %**% bind_draws mcse_mean ess_mean mcse_quantile ess_quantile
#' @importFrom utils capture.output modifyList combn head
#' @importFrom matrixStats logSumExp
#' @importFrom stats delete.response ftable median model.frame model.matrix model.response quantile rbinom sd terms model.matrix.default setNames update update.default .getXlevels as.formula na.fail qlogis dbinom uniroot qnorm
#' @import assertthat
#' @import checkmate
#' @import Formula
#' @import rstantools
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom tidyselect vars_select_helpers vars_select
#' @importFrom rlang env_bury
#' @importFrom scales number_format extended_breaks
#' @import abind
#' @export posterior_linpred posterior_predict posterior_interval
#' @export predictive_interval prior_summary nsamples
#'
#'
NULL
