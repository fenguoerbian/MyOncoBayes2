#'
#' # SBC Graphical Results
#' 
#'
#' ## Model 1: Single-agent logistic regression
#'
#' ### Component intercept/slopes
#'
#' #### Means
#'
print(pl_split$log2bayes_EXNEX.mu_log_beta$hist)
print(pl_split$log2bayes_EXNEX.mu_log_beta$ecdf_diff)
#'
#' #### Standard deviations (tau's)
#'
print(pl_split$log2bayes_EXNEX.tau_log_beta$hist)
print(pl_split$log2bayes_EXNEX.tau_log_beta$ecdf_diff)
#'
#' ### Component intercept/slopes: group estimates
#'
#'
#' #### Group estimates components
#'
print(pl_split$log2bayes_EXNEX.beta_group$hist)
print(pl_split$log2bayes_EXNEX.beta_group$ecdf_diff)

#'
#' ## Model 2: Double combination, fully exchangeable
#'
#' ### Component intercept/slopes: exchangeable mixture component
#'
#' #### Means
#'
print(pl_split$combo2_EX.mu_log_beta$hist)
print(pl_split$combo2_EX.mu_log_beta$ecdf_diff)
#'
#' #### Standard deviations (tau's)
#'
print(pl_split$combo2_EX.tau_log_beta$hist)
print(pl_split$combo2_EX.tau_log_beta$ecdf_diff)
#'
#' ### Interaction parameters (from exchangeable part)
#'
#' #### Mean
#'
print(pl_split$combo2_EX.mu_eta$hist)
print(pl_split$combo2_EX.mu_eta$ecdf_diff)
#'
#' #### Standard deviation
#'
print(pl_split$combo2_EX.tau_eta$hist)
print(pl_split$combo2_EX.tau_eta$ecdf_diff)
#'
#' ### Component intercept/slopes: group estimates
#'
#' #### Group estimates components
#'
print(pl_split$combo2_EX.beta_group$hist)
print(pl_split$combo2_EX.beta_group$ecdf_diff)
#'
#' #### Group estimates interaction(s)
#'
print(pl_split$combo2_EX.eta_group$hist)
print(pl_split$combo2_EX.eta_group$ecdf_diff)

#'
#'
#' ## Model 3: Double combination, EXchangeable/NonEXchangeable model
#'
#' ### Component intercept/slopes: exchangeable mixture component
#'
#' #### Means
#'
print(pl_split$combo2_EXNEX.mu_log_beta$hist)
print(pl_split$combo2_EXNEX.mu_log_beta$ecdf_diff)
#'
#' ### Standard deviations (tau's)
#'
print(pl_split$combo2_EXNEX.mu_log_beta$hist)
print(pl_split$combo2_EXNEX.mu_log_beta$ecdf_diff)
#'
#' ### Interaction parameters (from exchangeable part)
#'
#' #### Mean
#'
print(pl_split$combo2_EXNEX.mu_eta$hist)
print(pl_split$combo2_EXNEX.mu_eta$ecdf_diff)
#'
#' ### Standard deviation (tau)
#'
print(pl_split$combo2_EXNEX.tau_eta$hist)
print(pl_split$combo2_EXNEX.tau_eta$ecdf_diff)
#'
#' ### Component intercept/slopes: group estimates
#'
#'
#' #### Group estimates components
#'
print(pl_split$combo2_EXNEX.beta_group$hist)
print(pl_split$combo2_EXNEX.beta_group$ecdf_diff)
#'
#' #### Group estimates interaction(s)
#'
print(pl_split$combo2_EXNEX.eta_group$hist)
print(pl_split$combo2_EXNEX.eta_group$ecdf_diff)

#'
#' ## Model 4: Triple combination, EX/NEX model
#'
#' ### Component intercept/slopes: exchangeable mixture component
#'
#' #### Means
#'
print(pl_split$combo3_EXNEX.mu_log_beta$hist)
print(pl_split$combo3_EXNEX.mu_log_beta$ecdf_diff)
#'
#' ### Standard deviations (tau's)
#'
print(pl_split$combo3_EXNEX.mu_log_beta$hist)
print(pl_split$combo3_EXNEX.mu_log_beta$ecdf_diff)
#'
#' ### Interaction parameters (means from exchangeable part)
#'
#' #### Mean
#'
print(pl_split$combo3_EXNEX.mu_eta$hist)
print(pl_split$combo3_EXNEX.mu_eta$ecdf_diff)
#'
#' ### Standard deviation (tau)
#'
print(pl_split$combo3_EXNEX.tau_eta$hist)
print(pl_split$combo3_EXNEX.tau_eta$ecdf_diff)
#'
#' ### Component intercept/slopes: group estimates
#'
#'
#' #### Group estimates components
#'
print(pl_split$combo3_EXNEX.beta_group$hist)
print(pl_split$combo3_EXNEX.beta_group$ecdf_diff)
#'
#' #### Group estimates interaction(s)
#'
print(pl_split$combo3_EXNEX.eta_group$hist)
print(pl_split$combo3_EXNEX.eta_group$ecdf_diff)

