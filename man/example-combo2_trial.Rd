% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_combo2_trial.R
\name{example-combo2_trial}
\alias{example-combo2_trial}
\title{Two-drug combination example using BLRM Trial}
\description{
Example using \code{\link{blrm_trial}} to
guide the built-in two-drug combination study example.
}
\details{
\code{\link{blrm_trial}} is used to collect
and store all relevant design information for the example. Subsequent
use of the \code{\link{update.blrm_trial}} command
allows convenient model fitting via
\code{\link{blrm_exnex}}. The
\code{\link{summary.blrm_trial}} method allows
exploration of the design and modeling results.

To run this example, use \code{example_model("combo2_trial")}. See
\code{\link{example_model}}.
}
\examples{
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(OncoBayes2.MC.warmup=10, OncoBayes2.MC.iter=20, OncoBayes2.MC.chains=1,
                            OncoBayes2.MC.save_warmup=FALSE)

library(tibble)
library(dplyr)
library(tidyr)

# Combo2 example using blrm_trial functionality

# construct initial blrm_trial object from built-in example datasets
combo2_trial_setup <- blrm_trial(
  data = hist_combo2,
  dose_info = dose_info_combo2,
  drug_info = drug_info_combo2,
  simplified_prior = FALSE
)

dims <- summary(combo2_trial_setup, "dimensionality")

# Fit the initial model with the historical data and fully specified prior
combo2_trial_start <- update(
  combo2_trial_setup,
  prior_EX_mu_mean_comp = matrix(
    c(logit(0.2), 0,  # hyper-mean of (intercept, log-slope) for drug A
      logit(0.2), 0), # hyper-mean of (intercept, log-slope) for drug B
    nrow = dims$num_components,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_mu_sd_comp = matrix(
    c(2.0, 1,  # hyper-sd of mean mu for (intercept, log-slope) for drug A
      2.0, 1), # hyper-sd of mean mu for (intercept, log-slope) for drug B
    nrow = dims$num_components,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_tau_mean_comp = matrix(
    c(log(0.25), log(0.125),
      log(0.25), log(0.125)),
    nrow = dims$num_components,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_tau_sd_comp = matrix(
    c(log(4) / 1.96, log(4) / 1.96,
      log(4) / 1.96, log(4) / 1.96),
    nrow = dims$num_components,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_mu_mean_inter = 0,
  prior_EX_mu_sd_inter = 1.121,
  prior_EX_tau_mean_inter = matrix(log(0.125),
                                   nrow = dims$num_strata,
                                   ncol = dims$num_interaction_terms),
  prior_EX_tau_sd_inter = matrix(log(4) / 1.96,
                                 nrow = dims$num_strata,
                                 ncol = dims$num_interaction_terms),
  prior_is_EXNEX_comp = rep(FALSE, dims$num_components),
  prior_is_EXNEX_inter = rep(FALSE, dims$num_interaction_terms),
  prior_EX_prob_comp = matrix(1,
                              nrow = dims$num_groups,
                              ncol = dims$num_components),
  prior_EX_prob_inter = matrix(1,
                               nrow = dims$num_groups,
                               ncol = dims$num_interaction_terms),
  prior_tau_dist = 1
)


# print summary of prior specification
prior_summary(combo2_trial_start)

# summarize inference at observed dose levels
summary(combo2_trial_start, "data_prediction")

# summarize inference at specified dose levels
summary(combo2_trial_start, "dose_prediction")


# Update again with new data

# using update() with data argument supplied
# dem <- update(combo2_trial_start, data = codata_combo2)

# alternate way using update() with add_data argument for
# new observations only (those collected after the trial
# design stage).
new_data <- filter(codata_combo2, cohort_time > 0)

combo2_trial <- update(combo2_trial_start, add_data = new_data)

summary(combo2_trial, "data") # cohort_time is tracked
summary(combo2_trial, "data_prediction")
summary(combo2_trial, "dose_prediction")

rm(dims, new_data)

## Recover user set sampling defaults
options(.user_mc_options)

}
\seealso{
Other blrm_trial combo2 example: 
\code{\link{blrm_trial}()},
\code{\link{dose_info_combo2}},
\code{\link{drug_info_combo2}}
}
\concept{blrm_trial combo2 example}
