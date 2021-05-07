## ---- include=FALSE------------------------------------------------------
library(OncoBayes2)
library(knitr)
library(ggplot2)
  
theme_set(bayesplot::bayesplot_theme_get())
  
knitr::opts_chunk$set(
  fig.width = 1.62*4,
  fig.height = 4
)
## setup up fast sampling when run on CRAN
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
## NOTE: for running this vignette locally, please uncomment the
## following line:
## is_CRAN <- TRUE
.user_mc_options <- list()

if (is_CRAN) {
  .user_mc_options <- options(OncoBayes2.MC.warmup=10, OncoBayes2.MC.iter=20, OncoBayes2.MC.chains=1)
}

options(warn = -1)


## ---- message = FALSE----------------------------------------------------
## Load involved packages
library(dplyr)  ## for mutate
library(tidyr)  ## defines expand_grid
library(tibble) ## for tibbles

## ------------------------------------------------------------------------
kable(hist_combo2)

## ------------------------------------------------------------------------
levels(hist_combo2$group_id)

## ------------------------------------------------------------------------
kable(drug_info_combo2)

## ------------------------------------------------------------------------
dose_info <- filter(dose_info_combo2, group_id == "trial_AB", drug_A == 3)
kable(dose_info)

## ------------------------------------------------------------------------
combo2_trial_setup <- blrm_trial(
  data = hist_combo2,
  drug_info = drug_info_combo2,
  dose_info = dose_info,
  simplified_prior = FALSE
)

## ---- message = FALSE, echo = TRUE, results = "hide"---------------------
combo2_trial_start <- blrm_trial(
  data = hist_combo2,
  drug_info = drug_info_combo2,
  dose_info = dose_info,
  simplified_prior = TRUE,
  EX_prob_comp_hist = 0.8,
  EX_prob_comp_new = 1
)

## ---- eval = FALSE-------------------------------------------------------
#  prior_summary(combo2_trial_start) # not run here

## ------------------------------------------------------------------------
kable(summary(combo2_trial_start, "dose_prediction"), digits = 2)

## ------------------------------------------------------------------------
candidate_starting_dose <- summary(combo2_trial_start, "dose_info") %>%
  filter(drug_A == 3, drug_B == 400) %>%
  expand_grid(num_patients = 3:6) %>%
  mutate(num_toxicities = 0)

pp_summary <- summary(combo2_trial_start, interval_prob = c(-1, 0, 1, 6), predictive = TRUE,
                      newdata = candidate_starting_dose)

kable(bind_cols(select(candidate_starting_dose, num_patients),
                select(pp_summary, ends_with("]"))), digits = 3)

## ------------------------------------------------------------------------
new_cohort <- tibble(group_id = "trial_AB",
                     drug_A = 3,
                     drug_B = 400,
                     num_patients = 5,
                     num_toxicities = 1)

## ---- message = FALSE, echo = TRUE, results = "hide"---------------------
combo2_trial_update <- update(combo2_trial_start, add_data = new_cohort)

## ------------------------------------------------------------------------
kable(summary(combo2_trial_update, "newdata_prediction",
              newdata = tibble(group_id = "trial_AB",
                               drug_A = 4.5,
                               drug_B = c(400, 600, 800))), digits = 2)

## ---- message = FALSE, echo = TRUE, results = "hide"---------------------
# set up two scenarios at the starting dose level
# store them as data frames in a named list
scenarios <- expand_grid(
  group_id  = "trial_AB",
  drug_A = 3,
  drug_B = 800,
  num_patients = 3,
  num_toxicities = 0:2
) %>% split(1:3) %>% setNames(paste(0:2, "DLTs"))

candidate_doses <- expand_grid(
  group_id = "trial_AB",
  drug_A = c(3, 4.5),
  drug_B = c(600, 800)
)


scenario_inference <- lapply(scenarios, function(scenario_newdata) {
  
  # refit the model with each scenario's additional data
  scenario_fit <- update(combo2_trial_update, add_data = scenario_newdata)
  
  # summarize posterior at candidate doses
  summary(scenario_fit, "newdata_prediction", newdata = candidate_doses)
  
})


## ---- echo = FALSE-------------------------------------------------------
kable(scenario_inference[["0 DLTs"]], digits = 2,
      caption = "Model inference when 0 DLTs are observed in the next cohort")

## ---- echo = FALSE-------------------------------------------------------
kable(scenario_inference[["1 DLTs"]], digits = 2,
      caption = "Model inference when 1 DLT is observed in the next cohort")

## ---- echo = FALSE-------------------------------------------------------
kable(scenario_inference[["2 DLTs"]], digits = 2,
      caption = "Model inference when 2 DLTs are observed in the next cohort")

## ------------------------------------------------------------------------
kable(filter(codata_combo2, cohort_time > 0))

## ---- message = FALSE, echo = TRUE, results = "hide"---------------------
final_fit <- update(combo2_trial_start, data = codata_combo2)

## ---- fig.height = 1.05 * 4, fig.width=1.62 * 4--------------------------
plot_toxicity_curve(final_fit, x = vars(drug_B), group = vars(group_id, drug_A),
                    facet_args = list(ncol = 1))


## ------------------------------------------------------------------------
sessionInfo()

## ---- include=FALSE------------------------------------------------------
## restore previous global user options
options(.user_mc_options)

