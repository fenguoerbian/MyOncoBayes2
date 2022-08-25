## ---- include=FALSE-----------------------------------------------------------
library(OncoBayes2)
library(knitr)

ggplot2::theme_set(bayesplot::bayesplot_theme_get())

knitr::opts_chunk$set(
fig.width = 1.62*4,
fig.height = 4,
warning=FALSE,
message=FALSE
)
## setup up fast sampling when run on CRAN
is_CRAN <- !identical(Sys.getenv("NOT_CRAN", "false"), "true")
## NOTE: for running this vignette locally, please uncomment the
## following line:
##is_CRAN <- FALSE
.user_mc_options <- list()

if (is_CRAN) {
.user_mc_options <- options(OncoBayes2.MC.warmup=10, OncoBayes2.MC.iter=20, OncoBayes2.MC.chains=1, OncoBayes2.MC.save_warmup=FALSE, mc.cores=1)
} else {
.user_mc_options <- options(OncoBayes2.MC.warmup=500, OncoBayes2.MC.iter=1000, OncoBayes2.MC.chains=4, OncoBayes2.MC.save_warmup=FALSE, mc.cores=1)
}


## ---- message = FALSE---------------------------------------------------------
## Load involved packages
library(dplyr)   ## for mutate
library(tidyr)   ## defines expand_grid
library(tibble)  ## for tibbles
library(ggplot2) ## for plotting

## -----------------------------------------------------------------------------
kable(hist_combo2)

## ----echo=TRUE----------------------------------------------------------------
levels(hist_combo2$group_id)

## -----------------------------------------------------------------------------
kable(drug_info_combo2)

## -----------------------------------------------------------------------------
dose_info <- filter(dose_info_combo2, group_id == "trial_AB",
                    drug_A %in% c(3,6), drug_B %in% c(0,400, 800))
kable(dose_info)

## -----------------------------------------------------------------------------
combo2_trial_setup <- blrm_trial(
  data = hist_combo2,
  drug_info = drug_info_combo2,
  dose_info = dose_info
)

## ---- message = FALSE, echo = TRUE, results = "hide"--------------------------
combo2_trial_start <- blrm_trial(
  data = hist_combo2,
  drug_info = drug_info_combo2,
  dose_info = dose_info,
  simplified_prior = TRUE,
  EXNEX_comp=FALSE,
  EX_prob_comp_hist=1,
  EX_prob_comp_new=1
)

## ---- eval = FALSE------------------------------------------------------------
#  prior_summary(combo2_trial_start) # not run here

## -----------------------------------------------------------------------------
kable(summary(combo2_trial_start, "dose_prediction"), digits = 2)

## -----------------------------------------------------------------------------
kable(summary(combo2_trial_start, "ewoc_check"), digits = 3)

## ---- include=FALSE-----------------------------------------------------------
po <- summary(combo2_trial_start, "ewoc_check")$prob_overdose_stat
min_stat <- po[which.min(abs(po))]

## -----------------------------------------------------------------------------
candidate_starting_dose <- summary(combo2_trial_start, "dose_info") %>%
  filter(drug_A == 3, drug_B == 400) %>%
  crossing(num_toxicities = 0, num_patients = 3:6)

pp_summary <- summary(combo2_trial_start, interval_prob = c(-1, 0, 1, 6), predictive = TRUE,
                      newdata = candidate_starting_dose)

kable(bind_cols(select(candidate_starting_dose, num_patients),
                select(pp_summary, ends_with("]"))), digits = 3)

## -----------------------------------------------------------------------------
new_cohort <- tibble(group_id = "trial_AB",
                     drug_A = 3,
                     drug_B = 400,
                     num_patients = 5,
                     num_toxicities = 1)

## ---- message = FALSE, echo = TRUE, results = "hide"--------------------------
combo2_trial_update <- update(combo2_trial_start, add_data = new_cohort)

## -----------------------------------------------------------------------------
kable(summary(combo2_trial_update, "dose_prediction"), digits = 2)

## -----------------------------------------------------------------------------
kable(summary(combo2_trial_update, "newdata_prediction",
              newdata = tibble(group_id = "trial_AB",
                               drug_A = 4.5,
                               drug_B = c(400, 600, 800))), digits = 2)

## ---- message = FALSE, echo = TRUE, results = "hide"--------------------------
# set up two scenarios at the starting dose level
# store them as data frames in a named list
scenarios <- expand_grid(
  group_id  = "trial_AB",
  drug_A = 3,
  drug_B = 800,
  num_patients = 3,
  num_toxicities = 0:2
) %>% split(1:3) %>% setNames(paste0(0:2, "/3 DLTs"))

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
}) %>%
  bind_rows(.id="Scenario")


## ---- echo = FALSE------------------------------------------------------------
kable(select(scenario_inference, -group_id, -stratum_id, -dose_id), digits = 2,
      caption = "Model inference for trial AB when varying hypothetical DLT scenarios for a cohort of size 3")

## -----------------------------------------------------------------------------
trial_AB_data <- filter(codata_combo2, group_id == "trial_AB", cohort_time==1)
kable(trial_AB_data)

## ---- message = FALSE, echo = TRUE, results = "hide"--------------------------
combo2_trial_histdata <- update(combo2_trial_start, add_data = trial_AB_data)

## -----------------------------------------------------------------------------
trial_A_codata <- filter(codata_combo2, group_id == "trial_A", cohort_time==1)
kable(trial_A_codata)

## ---- message = FALSE, echo = TRUE, results = "hide"--------------------------
combo2_trial_codata <- update(combo2_trial_histdata, add_data = trial_A_codata)

## ---- fig.height = 1.2 * 4, fig.width=1.62 * 4--------------------------------

plot_toxicity_intervals_stacked(combo2_trial_histdata,
                                newdata=mutate(dose_info, dose_id=NULL, stratum_id="all"),
                                x = vars(drug_B),
                                group = vars(drug_A),
                                facet_args = list(ncol = 1)
) + ggtitle("Trial AB with historical data only")

plot_toxicity_intervals_stacked(combo2_trial_codata,
                                newdata=mutate(dose_info, dose_id=NULL, stratum_id="all"),
                                x = vars(drug_B),
                                group = vars(drug_A),
                                facet_args = list(ncol = 1)
) + ggtitle("Trial AB with historical and concurrent data on drug A")


## -----------------------------------------------------------------------------
trial_AB_stage_2_codata <- filter(codata_combo2, cohort_time==2)
kable(trial_AB_stage_2_codata)

## ---- message = FALSE, echo = TRUE, results = "hide"--------------------------
combo2_trial_final <- update(combo2_trial_start, data = codata_combo2)

## ---- fig.height = 1.05 * 4, fig.width=1.62 * 4-------------------------------

grid_length <- 25

dose_info_plot_grid <- expand_grid(stratum_id = "all",
                                   group_id = "trial_AB",
                                   drug_A=seq(min(dose_info_combo2$drug_A), max(dose_info_combo2$drug_A), length.out=grid_length),
                                   drug_B=seq(min(dose_info_combo2$drug_B), max(dose_info_combo2$drug_B), length.out=grid_length))


dose_info_plot_grid_sum <- summary(combo2_trial_final,
                                   newdata=dose_info_plot_grid,
                                   prob=0.5)

ggplot(dose_info_plot_grid_sum, aes(drug_A, drug_B, z = !!as.name("75%"))) +
  geom_contour_filled(breaks=c(0, 0.1, 0.16, 0.33, 1)) +
  scale_fill_brewer("Quantile Range", type="div", palette = "RdBu", direction=-1) +
  ggtitle("DLT Probability 75% Quantile")


## -----------------------------------------------------------------------------
sessionInfo()

## ---- include=FALSE-----------------------------------------------------------
## restore previous global user options
options(.user_mc_options)

