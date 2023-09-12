###### Step 0: Load the package
# library(BOIN)
# library(OncoBayes2)
# library(readxl)
# library(tidyverse)

####################################
###### At the design stage ######
####################################

###### Step 2: Step up initial design
set.seed(123)
.user_mc_options <- options(OncoBayes2.MC.warmup=5000, OncoBayes2.MC.iter=20000, OncoBayes2.MC.chains=10,
                            OncoBayes2.MC.save_warmup=FALSE,mc.cores = 10)

SA_trial_setup <- blrm_trial(
    data = tibble(group_id = as.factor("All"),
                  drug_A = 50,
                  num_patients = NA,
                  num_toxicities = NA),
    # dose-toxicity data available at design stage of trial.
    
    dose_info = tibble(group_id = as.factor("All"), 
                       drug_A = c(50,75,150,300,450,600,800,1000), 
                       dose_id = c(1, 2, 3, 4, 5, 6, 7, 8), stratum_id = "all"), 
    # specification of the dose levels as planned for the ongoing trial arms.
    
    drug_info = tibble(drug_name = "drug_A", 
                       dose_ref = 300, 
                       dose_unit = "ug/kg"
                       ,
                       reference_p_dlt = exp(-3.068) / (1 + exp(-3.068))
    )
    # specification of drugs used in trial arms
    
    ,simplified_prior = FALSE
    ,EXNEX_comp=TRUE
    ,EXNEX_inter=FALSE
    ,interval_prob = c(0,0.16,0.33,1)
    ,interval_max_mass = c(prob_underdose = 1 # The prob of under-dose is allowed maximum to 1.
                           , prob_target = 1 # The prob of target-dose is allowed maximum to 1.
                           , prob_overdose = 0.25) # The prob of over-dose is allowed maximum to 0.28.
    # ,cores = getOption("mc.cores", 4)
)


dims <- summary(SA_trial_setup,"dimensionality")
num_comp <- dims$num_components

SA_trial_start <- update(
    SA_trial_setup,
    ##component1 MAP prior
    prior_EX_mu_mean_comp = matrix(
        c(-3.068, # mean of intercept 
          0.564), # mean of log-slope 
        nrow = num_comp,
        ncol = 2
    ),
    prior_EX_mu_sd_comp = matrix(
        c(2.706, # sd of intercept
          0.728), # sd of log-slope
        nrow = num_comp,
        ncol = 2
    ),
    prior_EX_corr_mu_comp = -0.817, 
    prior_EX_tau_mean_comp = matrix(
        c(0, 0),
        nrow = num_comp,
        ncol = 2
    ),
    prior_EX_tau_sd_comp = matrix(
        c(1, 1),
        nrow = num_comp,
        ncol = 2
    ),
    # prior_EX_corr_eta_comp = 0.817,
    # prior_EX_corr_eta_inter = 0.817,
    # prior_EX_mu_mean_inter = -0.817,
    
    #component2: non-informative prior
    # prior_NEX_mu_mean_comp = matrix(
    #   c(logit(0.25),
    #     0),
    #   nrow=num_comp,
    #   ncol=2
    # ),
    # prior_NEX_mu_sd_comp = matrix(
    #   c(2,
    #     1),
    #   nrow=num_comp,
    #   ncol=2
    # ),
    prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
    prior_tau_dist = 0, #1=log-normal
    prior_PD = FALSE
)

prior_summary(SA_trial_start)    # Need to update summary to refect correlation information

#Some summary of posterior,only for the planned dose prediction not the same as the data_prediction.
#The data_prediction is for the emergent data. dose_prediction will return prob_underdose, prob_target, 
#prob_overdose and ewoc_ok
summary(SA_trial_start, "dose_prediction")

summary(SA_trial_start, "ewoc_check")

