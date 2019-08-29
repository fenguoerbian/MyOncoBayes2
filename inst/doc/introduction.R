## ---- include=FALSE------------------------------------------------------
library(OncoBayes2)
library(RBesT)
library(knitr)
library(ggplot2)
library(dplyr)
library(tidybayes)

theme_set(theme_bw())

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

## ------------------------------------------------------------------------
kable(hist_combo2)

## ------------------------------------------------------------------------
levels(hist_combo2$group_id)

## ---- eval = FALSE-------------------------------------------------------
#  # Design parameters ---------------------
#  
#  dref <- c(3, 960)
#  num_comp <- 2 # two investigational drugs
#  num_inter <- 1 # one drug-drug interaction needs to be modeled
#  num_groups <- nlevels(hist_combo2$group_id) # four groups of data
#  num_strata <- 1 # no stratification needed
#  
#  # Model fit -----------------------------
#  
#  blrmfit <- blrm_exnex(
#    cbind(Ntox, Npat - Ntox) ~
#      1 + I(log(DosesAdm1 / dref[1])) |
#      1 + I(log(DosesAdm2 / dref[2])) |
#      0 + I(DosesAdm1/dref[1] *DosesAdm2/dref[2]) |
#      group_id,
#    data = hist_combo2,
#    prior_EX_mu_mean_comp = matrix(
#      c(logit(0.1), 0, # hyper-mean of (intercept, log-slope) for drug A
#        logit(0.1), 0), # hyper-mean of (intercept, log-slope) for drug B
#      nrow = num_comp,
#      ncol = 2,
#      byrow = TRUE
#    ),
#    prior_EX_mu_sd_comp = matrix(
#      c(3.33, 1, # hyper-sd of mean mu for (intercept, log-slope) for drug B
#        3.33, 1), # hyper-sd of mean mu for (intercept, log-slope) for drug B
#      nrow = num_comp,
#      ncol = 2,
#      byrow = TRUE
#    ),
#    prior_EX_tau_mean_comp = matrix(
#      c(log(0.25), log(0.125),
#        log(0.25), log(0.125)),
#      nrow = num_comp,
#      ncol = 2,
#      byrow = TRUE
#    ),
#    prior_EX_tau_sd_comp = matrix(
#      c(log(4) / 1.96, log(4) / 1.96,
#        log(4) / 1.96, log(4) / 1.96),
#      nrow = num_comp,
#      ncol = 2,
#      byrow = TRUE
#    ),
#    prior_EX_mu_mean_inter = 0,
#    prior_EX_mu_sd_inter = 1.121,
#    prior_EX_tau_mean_inter = matrix(log(0.125), nrow = num_inter, ncol = num_strata),
#    prior_EX_tau_sd_inter = matrix(log(4) / 1.96, nrow = num_inter, ncol = num_strata),
#    prior_is_EXNEX_comp = rep(FALSE, num_comp),
#    prior_is_EXNEX_inter = rep(FALSE, num_inter),
#    prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
#    prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
#    prior_tau_dist = 1
#  )

## ---- include = FALSE----------------------------------------------------
# Design parameters ---------------------

dref <- c(3, 960)
num_comp <- 2 # two investigational drugs
num_inter <- 1 # one drug-drug interaction needs to be modeled
num_groups <- nlevels(hist_combo2$group_id) # four groups of data
num_strata <- 1 # no stratification needed

# Model fit -----------------------------

blrmfit <- blrm_exnex(
  cbind(Ntox, Npat - Ntox) ~
    1 + I(log(DosesAdm1 / dref[1])) |
    1 + I(log(DosesAdm2 / dref[2])) |
    0 + I(DosesAdm1/dref[1] *DosesAdm2/dref[2]) |
    group_id,
  data = hist_combo2,
  prior_EX_mu_mean_comp = matrix(
    c(logit(0.1), 0, # hyper-mean of (intercept, log-slope) for drug A
      logit(0.1), 0), # hyper-mean of (intercept, log-slope) for drug B
    nrow = num_comp,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_mu_sd_comp = matrix(
    c(3.33, 1, # hyper-sd of mean mu for (intercept, log-slope) for drug B
      3.33, 1), # hyper-sd of mean mu for (intercept, log-slope) for drug B
    nrow = num_comp,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_tau_mean_comp = matrix(
    c(log(0.25), log(0.125),
      log(0.25), log(0.125)),
    nrow = num_comp,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_tau_sd_comp = matrix(
    c(log(4) / 1.96, log(4) / 1.96,
      log(4) / 1.96, log(4) / 1.96),
    nrow = num_comp,
    ncol = 2,
    byrow = TRUE
  ),
  prior_EX_mu_mean_inter = 0,
  prior_EX_mu_sd_inter = 1.121,
  prior_EX_tau_mean_inter = matrix(log(0.125), nrow = num_inter, ncol = num_strata),
  prior_EX_tau_sd_inter = matrix(log(4) / 1.96, nrow = num_inter, ncol = num_strata),
  prior_is_EXNEX_comp = rep(FALSE, num_comp),
  prior_is_EXNEX_inter = rep(FALSE, num_inter),
  prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
  prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
  prior_tau_dist = 1
)

## ---- eval = FALSE-------------------------------------------------------
#  prior_summary(blrmfit) # not run here

## ------------------------------------------------------------------------
newdata <- expand.grid(
  group_id = c("trial_AB"),
  DosesAdm1 = c(0, 3, 4.5, 6),
  DosesAdm2 = c(0, 400, 600, 800),
  stringsAsFactors = FALSE
)
newdata$group_id <- factor(newdata$group_id, levels(hist_combo2$group_id))

## ------------------------------------------------------------------------
newdata <- expand.grid(
  group_id = factor(c("trial_AB"), levels(hist_combo2$group_id)),
  DosesAdm1 = c(0, 3, 4.5, 6),
  DosesAdm2 = c(0, 400, 600, 800)
)

## ------------------------------------------------------------------------
summ_stats <- summary(blrmfit, 
                      newdata = newdata,
                      prob = 0.95,
                      interval_prob = c(0, 0.16, 0.33, 1))

kable(cbind(newdata, summ_stats), digits = 3)

## ---- eval = FALSE-------------------------------------------------------
#  # set up two scenarios at the starting dose level
#  # store them as data frames in a named list
#  scenarios <- expand.grid(
#    group_id  = factor("trial_AB", levels(hist_combo2$group_id)),
#    DosesAdm1 = 3,
#    DosesAdm2 = 400,
#    Npat      = 3,
#    Ntox      = 1:2
#  ) %>% split(1:2) %>% setNames(paste(1:2, "DLTs"))
#  
#  candidate_doses = expand.grid(
#    group_id = factor("trial_AB", levels(hist_combo2$group_id)),
#    DosesAdm1 = c(3, 4.5),
#    DosesAdm2 = 400
#  )
#  
#  scenario_inference <- lapply(scenarios, function(scen_row){
#  
#    # refit the model with each scenario's additional data
#    scenario_data <- rbind(hist_combo2, scen_row)
#    scenario_fit <- update(blrmfit, data = scenario_data)
#  
#    # summarize posterior at candidate doses
#    scenario_summ <- summary(scenario_fit,
#                             newdata = candidate_doses,
#                             interval_prob = c(0, 0.16, 0.33, 1))
#  
#    cbind(candidate_doses, scenario_summ)
#  
#  })
#  

## ---- include = FALSE----------------------------------------------------
# set up two scenarios at the starting dose level
# store them as data frames in a named list
scenarios <- expand.grid(
  group_id  = factor("trial_AB", levels(hist_combo2$group_id)),
  DosesAdm1 = 3,
  DosesAdm2 = 400,
  Npat      = 3,
  Ntox      = 1:2
) %>% split(1:2) %>% setNames(paste(1:2, "DLTs"))

candidate_doses = expand.grid(
  group_id = factor("trial_AB", levels(hist_combo2$group_id)),
  DosesAdm1 = c(3, 4.5),
  DosesAdm2 = 400
)

scenario_inference <- lapply(scenarios, function(scen_row){
  
  # refit the model with each scenario's additional data
  scenario_data <- rbind(hist_combo2, scen_row)
  scenario_fit <- update(blrmfit, data = scenario_data)
  
  # summarize posterior at candidate doses
  scenario_summ <- summary(scenario_fit,
                           newdata = candidate_doses,
                           interval_prob = c(0, 0.16, 0.33, 1))
  
  cbind(candidate_doses, scenario_summ)
  
})


## ---- echo = FALSE-------------------------------------------------------
kable(scenario_inference[["1 DLTs"]], digits = 3,
      caption = "Model inference when 1 DLT is observed in first cohort")

## ---- echo = FALSE-------------------------------------------------------
kable(scenario_inference[["2 DLTs"]], digits = 3,
      caption = "Model inference when 2 DLTs are observed in first cohort")

## ------------------------------------------------------------------------
kable(codata_combo2)

## ---- include = FALSE----------------------------------------------------
final_fit <- update(blrmfit, data = codata_combo2)

summ <- summary(final_fit, newdata, prob = c(0.5, 0.95), interval_prob = c(0,0.33,1))

final_summ_stats <- cbind(newdata, summ) %>%
    mutate(EWOC=1*`(0.33,1]`<=0.25)

## ---- eval = FALSE-------------------------------------------------------
#  final_fit <- update(blrmfit, data = codata_combo2)
#  
#  summ <- summary(final_fit, newdata, prob = c(0.5, 0.95), interval_prob = c(0,0.33,1))
#  
#  final_summ_stats <- cbind(newdata, summ) %>%
#      mutate(EWOC=1*`(0.33,1]`<=0.25)

## ---- fig.height = 4 * 1.62----------------------------------------------

ggplot(final_summ_stats,
       aes(x=factor(DosesAdm2), colour=EWOC)) +
    facet_wrap(~DosesAdm1, labeller=label_both) +
    scale_y_continuous(breaks=c(0, 0.16, 0.33, 0.4, 0.6, 0.8, 1.0)) +
    coord_cartesian(ylim=c(0,0.8)) +
    geom_hline(yintercept = c(0.16, 0.33),
               linetype = "dotted") +
    geom_pointrange(aes(y=`50%`, ymin=`2.5%`, ymax=`97.5%`)) +
    geom_linerange(aes(ymin=`25%`, ymax=`75%`), size=1.5) +
    ggtitle("DLT Probability", "Shown is the median (dot), 50% CrI (thick line) and 95% CrI (thin line)") +
    ylab(NULL) + 
    xlab("Dose Drug B [mg]")


## ---- fig.height = 4 * 1.62----------------------------------------------
library(tidybayes)

expand.grid(group_id=factor("trial_AB", levels=levels(codata_combo2$group_id)),
            DosesAdm2=exp(seq(log(100),log(800),length=100)), DosesAdm1=c(0, 3, 4.5, 6)) %>%
    add_draws(posterior_linpred(final_fit, newdata = ., transform=TRUE)) %>%
    median_qi(.width=c(0.5, 0.95)) %>%
    ggplot(aes(y = .value, x = DosesAdm2)) +
    scale_x_log10(breaks=c(100,200,400,600,800)) +
    facet_wrap(~DosesAdm1, labeller=label_both) +
    scale_y_continuous(breaks=c(0, 0.16, 0.33, 0.4, 0.6, 0.8, 1.0)) +
    geom_lineribbon() +
    scale_fill_brewer() +
    coord_cartesian(ylim=c(0,0.8)) + 
    geom_hline(yintercept = c(0.16, 0.33),
               linetype = "dotted") +
    ggtitle("DLT Probability", "Shown is the median (line), 50% CrI (dark) and 95% CrI (light)") +
    ylab(NULL) + 
    xlab("Dose Drug B [mg]")


## ------------------------------------------------------------------------
sessionInfo()

## ---- include=FALSE------------------------------------------------------
## restore previous global user options
options(.user_mc_options)

