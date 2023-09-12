
<!-- README.md is generated from README.Rmd. Please edit that file -->

# My modified OncoBayes2

<!-- badges: start -->
<!-- badges: end -->

This is my modified version of the [original
OncoBayes2](https://CRAN.R-project.org/package=OncoBayes2) package.

The most significant difference is the support of providing correlation
parameter when sampling the component via the newly added parameter
`prior_EX_corr_mu_comp` in `blrm_exnex()`.

In the original package,
$\theta = (log \alpha, log \beta) \sim N(\mu, \Sigma)$ where the
**covariance** matrix $\Sigma$ is a diagonal matrix whose diagonals are
provided via `prior_EX_mu_sd_comp`(Note, this parameter provides the
**standard variance**).

Now, the off-diagonal part of $\Sigma$ can also be specified by
supplying the **correlation** via `prior_EX_corr_mu_comp`.

## NOTE

1.  This is just a **proof of concept** during my replication of Table
    14-7 in this [trial
    protocol](https://clinicaltrials.gov/ProvidedDocs/64/NCT02108964/Prot_001.pdf).
    The original OncoBayes2 misses the functionality to accept user
    defined correlation, which is also reported on the
    [internet](https://stats.stackexchange.com/questions/617633/correlation-for-priors-in-blrm-with-ewoc-in-r-package-oncobayes2/626132#626132)

2.  Currently this package is forked and modified from [CRAN’s github
    version of OncoBayes2](https://github.com/cran/OncoBayes2). After
    installation it will **overwrite** your existing OncoBayes2.

3.  This is currently WIP. Many aspect, like the `prior_summary` method
    should be updated accordingly to reflect the correlation
    information. Also many document template is missing hence the manual
    cannot be auto-updated properly.

## Installation

You can install the development version of OncoBayes2 like so:

``` r
remotes::install_github("fenguoerbian/MyOncoBayes2")
```

## Example

This is a basic example which reproduce the results in Table 14-7 in
this [trial
protocol](https://clinicaltrials.gov/ProvidedDocs/64/NCT02108964/Prot_001.pdf).

``` r
library(OncoBayes2)
#> This is OncoBayes2 version 0.8.123456 (released 2023-07-21, git-sha e5a8128)
## basic example code
```

### Setup MCMC options

``` r
set.seed(123)
.user_mc_options <- options(OncoBayes2.MC.warmup=500, OncoBayes2.MC.iter=20000, OncoBayes2.MC.chains=10,
                            OncoBayes2.MC.save_warmup=FALSE, mc.cores = 10)
```

### Setup basic trial infomation

``` r
SA_trial_setup <- blrm_trial(
    data = tibble::tibble(group_id = as.factor("All"),
                  drug_A = 50,
                  num_patients = NA,
                  num_toxicities = NA),
    # dose-toxicity data available at design stage of trial.
    
    dose_info = tibble::tibble(group_id = as.factor("All"), 
                       drug_A = c(50,75,150,300,450,600,800,1000), 
                       dose_id = c(1, 2, 3, 4, 5, 6, 7, 8), stratum_id = "all"), 
    # specification of the dose levels as planned for the ongoing trial arms.
    
    drug_info = tibble::tibble(drug_name = "drug_A", 
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
)
#> Please configure blrm_exnex using the update() function.


dims <- summary(SA_trial_setup,"dimensionality")
num_comp <- dims$num_components
```

### Update the prior

``` r
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
    prior_EX_corr_mu_comp = -0.817,    # HERE! supply the CORRELATION!
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
    prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
    prior_tau_dist = 0, # single-agent, without historical data, set this to 0.
                        # see single-agent example in this package for more details
    prior_PD = FALSE
)
```

**NOTE:** `prior_summary` needs to be updated to reflect user supplied
correlation

``` r
prior_summary(SA_trial_start)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Mixture configuration
#> ---------------------
#> EXNEX components : 0 
#> component
#> I(log(drug_A/300)) 
#>                  0 
#> 
#> Prior probability for exchangeability per group
#>      component
#> group I(log(drug_A/300))
#>   All                  1
#> 
#> EXchangable hyperparameter priors
#> ---------------------------------
#> Component parameters
#> Mean mu_log_beta
#>                    coefficient intercept       log_slope      
#>                    prior            mean    sd      mean    sd
#> component                                                     
#> I(log(drug_A/300))                 -3.07  2.71      0.56  0.73
#> 
#> Heterogeneity tau_log_beta (known)
#>                    coefficient intercept    log_slope   
#>                    stratum           all          all   
#>                    prior            mean sd      mean sd
#> component                                               
#> I(log(drug_A/300))                     0  1         0  1
#> 
#> Correlation LKJ
#> component
#> I(log(drug_A/300)) 
#>                  1 
#> 
#> Model has no interaction parameters.
#> 
#> NonEXchangable priors
#> ---------------------
#> Component parameters
#> Mean mu_log_beta
#>                    coefficient intercept       log_slope      
#>                    prior            mean    sd      mean    sd
#> component                                                     
#> I(log(drug_A/300))                 -3.07  2.71      0.56  0.73
#> 
#> Model has no interaction parameters.
```

### Check the results against Table 14-7

``` r
summary(SA_trial_start, "dose_prediction")
#> # A tibble: 8 × 13
#>   group_id drug_A dose_id stratum_id   mean    sd   `2.5%`   `50%` `97.5%`
#>   <fct>     <dbl>   <dbl> <fct>       <dbl> <dbl>    <dbl>   <dbl>   <dbl>
#> 1 All          50       1 all        0.0776 0.190 8.75e-10 0.00183   0.786
#> 2 All          75       2 all        0.0905 0.202 1.61e- 8 0.00372   0.819
#> 3 All         150       3 all        0.121  0.226 2.23e- 6 0.0128    0.868
#> 4 All         300       4 all        0.172  0.254 2.27e- 4 0.0447    0.907
#> 5 All         450       5 all        0.223  0.269 2.02e- 3 0.0961    0.926
#> 6 All         600       6 all        0.282  0.279 5.91e- 3 0.172     0.939
#> 7 All         800       7 all        0.366  0.291 1.37e- 2 0.289     0.959
#> 8 All        1000       8 all        0.441  0.298 2.37e- 2 0.400     0.975
#> # ℹ 4 more variables: prob_underdose <dbl>, prob_target <dbl>,
#> #   prob_overdose <dbl>, ewoc_ok <lgl>
```

``` r
summary(SA_trial_start, "ewoc_check")
#> # A tibble: 8 × 9
#>   group_id drug_A dose_id stratum_id prob_overdose_est prob_overdose_stat
#>   <fct>     <dbl>   <dbl> <fct>                  <dbl>              <dbl>
#> 1 All          50       1 all                   0.0336            -609.  
#> 2 All          75       2 all                   0.0525            -379.  
#> 3 All         150       3 all                   0.110             -172.  
#> 4 All         300       4 all                   0.225              -57.4 
#> 5 All         450       5 all                   0.335                2.33
#> 6 All         600       6 all                   0.447               58.5 
#> 7 All         800       7 all                   0.589              149.  
#> 8 All        1000       8 all                   0.696              258.  
#> # ℹ 3 more variables: prob_overdose_mcse <dbl>, prob_overdose_ess <dbl>,
#> #   prob_overdose_rhat <dbl>
```
