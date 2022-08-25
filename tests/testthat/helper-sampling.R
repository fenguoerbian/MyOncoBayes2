library(checkmate)
library(assertthat)
library(Formula)
library(abind)
library(rstan)

set_sampling_default  <- function(iter, warmup, chains, cores=getOption("mc.cores", 1), save_warmup=FALSE, backend="rstan") {
    options(OncoBayes2.MC.iter=iter,
            OncoBayes2.MC.warmup=warmup,
            OncoBayes2.MC.chains=chains,
            mc.cores=cores,
            OncoBayes2.MC.save_warmup=save_warmup,
            OncoBayes2.MC.backend=backend)
}

very_fast_sampling <- function() {
    message("Tests running with very fast sampling")
    ## note: 250 warmups are needed to get Stans NUTS adaptation to work
    set_sampling_default(500, 250, 1, 1)
}

fake_sampling <- function() {
    message("Tests running with fake sampling")
    set_sampling_default(4, 2, 1, 1)
}

default_sampling <- function() {
    set_sampling_default(NULL, NULL, NULL, NULL)
}

run_example <- function(example) {
    env <- new.env()
    suppressWarnings(example_model(example, env, silent=TRUE))
    invisible(env)
}



## set up slim sampling in case we are on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    very_fast_sampling()
} else {
    fake_sampling()
}

## takes a fitted blrmfit object and returns the blrmfit object with a
## posterior containing one draw with all slopes and intercept set to
## the prior mean (exactly).
sample_prior_mean <- function(blrmfit) {
    ps <- prior_summary(blrmfit)

    num_groups <- ps$num_groups
    num_strata <- ps$num_strata
    num_comp <- dim(ps$EX_mu_log_beta)[2]
    has_inter <- ps$has_inter
    num_inter <- 0
    if(has_inter)
        num_inter <- dim(ps$EX_mu_eta)[1]

    draw <- list(
        log_beta_raw=array(0, c(2*num_groups, num_comp, 2)),
        eta_raw=array(0, c(2*num_groups, num_inter)),
        mu_log_beta=array(0, c(num_comp, 2)),
        tau_log_beta_raw=array(1, c(num_strata, num_comp, 2)),
        L_corr_log_beta=abind(replicate(num_comp, diag(2), FALSE), along=-1),
        mu_eta=array(0, c(num_inter)),
        tau_eta_raw=array(1, c(num_strata, num_inter)),
        L_corr_eta=diag(num_inter)
    )

    draw$log_beta_raw <- abind(c(replicate(num_groups, aperm(ps$EX_mu_log_beta[,,"mean",drop=FALSE], c(3,2,1)), FALSE),
                                 replicate(num_groups, aperm(ps$NEX_mu_log_beta[,,"mean",drop=FALSE], c(3,2,1)), FALSE)), along=1)

    if(has_inter) {
        draw$eta_raw <- abind(c(replicate(num_groups, aperm(ps$EX_mu_eta[,"mean",drop=FALSE], c(2,1)), FALSE),
                                replicate(num_groups, aperm(ps$NEX_mu_eta[,"mean",drop=FALSE], c(2,1)), FALSE)), along=1)
    }

    msg <- capture.output(blrmfit$stanfit <- sampling(OncoBayes2:::stanmodels$blrm_exnex, data=blrmfit$standata, chains=1, iter=1, warmup=0, seed=23542, init=list(draw), algorithm="Fixed_param", open_progress=FALSE))
    blrmfit
}
