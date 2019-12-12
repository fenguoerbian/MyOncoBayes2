#'
#' Utilities for SBC validation
#'
library(devtools)
devtools::load_all("../..")
library(rstan)
library(RBesT)
library(mvtnorm)
library(checkmate) # only needed inside blrm_exnex
library(Formula) # only needed inside blrm_exnex
library(abind) # only needed inside blrm_exnex
library(dplyr)
library(tidyr)
library(assertthat)
source("lkj.R")

options(OncoBayes2.abbreviate.min = 0)

## Sample prior

## First sample EX hyperparameters
## prior_EX_mu_mean_comp + prior_EX_mu_sd_comp => EX_mu_comp (group means EX)
## prior_EX_tau_mean_comp + prior_EX_tau_sd_comp => EX_tau_comp (group taus EX)
## prior_EX_corr_eta_comp => rho (only 0 for now)

## prior_EX_mu_mean_inter + prior_EX_mu_sd_inter => EX_mu_inter (group means EX)
## prior_EX_tau_mean_inter + prior_EX_tau_sd_inter => EX_tau_comp (group taus EX)
## prior_EX_corr_eta_inter => rho (only 0 for now)
## => EX_eta (group specific means for inter (mvn normal))

## Then sample group parameters, EX
## => log_beta (group specific means for comp (mvn normal))
## => eta (group specific means for comp (mvn normal))

## and NEX parameters are in the same data structure with index starting at num_groups+1

## prior_NEX_mu_mean_comp + prior_NEX_mu_sd_comp => NEX_beta = NEX_comp (group means NEX iid groupwise)
## prior_NEX_mu_mean_inter + prior_NEX_mu_sd_inter => NEX_eta = NEX_inter (group means NEX iid groupwise)

## prior_is_EXNEX_comp + prior_EX_prob_comp => pick which one
## prior_is_EXNEX_inter + prior_EX_prob_inter => pick which one

sample_prior <- function(model) {

    num_strata <- model$num_strata
    num_groups <- model$num_groups
    num_comp <- model$num_comp
    num_inter <- model$num_inter
    blrm_args  <- model$blrm_args
    standata <- model$base_fit$standata
    group_stratum <- standata$group_stratum_cid

    ## group specific parameters: EX, then NEX
    log_beta <- array(NA, dim=c(2*num_groups, num_comp, 2))
    eta <- array(NA, dim=c(2*num_groups, num_inter))

    ## sample EX hyperparameters

    ## mu
    EX_mu_comp <- array(NA, dim=c(num_comp, 2))
    EX_mu_inter <- array(NA, dim=c(num_inter))

    for (j in 1:num_comp) {
        EX_mu_comp[j,1] <- rnorm(1, standata$prior_EX_mu_mean_comp[j,1], standata$prior_EX_mu_sd_comp[j,1])
        EX_mu_comp[j,2] <- rnorm(1, standata$prior_EX_mu_mean_comp[j,2], standata$prior_EX_mu_sd_comp[j,2])
    }

    if (num_inter > 0) {
        for (j in 1:num_inter) {
            EX_mu_inter[j] <- rnorm(1, standata$prior_EX_mu_mean_inter[j], standata$prior_EX_mu_sd_inter[j])
        }
    }

    ## tau
    EX_tau_comp <- array(NA, dim=c(num_strata, num_comp, 2))
    EX_tau_inter <- array(NA, dim=c(num_strata, num_inter))
    ## correlation matrix
    EX_corr_comp <- array(NA, dim=c(num_strata, num_comp, 2, 2))
    EX_corr_inter <- array(NA, dim=c(num_strata, num_inter, num_inter))
    EX_Sigma_comp <- array(NA, dim=c(num_strata, num_comp, 2, 2))
    EX_Sigma_inter <- array(NA, dim=c(num_strata, num_inter, num_inter))

    sample_tau_prior <- function(dist, a, b) {
        if(dist == 0)
            return(a)
        if(dist == 1)
            return(rlnorm(1, a, b))
        if(dist == 2)
            return(abs(rnorm(1, a, b)))
        stop("Unsupported tau prior density.")
    }

    for (s in 1:num_strata) {
        for (j in 1:num_comp) {
            EX_tau_comp[s,j,1] <- sample_tau_prior(standata$prior_tau_dist,
                                                   standata$prior_EX_tau_mean_comp[s,j,1],
                                                   standata$prior_EX_tau_sd_comp[s,j,1] )
            EX_tau_comp[s,j,2] <- sample_tau_prior(standata$prior_tau_dist,
                                                   standata$prior_EX_tau_mean_comp[s,j,2],
                                                   standata$prior_EX_tau_sd_comp[s,j,2] )
            EX_corr_comp[s,j,,] <- rcorvine(2, standata$prior_EX_corr_eta_comp[j], FALSE)
            EX_Sigma_comp[s,j,,] <- diag(as.vector(EX_tau_comp[s,j,]), 2, 2) %*% matrix(EX_corr_comp[s,j,,,drop=FALSE], 2, 2) %*% diag(as.vector(EX_tau_comp[s,j,]), 2, 2)
        }
        if (num_inter > 0) {
            for (j in 1:num_inter) {
                EX_tau_inter[s,j] <- sample_tau_prior(standata$prior_tau_dist,
                                                      standata$prior_EX_tau_mean_inter[s,j],
                                                      standata$prior_EX_tau_sd_inter[s,j] )
            }
            EX_corr_inter[s,,] <- diag(num_inter)
            if (num_inter > 1) {
                EX_corr_inter[s,,] <- rcorvine(num_inter, standata$prior_EX_corr_eta_inter, FALSE)
            }
            EX_Sigma_inter[s,,] <- diag(as.vector(EX_tau_inter[s,,drop=FALSE]), num_inter, num_inter) %*% matrix(EX_corr_inter[s,,,drop=FALSE], num_inter, num_inter) %*% diag(as.vector(EX_tau_inter[s,,drop=FALSE]), num_inter, num_inter)
        }
    }

  ## EX - group-specific parameters
  for (g in 1:num_groups) {
    s <- group_stratum[g]
    for (j in 1:num_comp) {
        log_beta[g,j,1:2] <- rmvnorm(1, EX_mu_comp[j,], EX_Sigma_comp[s,j,,])
        ##log_beta[g,j,1] <- rnorm(1, EX_mu_comp[j,1], EX_tau_comp[s,j,1] )
        ##log_beta[g,j,2] <- rnorm(1, EX_mu_comp[j,2], EX_tau_comp[s,j,2] )
        ##assert_that(standata$prior_EX_corr_eta_comp[j] == 1, msg="LKJ correlation == 1 is only supported.")
    }
    if (num_inter > 0) {
        if(num_inter > 1) {
            eta[g,] <- rmvnorm(1, EX_mu_inter, EX_Sigma_inter[s,,])
        } else {
            eta[g,1] <- rnorm(1, EX_mu_inter, EX_tau_inter[s,1])
        }
      ##for (j in 1:num_inter) {
      ##  eta[g,j] <- rnorm(1, EX_mu_inter[j], EX_tau_inter[s,j] )
      ##}
      ##assert_that(standata$prior_EX_corr_eta_inter == 1, msg="LKJ correlation == 1 is only supported.")
    }
  }
  ## NEX - group-specific parameters
  for (g in 1:num_groups) {
    for (j in 1:num_comp) {
      log_beta[num_groups+g,j,1] <- rnorm(1, standata$prior_NEX_mu_mean_comp[j,1], standata$prior_NEX_mu_sd_comp[j,1])
      log_beta[num_groups+g,j,2] <- rnorm(1, standata$prior_NEX_mu_mean_comp[j,2], standata$prior_NEX_mu_sd_comp[j,2])
    }
    if (num_inter > 0) {
      for (j in 1:num_inter) {
        eta[num_groups+g,j] <- rnorm(1, standata$prior_NEX_mu_mean_inter[j], standata$prior_NEX_mu_sd_inter[j])
      }
    }
  }

  ## convert slope to natural scale (enforced positivity)
  beta <- log_beta
  for (g in 1:(2*num_groups)) {
    for (j in 1:num_comp) {
      beta[g,j,2] <- exp(beta[g,j,2])
    }
  }

  ## sample EX / NEX membership
  is_EX_comp <- array(NA, dim=c(num_groups,num_comp))
  is_EX_inter <- array(NA, dim=c(num_groups,num_inter))
  draw_beta  <- array(NA, dim=c(1, num_groups, num_comp, 2))
  draw_eta  <- array(NA, dim=c(1, num_groups, num_inter))
  for(g in 1:num_groups) {
    for (j in 1:num_comp) {
      if(standata$prior_is_EXNEX_comp[j] == 1) {
        is_EX_comp[g,j] <- rbinom(1, 1, standata$prior_EX_prob_comp[g,j])
      } else {
        is_EX_comp[g,j] <- 1
      }
      gidx <- ifelse(is_EX_comp[g,j] == 1, g, num_groups + g)
      draw_beta[1,g,j,1]  <- beta[gidx,j,1]
      draw_beta[1,g,j,2]  <- beta[gidx,j,2]
    }
    if (num_inter > 0) {
      for (j in 1:num_inter) {
        if(standata$prior_is_EXNEX_inter[j] == 1) {
          is_EX_inter[g,j] <- rbinom(1, 1, standata$prior_EX_prob_inter[g,j])
        } else {
          is_EX_inter[g,j] <- 1
        }
        gidx <- ifelse(is_EX_inter[g,j] == 1, g, num_groups + g)
        draw_eta[1,g,j] <- eta[gidx,j]
      }
    }
  }
  list(draw_beta=draw_beta,
       draw_eta=draw_eta,
       EX_mu_comp=EX_mu_comp,
       EX_mu_inter=EX_mu_inter,
       EX_tau_comp=EX_tau_comp,
       EX_tau_inter=EX_tau_inter,
       EX_corr_comp=EX_corr_comp,
       EX_corr_inter=EX_corr_inter,
       log_beta=log_beta,
       eta=eta,
       is_EX_comp=is_EX_comp,
       is_EX_inter=is_EX_inter
  )
}

#'
#' Simulates a draw from the prior and fake data for it. This will be
#' the data generating step in the simulation. The function recieves
#' the problem data, job specifics and a blrmfit object which defines
#' the prior to sample and the design matrix.
#'
simulate_fake <- function(data, job, model, ...) {

  model  <- data$models[[model]]

  prior_draw  <- sample_prior(model)

  standata  <- model$base_fit$standata

  ## logit by data-row
  draw_mu <- with(standata, blrm_logit_grouped_vec(group, stratum, X_comp, X_inter, prior_draw$draw_beta, prior_draw$draw_eta))

  num_trials <- standata$r + standata$nr

  yrep <- rbinom(length(num_trials), num_trials, inv_logit(draw_mu))

  list(yrep = yrep, draw = prior_draw)

}

#'
#'
#' Procedure to fit each fake data set using our fitting
#' procedure. This method obtains the problem data, job details and an
#' **instance** of the scenario as generated by `simulate_fake`.
#'

fit_exnex <- function(data, job, instance, ...) {

  yrep <- instance$yrep
  draw <- instance$draw

  pars <- job$pars$prob.pars
  model  <- data$models[[pars$model]]

  dref <- model$dref
  # attach(model)
  sim_data <- model$base_fit$data
  sim_data$Ntox <- yrep

  fit <- update(model$base_fit,
                data = sim_data,
                iter = model$blrm_args$iter,
                warmup = model$blrm_args$warmup)

  # tmp <- summary(fit)
  # summ_fit <- cbind(sim_data, tmp)

  sampler_params <- rstan::get_sampler_params(fit$stanfit, inc_warmup=FALSE)
  n_divergent <- sum(sapply(sampler_params, function(x) sum(x[,'divergent__'])) )

  params <- c("mu_log_beta", "tau_log_beta", "mu_eta", "tau_eta")
  fit_sum <- rstan::summary(fit$stanfit)$summary
  min_Neff <- ceiling(min(fit_sum[apply(sapply(params, grepl, x = rownames(fit_sum)), 1, any),
                                  "n_eff"], na.rm=TRUE))


  post <- rstan::extract(fit$stanfit, pars = params, inc_warmup = FALSE)
  post_thin <- lapply(post, function(A){
    asub(A, idx = seq(1, dim(A)[1], length = 1024-1), dims = 1, drop = FALSE)
  })
  draw$EX_mu_inter <- matrix(draw$EX_mu_inter, nrow = 1)
  draw$EX_tau_inter <- matrix(draw$EX_tau_inter, nrow = 1)

  rank0 <- list(
    mu_log_beta = apply(
      sweep(post_thin$mu_log_beta, c(2, 3), draw$EX_mu_comp) < 0,
      c(2, 3),
      sum
    ),
    tau_log_beta = apply(
      sweep(post_thin$tau_log_beta, c(2, 3, 4), draw$EX_tau_comp) < 0,
      c(2, 3, 4),
      sum
    )
  )

  if(fit$has_inter){
    rank1 <- list(
      mu_log_beta = rank0$mu_log_beta,
      tau_log_beta = rank0$tau_log_beta,
      mu_eta = apply(
        sweep(post_thin$mu_eta, 2, draw$EX_mu_inter) < 0,
        2,
        sum
      ),
      tau_eta = apply(
        sweep(post_thin$tau_eta, c(2, 3), draw$EX_tau_inter) < 0,
        c(2, 3),
        sum
      )
    )
  } else{
    rank1 <- rank0
  }


  # relabeler <- setNames(
  #   as.matrix(fit$stanfit, pars = params) %>% colnames(),
  #   rstan::As.mcmc.list(fit$stanfit, pars = params) %>% '[['(1) %>% colnames()
  # )

  rank <- do.call(
    rbind,
    lapply(
      names(rank1),
      function(param){
        r0 <- rank1[[param]]
        if(is.null(dim(r0))){
          D <- 1
        } else{
          D <- length(dim(r0))
        }
        df0 <- as.data.frame(as.table(r0))
        for(v in paste0("Var", 1:D)){
          levels(df0[[v]]) <- 1:nlevels(df0[[v]])
        }
        idx <- apply(df0[paste0("Var", 1:D)], 1, paste, collapse = ",")
        df <- data.frame(
          param = paste0(param, "[", idx, "]"),
          rank = df0$Freq
        )
      }
    )
  )

  labels <- fit$labels
  relabeler <- setNames(as.character(rank$param), rank$param)
  relabeler <- .label_index(relabeler, "mu_log_beta", labels$component, labels$param_log_beta)
  relabeler <- .label_index(relabeler, "tau_log_beta", fit$strata_fct, labels$component, labels$param_log_beta)
  if(fit$has_inter){
    relabeler <- .label_index(relabeler, "mu_eta", labels$param_eta)
    relabeler <- .label_index(relabeler, "tau_eta", fit$strata_fct, labels$param_eta)
  }
  labeler <- setNames(names(relabeler), relabeler)
  rank$param <- labeler[rank$param]

  # rank_wide <- setNames(data.frame(t(rank$rank)), paste0("rank-", rank$param))
  rank_wide <- setNames(data.frame(t(rank$rank)), rank$param)

  list(rank = rank_wide,
       min_Neff = min_Neff,
       n_divergent = n_divergent)


}

# AB: not currently using this function from rbest SBC...
# scale_ranks <- function(Nbins, scale=1) {
#   ## scale must evenly divide the total number of bins
#   assert_that(round(Nbins/scale) == Nbins/scale)
#   breaks <- (0:(Nbins/scale))
#   Nbreaks <- length(breaks)
#   function(scen) {
#     vars <- grep("^rank.", names(scen), value=TRUE)
#     res <- lapply(vars, function(v) hist(ceiling((scen[[v]]+1)/scale), breaks=breaks, plot=FALSE, include.lowest=FALSE)$counts)
#     names(res) <- gsub("^rank", "count", vars)
#     res$rank <- breaks[-Nbreaks]
#     res <- as.data.frame(do.call(cbind, res))
#     res
#   }
# }


## Submits to batchtools cluster with fault tolerance, i.e.
## resubmitting failed jobs max_num_tries times
auto_submit <- function(jobs, registry, resources=list(), max_num_tries = 10) {
  all_unfinished_jobs <- jobs
  
  num_unfinished_jobs <- nrow(all_unfinished_jobs)
  num_all_jobs <- num_unfinished_jobs
  remaining_tries <- max_num_tries
  all_jobs_finished <- FALSE
  while (remaining_tries > 0 && !all_jobs_finished) {
    remaining_tries <- remaining_tries - 1
    
    print(paste("Submitting jobs at ", Sys.time()))
    # Once things run fine let's submit this work to the cluster.
    submitJobs(all_unfinished_jobs, resources=resources)
    # Wait for results.
    waitForJobs()
    print(paste("Finished waiting for jobs at ", Sys.time()))
    
    # Check status:
    print(getStatus())
    
    # Ensure that all jobs are done
    if (nrow(findNotDone()) != 0) {
      not_done_jobs <- findNotDone()
      print(getErrorMessages(not_done_jobs))
      ##browser()
      ##invisible(readline(prompt="Press [enter] to continue"))
      
      print(paste0("Some jobs did not complete. Please check the batchtools registry ", registry$file.dir))
      all_unfinished_jobs <- inner_join(not_done_jobs, all_unfinished_jobs)
      
      if (num_unfinished_jobs == nrow(all_unfinished_jobs) &&  nrow(all_unfinished_jobs) > 0.25 * num_all_jobs)
      {
        # Unfinished job count did not change -> retrying will probably not help. Abort!
        cat("Error: unfinished job count is not decreasing. Aborting job retries.")
        remaining_tries <- 0
      }
      
      if (num_unfinished_jobs == nrow(jobs))
      {
        # All jobs errored -> retrying will probably not help. Abort!
        cat("Error: all jobs errored. Aborting job retries.")
        remaining_tries <- 0
      }
      
      num_unfinished_jobs <- nrow(all_unfinished_jobs)
      print(paste0("Trying to resubmit jobs. Remaining tries: ", remaining_tries, " / ", max_num_tries))
    } else {
      all_jobs_finished <- TRUE
    }
  }
  
  all_jobs_finished
}
