#'
#' Utilities for SBC validation
#'

example_designs <- list(

  combo3_EXNEX = list(
    modeltype = "EXNEX",
    design = expand.grid(
      stratum = "STRAT",
      group_id = LETTERS[1:3],
      DosesAdm1 = c(0.25, 0.5, 1),
      DosesAdm2 = c(0.25, 0.5, 1),
      DosesAdm3 = c(0.25, 0.5, 1),
      Npat = 3
    ) %>% arrange(stratum, group_id, DosesAdm1, DosesAdm2, DosesAdm3),
    dref = c(1, 1, 1)
  ),

  combo2_EX = list(
    modeltype = "EX",
    design = expand.grid(
      stratum = "STRAT",
      group_id = LETTERS[1:3],
      DosesAdm1 = c(0.25, 0.5, 1),
      DosesAdm2 = c(0.25, 0.5, 1),
      Npat = 5
    ) %>% arrange(stratum, group_id, DosesAdm1, DosesAdm2),
    dref = c(1, 1)
  ),

  combo2_EXNEX = list(
    modeltype = "EXNEX",
    design = expand.grid(
      stratum = "STRAT",
      group_id = LETTERS[1:3],
      DosesAdm1 = c(0.25, 0.5, 1),
      DosesAdm2 = c(0.25, 0.5, 1),
      Npat = 5
    ) %>% arrange(stratum, group_id, DosesAdm1, DosesAdm2),
    dref = c(1, 1)
  ),

  log2bayes_EXNEX = list(
    modeltype = "EXNEX",
    design = expand.grid(
      stratum = "STRAT",
      group_id = LETTERS[1:3],
      DosesAdm1 = c(0.0625, 0.125, 0.25, 0.5, 1),
      Npat = 5
    ) %>% arrange(stratum, group_id, DosesAdm1),
    dref = c(1)
  )

)

example_models <- lapply(
  example_designs,
  function(example){

    design <- example$design
    dref <- example$dref
    modeltype <- example$modeltype
    is_exnex <- modeltype != "EX"
    p_exch <- switch(
      modeltype,
      "EX" = 1,
      "NEX" = 1e-06,
      "EXNEX" = 0.8
    )

    num_comp <- sum(grepl(names(design), pattern = "DosesAdm"))
    num_groups <- nlevels(design$group_id)
    num_strata <- nlevels(design$stratum)
    num_inter <- ifelse(num_comp > 1, sum(choose(num_comp, 2:num_comp)), 0)

    formula <- switch(
      as.character(num_comp),
      "1" = cbind(Ntox, Npat - Ntox) ~ 1 + I(log(DosesAdm1/dref[1])) | 0 | group_id,
      "2" = cbind(Ntox, Npat-Ntox) ~ 1 + I(log(DosesAdm1/dref[1])) | 1 + I(log(DosesAdm2/dref[2])) | 0 + I(DosesAdm1/dref[1] * DosesAdm2/dref[2]) | group_id,
      "3" = cbind(Ntox, Npat-Ntox) ~ 1 + I(log(DosesAdm1/dref[1])) | 1 + I(log(DosesAdm2/dref[2])) | 1 + I(log(DosesAdm3/dref[3])) | 0 + I(DosesAdm1/dref[1] * DosesAdm2/dref[2]) + I(DosesAdm1/dref[1] * DosesAdm3/dref[3]) + I(DosesAdm2/dref[2] * DosesAdm3/dref[3]) + I(DosesAdm1/dref[1] * DosesAdm2/dref[2] * DosesAdm3/dref[3]) | group_id
    )

    blrm_args <- list(

      formula = formula,

      data = design,

      prior_EX_mu_mean_comp = matrix(c(logit(1/3), 0), nrow=num_comp, ncol=2, TRUE),
      prior_EX_mu_sd_comp = matrix(c(2, 1), nrow=num_comp, ncol=2, TRUE),
      prior_EX_tau_mean_comp = matrix(log(  c(0.25, 0.125)), nrow=num_comp, ncol=2, TRUE),
      prior_EX_tau_sd_comp = matrix(log(4)/1.96, nrow=num_comp, ncol=2, TRUE),
      prior_EX_mu_mean_inter = rep(0, num_inter),
      prior_EX_mu_sd_inter = rep(log(2)/1.96, num_inter),
      prior_EX_tau_mean_inter = matrix(log(0.25)  , nrow=num_strata, ncol=num_inter),
      prior_EX_tau_sd_inter = matrix(log(2)/1.96, nrow=num_strata, ncol=num_inter),
      prior_EX_prob_comp = matrix(p_exch, nrow=num_groups, ncol=num_comp),
      prior_EX_prob_inter = matrix(1.0, nrow=num_groups, ncol=num_inter),
      prior_is_EXNEX_comp = rep(is_exnex, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_tau_dist=1,
      prior_NEX_mu_mean_comp = matrix(c(logit(1/3), 0), nrow=num_comp, ncol=2, TRUE),
      prior_NEX_mu_sd_comp = matrix(c(2, 1), nrow=num_comp, ncol=2, TRUE),

      prior_NEX_mu_mean_inter = rep(0, num_inter),
      prior_NEX_mu_sd_inter = rep(log(4)/1.96, num_inter),

      iter = 1500,
      warmup = 500,
      thin = 1,
      init = 0.5,
      chains = 2,
      cores = 1,
      control = list(),
      prior_PD = FALSE

    )

    base_args <- blrm_args
    base_args$data <- base_args$data %>% mutate(Ntox = 0)
    base_args$iter <- 2
    base_args$warmup <- 1
    base_fit <- do.call(blrm_exnex, base_args)

    return(
      list(dref = dref,
           num_strata = num_strata,
           num_groups = num_groups,
           num_comp = num_comp,
           num_inter = num_inter,
           blrm_args = blrm_args,
           base_fit = base_fit)
    )

  }
)
