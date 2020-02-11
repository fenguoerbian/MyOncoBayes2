#'
#' Utilities for SBC validation
#'

example_designs <- list(

  combo3_EXNEX = list(
    modeltype = "EXNEX",
    design = expand.grid(
      stratum_id = "STRAT",
      group_id = LETTERS[1:3],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      drug_C = c(0.25, 0.5, 1),
      num_patients = 3
    ) %>% arrange(stratum_id, group_id, drug_A, drug_B, drug_C),
    dref = c(1, 1, 1)
  ),

  combo2_EX = list(
    modeltype = "EX",
    design = expand.grid(
      stratum_id = "STRAT",
      group_id = LETTERS[1:3],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      num_patients = 5
    ) %>% arrange(stratum_id, group_id, drug_A, drug_B),
    dref = c(1, 1)
  ),

  combo2_EXNEX = list(
    modeltype = "EXNEX",
    design = expand.grid(
      stratum_id = "STRAT",
      group_id = LETTERS[1:3],
      drug_A = c(0.25, 0.5, 1),
      drug_B = c(0.25, 0.5, 1),
      num_patients = 5
    ) %>% arrange(stratum_id, group_id, drug_A, drug_B),
    dref = c(1, 1)
  ),

  log2bayes_EXNEX = list(
    modeltype = "EXNEX",
    design = expand.grid(
      stratum_id = "STRAT",
      group_id = LETTERS[1:3],
      drug_A = c(0.0625, 0.125, 0.25, 0.5, 1),
      num_patients = 5
    ) %>% arrange(stratum_id, group_id, drug_A),
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

      names(dref) <- grep(names(design), pattern = "drug_", value=TRUE)

      linear_formula  <- blrm_formula_linear(dref, 3)

      num_comp <- linear_formula$num_components
      num_inter <- linear_formula$num_interaction_terms
      formula <- as.formula(linear_formula$blrm_formula)

      num_groups <- nlevels(design$group_id)
      num_strata <- nlevels(design$stratum)

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
          ##iter = 150,
          ##warmup = 50,
          thin = 1,
          init = 0.5,
          chains = 2,
          ##cores = 1, ## control via mc.cores option
          control = list(),
          prior_PD = FALSE
      )

      base_args <- blrm_args
      base_args$data <- base_args$data %>% mutate(num_toxicities = 0)
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
