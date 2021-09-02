# load gold runs from cache if available

load_gold <- function(use_cache=FALSE) {
    # load gold runs
    gold_path <- test_path("_snaps/gold_runs.rds")
    if(use_cache && file.exists(gold_path)) {
        gold_runs <- readRDS(gold_path)
    } else {
        options(mc.cores=parallel::detectCores(logical=FALSE))

        set.seed(123144)

        single_agent  <- run_example("single_agent")
        combo2  <- run_example("combo2")
        # combo3  <- run_example("combo3")

        trial_examples <- lapply(examples[1:3], function(example){
            with(example, blrm_trial(histdata, dose_info, drug_info %>% mutate(reference_p_dlt = 0.1), simplified_prior = TRUE))
        })

        gold_runs <- list(
            single_agent = single_agent,
            combo2 = combo2,
            trial_examples = trial_examples,
            run_date=date()
        )
    }
    return(gold_runs)
}

