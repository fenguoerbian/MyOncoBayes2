
if(!exists("gold_runs")) {
    if(!testthat:::on_cran() & identical(Sys.getenv("TEST_VDIFFR"), "true")){

        # load gold runs, possibly from cache
        .mc_options <- default_sampling()
        gold_runs <- load_gold(TRUE)
        options(.mc_options)

    } else {

        # load gold runs, use cache if available, but do short
        # sampling if not
        if (identical(Sys.getenv("NOT_CRAN"), "true")) {
            .mc_options <- very_fast_sampling()
        } else {
            .mc_options <- fake_sampling()
        }
        gold_runs <- load_gold(TRUE)
        options(.mc_options)

    }
}

