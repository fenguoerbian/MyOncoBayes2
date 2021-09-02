## library(OncoBayes2)
devtools::load_all()

git_commit  <- system("git rev-parse HEAD", intern=TRUE)

set.seed(42)

source("tests/testthat/helper-sampling.R")
source("tests/testthat/helper-trial_examples.R")
source("tests/testthat/helper-load_gold.R")

gold_runs <- load_gold(FALSE)
gold_runs$commit <- git_commit

saveRDS(gold_runs, test_path("_snaps/gold_runs.rds"), compress = "xz")

cat("Consider updating reference plots using tools/update_ref_plots.R\n")
cat("... or rerun CI/CD with variable RESET_REF_PLOTS=true and then replace plots from artifacts.\n")
