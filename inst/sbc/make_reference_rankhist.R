library(pkgbuild)
## ensure that the current dev version of OncoBayes2 is loaded
pkgbuild::compile_dll("../..")

library(batchtools)
set.seed(453453)
## load utilities and current dev version of OncoBayes2
source("sbc_tools.R")

#' according to the docs this speeds up the reduce step
options(batchtools.progress = FALSE)

reg <- makeExperimentRegistry(file.dir = tempfile("sbc_"),
                              ## use the default configured batchtools configuration batchtools.conf
                              ## - found in the environemnt variable R_BATCHTOOLS_SEARCH_PATH
                              ## - current working directory
                              ## - $(HOME)/.batchtools.conf
                              # conf.file = NA,
                              seed = 47845854,
                              ## our worker functions and package loading
                              source="sbc_tools.R")

## resources of each job: Less than 55min, 2000MB RAM and 2 cores
job_resources <- list(walltime=55, memory=2000, ncpus=2)

if(FALSE) {
  ## for debugging here
  removeProblems("combo2_EX")
  removeProblems("combo2_NEX")
  removeProblems("combo2_EXNEX")
  removeProblems("base")
}

#' Evaluate dense and sparse data-scenario
source("sbc_example_models.R")
base_data  <- list(models = example_models)


addProblem("base",
           data = base_data,
           fun = simulate_fake,
           seed = 2345,
           ## caching speeds up the reduce step
           cache = TRUE
)


addAlgorithm("OncoBayes2", fit_exnex)

## family, mean_mu, sd_mu, sd_tau, samp_sd
scenarios <- data.frame(model = names(example_models),
                        stringsAsFactors=FALSE)

pdes <- list(base = scenarios)
ades <- list(OncoBayes2 = data.frame())

#' Add the defined problems and analysis methods to the registry and
#' set the number of replications:
S  <- 1E4
addExperiments(pdes, ades, repls = S)

summarizeExperiments()

if(FALSE) {
  ## used for debugging
  job1 <- testJob(7)
  job1
  job2 <- testJob(6)
  job3 <- testJob(11)
    job <- makeJob(1)
    attributes(job)
    names(job)
    example_model("combo2")
    sp  <- sample_prior(example_models$combo2_EXNEX)
    sp  <- sample_prior(example_models$log2bayes_EXNEX)
    out <- fit_exnex(data = base_data, job = job, instance = job$instance )
}


#'
#' Chunk the jobs into 4000 chunks to run
#'
ids <- getJobTable()
ids[, chunk:=chunk(job.id, 4000)]

num_jobs <- nrow(ids)

#' Once things run fine let's submit this work to the cluster.
##submitJobs(ids, job_resources)

#' This function deals with unstable nodes in the cluster
auto_submit(ids, reg, job_resources)

#' Check status:
print(getStatus())

#' Ensure that no error occured
assert_that(nrow(findErrors()) == 0)

#' Collect results.
calibration_data <- ijoin(
  ## grab job parameters
  unwrap(getJobPars()),
  unwrap(reduceResultsDataTable(fun=function(x) c(rank=as.list(x$rank),
                                                  list(n_divergent=x$n_divergent,
                                                       min_Neff=ceiling(x$min_Neff))) ))
)

# check that indeed all jobs have finished
assert_that(nrow(calibration_data) == num_jobs)

# there is only one algorithm. remove that column.
calibration_data <- calibration_data %>% select(-algorithm, -problem)

#' Bin raw data as used in the analysis.
B <- 1024L / 2^5

rank_params <- names(calibration_data)[grepl(names(calibration_data), pattern = "rank")]

# calibration_data_binned <- calibration_data[, scale64(.SD), by=c("problem", "model", params)]

calibration_data_binned <- calibration_data %>%
  mutate_at(.vars = rank_params, .funs = function(x) ceiling((x + 1) / (1024 / B) - 1))

names(calibration_data_binned) <- gsub(
  names(calibration_data_binned),
  pattern = "rank[.]",
  replacement = ""
)

params <- gsub(rank_params, pattern = "rank[.]", replacement = "")

calibration_binned <- calibration_data_binned %>%
  dplyr::select(-job.id, - n_divergent, - min_Neff) %>%
  tidyr::gather(key = "param", value = "bin", - model) %>%
  group_by(model, param, bin) %>%
  tally() %>%
  right_join(
    expand.grid(
      model = unique(calibration_data_binned$model),
      param = params,
      bin = 0:(B - 1),
      stringsAsFactors = FALSE
    ),
    c("model", "param", "bin")
  ) %>%
  replace_na(list(n = 0)) %>%
  arrange(model, param, bin) %>%
  spread(key = param, value = n)



#' Save as data.frame to avoid data.table dependency.
calibration_data <- as.data.frame(calibration_data)
calibration_binned <- as.data.frame(calibration_binned)

#' Further identification and verification data of run
git_hash <- system2("git", c("rev-parse", "HEAD"), stdout=TRUE)
created <- Sys.time()
created_str <- format(created, "%F %T %Z", tz="UTC")

calibration <- list(raw = calibration_data,
                    data = calibration_binned,
                    S = S,
                    B = B,
                    git_hash = git_hash,
                    created = created)

saveRDS(calibration, file = "calibration.rds")

library(tools)
md5 <- md5sum("calibration.rds")
cat(paste0("Created:  ", created_str, "\ngit hash: ", git_hash, "\nMD5:      ", md5, "\n"),
    file="calibration.md5")

#' Cleanup
removeRegistry(0)

#' Session info
sessionInfo()
