
context("blrm_exnex tests")


set.seed(123144)

eps <- 1E-4

## reference runs
single_agent  <- run_example("single_agent")
combo2  <- run_example("combo2")
combo3  <- run_example("combo3")

suppressPackageStartupMessages(library(dplyr))

## perform some basic checks on the shape of outputs
check_model_basic  <- function(fit, envir) {
    skip_on_cran()
    data <- fit$data
    num_obs  <- nrow(data)
    ss  <- summary(fit, interval_prob=c(0,0.16,0.33,1))

    expect_equal(nrow(ss), num_obs)
    expect_equal(ncol(ss), 8)
    expect_equal(sum(is.na(ss)), 0)
    ss2  <- summary(fit, interval_prob=c(0,1))
    expect_equal(sum((ss2[,"(0,1]"] - 1)<eps), nrow(data))

    ## check that the median and 50% interval prob align
    s_med  <- ss[,"50%"]
    for(i in 1:length(s_med)) {
        i <- 1
        q_med  <- s_med[i]
        s_q  <- summary(fit, interval_prob=c(0,q_med,1))
        expect_true(abs(s_q[i,6] - 0.5) < eps)
        expect_true(abs(s_q[i,7] - 0.5) < eps)
    }

    iter  <- getOption("OncoBayes2.MC.iter")
    warmup  <- getOption("OncoBayes2.MC.warmup")
    chains  <- getOption("OncoBayes2.MC.chains")
    num_sim <- chains * (iter - warmup)

    plin  <- posterior_linpred(fit)
    expect_equal(nrow(plin), num_sim)
    expect_equal(ncol(plin), nrow(data))
    expect_equal(sum(is.na(ss)), 0)

    envir$fit  <- fit
    envir$data1  <- data[1,,drop=FALSE]

    suppressMessages(capture.output(single_obs_fit <- with(envir, update(fit, data=data1))))
    ss1  <- summary(single_obs_fit)
    expect_equal(nrow(ss1), 1)
    expect_equal(ncol(ss1), 5)
    expect_equal(sum(is.na(ss1)), 0)

    plin1  <- posterior_linpred(single_obs_fit)
    expect_equal(nrow(plin1), num_sim)
    expect_equal(ncol(plin1), 1)
}


test_that("blrm_exnex data handling consistency single-agent",
          check_model_basic(single_agent$blrmfit, single_agent))

test_that("blrm_exnex data handling consistency combo2",
          check_model_basic(combo2$blrmfit, combo2))

test_that("blrm_exnex data handling consistency combo3",
          check_model_basic(combo3$blrmfit, combo3))

test_that("interval probabilites are consistent", {
    skip_on_cran()
    combo2_sens  <- combo2
    sens_low <- hist_combo2 %>%
        mutate(Ntox=0)
    combo2_sens$sens_low  <- sens_low
    outp <- suppressMessages(capture.output(sens_fit_low <- with(combo2_sens, update(blrmfit, data=sens_low))))
    s_low  <- summary(sens_fit_low, interval_prob=c(0,0.9,1))
    expect_equal(sum((s_low[,"(0,0.9]"] - 1)<eps), nrow(sens_low))
    expect_equal(sum((s_low[,"(0.9,1]"])<eps), nrow(sens_low))

    sens_high <- hist_combo2 %>%
        mutate(Ntox=Npat)
    combo2_sens$sens_high <- sens_high
    outp <- suppressMessages(capture.output(sens_fit_high <- with(combo2_sens, update(blrmfit, data=sens_high))))
    s_high  <- summary(sens_fit_high, interval_prob=c(0,0.1,1))
    expect_equal(sum((s_high[,"(0.1,1]"] - 1)<eps), nrow(sens_low))
    expect_equal(sum((s_high[,"(0,0.1]"])<eps), nrow(sens_low))
})


## TODO: Add proper runs which are compared against published results,
## e.g. from the book chapter

test_that("interval probabilites are not NaN", {
    ss1  <- summary(combo2$blrmfit, interval_prob=c(0.1, 0.4))
    expect_true(all(!is.na(ss1[,6])))
    ss2  <- summary(combo2$blrmfit, transform=FALSE, interval_prob=c(0.1, 0.4))
    expect_true(all(!is.na(ss2[,6])))
})

query_fit <- function(fit) {
    capture_output(print(fit))
    prior_summary(fit)
    summary(fit)
    summary(fit, interval_prob=c(0,0.16,0.33,1.0))
    pl <- posterior_linpred(fit)
    pi <- posterior_interval(fit)
    pp <- predictive_interval(fit)
}

test_that("blrmfit methods work for single agent models", {
    ## this test is successfull if methods run without errors as we do
    ## not have running examples given they are too costly.
    query_fit(single_agent$blrmfit)
})

test_that("blrmfit methods work for combo2 models", {
    ## this test is successfull if methods run without errors as we do
    ## not have running examples given they are too costly.
    query_fit(combo2$blrmfit)
})

test_that("blrmfit methods work for combo3 models", {
    ## this test is successfull if methods run without errors as we do
    ## not have running examples given they are too costly.
    query_fit(combo3$blrmfit)
})
