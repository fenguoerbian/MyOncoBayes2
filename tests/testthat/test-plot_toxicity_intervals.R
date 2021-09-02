context("plot intervals tests")

# use gold runs for these tests avoiding sampling uncertainty
single_agent <- gold_runs$single_agent
combo2 <- gold_runs$combo2
trial_examples <- gold_runs$trial_examples

# plot_toxicity_intervals.blrmfit ----------------------------------------------

test_that(
  "plot_toxicity_intervals.blrmfit works for single-agent example",
  {
    expect_gg(plot_toxicity_intervals(single_agent$blrmfit,
                                      x = vars(drug_A)))

    a <- plot_toxicity_intervals(single_agent$blrmfit,
                                 x = vars(drug_A),
                                 group = vars(group_id))
    expect_gg(a)

    expect_gg(plot_toxicity_intervals(single_agent$blrmfit,
                                      x = vars(drug_A),
                                      interval_prob = c(0, 0.1, 0.2, 0.5, 1),
                                      interval_max_mass = c(NA, NA, NA, 0.1)))


    # check nonsense variable
    expect_error(plot_toxicity_intervals(single_agent$blrmfit,
                                         x = vars(drug_B)))


    # check character
    expect_gg(plot_toxicity_intervals(single_agent$blrmfit,
                                      x = "drug_A"))

    # check nonsense character
    expect_error(plot_toxicity_intervals(single_agent$blrmfit,
                                         x = "drug_B"))

    # one var(), one formula
    b <- plot_toxicity_intervals(single_agent$blrmfit,
                                 x = vars(drug_A),
                                 group = ~ group_id)
    expect_gg(b)


    # one var(), one nonsense formula
    expect_error(plot_toxicity_intervals(single_agent$blrmfit,
                                         x = vars(drug_A),
                                         group = ~ drug_B))

    # duplicated arguments
    expect_error(plot_toxicity_intervals(single_agent$blrmfit,
                                         x = vars(drug_A),
                                         group = vars(drug_A)))

    # duplicated arguments
    expect_error(plot_toxicity_intervals(single_agent$blrmfit,
                                         x = vars(drug_A),
                                         group = ~drug_A))

  })


test_that(
  "plot_toxicity_intervals.blrmfit works for combo2 example",
  {
    a <- plot_toxicity_intervals(combo2$blrmfit,
                                 x = vars(drug_A),
                                 group = vars(group_id, drug_B))


    b <- plot_toxicity_intervals(combo2$blrmfit,
                                 x = "drug_A",
                                 group = c("group_id", "drug_B"))

    expect_gg(a)
    expect_gg(b)

  })


test_that(
  "plot_toxicity_intervals.blrm_trial works for blrm_trial examples",
  {
    for(trial in trial_examples){
      a <- plot_toxicity_intervals(trial)
      expect_gg(a)
    }

  })



# Visual regression testing ----------------------------------------------------

test_that("plot_toxicity_curve.blrmfit renders single-agent plots correctly", {

  testthat::skip_on_cran()
  testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
  testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

  # basic
  p <- plot_toxicity_intervals(single_agent$blrmfit, x = vars(drug_A))
  vdiffr::expect_doppelganger("single-agent blrmfit", p)

  # change limits
  p <- plot_toxicity_intervals(single_agent$blrmfit, x = vars(drug_A),
                               interval_prob = c(0, 0.25, 0.5, 1),
                               interval_max_mass = c(0.5, 0.1, 0.05))
  vdiffr::expect_doppelganger("single-agent blrmfit limits", p)

  # change width
  p <- plot_toxicity_intervals(single_agent$blrmfit, x = vars(drug_A),
                               ewoc_colors = c("blue", "purple"))
  vdiffr::expect_doppelganger("single-agent blrmfit colors", p)

})


test_that("plot_toxicity_curve.blrmfit renders combo2 plots correctly", {

  testthat::skip_on_cran()
  testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
  testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

  blrmfit <- combo2$blrmfit

  nd <- expand_grid(
    group_id = factor("trial_AB", levels(blrmfit$group_fct)),
    drug_A = c(3, 4.5, 6),
    drug_B = c(400, 600, 800)
  )

  # basic
  p <- plot_toxicity_intervals(blrmfit, x = vars(drug_A), group = ~ drug_B + group_id,
                           newdata = nd)
  vdiffr::expect_doppelganger("combo2 blrmfit", p)


})


for(ex in names(trial_examples)){

  test_that(paste("plot_toxicity_intervals.blrm_trial renders", ex, "plots correctly"), {

    testthat::skip_on_cran()
    testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
    testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

    example <- trial_examples[[ex]]

    # blrm_trial
    p <- plot_toxicity_intervals(example)
    vdiffr::expect_doppelganger(paste(ex, "blrm_trial"), p)

  })

}

