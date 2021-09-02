context("plot stacked tests")

# use gold runs for these tests avoiding sampling uncertainty
single_agent <- gold_runs$single_agent
combo2 <- gold_runs$combo2
trial_examples <- gold_runs$trial_examples

# see helper-load_gold.R and setup-gold.R for example runs: single_agent, combo2, and trial_examples


# plot_toxicity_intervals_stacked.blrmfit ----------------------------------------------
test_that(
  "plot_toxicity_intervals_stacked.blrmfit works for single-agent example",
  {
    expect_gg(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                              x = vars(drug_A)))

    a <- plot_toxicity_intervals(single_agent$blrmfit,
                                 x = vars(drug_A),
                                 group = vars(group_id))
    expect_gg(a)

    expect_gg(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                              x = vars(drug_A),
                                              interval_prob = c(0, 0.1, 0.2, 0.5, 1)))


    # check nonsense variable
    expect_error(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                                 x = vars(drug_B)))


    # check character
    expect_gg(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                              x = "drug_A"))

    # check nonsense character
    expect_error(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                                 x = "drug_B"))

    # one var(), one formula
    b <- plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                         x = vars(drug_A),
                                         group = ~ group_id)
    expect_gg(b)


    # one var(), one nonsense formula
    expect_error(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                                 x = vars(drug_A),
                                                 group = ~ drug_B))

    # duplicated arguments
    expect_error(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                                 x = vars(drug_A),
                                                 group = vars(drug_A)))

    # duplicated arguments
    expect_error(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                                 x = vars(drug_A),
                                                 group = ~drug_A))

    expect_warning(plot_toxicity_intervals_stacked(single_agent$blrmfit,
                                                   x = vars(drug_A),
                                                   transform = FALSE,
                                                   predictive = FALSE))

  })


test_that(
  "plot_toxicity_intervals_stacked.blrmfit works for combo2 example",
  {
    a <- plot_toxicity_intervals_stacked(combo2$blrmfit,
                                         x = vars(drug_A),
                                         group = vars(group_id, drug_B))


    b <- plot_toxicity_intervals_stacked(combo2$blrmfit,
                                         x = "drug_A",
                                         group = c("group_id", "drug_B"))

    di  <- dose_info_combo2 %>%
      filter(group_id=="trial_AB") %>%
      mutate(num_toxicities=0, num_patients=4, label=paste(group_id, num_patients, sep="/"))

    c <- plot_toxicity_intervals_stacked(combo2$blrmfit,
                                         x = vars(drug_A),
                                         newdata=di,
                                         group=~ label + drug_B,
                                         predictive=TRUE)

    expect_gg(a)
    expect_gg(b)
    expect_gg(c)

  })



# plot_toxicity_intervals_stacked.blrm_trial ----------------------------------------------


test_that(
  "plot_toxicity_intervals_stacked.blrm_trial works for blrm_trial examples",
  {
    for(trial in trial_examples){

      a <- plot_toxicity_intervals_stacked(trial)
      expect_gg(a)

      a <- plot_toxicity_intervals_stacked(trial, predictive = TRUE)
      expect_gg(a)

      expect_error({
        plot_toxicity_intervals_stacked(trial, predictive = TRUE,
                                        newdata = slice(summary(trial, "dose_info"), 1))
      })

      a <- plot_toxicity_intervals_stacked(trial, predictive = TRUE,
                                           newdata = slice(summary(trial, "dose_info"), 1) %>%
                                             mutate(num_patients = 10, num_toxicities = 0))
      expect_gg(a)

      a <- plot_toxicity_intervals_stacked(trial, predictive = TRUE,
                                           interval_prob = c(-1, 0, 1, 10),
                                           newdata = slice(summary(trial, "dose_info"), 1) %>%
                                             mutate(num_patients = 10, num_toxicities = 0))
      expect_gg(a)
    }

  })



# Visual regression testing ----------------------------------------------------

test_that("plot_toxicity_curve.blrmfit renders single-agent plots correctly", {

  testthat::skip_on_cran()
  testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
  testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

  blrmfit <- single_agent$blrmfit

  # basic
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A))
  vdiffr::expect_doppelganger("single-agent blrmfit", p)

  # change intervals
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), interval_prob = c((0:5)/10, 1))
  vdiffr::expect_doppelganger("single-agent blrmfit intervals", p)

  # warning when transform = FALSE
  expect_warning(
    p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), interval_prob = c((0:5)/10, 1),
                                         transform = FALSE)
  )
  vdiffr::expect_doppelganger("single-agent blrmfit transform", p)

  # predictive
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), predictive = TRUE, ylim = c(0, 1),
                                       group = ~ num_patients)
  vdiffr::expect_doppelganger("single-agent blrmfit predictive", p)

  #+ fig.height = 12
  # pass alternative newdata
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), predictive = TRUE, ylim = c(0, 1),
                                       newdata = select(blrmfit$data, -starts_with("num")) %>%
                                         expand_grid(num_patients = 3:6,
                                                     num_toxicities = 0),
                                       group = ~ num_patients,
                                       facet_args = list(ncol = 1))
  vdiffr::expect_doppelganger("single-agent blrmfit newdata", p)

  # change interval_prob
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), predictive = TRUE, ylim = c(0, 1),
                                       newdata = select(blrmfit$data, -starts_with("num")) %>%
                                         expand_grid(num_patients = 3:6,
                                                     num_toxicities = 0),
                                       group = ~num_patients,
                                       interval_prob = (-1):6,
                                       facet_args = list(ncol = 1))
  vdiffr::expect_doppelganger("single-agent blrmfit intervals newdata", p)

  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), predictive = TRUE, ylim = c(0, 1),
                                       newdata = select(blrmfit$data, -starts_with("num")) %>%
                                         expand_grid(num_patients = 3:6,
                                                     num_toxicities = 0),
                                       group = ~num_patients,
                                       transform = TRUE,
                                       facet_args = list(ncol = 1))
  vdiffr::expect_doppelganger("single-agent blrmfit predictive transform", p)

})


test_that("plot_toxicity_intervals_stacked.blrmfit renders combo2 plots correctly", {

  testthat::skip_on_cran()
  testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
  testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

  blrmfit <- combo2$blrmfit

  nd <- filter(select(blrmfit$data, -starts_with("num")), group_id == "trial_AB")

  # basic
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), group = ~ group_id + drug_B,
                                       newdata = nd)
  vdiffr::expect_doppelganger("combo2 blrmfit", p)

  # change intervals
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), group = ~ group_id + drug_B,
                                       newdata = nd, interval_prob = c((0:5)/10, 1))
  vdiffr::expect_doppelganger("combo2 blrmfit intervals", p)

  #+ fig.height = 1.62 * 8, fig.width = 8
  # predictive
  nd <- filter(blrmfit$data, group_id == "trial_AB")
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), group = ~ num_patients + group_id + drug_B,
                                       newdata = nd, predictive = TRUE)
  vdiffr::expect_doppelganger("combo2 blrmfit predictive", p)

  #+ fig.height = 12, fig.width = 8
  # pass alternative newdata
  nd <- expand_grid(
    select(nd, -starts_with("num")),
    num_patients = 3:6,
    num_toxicities = 0
  )
  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), group = ~ num_patients + group_id + drug_B,
                                       predictive = TRUE, ylim = c(0, 1),
                                       newdata = nd,
                                       facet_args = list(ncol = 3, dir = "v"))
  vdiffr::expect_doppelganger("combo2 blrmfit predictive newdata", p)


  p <- plot_toxicity_intervals_stacked(blrmfit, x = vars(drug_A), group = ~ num_patients + group_id + drug_B,
                                       predictive = TRUE, ylim = c(0, 1),
                                       newdata = nd,
                                       facet_args = list(ncol = 3, dir = "v"),
                                       transform = TRUE)

  vdiffr::expect_doppelganger("combo2 blrmfit predictive transform", p)

})


for(ex in names(trial_examples)){

  test_that(paste("plot_toxicity_intervals_stacked.blrm_trial renders", ex, "plots correctly"), {

    testthat::skip_on_cran()
    testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
    testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

    example <- trial_examples[[ex]]

    # blrm_trial
    p <- plot_toxicity_intervals_stacked(example)
    vdiffr::expect_doppelganger(paste(ex, "blrm_trial"), p)

    p <- plot_toxicity_intervals_stacked(example, ewoc_shading = FALSE)
    vdiffr::expect_doppelganger(paste(ex, "blrm_trial no ewoc"), p)

    p <- plot_toxicity_intervals_stacked(example, predictive = TRUE)
    vdiffr::expect_doppelganger(paste(ex, "blrm_trial predictive"), p)

    p <- plot_toxicity_intervals_stacked(example, predictive = TRUE, transform = TRUE)
    vdiffr::expect_doppelganger(paste(ex, "blrm_trial predictive transform"), p)

  })

}

