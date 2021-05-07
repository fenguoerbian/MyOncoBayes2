
context("plot_toxicity_curve tests")

set.seed(123144)

eps <- 1E-4

## reference runs
single_agent  <- run_example("single_agent")
combo2  <- run_example("combo2")
combo3  <- run_example("combo3")

# Basic plot with default options ----------------------------------------------

test_that("plot_toxicity_curve.blrmfit works for single-agent example",
{
  expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A)))

  expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    transform = FALSE))

  expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    transform = FALSE,
                    ylim = c(0, 0.5)))

  expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    transform = FALSE,
                    ylim = c(0, 0.5),
                    xlim = c(0.1, 25)) +
  scale_x_log10())

  expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    transform = FALSE,
                    ylim = c(0, 0.5),
                    xlim = c(0.1, 25)))


# check nonsense variable
expect_error(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_B)))


# check character
expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = "drug_A"))

# check nonsense character
expect_error(plot_toxicity_curve(single_agent$blrmfit,
                    x = "drug_B"))

# check passing two vars()
expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    group = vars(group_id)))

# one var(), one formula
expect_gg(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    group = ~ group_id))

# one var(), one nonsense formula
expect_error(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    group = ~ drug_B))

# duplicated arguments
expect_error(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    group = vars(drug_A)))

# duplicated arguments
expect_error(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    group = ~drug_A))

# manual ylim that obscures data
expect_warning(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    ylim = c(0, 0.1)))

# manual xlim and ylim that obscures data
expect_warning(plot_toxicity_curve(single_agent$blrmfit,
                    x = vars(drug_A),
                    ylim = c(0, 0.1),
                    xlim = c(10, 20)))
})


test_that("plot_toxicity_curve.blrmfit works for combo2 example",
{
  expect_gg(plot_toxicity_curve(combo2$blrmfit,
                    x = vars(drug_A),
                    group = vars(group_id, drug_B)))

 expect_gg(plot_toxicity_curve(combo2$blrmfit,
                    x = "drug_A",
                    group = c("group_id", "drug_B")))

  expect_gg(plot_toxicity_curve(combo2$blrmfit,
                    x = "drug_A",
                    group = ~ group_id + drug_B))

    expect_gg(plot_toxicity_curve(combo2$blrmfit,
                    x = "drug_A",
                    group = ~ group_id + drug_B,
                    transform = FALSE))

  expect_gg(plot_toxicity_curve(combo2$blrmfit,
                    x = "drug_A",
                    group = ~ group_id + drug_B,
                    transform = FALSE,
                    ylim = c(0, 0.5),
                    xlim = c(0.1, 25)))

})


test_that("plot_toxicity_curve.blrmfit works for combo3 example",
{
  expect_gg(plot_toxicity_curve(combo3$blrmfit,
                    x = vars(drug_A),
                    group = vars(group_id, drug_B, drug_C)))

})







# Basic plot with default options ----------------------------------------------

test_that("plot_toxicity_intervals.blrmfit works for single-agent example",
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


test_that("plot_toxicity_intervals.blrmfit works for combo2 example",
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


test_that("plot_toxicity_intervals.blrmfit works for combo3 example",
{
  a <- plot_toxicity_intervals(combo3$blrmfit,
                    x = vars(drug_A),
                    group = vars(group_id, drug_B, drug_C))

  expect_gg(a)

})




test_that("plot_toxicity_curve.blrm_trial works for blrm_trial examples",
{
  for(example in examples[1:3]){
    trial <- with(example,
    blrm_trial(histdata, dose_info, drug_info %>% mutate(reference_p_dlt = 0.1), simplified_prior = TRUE))
    expect_gg(plot_toxicity_curve(trial))
  }

})

test_that("plot_toxicity_intervals.blrm_trial works for blrm_trial examples",
{
  for(example in examples[1:3]){
    trial <- with(example,
    blrm_trial(histdata, dose_info, drug_info %>% mutate(reference_p_dlt = 0.1), simplified_prior = TRUE))
    a <- plot_toxicity_intervals(trial)
    expect_gg(a)
  }

})