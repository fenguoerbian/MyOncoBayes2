#' Plot a fitted model
#'
#' @description
#' **Warning**: these methods are at an experimental stage of development, and
#' may change with future releases.
#'
#' Plotting methods for \code{blrmfit} and \code{blrm_trial} objects.
#' \code{plot_toxicity_curve} plots continuous profiles of the dose-toxicity curve.
#' \code{plot_toxicity_intervals} plots the posterior probability mass in
#' subintervals of [0,1].
#'
#' @name plot_blrm
#'
#' @param object fitted model object
#' @param newdata optional data frame specifying for what to predict;
#'     if missing, then the data of the input model \code{object} is
#'     used. If \code{object} is a \code{blrmfit} object, \code{newdata} defaults to
#'     the \code{data} argument. If \code{object} is a \code{blrm_trial}, it defaults
#'     to \code{summary(object, "dose_info")}.
#' @template args-plot
#' @param xlim x-axis limits
#' @param ylim y-axis limits on the probability scale
#' @template args-transform
#' @param prob central probability mass to report for the inner ribbon, i.e.
#'     the quantiles 0.5-prob/2 and 0.5+prob/2 are displayed.
#' @param prob_outer central probability mass to report for the outer ribbon, i.e.
#'     the quantiles 0.5-prob/2 and 0.5+prob/2 are displayed.
#' @param alpha,size Arguments passed to geoms. For this plot, \code{alpha} is
#'   passed to [ggplot2::geom_ribbon()], and \code{size} is passed to
#'   [ggplot2::geom_line].
#' @param facet_args A named list of arguments (other than `facets`) passed
#'   to [ggplot2::facet_wrap()].
#' @param hline_at Location(s) of horizontal guide lines (passed to
#'   [bayesplot::hline_at]).
#' @param grid_length Number of grid points within \code{xlim} for plotting.
#' @param interval_prob defines the interval probabilities reported in
#'     the standard outputs. Defaults to \code{c(0, 0.16, 0.33, 1)}.
#' @param interval_max_mass vector defining for each interval of
#'     the \code{interval_prob} vector a maximal admissible
#'     probability mass for a given dose level. Whenever the posterior
#'     probability mass in a given interval exceeds the threshold,
#'     then the Escalation With Overdose Control (EWOC) criterion is
#'     considered to be not fulfilled. Dose levels not fulfilling
#'     EWOC are ineligible for the next cohort of patients. The
#'     default restricts the overdose probability to less than 0.25.
#' @param ewoc_colors Fill colors used for bars indicating EWOC OK or not.
#'     Vector of two characters, each of which must correspond to 
#'     \code{bayesplot} package color schemes 
#'     (see \code{?bayesplot::color_scheme_get})
#' @details Plots the dose-toxicity curve, with P(DLT) on the vertical axis
#' and dose of "drug_name_x" on the horizontal.
#'
#' @return A ggplot object that can be further
#'   customized using the **ggplot2** package.
#'
#' @seealso \code{\link{plot_toxicity_intervals}}
#'
#' @template start-example
#' @examples
#' library(dplyr) # for vars()
#'
#' example_model("combo2")
#' plot_toxicity_curve(blrmfit,
#'                     x = vars(drug_A),
#'                     group = ~ group_id * drug_B,
#'                     newdata = filter(dose_info_combo2, group_id == "trial_AB"),
#'                     facet_args = list(ncol = 4))
#' plot_toxicity_intervals(blrmfit,
#'                         x = vars(drug_A),
#'                         group = ~ group_id * drug_B,
#'                         newdata = filter(dose_info_combo2, group_id == "trial_AB"))
#' @template stop-example
NULL

#' @rdname plot_blrm
#' @export
plot_toxicity_curve <- function(object, newdata, x, group, xlim, ylim, transform,
                                prob, prob_outer, size, alpha, facet_args = list(),
                                hline_at, grid_length) UseMethod("plot_toxicity_curve")

#' @rdname plot_blrm
#' @export
plot_toxicity_intervals <- function(object, newdata, x, group, interval_prob, interval_max_mass,
                                    ewoc_colors) UseMethod("plot_toxicity_intervals")

#' @rdname plot_blrm
#' @method plot_toxicity_curve default
#' @export
plot_toxicity_curve.default <- function(object,
                                        newdata,
                                        x,
                                        group,
                                        xlim,
                                        ylim,
                                        transform = TRUE,
                                        prob = 0.5,
                                        prob_outer = 0.95,
                                        size = 0.75,
                                        alpha = 1,
                                        facet_args = list(),
                                        hline_at = c(0.16, 0.33),
                                        grid_length = 100)
{

  assert_that(inherits(object, "blrmfit"), msg = "object must be of class blrmfit or blrm_trial.")

  # check probabilities
  assert_numeric(prob, lower=0, upper=1, finite=TRUE, any.missing=FALSE, min.len=1)
  assert_numeric(prob_outer, lower=0, upper=1, finite=TRUE, any.missing=FALSE, min.len=1)
  probs <- sort(c(prob, prob_outer))
  assert_numeric(grid_length, lower = 2, upper = Inf, finite = TRUE, len = 1,
                 any.missing = FALSE)

  if(missing(newdata))
    newdata <- object$data

  variables <- check_plot_variables(x, group, newdata)
  x <- variables$x
  group_variables <- variables$group_variables
  group_formula <- variables$group_formula

  if(missing(xlim)){
    xlim <- c(0, max(newdata[[x]]))
  }

  xgrid <- setNames(list(seq(xlim[1], xlim[2], length.out = grid_length)), x)

  newdata_grid <- newdata[!names(newdata) %in% x] %>%
    distinct() %>%
    expand_grid(as.data.frame(xgrid))

  posterior_summary <- summary(object, newdata = newdata_grid, prob = probs,
                               transform = transform)

  ribbon_data <- lapply(probs,
                        function(p){
                          plot_data <- newdata_grid
                          lab <- paste0(100 * p, "%")
                          plot_data$prob <- lab
                          plot_data$lower <- posterior_summary[, paste0(100 * (1 - p) / 2, "%")]
                          plot_data$upper <- posterior_summary[, paste0(100 * (1 + p) / 2, "%")]
                          plot_data$middle <- posterior_summary[, "50%"]
                          plot_data
                        })
  plot_data <- rbind(ribbon_data[[1]], ribbon_data[[2]])
  plot_data$prob <- factor(plot_data$prob, rev(unique(plot_data$prob)))

  ymax <- max(plot_data$upper)
  xrange <- range(newdata[[x]])

  ylim_auto <- c(0, max(0.5, ceiling(10 * ymax) / 10))

  if(!transform){
    hline_at <- logit(hline_at)
    ylim_auto <- c(logit(0.02), ylim_auto[2])
    if(!missing(ylim)){
      ylim <- logit(c(max(0.02, ylim[1]), ylim[2]))
    }
  }

  # check if any data is obscured; truncate as needed
  if(missing(ylim)){
    ylim <- ylim_auto
  } else{
    if(transform && ylim[1] != 0){
      warning("Lower y-axis limit recommended to be set to 0.")
    }
    if(ylim[2] < ymax){
      warning(paste("Upper y-axis limit less than some data in the plot.",
                    "Modify ylim or some data",
                    "will be obscured."))
    }
  }

  plot_data$upper <- pmin(plot_data$upper, ylim[2])
  plot_data$lower <- pmax(plot_data$lower, ylim[1])

  if(xlim[2] < xrange[2]){
    warning(paste("Upper x-axis limit less than some data in the plot.",
                  "Modify xlim or some data",
                  "will be obscured."))
  }
  if(xlim[1] > xrange[1]){
    warning(paste("Lower x-axis limit less than some data in the plot.",
                  "Modify xlim or some data",
                  "will be obscured."))
  }

  scheme <- bayesplot::color_scheme_get()

  on.exit(NULL)
  pl <- ggplot(plot_data, aes_string(x = x)) +
    geom_ribbon(aes_(ymin = ~ lower, ymax = ~ upper, fill = ~ prob,
                     group = ~ prob),
                alpha = alpha) +
    geom_line(aes_(y = ~ middle, color = "Posterior median",
                   linetype = "Posterior median"),
              size = 0.75 * size, lineend = "round") +
    bayesplot::bayesplot_theme_get() + ## note: we want the currently active bayesplot theme (not always the default)
    bayesplot::hline_at(hline_at, color = "gray20", linetype = "dashed") +
    scale_fill_manual("Central posterior probability",
                      values = setNames(c(scheme[[2]], scheme[[1]]),
                                        unique(plot_data$prob))) +
    scale_x_continuous(breaks = function(u){
      breaks <- scales::extended_breaks(n = 4)
      sort(c(unique(newdata[[x]]), breaks(u)))
    }, limits = xlim) +
    scale_color_manual(values = setNames(scheme[[5]], "Posterior median")) +
    scale_linetype_manual(values = setNames(1, "Posterior median")) +
    guides(color = guide_legend(NULL),
           linetype = guide_legend(NULL))

  if(transform){
    pl <-  pl + scale_y_continuous(breaks = sort(c((0:10) / 10, hline_at)),
                                   minor_breaks = NULL,
                                   limits = ylim,
                                   labels = label_percent(accuracy = 1)) +
      labs(x = x,
           y = "P(DLT)")
  } else if(!transform){
    pl <-  pl + scale_y_continuous(
      breaks = function(u){
        breaks <- scales::extended_breaks(n = 6)
        round(sort(c(hline_at, breaks(u))), 2)
      },
      minor_breaks = NULL,
      limits = ylim) +
      labs(x = x,
           y = "logit(P(DLT))")
  }

  if(!missing(group)) {
    facet_args <- modifyList(list(facets = group_formula,
                                  scales = "free",
                                  strip.position = "top",
                                  labeller="label_both",
                                  nrow = NULL,
                                  ncol = NULL), facet_args)
    pl  <- pl + do.call(facet_wrap, facet_args)
  }

  pl

}

#' @rdname plot_blrm
#' @method plot_toxicity_curve blrm_trial
#' @export
plot_toxicity_curve.blrm_trial <- function(object,
                                           newdata,
                                           x,
                                           group,
                                           xlim,
                                           ylim,
                                           transform = TRUE,
                                           prob = 0.5,
                                           prob_outer = 0.95,
                                           size = 0.75,
                                           alpha = 1,
                                           facet_args = list(),
                                           hline_at = c(0.16, 0.33),
                                           grid_length = 100)
{

  .assert_is_blrm_trial_and_prior_is_set(object)
  drug_info <- summary(object, "drug_info")

  if(missing(x)) x <- drug_info$drug_name[1]
  if(missing(group)){
    if(nrow(drug_info) > 1){
      group <- c("group_id", drug_info$drug_name[2:nrow(drug_info)])
    } else{
      group <- "group_id"
    }
  }

  if(missing(hline_at)){
    hline_at <- object$interval_prob[object$interval_prob > 0 & object$interval_prob < 1]
  }

  if(missing(newdata)) newdata <- summary(object, "dose_info")

  plot_toxicity_curve(object$blrmfit,
                      newdata,
                      x,
                      group,
                      xlim,
                      ylim,
                      transform,
                      prob,
                      prob_outer,
                      size,
                      alpha,
                      facet_args,
                      hline_at,
                      grid_length)

}

#' @rdname plot_blrm
#' @method plot_toxicity_intervals default
#' @export
plot_toxicity_intervals.default <- function(object,
                                            newdata,
                                            x,
                                            group,
                                            interval_prob = c(0, 0.16, 0.33, 1),
                                            interval_max_mass = c(NA, NA, 0.25),
                                            ewoc_colors = c("green", "red"))
{

  assert_that(inherits(object, "blrmfit"), msg = "object must be of class blrmfit or blrm_trial.")

  # check probabilities
  assert_numeric(interval_prob, any.missing=FALSE, sorted=TRUE)

  # make R CMD CHECK happy
  .x  <- value <- ewoc_ok <- interval <- cutoff <- NULL

  if(missing(newdata))
    newdata <- object$data

  variables <- check_plot_variables(x, group, newdata)
  x <- variables$x
  group_variables <- variables$group_variables
  group_formula <- variables$group_formula


  posterior_summary <- summary(object, newdata = newdata, interval_prob = interval_prob)

  interval_labs <- paste0("(", paste(interval_prob[1:(length(interval_prob) - 1)],
                                     interval_prob[2:length(interval_prob)], sep = ","), "]")

  assert_that(length(interval_max_mass) == length(interval_labs),
              msg = paste("interval_prob and cutoffs have inconsistent lengths.",
                          "Need length(cutoffs) = length(interval_probs) - 1."))

  # cutoffs[is.na(cutoffs)] <- 1
  cutoffs <- data.frame(interval = factor(interval_labs, rev(interval_labs)),
                        cutoff = interval_max_mass)

  plot_data <- cbind(newdata, posterior_summary[,interval_labs]) %>%
    pivot_longer(ends_with("]"),
                 names_to = "interval",
                 values_to = "value") %>%
    mutate(interval = factor(interval, rev(interval_labs))) %>%
    left_join(cutoffs, "interval") %>%
    tidyr::replace_na(list(cutoff = 1)) %>%
    mutate(ewoc_ok = ifelse(value <= cutoff, "Yes", "No"))

  bayesplot::color_scheme_set(paste0("mix-", paste(ewoc_colors, collapse = "-")))
  scheme <- bayesplot::color_scheme_get()

  if(!missing(group)){
    plot_data$group <- apply(plot_data[, group_variables], 1, function(x){
      paste(paste(group_variables, x, sep = ": "), collapse = "\n")
    })
  }

  plot_data$.x <- factor(plot_data[[x]], sort(unique(plot_data[[x]])))

  on.exit(NULL)
  pl <- ggplot(data = plot_data,
               mapping = aes(x = .x, y = value, fill = ewoc_ok)) +
    geom_col() +
    (if(missing(group)) facet_grid(interval ~ .) else facet_grid(interval ~ group)) +
    scale_fill_manual(name = "Escalation\nWith\nOverdose\nControl\nOK?",
                      values = setNames(c(scheme[[4]], scheme[[3]]),
                                        c("Yes", "No")),
                      limits = c("Yes", "No"))+
    ylim(0, 1) +
    geom_hline(data = filter(cutoffs, !is.na(cutoff)),
               mapping = aes(yintercept = cutoff),
               linetype = "dashed") +
    bayesplot::bayesplot_theme_get() +
    ylab("Posterior probability mass") +
    xlab(x)

  pl

}

#' @rdname plot_blrm
#' @method plot_toxicity_intervals blrm_trial
#' @export
plot_toxicity_intervals.blrm_trial <- function(object,
                                               newdata,
                                               x,
                                               group,
                                               interval_prob,
                                               interval_max_mass,
                                               ewoc_colors = c("green", "red"))
{

  .assert_is_blrm_trial_and_prior_is_set(object)
  drug_info <- summary(object, "drug_info")

  if(missing(x)) x <- drug_info$drug_name[1]
  if(missing(group)){
    if(nrow(drug_info) > 1){
      group <- c("group_id", drug_info$drug_name[2:nrow(drug_info)])
    } else{
      group <- "group_id"
    }
  }

  if(missing(interval_prob)){
    interval_prob <- object$interval_prob
  }
  if(missing(interval_max_mass)){
    interval_max_mass <- object$interval_max_mass
  }

  if(missing(newdata)) newdata <- summary(object, "dose_info")

  plot_toxicity_intervals(object$blrmfit,
                          newdata,
                          x,
                          group,
                          interval_prob,
                          interval_max_mass,
                          ewoc_colors)

}

# internal ----------------------------------------------------------------

#' Internal function for tidy parameter selection. See bayesplot:::tidyselect_parameters
#'
#' @noRd
#' @param complete_pars A character vector of *all* parameter names.
#' @param pars_list A list of columns generated by `vars()`.
#' @return Character vector of selected parameter names.
#' @keywords internal
tidyselect_parameters <- function(complete_pars, pars_list) {
  # We use the list of helpers so that we don't have to keep track of any
  # changes to tidyselect. We use `env_bury()`` so that the definitions of
  # selection helpers are available. This pattern is taken from the example code
  # in `vars_select_helpers`.
  helpers <- tidyselect::vars_select_helpers
  pars_list <- lapply(pars_list, rlang::env_bury, !!! helpers)
  selected <- tidyselect::vars_select(.vars = complete_pars, !!! pars_list)
  if (!length(selected)) {
    stop("No parameters were found matching those names.")
  }
  unname(selected)
}

#' Internal function from scales package
#' @noRd
#' @keywords internal
label_percent <- function (accuracy = NULL, scale = 100, prefix = "", suffix = "%",
                           big.mark = " ", decimal.mark = ".", trim = TRUE,
                           ...)
{
  scales::number_format(accuracy = accuracy, scale = scale, prefix = prefix,
                        suffix = suffix, big.mark = big.mark, decimal.mark = decimal.mark,
                        trim = trim, ...)

}

#' Internal function for checking the variable specifications in the plotting
#' functions
#'
#' @noRd
#' @template args-plot
#' @param newdata data frame of covariate levels for plotting
#' @keywords internal
check_plot_variables <- function(x, group, newdata){
  # check x
  if (rlang::is_quosures(x)) {
    assert_that(length(x) == 1, msg = "x must have length 1.")
    x <- tidyselect_parameters(complete_pars = names(newdata), pars_list = x)
  } else{
    assert_character(x, len = 1, any.missing = FALSE)
    if(!x %in% names(newdata)){
      stop(paste0("Variable name x = '", x, "' doesn't match names(newdata)."))
    }
  }

  # check group
  if(!missing(group)){
    if(inherits(group, "formula")){
      group_formula <- group
      group_terms <- terms(group, data = newdata)
      group_variables <- attr(group_terms, "term.labels")
      group_variables <- group_variables[!grepl(":", group_variables)] # exclude interactions
    } else if(inherits(group, "character")){
      group_variables <- group
      group_formula <- as.formula(paste("~", paste(group_variables, collapse = "+")))
    } else if(rlang::is_quosures(group)){
      group_variables <- tidyselect_parameters(complete_pars = names(newdata), pars_list = group)
      group_formula <- as.formula(paste("~", paste(group_variables, collapse = "+")))
    } else{
      stop("Unrecognized class for group argument.")
    }
    assert_that(all(group_variables %in% names(newdata)),
                msg = paste("Variable names",
                            paste(group_variables[!group_variables %in% names(newdata)],
                                  sep = ","),
                            "do not match names(newdata)."))
    assert_that(!x %in% group_variables,
                msg = paste(x, "cannot be both the x-axis and a grouping variable."))
    return(list(x = x, group_variables = group_variables, group_formula = group_formula))
  } else{
    return(list(x = x))
  }

}
