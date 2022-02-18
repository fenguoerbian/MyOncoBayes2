#' Plot a fitted model
#'
#' @description
#' \strong{Warning}: these methods are at an experimental stage of development, and
#' may change with future releases.
#'
#' Plotting methods for \code{blrmfit} and \code{blrm_trial} objects.
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
#'     the quantiles \code{0.5-prob/2} and \code{0.5+prob/2} are displayed.
#' @param prob_outer central probability mass to report for the outer ribbon, i.e.
#'     the quantiles \code{0.5-prob/2} and \code{0.5+prob/2} are displayed.
#' @param alpha,size Arguments passed to geoms. For this plot, \code{alpha} is
#'   passed to \code{\link[ggplot2:geom_ribbon]{ggplot2::geom_ribbon()}}, and \code{size} is passed to
#'   \code{\link[ggplot2:geom_line]{ggplot2::geom_line()}}.
#' @param facet_args A named list of arguments (other than `facets`) passed
#'   to \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.
#' @param hline_at Location(s) of horizontal guide lines (passed to
#'   \code{\link[bayesplot:hline_at]{bayesplot::hline_at()}}).
#' @param predictive logical indicates if the posterior predictive is
#'     being summarized. Defaults to \code{FALSE}.
#' @param interval_prob defines the interval probabilities reported in
#'     the standard outputs. Defaults to \code{c(0, 0.16, 0.33, 1)},
#'     when \code{predictive = FALSE} and/or \code{transform = TRUE}, or to intervals
#'     giving 0, 1, or 2+ DLTs when \code{predictive = TRUE} and \code{transform = FALSE}.
#'     For \code{blrm_trial} methods, this is taken from
#'     \code{summary(blrm_trial, "interval_prob")} by default.
#' @param interval_max_mass vector defining for each interval of
#'     the \code{interval_prob} vector a maximal admissible
#'     probability mass for a given dose level. Whenever the posterior
#'     probability mass in a given interval exceeds the threshold,
#'     then the Escalation With Overdose Control (EWOC) criterion is
#'     considered to be not fulfilled. Dose levels not fulfilling
#'     EWOC are ineligible for the next cohort of patients. The
#'     default restricts the overdose probability to less than 0.25. For
#'     \code{blrm_trial} methods, this is taken from
#'     \code{summary(blrm_trial, "interval_max_mass")} by default.
#' @param grid_length Number of grid points within \code{xlim} for plotting.
#' @param ewoc_shading logical indicates if doses violating EWOC should be
#'     shaded in gray. Applies only to \code{blrm_trial} methods. Defaults to \code{TRUE}.
#' @param ewoc_colors Fill colors used for bars indicating EWOC OK or not.
#'     Vector of two characters, each of which must correspond to
#'     \code{\link[bayesplot]{bayesplot-package}} color schemes
#'     (see \code{?\link[bayesplot:color_scheme_get]{bayesplot::color_scheme_get()}})
#' @param ... currently unused
#' @details
#'
#' \code{plot_toxicity_curve} plots continuous profiles of the dose-toxicity curve.
#'
#' \code{plot_toxicity_intervals} plots the posterior probability mass in
#' subintervals of [0,1], at a discrete set of provisional doses.
#'
#' \code{plot_toxicity_intervals_stacked} is similar to
#' \code{plot_toxicity_intervals}, but over a continuous range of doses.
#'
#' @return A ggplot object that can be further
#'   customized using the \code{\link[ggplot2]{ggplot2}} package.
#'
#' @template start-example
#' @examples
#'
#' example_model("combo2", silent=TRUE)
#'
#' # Plot the dose-toxicity curve
#' plot_toxicity_curve(blrmfit,
#'                     x = "drug_A",
#'                     group = ~ group_id * drug_B,
#'                     newdata = subset(dose_info_combo2, group_id == "trial_AB"),
#'                     facet_args = list(ncol = 4))
#'
#' # Plot posterior DLT-rate-interval probabilities at discrete dose levels
#' plot_toxicity_intervals(blrmfit,
#'                         x = "drug_A",
#'                         group = ~ group_id * drug_B,
#'                         newdata = subset(dose_info_combo2, group_id == "trial_AB"))
#'
#' # Plot posterior DLT-rate-interval probabilities over continuous dose
#' plot_toxicity_intervals_stacked(blrmfit,
#'                                 x = "drug_A",
#'                                 group = ~ group_id * drug_B,
#'                                 newdata = subset(dose_info_combo2, group_id == "trial_AB"))
#'
#' # Plot predictive distribution probabilities over continuous dose
#' plot_toxicity_intervals_stacked(blrmfit,
#'                                 x = "drug_A",
#'                                 group = ~ group_id * drug_B,
#'                                 predictive = TRUE,
#'                                 interval_prob = c(-1, 0, 1, 6),
#'                                 newdata = transform(subset(dose_info_combo2,
#'                                                            group_id == "trial_AB"),
#'                                                     num_patients = 6,
#'                                                     num_toxicities = 0))
#' @template stop-example
NULL

#' @rdname plot_blrm
#' @export
plot_toxicity_curve <- function(object, ...) UseMethod("plot_toxicity_curve")

#' @rdname plot_blrm
#' @export
plot_toxicity_intervals <- function(object, ...) UseMethod("plot_toxicity_intervals")

#' @rdname plot_blrm
#' @export
plot_toxicity_intervals_stacked <- function(object, ...) UseMethod("plot_toxicity_intervals_stacked")

#' @method plot_toxicity_curve default
#' @noRd
plot_toxicity_curve.default <- function(object, ...){ stop("object must inherit blrmfit or blrm_trial class") }

#' @method plot_toxicity_intervals default
#' @noRd
plot_toxicity_intervals.default <- function(object, ...){ stop("object must inherit blrmfit or blrm_trial class") }

#' @method plot_toxicity_intervals_stacked default
#' @noRd
plot_toxicity_intervals_stacked.default <- function(object, ...){ stop("object must inherit blrmfit or blrm_trial class") }

#' @rdname plot_blrm
#' @method plot_toxicity_curve blrmfit
#' @export
plot_toxicity_curve.blrmfit <- function(object,
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
                                        grid_length = 100,
                                        ...)
{

  # make R CMD CHECK happy
  ewoc_lab <- lower <- upper <- middle <- NULL

  assert_that(inherits(object, "blrmfit"), msg = "object must be of class blrmfit or blrm_trial.")

  # check probabilities
  assert_numeric(prob, lower=0, upper=1, finite=TRUE, any.missing=FALSE, min.len=1)
  assert_numeric(prob_outer, lower=0, upper=1, finite=TRUE, any.missing=FALSE, min.len=1)
  probs <- sort(c(prob, prob_outer))
  assert_numeric(grid_length, lower = 2, upper = Inf, finite = TRUE, len = 1,
                 any.missing = FALSE)

  if(missing(newdata)){
    newdata <- object$data
  }

  variables <- check_plot_variables(x, group, newdata)
  x <- variables$x
  group_variables <- variables$group_variables
  group_formula <- variables$group_formula

  if(missing(xlim)){
    xlim <- c(0, max(newdata[[x]]))
  }

  newdata_grid <- expand_newdata(newdata, xlim, object, x, group_variables, grid_length)

  posterior_summary <- summary(object,
                               newdata = newdata_grid,
                               prob = probs,
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
  ymin <- min(plot_data$lower)
  xrange <- range(newdata[[x]])

  ylim_auto <- c(0, max(0.5, ceiling(10 * ymax) / 10))

  if(!transform){
    hline_at <- logit(hline_at)
    ylim_auto <- c(logit(0.02), ylim_auto[2])
    if(!missing(ylim)){
      ylim <- logit(c(max(0.02, ylim[1]), ylim[2]))
    }
  }

  if(missing(ylim)){
    ylim <- ylim_auto
  }

  scheme <- bayesplot::color_scheme_get()

  on.exit(NULL)
  plot_data <- filter(plot_data, !!sym(x) >= xlim[1], !!sym(x) <= xlim[2])
  pl <- ggplot(plot_data, aes(x = !!sym(x))) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = prob,
                    group = prob),
                alpha = alpha) +
    geom_line(aes(y = middle, color = "Posterior median",
                  linetype = "Posterior median"),
              size = 0.75 * size, lineend = "round") +
    bayesplot::bayesplot_theme_get() + ## note: we want the currently active bayesplot theme (not always the default)
    bayesplot::hline_at(hline_at, color = "gray20", linetype = "dashed") +
    scale_fill_manual("Central posterior probability",
                      values = setNames(c(scheme[[2]], scheme[[1]]),
                                        unique(plot_data$prob))) +
    scale_x_continuous(breaks = function(u){
      breaks <- scales::extended_breaks(n = 4)
      sort(c(unique(newdata[[x]][newdata[[x]] >= xlim[1] & newdata[[x]] <= xlim[2]]), breaks(u)))
    }, limits = xlim,
    oob = scales::squish) +
    scale_color_manual(values = setNames(scheme[[5]], "Posterior median")) +
    scale_linetype_manual(values = setNames(1, "Posterior median")) +
    guides(color = guide_legend(NULL),
           linetype = guide_legend(NULL))

  if(transform){
    pl <-  pl + scale_y_continuous(breaks = sort(c((0:10) / 10, hline_at)),
                                   minor_breaks = NULL,
                                   limits = ylim,
                                   labels = label_percent(accuracy = 1),
                                   oob = scales::squish) +
      labs(x = x,
           y = "P(DLT)")
  } else if(!transform){
    pl <-  pl + scale_y_continuous(
      breaks = function(u){
        breaks <- scales::extended_breaks(n = 6)
        round(sort(c(hline_at, breaks(u))), 2)
      },
      minor_breaks = NULL,
      limits = ylim,
      oob = scales::squish) +
      labs(x = x,
           y = "logit(P(DLT))")
  }

  if(!missing(group)) {
    facet_args <- modifyList(list(facets = group_formula,
                                  scales = "fixed",
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
                                           hline_at,
                                           grid_length = 100,
                                           ewoc_shading = TRUE,
                                           ...)
{

  .assert_is_blrm_trial_and_prior_is_set(object)
  drug_info <- summary(object, "drug_info")
  ewoc_lab <- NULL

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

  plot_base <- plot_toxicity_curve(object$blrmfit,
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

  if(ewoc_shading){

    variables <- check_plot_variables(x, group, newdata)
    x <- variables$x
    group_variables <- variables$group_variables
    group_formula <- variables$group_formula

    if(missing(xlim)){
      xlim <- c(0, max(newdata[[x]]))
    }

    newdata_grid <- as_tibble(expand_newdata(newdata, xlim, object$blrmfit, x, group_variables, grid_length))

    ewoc_data <- summary(object, "newdata_prediction", newdata = newdata_grid)
    ewoc_data$ewoc_lab <- factor(ewoc_data$ewoc_ok,
                                 c(TRUE, FALSE),
                                 c("OK", "Not OK"))

    ylim <- plot_base$scales$get_scales("y")$limits

    plot_out <- plot_base +
      scale_alpha_manual("EWOC", values=c(0, 0.5), breaks=c("OK", "Not OK"), drop=FALSE) +
      geom_ribbon(data = ewoc_data,
                  aes(ymin = ylim[1], ymax = ylim[2], alpha = ewoc_lab),
                  fill="black")

  } else{

    plot_out <- plot_base

  }

  on.exit(NULL)
  return(plot_out)

}

#' @rdname plot_blrm
#' @method plot_toxicity_intervals blrmfit
#' @export
plot_toxicity_intervals.blrmfit <- function(object,
                                            newdata,
                                            x,
                                            group,
                                            interval_prob = c(0, 0.16, 0.33, 1),
                                            interval_max_mass = c(NA, NA, 0.25),
                                            ewoc_colors = c("green", "red"),
                                            ...)
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

  scheme <- bayesplot::color_scheme_get(paste0("mix-", paste(ewoc_colors, collapse = "-")))

  posterior_summary <- summary(object, newdata = newdata, interval_prob = interval_prob)

  interval_labs <- paste0("(", paste(interval_prob[seq_len(length(interval_prob) - 1)],
                                     interval_prob[1 + seq_len(length(interval_prob) - 1)], sep = ","), "]")

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
    ylab("Posterior probability mass") +
    xlab(x)

  pl


  return(pl)

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
                                               ewoc_colors = c("green", "red"),
                                               ...)
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


#' @rdname plot_blrm
#' @method plot_toxicity_intervals_stacked blrmfit
#' @export
plot_toxicity_intervals_stacked.blrmfit <- function(object,
                                                    newdata,
                                                    x,
                                                    group,
                                                    xlim,
                                                    ylim = c(0, 0.5),
                                                    predictive = FALSE,
                                                    transform = !predictive,
                                                    interval_prob,
                                                    grid_length = 100,
                                                    facet_args = list(),
                                                    ...)
{



  # make R CMD CHECK happy
  prob <- cutoff <- ewoc_ok <- cprob <- interval <- NULL

  assert_numeric(grid_length, lower = 2, upper = Inf, finite = TRUE, len = 1,
                 any.missing = FALSE)

  if(!transform & !predictive){
    warning("transform = FALSE not meaningful when predictive = FALSE. Setting transform = TRUE.")
    transform <- TRUE
  }

  if(missing(newdata)){
    newdata <- object$data
  }

  variables <- check_plot_variables(x, group, newdata)
  x <- variables$x
  group_variables <- variables$group_variables
  group_formula <- variables$group_formula

  if(predictive){

    num_patients <- pp_binomial_trials(object, newdata)
    newdata$num_patients <- num_patients
    newdata$num_toxicities <- 0

    if(length(unique(num_patients)) > 1){

      covarnames <- all.vars(delete.response(terms(object$formula)))
      if(all(group_variables %in% covarnames)){
        message(paste("Data contains cohorts of different sizes. Suggest",
                      "including the cohort size as a grouping variable via",
                      "the group argument, or overplotting may occur."))
      }
    }

  }

  if(missing(xlim)){
    xlim <- c(0, max(newdata[[x]]))
  }

  newdata_grid <- expand_newdata(newdata, xlim, object, x, group_variables, grid_length, predictive)


  if(missing(interval_prob)){

    if(predictive & !transform){
      interval_prob <- unique(c(-1, 0, 1, max(num_patients)))
    } else if(predictive & transform){
      interval_prob <- c(-1, 0, 0.16, 0.33, 1)
    } else if(!predictive){
      interval_prob <- c(0, 0.16, 0.33, 1)
    }

  } else{

    if(predictive & !transform){
      assert_numeric(ylim, lower = 0, upper = Inf, any.missing = FALSE)
      validate_that(max(interval_prob) >= max(num_patients),
                    msg = paste("interval_prob does not cover the full support",
                                "of the predictive distribution: largest value",
                                "in interval_prob is smaller than",
                                "the largest cohort in newdata."))
      validate_that(min(interval_prob) < 0,
                    msg = paste("interval_prob does not cover the full support",
                                "of the predictive distribution: largest value",
                                "in interval_prob is smaller than",
                                "the largest cohort in newdata."))
    } else if(predictive & transform){
      assert_numeric(ylim, lower = 0, upper = 1, any.missing = FALSE)
      validate_that(max(interval_prob) >= 1,
                    msg = paste("interval_prob does not cover the full support",
                                "of the predictive distribution for the",
                                "proportion of DLTs: largest value",
                                "in interval_prob is less than 1."))
      validate_that(min(interval_prob) < 0,
                    msg = paste("interval_prob does not cover the full support",
                                "of the predictive distribution for the",
                                "proportion of DLTs: smallest",
                                "in interval_prob not less than 0."))

    } else if(!predictive){
      assert_numeric(ylim, lower = 0, upper = 1, any.missing = FALSE)
      validate_that(max(interval_prob) >= 1,
                    msg = paste("interval_prob does not cover the full support",
                                "of the posterior distribution for the",
                                "DLT rate: largest value",
                                "in interval_prob is less than 1."))
      validate_that(min(interval_prob) <= 0,
                    msg = paste("interval_prob does not cover the full support",
                                "of the posterior distribution for the",
                                "DLT rate: smallest value",
                                "in interval_prob is not 0."))

    }

  }

  # check probabilities
  assert_numeric(interval_prob, any.missing = FALSE, sorted = TRUE)

  interval_labs <- paste0("(", paste(interval_prob[seq_len(length(interval_prob) - 1)],
                                     interval_prob[1 + seq_len(length(interval_prob) - 1)], sep = ","), "]")

  sum_data <- summary(object,
                      newdata = newdata_grid,
                      interval_prob = interval_prob,
                      transform = transform,
                      predictive = predictive)

  x_name <- as.name(x)
  sum_data <- bind_cols(newdata_grid, sum_data)

  stacked0 <- select(sum_data, group_variables, x, starts_with("("))

  if(!is.null(group_variables)){
    stacked0$group <- apply(stacked0[c(group_variables)], 1, paste, collapse = "_")
  } else{
    stacked0$group <- "one_group"
  }

  stacked1 <- bind_rows(lapply(split(stacked0, stacked0$group),
                               function(stacked_group){
                                 out0 <- stacked_group[order(stacked_group[[x]]),]
                                 out0 %>%
                                   pivot_longer(starts_with("("), names_to="interval", values_to="prob")
                               }))

  stacked1$group <- apply(stacked1[c(group_variables, x)], 1, paste, collapse = "_")
  stacked <- bind_rows(lapply(split(stacked1, stacked1$group),
                              function(stacked_group){
                                stacked_group$cprob <- cumsum(stacked_group$prob)
                                stacked_group
                              }))

  legend_lab <- "Toxicity Interval"
  if(predictive){

    if(!transform){

      support <- tibble(
        support = 0:max(num_patients),
        interval = cut(support, interval_prob)
      ) %>%
        group_by(interval) %>%
        summarize(from = min(support),
                  to = max(support)) %>%
        mutate(label = case_when(
          from == to ~ as.character(from),
          to == max(num_patients) ~ paste0(from, "+"),
          from < to ~ paste(from, to, sep = "-")
        ))
      relabel <- setNames(split(support$interval, 1:nrow(support)),
                          support$label)

      stacked$interval <- factor(stacked$interval)
      levels(stacked$interval) <- relabel

      legend_lab <- "Number of DLTs"



    } else if(transform){

      legend_lab <- "Toxicity Interval\n(Proportion of DLTs)"

    }


  }

  on.exit(NULL)
  pl <- ggplot(mapping = aes(x=!!sym(x)), data = stacked) +
    geom_ribbon(aes(ymin=1-(cprob-prob), ymax=1-cprob, fill=factor(interval))) +
    scale_y_continuous(limits = ylim, oob = scales::squish) +
    scale_fill_brewer(legend_lab, type="div", palette="RdYlBu", direction=-1,
                      drop = FALSE) +
    ylab(paste0(ifelse(predictive, "Predictive\n", "Toxicity Interval\n"), "Probability")) +
    scale_x_continuous(x, breaks = function(u){
      breaks <- scales::extended_breaks(n = 4)
      sort(c(unique(newdata[[x]]), breaks(u)))
    }, limits = xlim) +
    bayesplot::bayesplot_theme_get()

  if(!missing(group)) {
    facet_args <- modifyList(list(facets = group_formula,
                                  scales = "fixed",
                                  strip.position = "top",
                                  labeller="label_both",
                                  nrow = NULL,
                                  ncol = NULL), facet_args)
    pl  <- pl + do.call(facet_wrap, facet_args)
  }

  return(pl)
}


#' @rdname plot_blrm
#' @method plot_toxicity_intervals_stacked blrm_trial
#' @export
plot_toxicity_intervals_stacked.blrm_trial <- function(object,
                                                       newdata,
                                                       x,
                                                       group,
                                                       xlim,
                                                       ylim = c(0, 0.5),
                                                       predictive = FALSE,
                                                       transform = !predictive,
                                                       interval_prob,
                                                       grid_length = 100,
                                                       ewoc_shading = TRUE,
                                                       facet_args = list(),
                                                       ...)
{

  .assert_is_blrm_trial_and_prior_is_set(object)
  drug_info <- summary(object, "drug_info")
  ewoc_lab <- NULL

  if(missing(newdata)){
    newdata <- summary(object, "dose_info")
    if(predictive){
      message('By default, predictive distributions for cohorts of size 6 at summary(object, "dose_info") will be used.')
      newdata <- mutate(newdata, num_patients = 6, num_toxicities = 0)
    }
  }

  if(predictive){
    assert_that(has_name(newdata, c("num_patients", "num_toxicities")),
                msg = paste("For predictive plot, newdata must contain",
                            "num_patients and num_toxicities columns"))
  }

  if(missing(x)) x <- drug_info$drug_name[1]
  if(missing(group)){
    num_patients_name <- NULL
    if(predictive && length(unique(newdata$num_patients)) > 1){
      num_patients_name <- "num_patients"
    }
    if(nrow(drug_info) > 1){
      group <- c("group_id", drug_info$drug_name[2:nrow(drug_info)], num_patients_name)
    } else{
      group <- c("group_id", num_patients_name)
    }
  }

  if(missing(interval_prob)){
    if(predictive){
      if(!transform) interval_prob <- unique(c(-1, 0, 1, max(newdata$num_patients)))
      if(transform) interval_prob <- c(-1, object$interval_prob)
    } else{
      interval_prob <- object$interval_prob
    }
  }

  plot_base <- plot_toxicity_intervals_stacked(object$blrmfit,
                                               newdata,
                                               x,
                                               group,
                                               xlim,
                                               ylim,
                                               predictive,
                                               transform,
                                               interval_prob,
                                               grid_length,
                                               facet_args)

  on.exit(NULL)

  if(ewoc_shading){

    variables <- check_plot_variables(x, group, newdata)
    x <- variables$x
    group_variables <- variables$group_variables
    group_formula <- variables$group_formula

    if(missing(xlim)){
      xlim <- c(0, max(newdata[[x]]))
    }


    newdata_grid <- as_tibble(expand_newdata(newdata, xlim, object$blrmfit, x, group_variables, grid_length, predictive))

    ewoc_data <- summary(object, "newdata_prediction", newdata = newdata_grid)

    ewoc_data$ewoc_lab <- factor(ewoc_data$ewoc_ok,
                                 c(TRUE, FALSE),
                                 c("OK", "Not OK"))

    plot_out <- plot_base +
      scale_alpha_manual("EWOC", values=c(0, 0.5), breaks=c("OK", "Not OK"), drop=FALSE) +
      geom_ribbon(data = ewoc_data,
                  aes(ymin = ylim[1], ymax = ylim[2], alpha = ewoc_lab),
                  fill="black")

  } else{

    plot_out <- plot_base

  }


  return(plot_out)

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


#' Internal function for checking the variable specifications in the plotting
#' functions
#'
#' @noRd
#' @template args-plot
#' @param newdata data frame of covariate levels for plotting
#' @param xlim x-axis limits
#' @param object blrmfit object
#' @keywords internal
expand_newdata <- function(newdata,
                           xlim,
                           object,
                           x,
                           group_variables,
                           grid_length,
                           predictive = FALSE){


  formula <- object$formula
  env <- environment(formula)

  covarnames <- all.vars(delete.response(terms(formula)))
  covarnames <- unique(c(covarnames, x, group_variables))
  covarnames <- covarnames[covarnames %in% names(newdata)] # modified from base::get_all_vars()

  if(predictive){

    covarnames <- unique(c(covarnames, "num_patients", "num_toxicities"))

  }

  inp <- parse(text = paste0("list(", paste(covarnames, collapse = ","),
                             ")"), keep.source = FALSE)
  covariates <- as_tibble(setNames(eval(inp, newdata, env),
                                   covarnames))

  xseq <- sort(c(unique(covariates[[x]])[unique(covariates[[x]]) >= xlim[1] & unique(covariates[[x]]) <= xlim[2]],
                 seq(xlim[1], xlim[2], length.out = grid_length)))
  xgrid <- setNames(list(xseq), x)
  len <- length(xseq)

  newdata_unique <- unique(covariates[names(covariates) != x])

  newdata_grid <- newdata_unique[0,]
  for(j in 1:nrow(newdata_unique)){
    nd_j <- cbind(
      bind_rows(lapply(seq_len(len), function(i) newdata_unique[j,, drop = FALSE])),
      xgrid
    )
    newdata_grid <- rbind(newdata_grid,
                          nd_j)
  }

  newdata_grid
}

