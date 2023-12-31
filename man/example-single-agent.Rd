% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example_single_agent.R
\name{example-single-agent}
\alias{example-single-agent}
\title{Single Agent Example}
\description{
Example using a single experimental drug.
}
\details{
The single agent example is described in the reference
Neuenschwander, B. et al (2008). The data are described
in the help page for \code{hist_SA}. In this case, the data
come from only one study, with the treatment being only single
agent. Hence the model specified does not involve a hierarchical
prior for the intercept and log-slope parameters. The model
described in Neuenschwander, et al (2008) is adapted as follows:
\deqn{\mbox{logit}\, \pi(d) = \log\, \alpha + \beta \, \log\, \Bigl(\frac{d}{d^*}\Bigr),}
where \eqn{d^* = 250}, and the prior for
\eqn{\boldsymbol\theta = (\log\, \alpha, \log\, \beta)} is
\deqn{\boldsymbol\theta \sim \mbox{N}(\boldsymbol m, \boldsymbol S),}
and \eqn{\boldsymbol m = (\mbox{logit}\, 0.5, \log\, 1)} and
\eqn{\boldsymbol S = \mbox{diag}(2^2, 1^2)} are constants.

In the \code{blrm_exnex} framework, in which
the prior must be specified as a hierarchical model
\eqn{\boldsymbol\theta \sim \mbox{N}(\boldsymbol \mu, \boldsymbol \Sigma)}
with additional priors on \eqn{\boldsymbol\mu} and \eqn{\boldsymbol\Sigma},
the simple prior distribution above is accomplished
by fixing the diagonal elements \eqn{\tau^2_\alpha} and \eqn{\tau^2_\beta}
of \eqn{\boldsymbol\Sigma} to zero, and taking
\deqn{\boldsymbol\mu \sim \mbox{N}(\boldsymbol m, \boldsymbol S).}

The arguments \code{prior_tau_dist} and \code{prior_EX_tau_mean_comp}
as specified below ensure that the \eqn{\tau}'s are fixed at zero.
}
\examples{
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(OncoBayes2.MC.warmup=10, OncoBayes2.MC.iter=20, OncoBayes2.MC.chains=1,
                            OncoBayes2.MC.save_warmup=FALSE)

## Example from Neuenschwander, B., et al. (2009). Stats in Medicine

num_comp <- 1 # one investigational drug
num_inter <- 0 # no drug-drug interactions need to be modeled
num_groups <- nlevels(hist_SA$group_id) # no stratification needed
num_strata <- 1 # no stratification needed


dref <- 50

## Since there is no prior information the hierarchical model
## is not used in this example by setting tau to (almost) 0.
blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + log(drug_A / dref) |
    0 |
    group_id,
  data = hist_SA,
  prior_EX_mu_mean_comp = matrix(
    c(logit(1/2), # mean of intercept on logit scale
      log(1)),    # mean of log-slope on logit scale
    nrow = num_comp,
    ncol = 2
  ),
  prior_EX_mu_sd_comp = matrix(
    c(2,  # sd of intercept
      1), # sd of log-slope
    nrow = num_comp,
    ncol = 2
  ),
  ## Here we take tau as known and as zero.
  ## This disables the hierarchical prior which is
  ## not required in this example as we analyze a
  ## single trial.
  prior_EX_tau_mean_comp = matrix(
    c(0, 0),
    nrow = num_comp,
    ncol = 2
  ),
  prior_EX_tau_sd_comp = matrix(
    c(1, 1),
    nrow = num_comp,
    ncol = 2
  ),
  prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
  prior_tau_dist = 0,
  prior_PD = FALSE
)
## Recover user set sampling defaults
options(.user_mc_options)

}
\references{
Neuenschwander, B., Branson, M., & Gsponer, T. (2008). 
 Critical aspects of the Bayesian approach to phase I cancer trials. 
 \emph{Statistics in medicine}, 27(13), 2420-2439.
}
