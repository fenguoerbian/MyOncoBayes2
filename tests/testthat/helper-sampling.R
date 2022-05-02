library(checkmate)
library(assertthat)
library(Formula)
library(abind)
library(rstan)

set_sampling_default  <- function(iter, warmup, chains, cores=getOption("mc.cores", 1), save_warmup=FALSE, backend="rstan") {
    options(OncoBayes2.MC.iter=iter,
            OncoBayes2.MC.warmup=warmup,
            OncoBayes2.MC.chains=chains,
            mc.cores=cores,
            OncoBayes2.MC.save_warmup=save_warmup,
            OncoBayes2.MC.backend=backend)
}

very_fast_sampling <- function() {
    message("Tests running with very fast sampling")
    set_sampling_default(300, 150, 1, 1)
}

fake_sampling <- function() {
    message("Tests running with fake sampling")
    set_sampling_default(4, 2, 1, 1)
}

default_sampling <- function() {
    set_sampling_default(NULL, NULL, NULL, NULL)
}

run_example <- function(example) {
    env <- new.env()
    suppressWarnings(example_model(example, env, silent=TRUE))
    invisible(env)
}



## set up slim sampling in case we are on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    very_fast_sampling()
} else {
    fake_sampling()
}

