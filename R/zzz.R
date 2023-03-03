.onLoad <- function(libname, pkgname) {
    if (!("methods" %in% .packages())) attachNamespace("methods")
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) loadModule(m, what = TRUE)
}

## define globally cmdstanr model object as being undefined, it will
## be instantiated upon the first use of cmdstanr as backend
pkg_env <- new.env()
pkg_env$.cmdstanr_blrm_exnex_model <- NULL

.onAttach <- function(...) {
    ver <- utils::packageVersion("OncoBayes2")
    packageStartupMessage("This is OncoBayes2 version ", ver, " (released ", format(pkg_create_date, "%F"), ", git-sha ", pkg_sha, ")")
}
