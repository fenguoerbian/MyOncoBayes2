#' Internal function to simulate from the posterior new parameter draws
#'
#' @keywords internal
pp_data <- function(object, newdata, draws, re.form) {
    if(!missing(re.form))
        stop("ERROR: re.form not yet supported")

    idx_group_term <- object$idx_group_term
    idx_inter_term <- object$idx_inter_term
    has_inter <- object$has_inter
    num_comp <- object$standata$num_comp

    ## setup design matrices which must have intercept and slope

    if(!missing(newdata)) {
        data <- newdata
        f <- object$formula
        orig_mf <- object$model
        tt <- terms(f, data=orig_mf, lhs=1, rhs=seq_len(idx_group_term))
        Terms <- delete.response(tt)
        mf <- model.frame(Terms, data)

        X <- .get_X(f, mf, num_comp, has_inter, idx_inter_term)
        X_comp  <- X$comp
        X_inter <- X$inter
    } else {
        data <- object$data
        X_comp <- object$standata$X_comp
        X_inter <- object$standata$X_inter
    }

    strata_group_fct <- .get_strata_group_fct(object, data)

    num_obs <- dim(X_comp)[2]

    group_idx  <- as.integer(unclass(strata_group_fct$group_fct))
    stratum_idx  <- as.integer(unclass(strata_group_fct$strata_fct))

    if(has_inter) {
        post_rv <- as_draws_rvars(as.array(object$stanfit, pars=c("beta_group", "eta_group")))
    } else {
        post_rv <- as_draws_rvars(as.array(object$stanfit, pars=c("beta_group")))
        post_rv$eta_group <- rvar(0)
    }

    if(!missing(draws))
        post_rv <- resample_draws(post_rv, method="deterministic", ndraws=draws)

    pp <- blrm_logit_grouped_rv(group_idx, X_comp, X_inter, post_rv$beta_group, post_rv$eta_group)
    colnames(pp) <- rownames(data)
    pp
}

#' @keywords internal
.validate_factor <- function(test, expected, name) {
    expected_levels  <- levels(expected)
    if(is.factor(test)) {
        test_levels  <- levels(test)
        unseen_levels <- setdiff(test_levels, expected_levels)
        assert_that(length(unseen_levels) == 0,
                    msg=paste0("Found unkown factor levels in ", name, ": ", paste(unseen_levels, collapse=", ")))
        assert_that(all(expected_levels == test_levels), msg=paste0("Mismatch in factor level defintion of ", name))
        return(test)
    }

    unseen_levels <- setdiff(unique(test), expected_levels)
    assert_that(length(unseen_levels) == 0,
                msg=paste0("Found unkown factor levels in ", name, ": ", paste(unseen_levels, collapse=", ")))
    factor(test, levels=expected_levels)
}



#' Numerically stable summation of log inv logit
#' @keywords internal
log_inv_logit <- function(mat) {
    ##- ifelse(is.finite(mat) & (mat < 0), log1p(exp(mat)) - mat, log1p(exp(-mat)))
    ##idx <- is.finite(mat) & (mat < 0)
    idx <- mat < 0
    mat[idx] <- mat[idx] - log1p(exp(mat[idx]))
    mat[!idx] <- -1*log1p(exp(-mat[!idx]))
    mat
}

## numerically stable version of log1m_exp(x) = log(1-exp(x)) for x < 0
log1m_exp_max0  <- function(x) {
    ## qlogis(x) = logit(exp(x)) = x - log(1-exp(x))
    x - qlogis(x, log.p=TRUE)
}

blrm_logit_grouped_rv <- function(group, X_comp, X_inter, beta, eta) {
    ## respective Stan declarations:
    ## vector[2] beta_group[num_groups,num_comp];
    ## vector[num_inter] eta_group[num_groups];
    ## matrix[num_obs,2] X_comp[num_comp];
    num_comp <- dim(X_comp)[1]
    num_inter <- dim(X_inter)[2]
    num_obs <- length(group)
    abeta <- draws_of(beta)
    aeta <- draws_of(eta)
    ## dropping dimnames speeds up subsetting below
    dimnames(abeta) <- NULL
    dimnames(aeta) <- NULL
    dimnames(X_comp) <- NULL
    dimnames(X_inter) <- NULL
    S <- dim(abeta)[1]
    mu <- matrix(0.0*NA, S, num_obs)
    ##log_p0_nr_comp  <- matrix(0.0*NA, S, num_comp) ## obsolete definition
    for (i in seq_len(num_obs)) { # LW: patched to NOT run if num_obs
        g <- group[i]
        log_p0_nr <- numeric(S)
        for(j in seq_len(num_comp)) {
            ##log_p0_nr_comp[,j] <- log_inv_logit(-1 * tcrossprod(adrop(X_comp[j,i,,drop=FALSE], drop=1), adrop(abeta[,g,j,,drop=FALSE], c(2,3))))
            log_p0_nr <- log_p0_nr + log_inv_logit_fast(adrop(abeta[,g,j,,drop=FALSE], c(2,3)) %*% (-1*X_comp[j,i,,drop=FALSE]) )
        }
        if(num_inter > 0) {
            ##mu[,i] <- log1m_exp_max0(log_p0_nr) - log_p0_nr + tcrossprod(X_inter[i,,drop=FALSE], adrop(aeta[,g,,drop=FALSE], 2))
            ##mu[,i] <- -qlogis(log_p0_nr, log.p=TRUE) + adrop(aeta[,g,,drop=FALSE], 2) %*% X_inter[i,,drop=TRUE]
            mu[,i] <- log1m_exp_max0_fast(log_p0_nr) - log_p0_nr + adrop(aeta[,g,,drop=FALSE], 2) %*% X_inter[i,,drop=TRUE]
        } else {
            ##mu[,i] <- log1m_exp_max0(log_p0_nr) - log_p0_nr
            ##mu[,i] <- -qlogis(log_p0_nr, log.p=TRUE)
            mu[,i] <- log1m_exp_max0_fast(log_p0_nr) - log_p0_nr
        }
    }
    mu
}

pp_binomial_trials <- function(object, newdata) {
    data <- object$data
    if(!missing(newdata))
        data <- newdata

    f <- object$formula
    orig_mf <- object$model
    idx_group_term <- object$idx_group_term
    tt <- terms(f, data=orig_mf, lhs=1, rhs=seq_len(idx_group_term))
    mf <- model.frame(tt, data)
    y <- model.response(mf)
    return(rowSums(y))
}


#' extracts from a blrmfit object and a given data-set the group and
#' stratum factor
#'
#' @keywords internal
.get_strata_group_fct  <- function(object, data) {
    f <- object$formula
    orig_mf <- object$model
    idx_group_term <- object$idx_group_term
    tt <- terms(f, data=orig_mf, lhs=1, rhs=seq_len(idx_group_term))
    Terms <- delete.response(tt)
    mf <- model.frame(Terms, data)

    group_index_term <- model.part(f, data = mf, rhs = idx_group_term)
    if(ncol(group_index_term) == 2) {
        idx_group_index <- 2
        idx_strata_index <- 1
    } else {
        idx_group_index <- 1
        idx_strata_index <- NA
    }

    model_group_fct <- object$group_fct
    group_fct <- model.part(f, data = mf, rhs = idx_group_term)
    group_fct <- group_fct[,idx_group_index]
    group_fct <- .validate_factor(group_fct, model_group_fct, "grouping")

    ## enforce that all strata labels are known and match what has
    ## been defined at model fit
    if(!is.na(idx_strata_index)) {
        model_strata_fct <- object$strata_fct
        strata_fct <- model.part(f, data = mf, rhs = idx_group_term)[,idx_strata_index]
        strata_fct <- .validate_factor(strata_fct, model_strata_fct, "stratum")
        strata_group <- data.frame(strata_fct=strata_fct, group_fct=group_fct)
    } else {
        strata_group <- data.frame(strata_fct=1, group_fct=group_fct)
    }

    strata_group
}
