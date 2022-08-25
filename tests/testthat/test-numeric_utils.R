context("numeric utils")

test_that("log_inv_logit_fast evaluates correctly", {
    l <- seq(-20,20,length=41)
    ob2 <- OncoBayes2:::log_inv_logit_fast(l)

    bfam <- binomial()
    inv_logit <- bfam$linkinv

    expect_equal(ob2, log(inv_logit(l)), tolerance=1E-8)
})

test_that("log1m_exp_max0_fast evaluates correctly", {
    l <- seq(-40,1,length=41)
    ob2 <- OncoBayes2:::log1m_exp_max0_fast(l)

    suppressWarnings(ref <- log1p(-exp(l)))

    expect_equal(ob2, ref, tolerance=1E-8)
})

