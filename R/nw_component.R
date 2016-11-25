library(purrr)
library(MCMCpack)
# library(mvtnorm)
library(mvnfast)
## x = y[i, , drop = FALSE]

## use a closure to represent and update the parameters of the component
## distributions
normal_component <- function(prior_mean, prior_n_mean, prior_n_precision, prior_precision) {
    prior <- list(m = prior_mean, kappa = prior_n_mean,
                  nu = prior_n_precision, precision = prior_precision)
    m <- prior_mean
    precision <- prior_precision
    n_mean <- prior_n_mean
    n_precision <- prior_n_precision
    obj <- list(
        ## update parameters of distribution
        update = function(x) {
            ## cache some values
            n <- nrow(x); d = ncol(x)
            prior_mean <- m
            prior_precision <- precision
            prior_n_mean <- n_mean
            prior_n_precision <- n_precision
            x_mean <- colMeans(x)

            ## update parameters with data
            m <<-
                prior_n_mean /  (prior_n_mean + n) * prior_mean +
                n / (prior_n_mean + n) * x_mean
            n_mean <<-
                prior_n_mean + n
            n_precision <<-
                prior_n_precision + n
            precision <<-
                prior_precision +
                crossprod(x-x_mean, x-x_mean) +
                prior_n_mean * n / n_mean *
                tcrossprod(matrix(x_mean - prior_mean), matrix(x_mean - prior_mean, byrow = TRUE))
        },
        ## evaluate likelihood of observations x
        ## using the current parameter values
        loglik = function(x) {
          sigma = solve(precision)
          return(mvtnorm::dmvnorm(x, m,    sigma, log = TRUE))
        },
        ## return parameters
        params = function(what = NULL) {
            ## return(dnorm(x, m, S))

            if (!is.null(what))
                return(list(m = m, precision = precision,
                            n_mean = n_mean, n_precision = n_precision)[[what]])
            list(m = m, precision = precision,
                 n_mean = n_mean, n_precision = n_precision)
        },
        prior = prior,
        ## evaluate posterior predictive for x
        log_post_pred = function(x) {
            n <- nrow(x); d = ncol(x)
            sum(mvnfast::dmvt(x,
                              m,
                              (n_mean + n) / (n_mean * (n_precision - d + n)) * precision,
                              n_precision - d + n, log = TRUE))
            },
        ## evaluate posterior predictive, but with prior parameters
        ## not updated with any data
        log_prior_pred = function(x) {
            n <- nrow(x); d = ncol(x)
            sum(mvnfast::dmvt(x,
                              prior$m,
                              (prior$kappa + n) / (prior$kappa * (prior$nu - d + n)) * prior$precision,
                              prior$nu - d + n, log = TRUE))
        }
    )
    class(obj) <- c(class(obj), "component")
    obj
}

params = function(x, ...) UseMethod("params")
params.component = function(obj, what = NULL) {
    obj$params(what)
}

####################################
##============ tests =============##
####################################
library(testthat)
test_that("prior_predictive works", {
  y = mvtnorm::rmvnorm(400, mean = c(0, 0), sigma = matrix(c(1, 0, 0, 1), 2, 2))
  d = 2
  m_n = c(2, 2)
  S_n = matrix(c(1, .3, .3, 1), 2, 2)
  nu_n = 2; kappa_n = 2
  nv1 = normal_component(m_n, kappa_n, nu_n, S_n)
  expect_equal(
    sapply(1:nrow(y), function(i) nv1$log_prior_pred(y[i, , drop = FALSE])),
    apply(y, 1, function(xi) mvtnorm::dmvt(xi, m_n, (kappa_n + 1) / (kappa_n * (nu_n - d + 1)) * S_n, nu_n - d + 1))
  )
})

test_that("nw distribution operations work", {
  prior_mean <- c(2,2)
  prior_precision = matrix(c(1, 0, 0, 1), 2, 2)
  prior_n_mean = prior_n_precision = 2
  nv1 = normal_component(c(2, 2), 2, 2, prior_precision)
  y = mvtnorm::rmvnorm(400, mean = c(0, 0), sigma = matrix(c(1, 0, 0, 1), 2, 2))
  nv1$update(y)
  nv1$log_prior_pred(y)
  nv1$log_prior_pred(y)
  nv1$log_post_pred(y)
  nv1$params("precision")
  nv1$update(y)
})

test_that("marginals", {
prior_mean <- c(2,2)
prior_precision = matrix(c(1, 0, 0, 1), 2, 2)
prior_n_mean = prior_n_precision = 2
x = seq(-5, 5, len = 100)
# plot(mvtnorm::dmvnorm(x,  prior_precision))
nv1 = normal_component(c(2, 2), 2, 2, prior_precision)
y = mvtnorm::rmvnorm(400, mean = c(0, 0), sigma = matrix(c(1, 0, 0, 1), 2, 2))
})


