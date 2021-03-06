# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Mahalanobis <- function(x, center, cov) {
    .Call('PDP_Mahalanobis', PACKAGE = 'PDP', x, center, cov)
}

dmvt_arma <- function(x, mean, sigma, nu) {
    .Call('PDP_dmvt_arma', PACKAGE = 'PDP', x, mean, sigma, nu)
}

log_pred <- function(x, mean, kappa, nu, sigma) {
    .Call('PDP_log_pred', PACKAGE = 'PDP', x, mean, kappa, nu, sigma)
}

logsumexp <- function(x) {
    .Call('PDP_logsumexp', PACKAGE = 'PDP', x)
}

one_pass_cpp <- function(z, y, components, prior, P, alpha, lambda, debug) {
    invisible(.Call('PDP_one_pass_cpp', PACKAGE = 'PDP', z, y, components, prior, P, alpha, lambda, debug))
}

qik_cpp <- function(x, k, z, components, Pi, lambda) {
    .Call('PDP_qik_cpp', PACKAGE = 'PDP', x, k, z, components, Pi, lambda)
}

qi0_cpp <- function(x, prior, alpha) {
    .Call('PDP_qi0_cpp', PACKAGE = 'PDP', x, prior, alpha)
}

get_assigment_prob <- function(x, z, components, prior, Pi, lambda, alpha) {
    .Call('PDP_get_assigment_prob', PACKAGE = 'PDP', x, z, components, prior, Pi, lambda, alpha)
}

update <- function(components, k, x, debug) {
    invisible(.Call('PDP_update', PACKAGE = 'PDP', components, k, x, debug))
}

new_component <- function(components, k, x, debug) {
    .Call('PDP_new_component', PACKAGE = 'PDP', components, k, x, debug)
}

