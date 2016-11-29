samplefun = function(K, q){
  sample.int(K+1, 1, prob = q) - 1
}
samplefun2 = function(prior) {
  prior_mean = prior[[1]]
  prior_kappa = prior[[2]]
  prior_nu = prior[[3]]
  prior_sigma = prior[[4]]
  S0 = MCMCpack::rwish(prior_nu, prior_sigma)
  m0 = mvtnorm::rmvnorm(1, prior_mean, 1 / prior_kappa * S0)
  list(m = m0, kappa = prior_kappa, nu = prior_nu, sigma = S0)
}
