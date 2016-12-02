library(SamplerCompare)
# set.seed(1)

## sample a random covariance matrix
sigma = matrix(c(1, .3, .3, 1), 2, 2)
S0 = MCMCpack::rwish(2, sigma)

## compute the cholesky decomposition
c0 = chol(S0)
testthat::expect_equal(t(c0) %*% c0, S0)

## now update with a data vector
x = rnorm(2)
c1 = chud(c0, x)
S1 = S0 + x %*% t(x)
testthat::expect_equal(t(c1) %*% c1, S1)

## now downdate again
c2 = chdd(c1, x)

## in some cases, c2[l,l] will be negative, 
## in these cases multiplit c2[l,] by -1:
if (c2[1,1] < 0)
  c2[1,] = c2[1,] * -1
if (c2[2,2] < 0)
  c2[2,] = c2[2,] * -1

S2 = S1 - x %*% t(x)
testthat::expect_equal(c0, c2)
testthat::expect_equal(S0, S2)
testthat::expect_equal(t(c2) %*% c2, S2)
