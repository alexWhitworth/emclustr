

library(testthat)
library(emclustr)

context("exponential clustering")

test_that("can catch errors", {
  expect_error(em_clust_exp(data= letters[1:20], nclust= 1))
  expect_error(em_clust_exp(data= rexp(100, 1), nclust= 0))
  expect_error(em_clust_exp(data= rexp(100, 1), nclust= 0.5))
  expect_error(em_clust_exp(data= letters[1:20], nclust= 1, itmax= 0))
  expect_error(em_clust_exp(data= letters[1:20], nclust= 1, itmax= 0.5))
  expect_error(em_clust_exp(data= letters[1:20], nclust= 1, itmax= 10.5))
})

test_that("em_clust_exp returns approrpiate output", {
  ## set up tests
  c1 <- rexp(100, 1); c2 <- rexp(100, 50); c3 <- rexp(100, 100)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  nclust <- 3
  test_exp <- em_clust_exp(c_tot, nclust= nclust, itmax= 10)
  
  ## test
  expect_true(test_exp$it <= 10)
  expect_equal(sum(test_exp$clust_prop), 1)
  expect_equal(length(test_exp$clust_prop), nclust)
  expect_equal(length(test_exp$clust_params), nclust)
  expect_equal(length(test_exp$mix_est), length(c_tot))
  expect_true(is.numeric(test_exp$bic))
  expect_true(is.numeric(test_exp$log_lik))
  expect_true(is.numeric(test_exp$mix_est))
  expect_true(is.numeric(test_exp$clust_prop))
  expect_true(is.numeric(test_exp$it))
})

test_that("em_clust_exp is converging", {
  ## set up tests
  c1 <- rexp(100, 1); c2 <- rexp(100, 50); c3 <- rexp(100, 100)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  nclust <- 3
  iter2 <- em_clust_exp(c_tot, nclust= nclust, itmax= 2)
  iter4 <- em_clust_exp(c_tot, nclust= nclust, itmax= 4)
  iter6 <- em_clust_exp(c_tot, nclust= nclust, itmax= 6)
  iter8 <- em_clust_exp(c_tot, nclust= nclust, itmax= 8)
  
  ## test
  expect_lt(iter2$log_lik, iter4$log_lik)
  expect_lt(iter4$log_lik, iter6$log_lik)
  expect_lt(iter6$log_lik, iter8$log_lik)
})

test_that("... does converge before itmax", {
  ## set up tests
  c1 <- rexp(100, 1); c2 <- rexp(100, 50); c3 <- rexp(100, 100)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  itmax <- 10000
  test_exp <- em_clust_exp(c_tot, nclust= 3, itmax= itmax)
  
  expect_lt(test_exp$it, itmax)
})