
library(testthat)
library(emclustr)

context("poisson clustering")

test_that("can catch errors", {
  expect_error(em_clust_pois(data= letters[1:20], nclust= 1))
  expect_error(em_clust_pois(data= rexp(100, 1), nclust= 0))
  expect_error(em_clust_pois(data= rexp(100, 1), nclust= 0.5))
  expect_error(em_clust_pois(data= letters[1:20], nclust= 1, itmax= 0))
  expect_error(em_clust_pois(data= letters[1:20], nclust= 1, itmax= 0.5))
  expect_error(em_clust_pois(data= letters[1:20], nclust= 1, itmax= 10.5))
})

test_that("... returns approrpiate output", {
  ## set up tests
  c1 <- rpois(100, 3); c2 <- rpois(100, 20); c3 <- rpois(100, 100)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  nclust <- 3
  test <- em_clust_pois(c_tot, nclust= nclust, itmax= 5)
  
  ## test
  expect_true(test$it <= 10)
  expect_equal(sum(test$clust_prop), 1)
  expect_equal(length(test$clust_prop), nclust)
  expect_equal(length(test$clust_params), nclust)
  expect_equal(length(test$mix_est), length(c_tot))
  expect_true(is.numeric(test$bic))
  expect_true(is.numeric(test$log_lik))
  expect_true(is.numeric(test$mix_est))
  expect_true(is.numeric(test$clust_prop))
  expect_true(is.numeric(test$it))
  expect_true(is.numeric(test$clust_params))
})

test_that("... is converging", {
  ## set up tests
  c1 <- rpois(100, 3); c2 <- rpois(100, 20); c3 <- rpois(100, 100)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  nclust <- 3
  
  iter2 <- em_clust_pois(c_tot, nclust= nclust, itmax= 2)
  iter4 <- em_clust_pois(c_tot, nclust= nclust, itmax= 4)
  
  ## test
  expect_lte(iter2$log_lik, iter4$log_lik)
})

test_that("... does converge before itmax", {
  ## set up tests
  c1 <- rpois(100, 3); c2 <- rpois(100, 20); c3 <- rpois(100, 100)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  itmax <- 10000
  test <- em_clust_pois(c_tot, nclust= 3, itmax= itmax)
  
  expect_lt(test$it, itmax)
})