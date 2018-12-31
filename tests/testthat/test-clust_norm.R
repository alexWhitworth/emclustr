
library(testthat)
library(emclustr)

context("univariate normal clustering")

test_that("can catch errors", {
  expect_error(em_clust_norm(data= letters[1:20], nclust= 1))
  expect_error(em_clust_norm(data= rexp(100, 1), nclust= 0))
  expect_error(em_clust_norm(data= rexp(100, 1), nclust= 0.5))
  expect_error(em_clust_norm(data= letters[1:20], nclust= 1, itmax= 0))
  expect_error(em_clust_norm(data= letters[1:20], nclust= 1, itmax= 0.5))
  expect_error(em_clust_norm(data= letters[1:20], nclust= 1, itmax= 10.5))
})

test_that("... returns approrpiate output", {
  ## set up tests
  c1 <- rnorm(100, 5, 1); c2 <- rnorm(100, 15, 1); c3 <- rnorm(100, 20, 1)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  nclust <- 3
  test <- em_clust_norm(c_tot, nclust= nclust, itmax= 5)
  
  ## test
  expect_true(test$it <= 10)
  expect_equal(sum(test$clust_prop), 1)
  expect_equal(length(test$clust_prop), nclust)
  expect_equal(nrow(test$clust_params), nclust)
  expect_equal(length(test$mix_est), length(c_tot))
  expect_true(is.numeric(test$bic))
  expect_true(is.numeric(test$log_lik))
  expect_true(is.numeric(test$mix_est))
  expect_true(is.numeric(test$clust_prop))
  expect_true(is.numeric(test$it))
  expect_true(is.matrix(test$clust_params))
  expect_true(is.numeric(test$clust_params))
})

test_that("... is converging", {
  ## set up tests
  c1 <- rnorm(100, 5, 1); c2 <- rnorm(100, 15, 1); c3 <- rnorm(100, 20, 1)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  nclust <- 3
  
  iter1 <- em_clust_norm(c_tot, nclust= nclust, itmax= 1)
  iter3 <- em_clust_norm(c_tot, nclust= nclust, itmax= 3)
  iter5 <- em_clust_norm(c_tot, nclust= nclust, itmax= 5)
  
  ## test
  expect_lte(iter1$log_lik, iter3$log_lik)
  expect_lte(iter3$log_lik, iter5$log_lik)
})

test_that("... does converge before itmax", {
  ## set up tests
  c1 <- rnorm(100, 5, 1); c2 <- rnorm(100, 15, 1); c3 <- rnorm(100, 20, 1)
  c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
  itmax <- 10000
  test <- em_clust_norm(c_tot, nclust= 3, itmax= itmax)
  
  expect_lt(test$it, itmax)
})