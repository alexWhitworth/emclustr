
library(testthat)
library(emclustr)

context("mvn missing clustering")

test_that("can catch errors", {
  expect_error(em_clust_mvn_miss(data= letters[1:20], nclust= 1))
  expect_error(em_clust_mvn_miss(data= rexp(100, 1), nclust= 0))
  expect_error(em_clust_mvn_miss(data= rexp(100, 1), nclust= 0.5))
  expect_error(em_clust_mvn_miss(data= letters[1:20], nclust= 1, itmax= 0))
  expect_error(em_clust_mvn_miss(data= letters[1:20], nclust= 1, itmax= 0.5))
  expect_error(em_clust_mvn_miss(data= letters[1:20], nclust= 1, itmax= 10.5))
})

test_that("... returns approrpiate output", {
  ## set up tests
  c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
  c2 <- gen_clust(100, 10, mean= rep(0, 10), sd= rep(2, 10))
  c3 <- gen_clust(100, 10, mean= rep(10, 10), sd= rep(1, 10))
  c_tot <- rbind(c1,c2,c3); rm(c1,c2,c3)
  c_tot <- apply(c_tot, 2, function(x) {
    samp <- sample(1:length(x), floor(length(x) * .2), replace=FALSE)
    x[samp] <- NA
    return(x)
  })
  
  nclust <- 3
  test <- em_clust_mvn_miss(c_tot, nclust= nclust, itmax= 5)
  
  ## test
  expect_true(test$it <= 10)
  expect_equal(sum(test$clust_prop), 1)
  expect_equal(length(test$clust_prop), nclust)
  expect_equal(length(test$clust_params), nclust)
  expect_equal(length(test$mix_est), nrow(c_tot))
  expect_true(is.numeric(test$bic))
  expect_true(is.numeric(test$pseudo_log_lik))
  expect_true(is.numeric(test$mix_est))
  expect_true(is.numeric(test$clust_prop))
  expect_true(is.numeric(test$it))
  expect_true(is.list(test$clust_params))
  expect_true(is.numeric(test$clust_params[[1]]$mu))
  expect_true(is.numeric(test$clust_params[[1]]$sigma))
})

test_that("... is converging", {
  ## set up tests
  c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
  c2 <- gen_clust(100, 10, mean= rep(0, 10), sd= rep(2, 10))
  c3 <- gen_clust(100, 10, mean= rep(10, 10), sd= rep(1, 10))
  c_tot <- rbind(c1,c2,c3); rm(c1,c2,c3)
  c_tot <- apply(c_tot, 2, function(x) {
    samp <- sample(1:length(x), floor(length(x) * .2), replace=FALSE)
    x[samp] <- NA
    return(x)
  })
  
  nclust <- 3
  
  iter1 <- em_clust_mvn_miss(c_tot, nclust= nclust, itmax= 1)
  iter2 <- em_clust_mvn_miss(c_tot, nclust= nclust, itmax= 2)
  iter3 <- em_clust_mvn_miss(c_tot, nclust= nclust, itmax= 3)
  
  ## test
  expect_lte(iter1$pseudo_log_lik, iter2$pseudo_log_lik)
  expect_lte(iter2$pseudo_log_lik, iter3$pseudo_log_lik)
})

test_that("... does converge before itmax", {
  ## set up tests
  c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
  c2 <- gen_clust(100, 10, mean= rep(0, 10), sd= rep(2, 10))
  c3 <- gen_clust(100, 10, mean= rep(10, 10), sd= rep(1, 10))
  c_tot <- rbind(c1,c2,c3); rm(c1,c2,c3)
  c_tot <- apply(c_tot, 2, function(x) {
    samp <- sample(1:length(x), floor(length(x) * .2), replace=FALSE)
    x[samp] <- NA
    return(x)
  })
  itmax <- 10000
  
  test <- em_clust_mvn_miss(c_tot, nclust= 3, itmax= itmax)
  
  expect_lt(test$it, itmax)
})