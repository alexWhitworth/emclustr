

# test em_clust_exp
c1 <- rexp(100, 1); c2 <- rexp(100, 50); c3 <- rexp(100, 100)
c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
test_exp <- em_clust_exp(c_tot, nclust= 3) 

# test em_clust_norm
c1 <- rnorm(100, 5, 1); c2 <- rnorm(100, 15, 1); c3 <- rnorm(100, 20, 1)
c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
text_norm <- em_clust_norm(c_tot, nclust= 3)
 
# test em_clust_pois
c1 <- rpois(100, 3); c2 <- rpois(100, 20); c3 <- rpois(100, 100)
c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
test_pois <- em_clust_pois(c_tot, nclust= 3)

# test em_clust_mvn
c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
c2 <- gen_clust(100, 10, mean= rep(0, 10), sd= rep(2, 10))
c3 <- gen_clust(100, 10, mean= rep(10, 10), sd= rep(1, 10))
c_tot <- rbind(c1,c2,c3); rm(c1,c2,c3)
test_mvn <- em_clust_mvn(c_tot, nclust= 3)

# test em_clust_mvn_miss
c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
c2 <- gen_clust(100, 10, mean= rep(0, 10), sd= rep(2, 10))
c3 <- gen_clust(100, 10, mean= rep(10, 10), sd= rep(1, 10))
c_tot <- rbind(c1,c2,c3); rm(c1,c2,c3)
c_tot <- apply(c_tot, 2, function(x) {
  samp <- sample(1:length(x), floor(length(x) * .2), replace=FALSE)
  x[samp] <- NA
  return(x)
})

test_mvn_miss <- em_clust_mvn_miss(c_tot, nclust= 3)
