
#' @title
#' Clustering for Multivariate Normal with missing data via EM
#' @description
#' This function uses the EM algorithm to do clustering in P-dimensions.
#' It assumes all clusters are spherically \eqn{N(\mu_m, \Sigma_m I)} and allows for
#' missing data. Observations may have missing data elements but entire rows may
#' not be missing.
#' @param data An `n x p` data matrix or data frame.
#' @param nclust The number of clusters.
#' @param itmax The maximum number of iterations allowed. Defaults to 10000.
#' @param tol Tuning parameter for convergence. Defaults to 10^-6.
#' @return A list containing: \code{it} the number of iterations; \code{clust_prop}
#' the estimated mixture proportions; \code{clust_params} the estimated mixture parameters; 
#' \code{mix_est} a vector of the estimated mixture for each data point; \code{pseudo_log_lik} 
#' the pseudo log likelihood of the data; \code{bic} the modeled BIC.
#' and \code{mix_est} a vector of the estimated mixture for each data point.
#' @seealso \code{\link{em_clust_mvn}}, \code{\link{em_clust_norm}}, \code{\link{gen_clust}}
#' @export
#' @examples
#' # generate test data
#' c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
#' c2 <- gen_clust(100, 10, mean= rep(0, 10), sd= rep(2, 10))
#' c3 <- gen_clust(100, 10, mean= rep(10, 10), sd= rep(1, 10))
#' c_tot <- rbind(c1,c2,c3); rm(c1,c2,c3)
#' c_tot <- apply(c_tot, 2, function(x) {
#'   samp <- sample(1:length(x), floor(length(x) * .2), replace=FALSE)
#'   x[samp] <- NA
#'   return(x)
#' })
#' # run example
#' mvn_miss <- em_clust_mvn_miss(c_tot, nclust= 3)

em_clust_mvn_miss <- function(data, nclust, itmax= 10000, tol= 10^-6) {
  # 01. initiate values
  data <- data.matrix(data) # coerce
  n      <- nrow(data); p <- ncol(data); it <- 1
  mu_sq <- lam_1  <- cdist <- vector(mode= "numeric", length= nclust)
  data2 <- data
  # cluster proportions
  lam_0  <- c(rep(1/nclust, nclust))
  # cluster variances
  sig_1  <- sig_0 <- c(rep(1, nclust))
  # cluster means
  init_m <- apply(data, 2, quantile, seq(0, 1, 1/nclust), na.rm= TRUE)
  mu1 <- mu0 <- matrix(NA, nrow= nclust, ncol= p)
  for (i in 1:nclust) {
    mu0[i,] <- rclust_mean(vector(mode= "numeric", length= p), 
                           xmin= init_m[i,], xmax= init_m[i+1,])
  }
  # initialize values needed for EM algorithm
  exp_mat1 <- n_mat <- matrix(NA, ncol= nclust, nrow= n)
  internal_vec1 <- c(rep(NA, nclust))
  d_vec <- vector(mode= "numeric", length= n)
  w_mat0 <- matrix(data= 1/nclust, nrow= n, ncol= nclust)
    
  # 02. iterate to convergence:
  repeat{
    # Step 1 --- EM to impute missing data
    # (1.a) M-Step
    miss_ind <- which(apply(data, 1, anyNA))
    for (i in miss_ind) {
      r_work <- data[i,]; r_ind <- which(is.na(r_work))
      for (k in r_ind) {
        data2[i,k] <- sum(w_mat0[i,] * mu0[,k]) # impute missing values via EM
      }
    }
    data_ssq <- apply(data2^2, 1, sum) # internal- cluster var
    
    # Step 2 --- EM algorithm on imputed data
    # (2.a) E-Step
    # calculate E(clust_m|...) for all data points -- ie w_mat
    for (i in 1:nclust) {
      n_mat[,i] <- lam_0[i] * dmvnorm(data2, mean= mu0[i,], sigma= diag(rep(sig_0[i], p)))
    }
    d_vec <- apply(n_mat, 1, sum)
    w_mat <- n_mat / d_vec
    
    # (2.b) M-Step
    # update mixture proportions
    w_dotm <- colSums(w_mat)
    lam_1 <- w_dotm / sum(w_dotm)
    # update mean/var for each cluster
    for (i in 1:nclust) { 
      mu1[i,]  <-  colSums(w_mat[, i] * data2 / w_dotm[i]) # cluster mean-vector
      mu_sq[i] <- sum(mu1[i,]^2) / p # internal- cluster var
      sig_1[i] <- sum(w_mat[,i] * data_ssq / (p * w_dotm[i])) - mu_sq[i] # cluster var-vector
    }
    
    # 03. check to exit
    if ((supDist(mu0, mu1) < tol & supDist(sig_0, sig_1) < tol &
           supDist(lam_0, lam_1) < tol) || it == itmax) {
      # build list of distributional parameters for return
      out <- list()
      for (i in 1:nclust) { 
        assign(paste0("N_", i), list(mu= mu1[i,], Sigma= sig_1[i]))
        out[[i]] <- get(paste0("N_", i)) 
        n_mat[,i] <- lam_1[i] * dmvnorm(data2, mean= mu1[i,], sigma= diag(rep(sig_1[i], p)))
      }
      m_max <- apply(w_mat, 1, which.max)
      log_lik <- sum(apply(n_mat,1, function(x) {log(sum(x))}))
      bic <- -2 * log_lik + (log(n) * (nclust + length(mu1)))
      
      return(list(it= it, clust_prop= lam_1, clust_params= out, mix_est= m_max,
                  pseudo_log_lik= log_lik, bic= bic))
    }
    # 04. update for next iteration
    it <- it + 1
    sig_0 <- sig_1; lam_0 <- lam_1
    w_mat0 <- w_mat
    mu0 <- mu1
  }
}