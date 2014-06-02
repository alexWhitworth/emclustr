
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
#' @return A list with four elements: \code{it} - the number of iterations, \code{clust_prop} 
#' the estimated mixture proportions, \code{clust_params} a list of the estimated mixture parameters, 
#' and \code{mix_est} a vector of the estimated mixture for each data point.
#' @export

# 02. Update 201C HW2 algorithm to P-dim
#----------------------------------
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
      exp_mat1[,i] <- apply(data2, 1, diff_func, mean_vec= mu0[i,]) / (2 * sig_0[i])
      n_mat[, i] <- log(lam_0[i]) - p/2 * log(sig_0[i]) - exp_mat1[,i]
    }
    for (i in 1:n) {
      d_vec[i] <- log_plus2(lam_0, sig_0^(p/2), -exp_mat1[i,])
    }
    w_mat <- exp(n_mat - d_vec) 
    
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
      }
      m_max <- apply(w_mat, 1, which.max)
      return(list(it= it, clust_prop= lam_1, clust_params= out, mix_est= m_max))
    }
    # 04. update for next iteration
    it <- it + 1
    sig_0 <- sig_1; lam_0 <- lam_1
    w_mat0 <- w_mat
    mu0 <- mu1
  }
}
