#' @title
#' Clustering for Multivariate Normal data via EM
#' @description
#' This function uses the EM algorithm to do clustering in P-dimensions.
#' It assumes all clusters are spherically \eqn{N(\mu_m, \Sigma_m I)}.
#' @param data An `n x p` data matrix or data frame.
#' @param nclust The number of clusters.
#' @param itmax The maximum number of iterations allowed. Defaults to 10000.
#' @param tol Tuning parameter for convergence. Defaults to 10^-6.
#' @return A list with four elements \code{it} - the number of iterations, \code{clust_prop} 
#' the estimated mixture proportions, \code{clust_params} a list of the estimated mixture parameters, 
#' and \code{mix_est} a vector of the estimated mixture for each data point.
#' @export

# 02. Update 201C HW2 algorithm to P-dim
#----------------------------------
em_clust_mvn <- function(data, nclust, itmax= 10000, tol= 10^-6) {
  # 01. initiate values
  data <- data.matrix(data) # coerce
  n      <- nrow(data); p <- ncol(data); it <- 1
  mu_sq <- lam_1  <- cdist <- vector(mode= "numeric", length= nclust)
  # cluster proportions
  lam_0  <- c(rep(1/nclust, nclust))
  # cluster variances
  sig_1  <- sig_0 <- c(rep(1, nclust))
  # cluster means
  mu1 <- mu0 <- matrix(NA, nrow= nclust, ncol= p)
  init_m <- apply(data, 2, quantile, seq(0, 1, 1/nclust))
  for (i in 1:nclust) {
    mu0[i,] <- rclust_mean(vector(mode= "numeric", length= p), 
                           xmin= init_m[i,], xmax= init_m[i+1,])
  }
  # initialize values needed for EM algorithm
  data_ssq <- apply(data^2, 1, sum) # internal- cluster var
  exp_mat1 <- n_mat <- matrix(NA, ncol= nclust, nrow= n)
  internal_vec1 <- c(rep(NA, nclust))
  d_vec <- vector(mode= "numeric", length= n)
  
  # 02. iterate to convergence:
  repeat{
    # (a) E-Step
    # calculate E(clust_m|...) for all data points -- ie w_mat
    for (i in 1:nclust) {
      exp_mat1[,i] <- apply(data, 1, diff_func, mean_vec= mu0[i,]) / (2 * sig_0[i])
      n_mat[, i] <- log(lam_0[i]) - p/2 * log(sig_0[i]) - exp_mat1[,i]
      }
    for (i in 1:n) {
      d_vec[i] <- log_plus2(lam_0, sig_0^(p/2), -exp_mat1[i,])
    }
    w_mat <- exp(n_mat - d_vec) 
    
    # (b) M-Step
    # update mixture proportions
    w_dotm <- colSums(w_mat)
    lam_1 <- w_dotm / sum(w_dotm)
    # update mean/var for each cluster
    for (i in 1:nclust) { 
      mu1[i,]  <-  colSums(w_mat[, i] * data / w_dotm[i]) # cluster mean-vector
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
    mu0 <- mu1
  }
}