
#' @title Mixture modeling of normally distributed univariate data.
#' @description This function uses the EM algorithm to do clustering of k-mixture components
#' where each component is one-dimensional \eqn{N(\mu, \sigma^2)}.
#' @param data An n-length vector. Must not be character.
#' @param nclust The number of clusters.
#' @param itmax The maximum number of iterations allowed. Defaults to 10000.
#' @param tol Tuning parameter for convergence. Defaults to 10^-8.
#' @return A list containing: \code{it} the number of iterations; \code{clust_prop}
#' the estimated mixture proportions; \code{clust_params} the estimated mixture parameters; 
#' \code{mix_est} a vector of the estimated mixture for each data point; \code{log_lik} the 
#' log likelihood of the data; \code{bic} the modeled BIC.
#' @seealso \code{\link{em_clust_mvn}}, \code{\link{em_clust_mvn_miss}}, \code{\link{gen_clust}}
#' @export
#' @examples
#' # generate test data
#' c1 <- rnorm(100, 5, 1); c2 <- rnorm(100, 15, 1); c3 <- rnorm(100, 20, 1)
#' c_tot <- c(c1, c2, c3); rm(c1,c2,c3)
#' # run example
#' norm_ex <- em_clust_norm(c_tot, nclust= 3) 

# 02. clustering with poisson distribution
#----------------------------------
em_clust_norm <- function(data, nclust, itmax= 10000, tol= 10^-8) {  
  if (!is.numeric(data)) {
    stop("Please input numeric data.")
  }
  if (itmax < 1 | itmax %% 1 != 0) stop("it_max must be a positive integer.")
  if (nclust < 1 | nclust %% 1 != 0) stop("nclust must be a positive integer.")
  
  # 01. initiate values
  n    <- length(data); it <- 1
  data <- as.vector(data) # coerce
  xmin <- min(data); xmax <- max(data)
  mu_1 <- sig_1 <- lam_1 <- vector(mode= "numeric", length= nclust)
  sig_0 <- c(rep(1, nclust))
  # cluster proportions
  lam_0  <- c(rep(1/nclust, nclust));
  # cluster means
  mu_0 <- vector(mode= "numeric", length= nclust)
  init_m <- quantile(data, seq(0,1, 1/nclust))
  for (i in 1:nclust) {
    mu_0[i] <- runif(1, min= init_m[i], max= init_m[i+1])
  }  
  
  # 02. iterate to convergence:
  repeat{
    # calculate E(clust_m|...) for all data points
    w_mat <- n_mat <- matrix(NA, ncol= nclust, nrow= n)
    d_vec <- vector(mode= "numeric", length= n)
    for (i in 1:nclust) {
      n_mat[, i] <- lam_0[i] * dnorm(data, mean= mu_0[i],, sd= sqrt(sig_0[i]))
    }
    d_vec <- apply(n_mat, 1, sum)
    w_mat <- n_mat / d_vec
    
    # update mixture proportions
    w_dotm <- colSums(w_mat)
    lam_1 <- w_dotm / sum(w_dotm) 
    # update mixture parameters
    mu_1 <- colSums(w_mat * data) / w_dotm
    for (i in 1:nclust) {
      sig_1[i] <- sum(w_mat[,i] * (data - mu_1[i])^2) / w_dotm[i]
    }
    
    # 03. check to exit
    if ((supDist(mu_0, mu_1) < tol & supDist(lam_0, lam_1) < tol) || it == itmax) {
      # calculate final parameters for return
      m_max <- apply(w_mat, 1, which.max)
      # log-lik 
      for (i in 1:nclust) {
        n_mat[,i] <- lam_1[i] * dnorm(data, mean=mu_1[i], sd= sqrt(sig_1[i]), log= FALSE)
      }
      log_lik <- sum(apply(n_mat,1, function(x) {log(sum(x))}))
      bic <- -2 * log_lik + (log(n) * 2 * nclust)
      
      return(list(it= it, clust_prop= lam_1, clust_params= cbind(mean=mu_1, sigma_sq=sig_1), mix_est= m_max,
                  log_lik= log_lik, bic= bic))
    }
    # 04. update for next iteration
    it <- it + 1
    mu_0 <- mu_1; lam_0 <- lam_1
  }
}
