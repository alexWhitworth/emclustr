
#' @title Mixture modeling of poisson distributed data.
#' @description This function uses the EM algorithm to do clustering of k-mixture components
#' where each component is \eqn{pois(\lambda)}.
#' @param data An n-length vector. Must not be character.
#' @param nclust The number of clusters.
#' @param itmax The maximum number of iterations allowed. Defaults to 10000.
#' @param tol Tuning parameter for convergence. Defaults to 10^-6.
#' @return A list containing: \code{it} the number of iterations; \code{clust_prop}
#' the estimated mixture proportions; \code{clust_params} the estimated mixture parameters; 
#' \code{mix_est} a vector of the estimated mixture for each data point; \code{log_lik} the 
#' log likelihood of the data; \code{bic} the modeled BIC.
#' @export
#' @examples
#' \dontshow{c1 <- rpois(100, 3); c2 <- rpois(100, 20); c3 <- rpois(100, 100)
#' c_tot <- c(c1, c2, c3); rm(c1,c2,c3)}
#' pois <- em_clust_pois(c_tot, nclust= 3) 

em_clust_pois <- function(data, nclust, itmax= 10000, tol= 10^-6) {  
  if (typeof(data) == "character") {
    stop("Please input numeric data.")
  }
  
  # 01. initiate values
  n      <- length(data); it <- 1
  data <- as.vector(data) # coerce
  xmin   <- min(data); xmax <- max(data)
  mu_1 <- lam_1 <- vector(mode= "numeric", length= nclust)
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
      n_mat[,i] <- lam_0[i] * dpois(data, mu_0[i])
    }
    d_vec <- apply(n_mat, 1, sum)
    w_mat <- n_mat / d_vec
    
    # update mixture proportions
    w_dotm <- colSums(w_mat)
    lam_1 <- w_dotm / sum(w_dotm) 
    # update mixture parameters
    mu_1 <- colSums(w_mat * data) / w_dotm
    
    # 03. check to exit
    if ((supDist(mu_0, mu_1) < tol & supDist(lam_0, lam_1) < tol) || it == itmax) {
      # calculate final parameters for return
      m_max <- apply(w_mat, 1, which.max)
      # log-lik 
      for (i in 1:nclust) {
        n_mat[,i] <- lam_1[i] * dpois(data, mu_1[i], log= FALSE)
      }
      log_lik <- sum(apply(n_mat,1, function(x) {log(sum(x))}))
      bic <- -2 * log_lik + log(n) * nclust
      
      return(list(it= it, clust_prop= lam_1, clust_params= mu_1, mix_est= m_max,
                  log_lik= log_lik, bic= bic))
    }
    # 04. update for next iteration
    it <- it + 1
    mu_0 <- mu_1; lam_0 <- lam_1
  }
}
