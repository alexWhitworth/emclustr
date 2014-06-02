
# 01. Helper functions
#----------------------------------
supDist <- function (x, y) return (max (abs (x - y)))

diff_func <- function(data_vec, mean_vec) {
  sum((data_vec - mean_vec)^2)
}

# This function calculates the log of the sum of a vector ie- \eqn{log(\sum_{i=1}n x_i)}
log_plus <- function(xvec) {
  m <- length(xvec)
  x <- log(xvec[1])
  for (j in 2:m) {
    sum_j <- sum(xvec[1:j-1])
    x <- x + log(1 + xvec[j]/sum_j)
  }
  return(x)
}

# This function calculates the log of the sum of a more complicated vector, where
# each element is \eqn{\frac{a}{b}e^c} ie- \eqn{log\left(\sum_{i=1}n \frac{a_i}{b_i}e^{c_i}\right)}.
log_plus2 <- function(a, b, c) {
  if ((length(a) != length(b)) || (length(a) != length(c))) {
    stop("Input equal length vectors.")
  }
  if (!(all(c > 0) || all(c < 0))) {
    stop("All values of c must be either > 0 or < 0.")
  }
  m <- length(a)
  # initilialize log sum
  x <- log(a[1]) - log(b[1]) + c[1]
  # aggregate / loop log sum
  for (j in 2:m) { 
    # build denominator
    b2 <- b[1:j-1]
    for (i in 1:j-1) {
      d1 <- 0
      c2 <- c[1:i]
      if (all(c2 > 0)) {
        c_min <- min(c2[1:j-1])
        c2 <- c2 - c_min 
      } else if (all(c2 < 0)) {
        c_min <- max(c2[1:j-1])
        c2 <- c2 - c_min
      }
      d1 <- d1 + a[i] * prod(b2[-i]) * exp(c2[i])
    }
    den <- b[j] * (d1)
    num <- a[j] * prod(b[1:j-1]) * exp(c[j] - c_min)
    x <- x + log(1 + num / den)
  }
  return(x)
}


# This function generates a random vector using the uniform distribution. 
# The elements of the vector may be generated from different unif(a,b).
# @param x A vector to be generated
# @param xmin A vector of minimum values to be used as the min of the uniform distribution
# for each element of \code{x}.
# @param xmax A vector of maximum values to be used as the max of the uniform distribution
# for each element of \code{x}.
rclust_mean <- function(x, xmin, xmax) {
  p <- length(x)
  for (j in 1:p) {
    x[j] <- runif(1, min= xmin[j], max= xmax[j])
  }
  return(x)
}

#' @title Random multivariate normal cluster generation
#' @description
#' This function generates a random p-dimensional multivariate normal cluster.
#' @param obs The number of observations you wish to generate.
#' @param p The number of dimensions for the multivariate normal cluster to be generated.
#' @param mean A vector of means, one for each dimension of the cluster.
#' @param sd A vector of standard deviations, one for each dimension of the cluster.
#' @return Returns an `obs x p` matrix where each column $p_i$ is an obs-length Gaussian 
#' distributed vector with \eqn{N(\mu= mean[i], \sigma= sd[i])}.
#' @export
#' @examples
#' c1 <- gen_clust(100, 10, mean= c(seq(-8, 10, 2)), sd= rep(1, 10))
gen_clust <- function(obs, p, mean, sd) {
  clust <- matrix(nrow= obs, ncol= p)
  for (i in 1:p) {
    clust[,i] <- rnorm(obs, mean= mean[i], sd= sd[i])
  }
  return(clust)
}
