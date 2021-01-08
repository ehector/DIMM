posdef.matrix <- function(n, ev = runif(n, 0, 10), seed){
  set.seed(seed)
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

#' @name dataset.normal.X1
#' @title simulate data
#' @description Simulate data for illustration of DIMM
#' @param family "normal" for continuous response
#' @param N sample size
#' @param M dimension of response
#' @param theta vector of true parameter values 
#' @param response_indicator a vector of integers from 1 up to M, indicating which block the response belongs to
#' @param seed set random seed for data generation
#' @keywords data simulator
#' @details Simulates an artificial dataset with one covariate and response with block AR-1 covariance structure
#' @return returns a simulated dataset with subject id, response, and one covariate X1
#' @export
#' 
dataset.normal.X1 <- function(family, N, M, theta, response_indicator, seed){
  set.seed(seed)
  sd <- 2.0
  r <- 0.5
  A <- sd^2 * r^abs(outer(1:max(table(response_indicator)), 1:max(table(response_indicator)) , "-"))
  S <- posdef.matrix(length(unique(response_indicator)), seed=seed)
  
  if (length(unique(table(response_indicator)))==1){
    Sigma <- kronecker(S,A)
  } else{
    Sigma <- kronecker(S,A)
    for(i in 1:length(unique(response_indicator))){
      if(table(response_indicator)[i] < max(table(response_indicator))){
        start <- sum(table(response_indicator)[1 : i])+1
        end <- start + max(table(response_indicator))-table(response_indicator)[i]-1
        Sigma <- Sigma[-c(start:end), -c(start:end)]
      }
    }
  }
  
  ## Generate N replicates of an M-dimensional Normal vector with mean beta_0+beta_1*x
  X <- rnorm(N,0,1)
  epsilon <- MASS::mvrnorm(N, rep(0,M) , Sigma, tol = 1e-6)
  Y <- t(epsilon + matrix(theta[1]+theta[2]*X, N, M))
  data_short <- cbind(id=seq(1,N), X1=X, response=t(Y))
  colnames(data_short)[3:(M+2)] <- seq(1,M,1)
  data <- reshape(as.data.frame(data_short), direction="long", varying=c(3:(M+2)), v.names="response")
  data <- data[order(data$id),-which(names(data)=="time")]
  return(data)
}
