dimm <- function(x, ...) UseMethod("dimm")

print.dimm <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nVariance:\n")
  print(x$vcov)
}

summary.dimm <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.dimm"
  return(res) 
}

print.summary.dimm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

objective.dimm <- function(object, ...){
  if (object$family == "normal") {
    func1 <- logCLnormal_par
    func2 <- eenormalmean_par
    func3 <- eenormalderivmean_par
    if (corstr=="CS")
      func4 <- eenormalCSvar_par
    if (corstr=="AR1")
      func4 <- eenormalAR1var_par
    if (corstr=="independence")
      func4 <- eenormalindvar_par 
  }
  obj <- list()
  obj$call <- object$call
  N <- length(object$subject_indicator)
  response <- cbind(matrix(object$response[,-which(colnames(object$response) == object$id)], 
                           length(object$response_indicator), N), object$response_indicator)
  colnames(response) <- c(object$subject_indicator,"response_indicator")
  obj$Q <- as.numeric(combined.estimation.mean(object$coefficients, object$MCLE$vcov, object$V_psi.inv, 
                                               response, object$covariates, object$J, object$K, N, div=1, func2, func3))
  obj$df <- (object$J-1)*(dim(object$covariates)[2]-1)
  obj$p.value <- as.numeric(pchisq(obj$Q, obj$df))
  class(obj) <- "objective.dimm"
  return(obj)
}

print.objective.dimm <- function(x, ...){
  res <- cbind(Q = x$Q,
               df = x$df,
               p.value = x$p.value )
  cat("Call:\n")
  print(x$call)
  cat("Testing fit of model:\n")
  print(res)
}

CS.CL.estimation <- function(par, block_y, block_x, m, n, div, func){
  ## This is the function to be minimized over the set of parameters within each block
  -func(beta=par[1:dim(block_x)[2]], cov=par[dim(block_x)[2]+1]^2*(matrix(par[dim(block_x)[2]+2],m,m)-diag(par[dim(block_x)[2]+2]-1,m,m)), 
        block_y, block_x, m, n)/div
}

AR1.CL.estimation <- function(par, block_y, block_x, m, n, div, func){
  ## This is the function to be minimized over the set of parameters within each block
  -func(beta=par[1:(dim(block_x)[2])], cov=par[dim(block_x)[2]+1]^2*par[dim(block_x)[2]+2]^abs(outer(1:m, 1:m , "-")), 
        block_y, block_x, m, n)/div
}

independence.CL.estimation <- function(par, block_y, block_x, m, n, div, func){
  ## This is the function to be minimized over the set of parameters within each block
  -func(beta=par[1:(dim(block_x)[2])], cov=diag(par[(dim(block_x)[2]+1)],m), block_y, block_x, m, n)/div
}

psi.mean <- function(par, MCLE, response, covariates, J, K, N, func){
  psi <- c()
  beta <- par[1:(dim(covariates)[2]-1)]
  for (k in 1:K){
    for (j in 1:J){
      block_y <- response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)]
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[colnames(response)==k],]
      block_x <- as.matrix(ids_block_x[rep(response[,which(colnames(response)=="response_indicator")], 
                                           length(unique(ids_block_x[,"id"])))==j,
                                       -which(colnames(ids_block_x)=="id")])
      
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      psi <- c(psi, Reduce("+", func(beta, cov=matrix(MCLE[[k]][[j]],m,m), 
                                     block_y, block_x, m=m, n=n))/N)
    }
  }
  return(psi)
}

psi.deriv.mean <- function(par, MCLE, response, covariates, J, K, N, funcderiv){
  psi <- c()
  beta <- par[1:(dim(covariates)[2]-1)]
  for (k in 1:K){
    for (j in 1:J){
      block_y <- response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)]
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[colnames(response)==k],]
      block_x <- as.matrix(ids_block_x[rep(response[,which(colnames(response)=="response_indicator")], 
                                           length(unique(ids_block_x[,"id"])))==j,
                                       -which(colnames(ids_block_x)=="id")])
      
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      psi <- cbind(psi, Reduce("+", funcderiv(cov=matrix(MCLE[[k]][[j]],m,m), 
                                              block_y, block_x, m=m, n=n))/N)
    }
  }
  return(psi)
}

combined.estimation.mean <- function(par, MCLE, W, response, covariates, J, K, N, div, func, funcderiv){
  N * psi.mean(par, MCLE, response, covariates, J, K, N, func) %*% W %*% psi.mean(par, MCLE, response, covariates, J, K, N, func)/div
}

estimation.deriv.mean <- function(par, MCLE, W, response, covariates, J, K, N, div, func, funcderiv){
  2 * N * psi.mean(par, MCLE, response, covariates, J, K, N, func) %*% W %*% 
    t(psi.deriv.mean(par, MCLE, response, covariates, J, K, N, funcderiv))/div
}

ridge.estimate.mean <- function(psi_list, MCLE, MCLE_mean, response, covariates, J, K, N, funcderiv, folds, lam){
  S <- list()
  W <- list()
  lam_list <- vector("numeric", K)
  scale <- matrix(0, dim(covariates)[2]-1, dim(covariates)[2]-1)
  main <- matrix(0, dim(covariates)[2]-1, 1)
  for (k in 1:K){
    S[[k]] <- list()
    for (j in 1:J){
      block_y <- response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)]
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[colnames(response)==k],]
      block_x <- as.matrix(ids_block_x[rep(response[,which(colnames(response)=="response_indicator")], 
                                           length(unique(ids_block_x[,"id"])))==j,
                                       -which(colnames(ids_block_x)=="id")])
      
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      S[[k]][[j]] <- Reduce("+", funcderiv(matrix(MCLE[[k]][[j]],m,m), 
                                           block_y=block_y, block_x=block_x, 
                                           m=m, n=n))/n
    }
    S_k <- do.call(cbind, S[[k]])
    if (length(lam) != 1){
      CV <- vector("numeric", length(lam))
      for(f in 1:length(lam))
        CV[f] <- modified.cholesky.cv(psi_list[[k]], J, dim(covariates)[2]-1, N, lam[f], folds)
      lam_list[k] <- lam[which.min(CV)] 
    } else {
      lam_list[k] <- lam
    }
    W[[k]] <- modified.cholesky(psi_list[[k]], J, dim(covariates)[2]-1, N, lam_list[k])
    for(j in 1:J){
      scale <- scale + n*S_k %*% W[[k]][,((j-1)*(dim(covariates)[2]-1)+1):(j*(dim(covariates)[2]-1))] %*% S[[k]][[j]]
      main <- main + n*S_k %*% W[[k]][,((j-1)*(dim(covariates)[2]-1)+1):(j*(dim(covariates)[2]-1))] %*% S[[k]][[j]] %*%
        MCLE_mean[[k]][[j]]
    }
  }
  return(list(coefficients=solve(scale)%*%main, 
              W=as.matrix(Matrix::bdiag(lapply(W, function(x) matrix(unlist(x), ncol=J*(dim(covariates)[2]-1), byrow=TRUE)))), 
              lam=lam_list))
}

modified.cholesky <- function(psi_list, J, p, N, lam){
  library(MASS)
  psi_mat <- c()
  D <- matrix(0, p*J, p*J)
  Gamma <- matrix(0, p*J, p*J)
  diag(Gamma) <- 1
  for(j in 1:J){
    psi_mat <- rbind(psi_mat, matrix(unlist(psi_list[[j]]), p, N))
  }
  psi_mat <- t(psi_mat)
  for(r in 2:dim(psi_mat)[2]){
    ridge <- lm.ridge(psi_mat[,r] ~ 0 + psi_mat[,1:(r-1)], lambda=lam)
    ## beta <- solve(t(psi_mat[,1:(r-1)]) %*% psi_mat[,1:(r-1)], t(psi_mat[,1:(r-1)]) %*% psi_mat[,r])
    Gamma[r, 1:(r-1)] <- -as.matrix(ridge$coef)
    D[r,r] <- var(psi_mat[,r] - psi_mat[,1:(r-1)]%*%as.matrix(ridge$coef))
  }
  D[1,1] <- var(psi_mat[,1])
  return((t(Gamma)%*%solve(D)%*%Gamma))
}

modified.cholesky.cv <- function(psi_list, J, p, N, lam, folds){
  library(MASS)
  psi_mat <- c()
  for(j in 1:J){
    psi_mat <- rbind(psi_mat, matrix(unlist(psi_list[[j]]), p, N))
  }
  psi_mat <- t(psi_mat)
  s_v <- c()
  CV <- vector("numeric", folds)
  for(f in 1:folds){
    set.seed(f)
    indices <- sample(setdiff(seq(1:dim(psi_mat)[1]), s_v), floor(dim(psi_mat)[1]/folds), replace=FALSE)
    s_v <- c(s_v, indices)
    psi_mat_sub <- psi_mat[setdiff(seq(1:dim(psi_mat)[1]), indices),]
    D <- matrix(0, p*J, p*J)
    Gamma <- matrix(0, p*J, p*J)
    diag(Gamma) <- 1
    for(r in 2:dim(psi_mat_sub)[2]){
      ridge <- lm.ridge(psi_mat_sub[,r] ~ 0 + psi_mat_sub[,1:(r-1)], lambda=lam)
      Gamma[r, 1:(r-1)] <- -as.matrix(ridge$coef)
      D[r,r] <- var(psi_mat_sub[,r] - psi_mat_sub[,1:(r-1)]%*%as.matrix(ridge$coef))
    }
    D[1,1] <- var(psi_mat_sub[,1])
    W <- (t(Gamma)%*%solve(D)%*%Gamma)
    V <- cov(psi_mat_sub)
    CV[f] <- length(indices)*log(det(V)) + sum(apply(psi_mat[indices,],1,function(x) t(x) %*% W %*% x))
  }
  return(sum(CV, na.rm=T)/folds)
}

dimm <- function(formula, data, id=id, response_indicator=NULL, subject_indicator=NULL, family, corstr, method=NULL,
                  lam=NULL, folds=NULL, cluster=NULL){
  cl <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mf[,"(id)"] <- data[,id]
  mt <- attr(mf, "terms")
  response <- cbind(id=mf[,"(id)"], model.response(mf, "numeric"))
  covariates <- cbind(id=mf[,"(id)"], model.matrix(mt, mf))
  
  if (is.null(response_indicator)){
    warning(gettextf("Response indicator is null. Using one response block"), domain = NA)
    response_indicator <- as.vector(rep(1, table(response[,"id"])[1]))
  }
  if (is.null(subject_indicator)){
    warning(gettextf("Subject indicator is null. Using one subject block"), domain=NA)
    subject_indicator <- as.vector(rep(1, length(unique(response[,"id"]))))
  }
  if (is.empty.model(mt)) 
    stop("Model is empty")
  if (family != "normal") 
    warning(gettextf("family = '%s' is not supported. Using 'normal'", family), domain = NA)
  if (corstr != "AR1" & corstr != "CS" & corstr != "independence")
    stop("corstr must be supported type")
  
  response <- cbind(matrix(response[,-which(colnames(response) == "id")], length(response_indicator), length(subject_indicator)), response_indicator)
  
  if (is.null(method))
    method <- "iterative"
  
  if (is.null(cluster)){
    output <- dimm.compute.mean(response, covariates, response_indicator, subject_indicator, family, corstr, method, lam, folds)
  }
  if (!is.null(cluster))
    if (cluster == 1){
      output <- dimm.compute.mean(response, covariates, response_indicator, subject_indicator, family, corstr, method, lam, folds)
    }
  else {
    output <- dimm.compute.mean.parallel(response, covariates, response_indicator, subject_indicator, 
                                            family, corstr, method, lam, folds, cluster)
  }
  names(output$coefficients) <- c(colnames(covariates)[-1])
  colnames(output$vcov) <- names(output$coefficients)
  rownames(output$vcov) <-names(output$coefficients)
  output$response <- cbind(id=mf[,"(id)"], model.response(mf, "numeric"))
  output$covariates <- cbind(id=mf[,"(id)"], model.matrix(mt, mf))
  if (family=="normal")
    output$fitted.values <- cbind(output$response[,which(colnames(output$response)=="id")],
                                  output$response[,-which(colnames(output$response)=="id")] - 
                                    output$covariates[,-which(colnames(output$covariates)=="id")] %*% output$coefficients)
  output$id <- "id"
  
  output <- c(output, list(call=cl, formula=formula))
  class(output) <- "dimm"
  return(output)
}

dimm.compute.mean <- function(response, covariates, response_indicator, subject_indicator, family, corstr, method, lam, folds){
  output <- list()
  library(nlme)
  
  if (family == "normal") {
    func1 <- logCLnormal_par
    func2 <- eenormalmean_par
    func3 <- eenormalderivmean_par
    if (corstr=="CS")
      func4 <- eenormalCSvar_par
    if (corstr=="AR1")
      func4 <- eenormalAR1var_par
    if (corstr=="independence")
      func4 <- eenormalindvar_par 
  }
  
  colnames(response) <- c(subject_indicator,"response_indicator")
  
  J <- length(unique(response_indicator))
  K <- length(unique(subject_indicator))
  N <- length(subject_indicator)
  p <- dim(covariates)[2]-1
  
  V <- list()
  S <- list()
  psi_list <- list()
  MCLE <- list()
  MCLE_mean <- list()
  n_k <- c()
  print("Computing block coefficients.", quote=FALSE)
  for (k in 1:K) {
    S[[k]] <- list()
    psi_list[[k]] <- list()
    MCLE[[k]] <- list()
    MCLE_mean[[k]] <- list()
    for (j in 1:J) {
      block_y <- response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)]
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.matrix(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,
                                       -which(colnames(ids_block_x)=="id")])
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      new_data <- data.frame(cbind(resp=c(block_y), block_x[,-1], newid=sort(rep(seq(1,n),m))))
      
      if (corstr == "CS") {
        block_analysis <- gls(resp ~ . -newid, data=new_data, correlation=corCompSymm(form= ~ 1 | newid))
        MCLE_jk <- block_analysis$coefficients
        MCLE[[k]][[j]] <- getVarCov(block_analysis)
        psi_list[[k]][[j]] <- func2(MCLE_jk[1:dim(block_x)[2]], MCLE[[k]][[j]], block_y=block_y, block_x=block_x, m=m, n=n)
      }
      if (corstr == "AR1") {
        block_analysis <- gls(resp ~ . -newid, data=new_data, correlation=corAR1(form= ~ 1 | newid))
        MCLE_jk <- block_analysis$coefficients
        MCLE[[k]][[j]] <- getVarCov(block_analysis)
        psi_list[[k]][[j]] <- func2(MCLE_jk[1:dim(block_x)[2]], MCLE[[k]][[j]], block_y=block_y, block_x=block_x, m=m, n=n)
      }
      if (corstr == "independence"){
        block_analysis <- lm(resp ~ . -newid, data=new_data)
        MCLE_jk <- block_analysis$coefficients
        MCLE[[k]][[j]] <- diag(summary(block_analysis)$sigma^2, m)
        psi_list[[k]][[j]] <- func2(MCLE_jk[1:dim(block_x)[2]], MCLE[[k]][[j]], block_y=block_y, block_x=block_x, m=m, n=n)
      }
      
      n_k[k] <- n
      MCLE_mean[[k]][[j]] <- MCLE_jk[1:dim(block_x)[2]]
      S[[k]][[j]] <- Reduce("+", func3(matrix(MCLE[[k]][[j]],m,m), block_y=block_y, block_x=block_x, m=m, n=n))/n
    }
    sample_cov <- list()
    for (j in 1:J) {
      cov <- list()
      for (i in 1:J) {
        cov <- cbind(cov, Reduce("+", mapply(function(x,y) x %*% t(y), psi_list[[k]][[j]], psi_list[[k]][[i]], SIMPLIFY=FALSE))/N)
      }
      sample_cov[[j]] <- cov
    }
    V[[k]] <- do.call(rbind,sample_cov) 
  }

  V_psi <- as.matrix(Matrix::bdiag(lapply(V, function(x) matrix(unlist(x), ncol=J*p, byrow=TRUE))))
  W <- solve(V_psi)
  
  print("Computing combined estimate.", quote=FALSE)
  scale <- matrix(0, p, p)
  main <- matrix(0, p, 1)
  for (k in 1:K){
    for(i in 1:J){
      for (j in 1:J){
        s_i <- matrix(unlist(S[[k]]), ncol=J*p, byrow=FALSE)[,((i-1)*p+1):(i*p)]
        s_j <- matrix(unlist(S[[k]]), ncol=J*p, byrow=FALSE)[,((j-1)*p+1):(j*p)]
        W_ij <- W[((k-1)*(J*p)+1):(k*J*p),((k-1)*(J*p)+1):(k*J*p)][((i-1)*p+1):(i*p), ((j-1)*p+1):(j*p)]
        main <- main + n_k[k]*s_i %*% W_ij %*% s_j %*% MCLE_mean[[k]][[j]]
        scale <- scale + n_k[k]*s_i %*% W_ij %*% s_j
      }
    }
  }
  
  if (method=="exact")
    output$coefficients <- as.vector(solve(scale)%*%main)
  if (method=="iterative"){
    init_betas <- colMeans(do.call(rbind,lapply(MCLE_mean,function(x){do.call(rbind,x)})))
    optimization <- optim(par=init_betas, 
                          combined.estimation.mean, gr=estimation.deriv.mean, 
                          MCLE=MCLE, W=W, response=response, 
                          covariates=covariates, J=J, K=K, N=N, func=func2, funcderiv=func3, method="Nelder-Mead", control=list(maxit=500))
    init_betas <- optimization$par
    output$coefficients <- optim(par=init_betas, 
                              combined.estimation.mean, gr=estimation.deriv.mean, 
                              MCLE=MCLE, W=W, response=response, 
                              covariates=covariates, J=J, K=K, N=N, func=func2, funcderiv=func3, method="BFGS")$par
  }
  
  output$vcov <- solve(scale)
  
  if (method=="ridge"){
    if (is.null(folds)) 
     folds <- 1
    if (is.null(lam))
     lam <- 0
    ridge_optim <- ridge.estimate.mean(psi_list, MCLE, MCLE_mean, response, covariates, J, K, N, funcderiv=func3, folds, lam)
    output$coefficients <- as.vector(ridge_optim$coefficients)
    output$W <- ridge_optim$W
    output$lam <- ridge_optim$lam
    sensitivity_psi <- matrix(unlist(S[[k]]), ncol=J*p, byrow=FALSE)
    output$vcov <- solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))%*%
      sensitivity_psi%*%output$W%*%solve(W)%*%output$W%*%t(sensitivity_psi)%*%
      solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))/N
  } 
  
  output$family <- family
  output$response_indicator <- response_indicator
  output$subject_indicator <- subject_indicator
  output$J <- J
  output$K <- K
  output$MCLE <- list()
  output$MCLE$beta <- MCLE_mean
  output$MCLE$vcov <- MCLE
  output$V_psi.inv <- W
  
  return(output)
}

dimm.compute.mean.parallel <- function(response, covariates, response_indicator, subject_indicator, 
                                        family, corstr, method, lam, folds, cluster){
  library(foreach)
  library(doSNOW)
  time1 <- proc.time()
  output <- list()
  
  if (family == "normal") {
    func1 <- logCLnormal_par
    func2 <- eenormalmean_par
    func3 <- eenormalderivmean_par
    if (corstr=="CS")
      func4 <- eenormalCSvar_par
    if (corstr=="AR1")
      func4 <- eenormalAR1var_par
    if (corstr=="independence")
      func4 <- eenormalindvar_par 
  }
  
  colnames(response) <- c(subject_indicator,"response_indicator")
  
  J <- length(unique(response_indicator))
  K <- length(unique(subject_indicator))
  N <- length(subject_indicator)
  p <- dim(covariates)[2]-1
  
  V <- list()
  S <- list()
  psi_list <- list()
  MCLE <- list()
  MCLE_mean <- list()
  n_k <- c()
  psi_list <- list()
  n_k <- c()
  time_list <- list()
  init_betas <- rep(1, p)
  if (corstr=="AR1" | corstr=="CS")
    init_cov_parameters <- c(1,0.5)
  if(corstr=="independence")
    init_cov_parameters <- 1
  time1 <- proc.time()-time1
  
  sock <- parallel::makeCluster(rep("localhost", cluster), type = "SOCK")
  doParallel::registerDoParallel(sock)
  
  print("Computing block coefficients.", quote=FALSE)
  MCLE_psilist_n <- foreach(k=1:K, .packages="nlme") %:% foreach(j=1:J, .packages="nlme") %dopar% {
    time_j <- proc.time()
    block_y <- response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)]
    ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
    block_x <- as.matrix(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,
                                     -which(colnames(ids_block_x)=="id")])
    m <- dim(block_y)[1]
    n <- dim(block_y)[2]
    
    new_data <- data.frame(cbind(resp=c(block_y), block_x[,-1], newid=sort(rep(seq(1,n),m))))
    
    if (corstr == "CS") {
      block_analysis <- gls(resp ~ . -newid, data=new_data, correlation=corCompSymm(form= ~ 1 | newid))
      MCLE_jk <- block_analysis$coefficients
      MCLE_matrix <- getVarCov(block_analysis)
      return(list(MCLE_matrix,
                  MCLE_jk[1:dim(block_x)[2]],
                  func2(MCLE_jk[1:dim(block_x)[2]], MCLE_matrix, block_y=block_y, block_x=block_x, m=m, n=n),
                  Reduce("+", func3(matrix(MCLE_matrix,m,m), block_y=block_y, block_x=block_x, m=m, n=n))/n,
                  n,
                  proc.time()-time_j))
    }
    if (corstr == "AR1") {
      block_analysis <- gls(resp ~ . -newid, data=new_data, correlation=corAR1(form= ~ 1 | newid))
      MCLE_jk <- block_analysis$coefficients
      MCLE_matrix <- getVarCov(block_analysis)
      return(list(MCLE_matrix, 
                  MCLE_jk[1:dim(block_x)[2]],
                  func2(MCLE_jk[1:dim(block_x)[2]], MCLE_matrix, block_y=block_y, block_x=block_x, m=m, n=n),
                  Reduce("+", func3(matrix(MCLE_matrix,m,m), block_y=block_y, block_x=block_x, m=m, n=n))/n,
                  n,
                  proc.time()-time_j))
    }
    if (corstr == "independence") {
      block_analysis <- lm(resp ~ . -newid, data=new_data)
      MCLE_jk <- block_analysis$coefficients
      MCLE_matrix <- diag(summary(block_analysis)$sigma^2, m)
      return(list(MCLE_matrix, 
                  MCLE_jk[1:dim(block_x)[2]],
                  func2(MCLE_jk[1:dim(block_x)[2]], MCLE_matrix, block_y=block_y, block_x=block_x, m=m, n=n),
                  Reduce("+", func3(matrix(MCLE_matrix,m,m), block_y=block_y, block_x=block_x, m=m, n=n))/n,
                  n,
                  proc.time()-time_j))
    }
  }

  for(k in 1:K){
    MCLE[[k]] <- list()
    MCLE_mean[[k]] <- list()
    psi_list[[k]] <- list()
    S[[k]] <- list()
    n_k[k] <- MCLE_psilist_n[[k]][[1]][[5]]
    time_list[[k]] <- list()
    for(j in 1:J){
      MCLE[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[1]]
      MCLE_mean[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[2]]
      psi_list[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[3]]
      S[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[4]]
      time_list[[k]][[j]] <- MCLE_psilist_n[[k]][[j]][[6]]
    }
    time_list[[k]] <- do.call(rbind, time_list[[k]])
  }
  time_list <- do.call(rbind, time_list)
  time_max <- time_list[which.max(time_list[,1]),]
  time2 <- proc.time()
  V <- foreach(k = 1:K) %dopar% {
    sample_cov <- list()
    for (j in 1:J) {
      cov <- list()
      for (i in 1:J) {
        cov <- cbind(cov, Reduce("+", mapply(function(x,y) x %*% t(y), psi_list[[k]][[j]], psi_list[[k]][[i]], SIMPLIFY=FALSE))/N)
      }
      sample_cov[[j]] <- cov
    }
    do.call(rbind, sample_cov)
  }
  V_psi <- as.matrix(Matrix::bdiag(lapply(V, function(x) matrix(unlist(x), ncol=J*p, byrow=TRUE))))
  W <- solve(V_psi)
  
  print("Computing combined estimate.", quote=FALSE)
  scale <- matrix(0, p, p)
  main <- matrix(0, p, 1)
  for (k in 1:K){
    for(i in 1:J){
      for (j in 1:J){
        s_i <- matrix(unlist(S[[k]]), ncol=J*p, byrow=FALSE)[,((i-1)*p+1):(i*p)]
        s_j <- matrix(unlist(S[[k]]), ncol=J*p, byrow=FALSE)[,((j-1)*p+1):(j*p)]
        W_ij <- W[((k-1)*(J*p)+1):(k*J*p),((k-1)*(J*p)+1):(k*J*p)][((i-1)*p+1):(i*p), ((j-1)*p+1):(j*p)]
        main <- main + n_k[k]*s_i %*% W_ij %*% s_j %*% MCLE_mean[[k]][[j]]
        scale <- scale + n_k[k]*s_i %*% W_ij %*% s_j
      }
    }
  }
  
  if (method=="exact")
    output$coefficients <- as.vector(solve(scale)%*%main)
  time2 <- proc.time()-time2
  if (method=="iterative"){
    init_betas <- colMeans(do.call(rbind,lapply(MCLE_mean,function(x){do.call(rbind,x)})))
    optimization <- optim(par=init_betas, 
                          combined.estimation.mean, gr=estimation.deriv.mean, 
                          MCLE=MCLE, W=W, response=response, 
                          covariates=covariates, J=J, K=K, N=N, func=func2, funcderiv=func3, method="Nelder-Mead", control=list(maxit=500))
    init_betas <- optimization$par
    output$coefficients <- optim(par=init_betas, 
                              combined.estimation.mean, gr=estimation.deriv.mean, 
                              MCLE=MCLE, W=W, response=response, 
                              covariates=covariates, J=J, K=K, N=N, func=func2, funcderiv=func3, method="BFGS")$par
  }
  
  output$vcov <- solve(scale)
  
  if (method=="ridge"){
    if (is.null(folds)) 
      folds <- 1
    if (is.null(lam))
      lam <- 0
    if (is.null(folds)) 
      folds <- 1
    if (is.null(lam))
      lam <- 0
    ridge_optim <- ridge.estimate.mean(psi_list, MCLE, MCLE_mean, response, covariates, J, K, N, funcderiv=func3, folds, lam)
    output$coefficients <- as.vector(ridge_optim$coefficients)
    output$W <- ridge_optim$W
    output$lam <- ridge_optim$lam
    sensitivity_psi <- matrix(unlist(S[[k]]), ncol=J*p, byrow=FALSE)
    output$vcov <- solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))%*%
    sensitivity_psi%*%output$W%*%solve(W)%*%output$W%*%t(sensitivity_psi)%*%
    solve(sensitivity_psi%*%output$W%*%t(sensitivity_psi))/N
  }
  
  stopCluster(sock)
  print("Cluster stopped.", quote=FALSE)
  registerDoSEQ()
  
  output$family <- family
  output$response_indicator <- response_indicator
  output$subject_indicator <- subject_indicator
  output$J <- J
  output$K <- K
  output$MCLE <- list()
  output$MCLE$beta <- MCLE_mean
  output$MCLE$vcov <- MCLE
  output$V_psi.inv <- W
  output$time <- time1+time2+time_max

  return(output)
}
