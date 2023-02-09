################################### CV-MSE Solution with nugget effect ###############################
# This method allows to compute the CV-MSE solution as described in the paper with a nugget effect
library(parallel)

nc <- detectCores()
cl <- makeCluster(rep("localhost", nc))

matprod.par <- function(cl, A, B){
  if (ncol(A) != nrow(B)) stop("Matrices do not conforme")
  idx <- splitIndices(nrow(A), length(cl))
  Alist <- lapply(idx, function(ii) A[ii,,drop=FALSE])
  ans <- clusterApply(cl, Alist, get("%*%"), B)
  do.call(rbind, ans)
}

CVMSE_nugget_quadr <- function(covObject, GP, X, y, multistart = 5, sigma_eps){
  n <- nrow(X)
  d <- ncol(X)
  LOO.MSE.Par <- function(param){
    coef(covObject)[1:d] <- sqrt(param[1:d]^2)
    coef(covObject)[d+1] <- param[d+1]^2
    F <- GP$F
    K <- covMat(covObject, X) + sigma_eps^2*diag(n)
    Kinv <- chol2inv(chol(K))
    KF <- matprod.par(cl, K, F)
    tFKF <- matprod.par(cl, t(F), KF)
    tFKinvF_inv <- chol2inv(chol(tFKF))
    
    KinvF <- matprod.par(cl, Kinv, F)
    tFKinv <-matprod.par(cl, t(F), Kinv)
    KinvF_FtKinvF_inv <- matprod.par(cl, KinvF, tFKinvF_inv)
    Kbar_tilde <- matprod.par(cl, KinvF_FtKinvF_inv, tFKinv)
    
    Kbar <- Kinv - Kbar_tilde
    diagKbarinv2 <- diag(1/(diag(Kbar)^2), nrow = nrow(K), ncol = ncol(K))
    
    Kbary <- matprod.par(cl, Kbar, as.matrix(y))
    diagKbarinv2Kbary <- matprod.par(cl, diagKbarinv2, Kbary) 
    
    f <- 1/n * matprod.par(cl, t(Kbary), diagKbarinv2Kbary)
    return(f)
  }
  
  
  parIni <- matrix(runif(multistart*(d+1), min = 0, max = 20), multistart, d+1) 
  colnames(parIni) <- c(colnames(X), "sigma2")
  
  fitList <- list()
  fitList <- foreach::"%dopar%"(foreach::foreach(
    i = 1:multistart, 
    .errorhandling='remove', .packages = "kergp"), {
      args <- list(p = parIni[i, ], f = LOO.MSE.Par)
      do.call(nlm, args = args)  # do.call(stats::optim, args = args)
    })
  
  nlist <- length(fitList)
  
  optValueVec <- sapply(fitList, function(x) x$minimum)
  bestIndex <- which.min(optValueVec)
  report <- list(parIni = parIni,
                 par = t(sapply(fitList, function(x) x$estimate)),         
                 MSE_val = - sapply(fitList, function(x) x$minimum),
                 nIter  = sapply(fitList, function(x) x$iterations))
  colnames(report$par) <- colnames(report$parIni)
  opt <- fitList[[bestIndex]] 
  SolOpt <- opt$estimate
  
  result <- list(sigmaOpt = SolOpt[d+1]^2, thetaOpt = sqrt(SolOpt[1:d]^2))
  return(result)
}

CVMSE_nugget_LOO <-function(covObject, X, y, multistart = 5, var_noise = 0){
  n <- nrow(X)
  d <- ncol(X)
  CV.func <- function(param){
    coef(covObject)[1:d] <- sqrt(param[1:d]^2)
    coef(covObject)[d+1] <- param[d+1]^2
    GP <- gp(formula = y ~ 1, data = cbind(X, y), cov = covObject, estim = FALSE, varNoise = var_noise)
    m <- influence(GP)$mean
    MSE <- 1/n*sum((y-m)^2)
    return(MSE)
  }
  
  parIni <- matrix(runif(multistart*(d+1), min = 0, max = 50), multistart, d+1) 
  colnames(parIni) <- c(colnames(X), "sigma")
  
  fitList <- list()
  fitList <- foreach::"%dopar%"(foreach::foreach(
    i = 1:multistart, 
    .errorhandling='remove'), {
      args <- list(p = parIni[i, ], f = CV.func)
      do.call(nlm, args = args)
    })
  
  nlist <- length(fitList)
  
  optValueVec <- sapply(fitList, function(x) x$minimum)
  bestIndex <- which.min(optValueVec)
  report <- list(parIni = parIni,
                 par = t(sapply(fitList, function(x) x$estimate)),         
                 MSE_val = - sapply(fitList, function(x) x$minimum),
                 nIter  = sapply(fitList, function(x) x$iterations))
  colnames(report$par) <- colnames(report$parIni)
  opt <- fitList[[bestIndex]] 
  SolOpt <- opt$estimate
  
  # coef(covObject)[1:d] <- sqrt(thetaOpt^2)
  # GPOpt <- gp(formula = y ~ 1, data = cbind(X, y), cov = covObject, estim = FALSE,  varNoise = noise)
  # m <- influence(GPOpt)$mean
  # sd <- influence(GPOpt)$mea
  # sigmaOpt <- 1/n*sum((y-m)^2/sd^2)
  
  result <- list(sigmaOpt2 = SolOpt[1+d]^2, thetaOpt = sqrt(SolOpt[1:d]^2))
  return(result)
}
