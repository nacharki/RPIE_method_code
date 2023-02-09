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

# A modifier
CVMSE_multistart_par <- function(covObject, GP, X, y, multistart = 5, noise){
  n <- nrow(X)
  d <- ncol(X)
  LOO.MSE.Par <- function(theta){
    coef(covObject)[1:d] <- sqrt(theta^2)
    F <- GP$F
    R <- covMat(covObject, X)
    Rinv <- chol2inv(chol(R))
    RF <- matprod.par(cl, R, F)
    tFRF <- matprod.par(cl, t(F), RF)
    tFRinvF_inv <- chol2inv(chol(tFRF))
  
    RinvF <- matprod.par(cl, Rinv, F)
    tFRinv <-matprod.par(cl, t(F), Rinv)
    RinvF_FtRinvF_inv <- matprod.par(cl, RinvF, tFRinvF_inv)
    Rbar_tilde <- matprod.par(cl, RinvF_FtRinvF_inv, tFRinv)
  
    Rbar <- Rinv - Rbar_tilde
    diagRbarinv2 <- diag(1/(diag(Rbar)^2), nrow = nrow(R), ncol = ncol(R))
    
    Rbary <- matprod.par(cl, Rbar, as.matrix(y))
    diagRbarinv2Rbary <- matprod.par(cl, diagRbarinv2, Rbary) 
    
    f <- 1/n* matprod.par(cl, t(Rbary), diagRbarinv2Rbary)
    return(f)
  }


  parIni <- matrix(runif(multistart*d, min = 0, max = 20), multistart, d) 
  colnames(parIni) <- colnames(X)
  
  fitList <- list()
  fitList <- foreach::"%dopar%"(foreach::foreach(
    i = 1:multistart, 
    .errorhandling='remove', .packages = "kergp"), {
      args <- list(p = parIni[i, ], f = LOO.MSE.Par)
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
  thetaOpt <- opt$estimate
  
  coef(covObject)[1:d] <- sqrt(thetaOpt^2)
  F <- GP$F
  R <- covMat(covObject, X)
  Rinv <- chol2inv(chol(R))
  RF <- matprod.par(cl, R, F)
  tFRF <- matprod.par(cl, t(F), RF)
  tFRinvF_inv <- chol2inv(chol(tFRF))
  
  RinvF <- matprod.par(cl, Rinv, F)
  tFRinv <-matprod.par(cl, t(F), Rinv)
  RinvF_FtRinvF_inv <- matprod.par(cl, RinvF, tFRinvF_inv)
  Rbar_tilde <- matprod.par(cl, RinvF_FtRinvF_inv, tFRinv)
  
  Rbar <- Rinv - Rbar_tilde
  
  diagRbarinv1 <- diag(1/(diag(Rbar)^1), nrow = nrow(R), ncol = ncol(R))
  
  Rbary <- matprod.par(cl, Rbar, as.matrix(y))
  diagRbarinv1Rbary <- matprod.par(cl, diagRbarinv1, Rbary) 
  
  sigmaOpt <- 1/n* matprod.par(cl, t(Rbary), diagRbarinv1Rbary)
  coef(covObject)[d+1] <- sigmaOpt
  
  result <- list(sigmaOpt = sigmaOpt, thetaOpt = sqrt(thetaOpt^2))
  return(result)
}






