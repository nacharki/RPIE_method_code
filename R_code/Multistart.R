# Fitting GP model with appropriate nugget effect
GP.joint.model <- function(data, CovModel, d = ncol(data)-1){
  nugget = 0
  
  GPm1 <- gp(formula = y ~ 1, data, cov = CovModel, estim = TRUE)
  ym1 <- plot(GPm1, trend.reestim = TRUE, kriging.type = "UK")$mean
  
  dataFitv1 <- data.frame(data[,1:d], y = (data$y-ym1)^2)
  GPv1 <- gp(formula = y ~ 1, data = dataFitv1, cov = CovModel, estim = TRUE)
  
  nugget <- predict(GPv1, data[,1:d], "UK")$mean 
  nuggetMat <- diag(as.vector(nugget), nrow = length(nugget))
  noise <- mean(nugget)
  
  GPm2 <- gp(formula = y ~ 1, data = data, cov = CovModel, estim = TRUE, varNoiseLower = noise, varNoiseUpper = noise)
  ym2 <- plot(GPm2, trend.reestim = TRUE, kriging.type = "UK")$mean
  
  dataFitv2 <- data.frame(data[,1:d], y = (data$y-ym2)^2)
  GPv2 <- gp(formula = y ~ 1, data = dataFitv2, cov = CovModel, estim = TRUE)
  
  nugget <- predict(GPv2, data[,1:d], "UK")$mean 
  nuggetMat <- diag(as.vector(nugget), nrow = length(nugget))
  noise <- mean(nugget)
  
  GPm <- gp(formula = y ~ 1, data = dataFit, cov = CovModel, estim = TRUE, varNoiseLower = noise, varNoiseUpper = noise)
  return(GPm)
}


# Fitting GP variance noise by Cross-Validation
Optim.Noise <- function(noise){
  GP <- gp(formula = y ~ 1, data = dataFit, cov = CovModel, estim = FALSE, varNoise = noise)
  m <- influence(GP)$mean
  MSE <- 1/n*sum((dataFit$y-m)^2)
  return(MSE)
}

# Optimization with multistart
CVMSE_multistart <- function(covObject, X, y, multistart = 5, noise = 0){
  n <- nrow(X)
  d <- ncol(X)
  CV.func <- function(theta){
    coef(covObject)[1:d] <- sqrt(theta^2)
    GP <- gp(formula = y ~ 1, data = cbind(X, y), cov = covObject, estim = FALSE, varNoise = noise)
    m <- influence(GP)$mean
    MSE <- 1/n*sum((y-m)^2)
    return(MSE)
  }
  
  parIni <- matrix(runif(multistart*d, min = 0, max = 50), multistart, d) 
  colnames(parIni) <- colnames(X)
  
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
  thetaOpt <- opt$estimate
  
  coef(covObject)[1:d] <- sqrt(thetaOpt^2)
  GPOpt <- gp(formula = y ~ 1, data = cbind(X, y), cov = covObject, estim = FALSE,  varNoise = noise)
  m <- influence(GPOpt)$mean
  sd <- influence(GPOpt)$mean
  sigmaOpt <- 1/n*sum((y-m)^2/sd^2)
  result <- list(sigmaOpt = sigmaOpt, thetaOpt = sqrt(thetaOpt^2))
  return(result)
}
