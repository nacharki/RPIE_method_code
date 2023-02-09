# Importing all required packages and methods 
source("test_func.R")
source("Multistart.R")
source("Sigma_LOO_alpha.R")
source("CVMSE_nugget.R")
library(kergp)
library("mvtnorm")

# Setting seed for numerical experiments
set.seed(1000) 

# Initializing parameters and design of experiment
n <- 600
d <- 10
inputs <- paste("x", 1L:d, sep = "")

X <- matrix(runif(n*d), ncol = d)
colnames(X) <- inputs

# Applying Bratleyetal function
testfun <- zhou98 
data <- data.frame(X, y = apply(X, 1, testfun) ) # + rnorm(n, sd=0.08)
data$y <- log(data$y)/(d*log(10))

# Dividing data into training and validation set
ind.sample <- sample(1:n, size = 0.75*n, replace = FALSE)
dataFit <- data[ind.sample, ]
y = dataFit$y
X = dataFit[,1:d]
dataVal <- data[-ind.sample, ]

# Creation of "CovRadial" object with one range by input
CovModel <- covRadial(inputs = inputs, d = d, k1Fun1 = k1Fun1Matern3_2, cov = "homo")
fitMLE <- gp(formula = y ~ 1, dataFit, cov = CovModel, estim = TRUE)

# Predicting the output for the validation data set : mean and std
predMLE <- predict(fitMLE, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = TRUE)$mean
sdMLE <- predict(fitMLE, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = TRUE)$sd

QMLE <- 1 - sum((dataVal$y - predMLE)^2)/sum((dataVal$y - mean(dataVal$y))^2)
print(paste("Accuracy :", QMLE))

# Printing the quasi-Gaussian Percentile on training set
alpha = 0.10
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitMLE)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitMLE)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitMLE, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitMLE, nugget = TRUE)
P_upper - P_lower


# MPIW and SdPIW MLE model
mean(abs(qnorm(1-alpha/2)*sdMLE - qnorm(alpha/2)*sdMLE))
sd(abs(qnorm(1-alpha/2)*sdMLE - qnorm(alpha/2)*sdMLE))

timeMLE


################################### MSE Cross-Validation approach ###################################

# Initializing a GP model by Cross-Validation method
var_nugget = fitMLE$varNoise
covFitCV = covRadial(inputs = inputs, d = d, k1Fun1 = k1Fun1Matern3_2, cov = "homo")
fitCV <- gp(formula = y ~ 1, data = dataFit, cov = covFitCV, estim = FALSE, varNoise = var_nugget)

# Searching for optimum hyper-parameters with MSE criterion : this operation may take time
param_LOO = CVMSE_nugget_LOO(covFitCV, X = dataFit[ ,1:d], y = dataFit$y, 
                             multistart = 5, var_noise = var_nugget)
param_LOO


source("Multistart_par.R")
system.time(param_quadr <- CVMSE_nugget_quadr(covObject = covFitCV, GP = fitCV, X = dataFit[ ,1:d], 
                                              y = dataFit$y, multistart = 8, sigma_eps = sqrt(var_nugget)))

stopCluster(cl)
param_quadr

# TimeCV = 1882.22 

# Creation of a GP model by Cross-Validation method with MSE criterion 
coef(covFitCV)[1:d] <- param_quadr$thetaOpt
coef(covFitCV)[d+1] <- param_quadr$sigmaOpt

fitCV <- gp(formula = y ~ 1, data = dataFit, cov = covFitCV, estim = FALSE,
            noise = TRUE, varNoise = var_nugget)
predCV <- predict(fitCV, newdata = as.matrix(dataVal[ , inputs]), type = "UK")$mean
sdCV <- predict(fitCV, newdata = as.matrix(dataVal[ , inputs]), type = "UK")$sd

Q <- data.frame(MLE = 1 - sum((dataVal$y - predMLE)^2)/sum((dataVal$y - mean(dataVal$y ))^2),
                CV = 1 - sum((dataVal$y - predCV)^2)/sum((dataVal$y - mean(dataVal$y ))^2),
                row.names = "Accuracy")
print(Q)


# Printing the quasi-Gaussian Percentile on training set
alpha = 0.10
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitCV)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitCV)
P_upper_CV - P_lower_CV


# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitCV, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitCV, nugget = TRUE)
P_upper - P_lower

# MPIW and SdPIW MLE model
mean(abs(qnorm(1-alpha/2)*sdCV - qnorm(alpha/2)*sdCV))
sd(abs(qnorm(1-alpha/2)*sdCV - qnorm(alpha/2)*sdCV))


######################################## Full-Bayesian Method ######################################

# Define prior on both sigma2 and vtheta
library(mcmc)

covObject = CovModel
sigma2_eps = fitMLE$varNoise
GP <- gp(formula = y ~ 1, dataFit, cov = covObject, estim = FALSE, varNoise = sigma2_eps)
par =  fitMLE$covariance@par

logLikFun0 <- function(par) {
  
  n <- length(y)
  lpar <- length(par)
  
  lparNN <- lpar
  coef(covObject) <- sqrt(par^2) ## MMM METHOD HERE
  C <- covMat(covObject, X, checkNames = FALSE, index = 1L) 
  diag(C) <- diag(C) + sigma2_eps
  
  L <- chol(C)
  L <- t(L)
  
  ## yStar <- backsolve(L, y, upper.tri = FALSE)
  
  if (!is.null(GP$F)) {
    pF <- NCOL(GP$F)
    p1 <- pF + 1L
    FStarPlus <- forwardsolve(L, cbind(GP$F, y))
    qrFStarPlus <- qr(FStarPlus)
    dStar <- qr.R(qrFStarPlus)[p1, p1]
    eStar  <- qr.Q(qrFStarPlus)[ , p1] * dStar
    sseStar <- dStar^2
  } else {
    eStar <- forwardsolve(L, y)
    sseStar <- crossprod(eStar)
  }
  
  logLik <-  -0.5 * ( 2 * sum(log(diag(L))) + sseStar )
  
  logLik <- as.numeric(logLik)
  
  
  return(logLik)
  
}

out = metrop(logLikFun0, initial = fitMLE$covariance@par, nbatch = 1e5,  scale = 0.3)
# out = metrop(out, scale = 0.0053, nbatch = 1e2)
out$accept
out$time
out$final
plot(out$batch[ , 1])


# Build a sample of Y(x_new) given X,y for all MCMC realizations
burn.in = 5e4
batch = out$batch[-(1:burn.in ),]
ypred_Val = matrix( nrow = nrow(batch), ncol = nrow(dataVal))
ypred_Train = matrix( nrow = nrow(batch), ncol = nrow(dataFit))
for(N in 1:nrow(batch)){
  coef(covObject) <- sqrt(batch[N,]^2)
  
  GP <- gp(formula = y ~ 1, dataFit, cov = covObject, estim = FALSE, varNoise = sigma2_eps)
  
  predCV_yTrain <- influence(GP)$mean
  sdCV_yTrain <- influence(GP)$sd
  ypred_Train[N,] = sapply(X = 1:nrow(dataFit),  FUN = function(xnew) rnorm(1, mean = predCV_yTrain[xnew], sd = sdCV_yTrain[xnew] ) )
  
  
  predVal_yN <- predict(GP, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = TRUE)$mean
  sdVal_yN <- predict(GP, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = TRUE)$sd
  ypred_Val[N,] = sapply(X = 1:nrow(dataVal),  FUN = function(xnew) rnorm(1, mean = predVal_yN[xnew], sd = sdVal_yN[xnew] ) )
}

# Predict the mean prediction on training and validation set
predFullBayes = as.vector(apply(ypred_Val, 2, mean))
Q_FullBayes <- 1 - sum((dataVal$y - predFullBayes)^2)/sum((dataVal$y - mean(dataVal$y))^2)
print(paste("Accuracy Q2:", Q_FullBayes))

# Deduce the quantile Y_alpha/2 and Y_1-alpha/2 
PI_boundsVal = t(apply(ypred_Val, 2, quantile, probs = c(0.05,0.95)))
PI_boundsTrain = t(apply(ypred_Train, 2, quantile, probs = c(0.05,0.95)))

CP_90 <- c()
n_Val <- length(dataVal$y)
for (k in 1:n_Val) {
  if (dataVal$y[k] <= PI_boundsVal[k, 2] && PI_boundsVal[k, 1] <= dataVal$y[k] ){
    CP_90 <- c(CP_90, 1)
  }else {
    CP_90 <- c(CP_90, 0)
  }
  
}
CP_90 <- mean(CP_90)
print(CP_90)

Ptile_90 <- c()
n_Val <- length(dataFit$y)
for (k in 1:n_Val) {
  if ( dataFit$y[k] <= PI_boundsTrain[k, 2] && PI_boundsTrain[k, 1] <= dataFit$y[k] ){
    Ptile_90 <- c(Ptile_90, 1)
  }else {
    Ptile_90 <- c(Ptile_90, 0)
  }
}
Ptile_90 <- mean(Ptile_90)
print(Ptile_90)

# MPIW CV model
mean(abs(PI_boundsVal[,2] - PI_boundsVal[,1]))
# SdPIW CV model
sd(abs(PI_boundsVal[,2] - PI_boundsVal[,1]))





################################## Apply the RPIE method on the MLE solution ################################

# Define a function that update the GLS estimator of hat_beta for a given Covariance function CovModel and nugget effect var_nugget
betaHat <- function(fit, CovModel, var_nugget, data){
  CovRef <-  CovModel
  coef(CovRef) <- fit$covariance@par
  C <- covMat(CovModel, data)
  diag(C) <- diag(C)  + var_nugget
  
  p <- NCOL(fit$F)
  p1 <- p + 1L
  L <- t(chol(C))
  FStarPlus <- forwardsolve(L, cbind(fit$F, dataFit$y))
  qrFStarPlus <- qr(FStarPlus)
  RStarPlus <- qr.R(qrFStarPlus)
  beta <- backsolve(RStarPlus[1L:p, 1L:p], RStarPlus[1L:p,  p1]) 
  return(beta)
}

## Import the update version of Wasserstein distance, the old version has some computing issues
source("Wassersteinpar.R")

## Initiate the RPIE method
CovRef =  CovModel
fit = fitMLE
coef(CovRef) <- fit$covariance@par
S1 <- covMat(CovRef, dataFit[, 1:d])
m1 <- fit$F %*% betaHat(fit, CovRef, var_nugget, data = dataFit[, 1:d])
max <- 1e6; min <- 1e-6
alpha = 0.95

wasserstein_shift <- function(lambda){
  CovDist <- CovModel
  kergp::coef(CovDist) <- lambda*fit$covariance@par
  GP <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)
  
  if(alpha>=1/2){
    sigma_min <- MinSearch(func = Sigma2_Percent, low = min, high = max, fit = GP, 
                           data = dataFit, CovModel = CovDist, tol = 1e-3)
  }else{
    sigma_min <- MinSearchLow(func = Sigma2_Percent, low = min, high = max, fit = GP, 
                              data = dataFit, CovModel = CovDist, tol = 1e-3)
  }
  kergp::coef(CovDist)[d+1] <- sigma_min
  S2 <- covMat(CovDist, dataFit[, 1:d])
  m2 <- fit$F %*% betaHat(fit, CovDist, var_nugget, data = dataFit[, 1:d])
  # dist_lambda <- wasserstein.tol(rep(0, nrow(dataFit)), S1, 
  #                               rep(0, nrow(dataFit)), S2, check = FALSE, tol = 1e-4)
  
  dist_lambda <- wasserstein.tol(m1, S1, m2, S2, check = FALSE, tol = 1e-4)
  return(dist_lambda)
}

## Searching for the best shifting parameter
fit = fitMLE
result <- optimize(wasserstein_shift, interval = c(0.1,4))
lambda_opt_MLE <- result$minimum
lambda_opt_MLE

# Validating result and comparison with initial model
alpha = 0.95
fit = fitMLE
CovDist <- CovModel
kergp::coef(CovDist) <- lambda_opt_MLE*fit$covariance@par
GP_quantile_MLE <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

sigma <- MinSearch(func = Sigma2_Percent, low = min, high = max, fit = GP_quantile_MLE, 
                   data = dataFit, CovModel = CovDist, tol = 1e-6)

kergp::coef(CovDist)[d+1] <- sigma
GP_quantile_MLE95 <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

predUpperMLE <- predict(GP_quantile_MLE95, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterpert = TRUE)$mean
sdUpperMLE <- predict(GP_quantile_MLE95, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterpert = TRUE)$sd


## Shifting lambda for wasserstein distance
alpha = 0.05
result <- optimize(wasserstein_shift, interval = c(0.1,10))
lambda_opt_MLE_P5 <- result$minimum
lambda_opt_MLE_P5

# Validating result and comparison with initial model
CovDist <- CovModel
kergp::coef(CovDist) <- lambda_opt_MLE_P5*fit$covariance@par
GP_quantile_MLE5 <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

sigma <- MinSearchLow(func = Sigma2_Percent, low = min, high = max, fit = GP_quantile_MLE5, 
                      data = dataFit, CovModel = CovDist, tol = 1e-4)

kergp::coef(CovDist)[d+1] <- sigma
GP_quantile_MLE5 <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

predLowerMLE <- predict(GP_quantile_MLE5, newdata = as.matrix(dataVal[ , inputs]), type = "UK")$mean
sdLowerMLE <- predict(GP_quantile_MLE5, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterpert = TRUE)$sd

# Printing the quasi-Gaussian Percentile on training set
alpha = 0.10
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, GP_quantile_MLE95)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, GP_quantile_MLE5)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, GP_quantile_MLE95, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, GP_quantile_MLE5, nugget = TRUE)
P_upper - P_lower

# MPIW and SdPIW MLE model
mean(abs(qnorm(1-alpha/2)*sdUpperMLE - qnorm(alpha/2)*sdLowerMLE))
sd(abs(qnorm(1-alpha/2)*sdUpperMLE - qnorm(alpha/2)*sdLowerMLE))

################################## Apply the RPIE method on the CV-MSE solution ################################

## Shifting lambda for wasserstein distance
fit = fitCV
alpha = 0.95
result95 <- optimize(wasserstein_shift, interval = c(0.1,2))
lambda_opt_CV <- result95$minimum
lambda_opt_CV

# Validating result and comparison with initial model
alpha = 0.95
CovDist <- CovModel
kergp::coef(CovDist) <- lambda_opt_CV*fit$covariance@par
GP_quantile_CV <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

sigma <- MinSearch(func = Sigma2_Percent, low = min, high = max, fit = GP_quantile_CV, 
                   data = dataFit, CovModel = CovDist, tol = 1e-6)

kergp::coef(CovDist)[d+1] <- sigma

GP_quantile_CV95 <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

# Evaluating the performance of the P95 model 
predUpperCV <- predict(GP_quantile_CV95, newdata = as.matrix(dataVal[ , inputs]), type = "UK")$mean
sdUpperCV <- predict(GP_quantile_CV95, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterpert = TRUE)$sd


## Shifting lambda for wasserstein distance
alpha = 0.05
coef(CovRef) <- fit$covariance@par

## Searching for the best shifting parameter
result5 <- optimize(wasserstein_shift, interval = c(0.1,5))
lambda_opt_CV_P5 <- result5$minimum
lambda_opt_CV_P5

# Validating result and comparison with initial model
alpha = 0.05
fit = fitCV
CovDist <- CovModel
kergp::coef(CovDist) <- lambda_opt_CV_P5*fit$covariance@par
GP_quantile_CV5 <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

sigma <- MinSearchLow(func = Sigma2_Percent, low = min, high = max, fit = GP_quantile_CV5, 
                      data = dataFit, CovModel = CovDist, tol = 1e-4)

kergp::coef(CovDist)[d+1] <- sigma

GP_quantile_CV5 <- kergp::gp(formula = y ~ 1, data = dataFit, cov = CovDist, estim = FALSE, varNoise = var_nugget)

# Evaluating the performance of the P05 model 
predLowerCV <- predict(GP_quantile_CV5, newdata = as.matrix(dataVal[ , inputs]), type = "UK")$mean
sdLowerCV <- predict(GP_quantile_CV5, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterpert = TRUE)$sd


# Evaluating the performance of the CV P95+P5 model 
predP95CV <- predUpperCV + qnorm(1-alpha/2)*sdUpperCV
predP05CV <- predLowerCV + qnorm(alpha/2)*sdLowerCV
MeanP50CV <- (predP05CV + predP95CV)/2

Q <- data.frame(CV = 1 - sum((dataVal$y - predCV)^2)/sum((dataVal$y - mean(dataVal$y ))^2),
                GP_Mean_CV = 1 - sum((dataVal$y - MeanP50CV)^2)/sum((dataVal$y - mean(dataVal$y ))^2),
                row.names = "Accuracy")
print(Q)

# Printing the quasi-Gaussian Percentile on training set
alpha = 0.10
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, GP_quantile_CV95)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, GP_quantile_CV5)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, GP_quantile_CV95, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, GP_quantile_CV5, nugget = TRUE)
P_upper - P_lower

# MPIW and SdPIW CV model
mean(abs(qnorm(1-alpha/2)*sdUpperCV - qnorm(alpha/2)*sdLowerCV))
sd(abs(qnorm(1-alpha/2)*sdUpperCV - qnorm(alpha/2)*sdLowerCV))


