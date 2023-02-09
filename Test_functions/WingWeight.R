# Importing all required packages and methods 
source("test_func.R")
source("Multistart.R")
source("Sigma_LOO_alpha.R")
source("CVMSE_nugget.R")
library(kergp)
library("mvtnorm")

# Setting seed for numerical experiments
set.seed(123456)

# Initializing parameters and design of experiment
normalize <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
}

n <- 600
d <- 10

x1 <- runif(n, min = 150, max = 200)
x2 <- runif(n, min = 220, max = 300)
x3 <- runif(n, min = 6, max = 10)
x4 <- runif(n, min = -10, max = 10)
x5 <- runif(n, min = 16, max = 45)
x6 <- runif(n, min = 0.5, max = 1)
x7 <- runif(n, min = 0.08, max = 0.18)
x8 <- runif(n, min = 2.5, max = 6)
x9 <- runif(n, min = 1700, max = 2500)
x10 <- runif(n, min = 0.025, max = 0.08)

X <- cbind(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6,
           x7 = x7, x8 = x8, x9 = x9, x10 = x10)
inputs <- paste("x", 1L:d, sep = "")


# Applying Wing Weight function
testfun <- wingweight
data <- data.frame(X = apply(X, 2, normalize), y = apply(X, 1, testfun) ) + rnorm(n, sd = 5)
colnames(data)[1:d] = inputs
hist(data$y)

# Dividing data into training and validation set
ind.sample <- sample(n, size = 0.75*n, replace = FALSE)
dataFit <- data[ind.sample, ]
y = dataFit$y
X = dataFit[,1:d]
dataVal <- data[-ind.sample, ]

# Creation of "CovRadial" object with one range by input
CovModel <- covRadial(inputs = inputs, d = d, k1Fun1 = k1Fun1Matern3_2, cov = "homo")

######################################## MLE Estimator ######################################

# Creation of a GP model by Maximum Likelihood Estimator method, 
# fitMLE <- gp(formula = y ~ 1, dataFit, cov = CovModel, estim = TRUE)

# Creation of a GP model by Maximum Likelihood Estimator method, using the sequential approach
fitMLE <- GP.joint.model(dataFit, CovModel, d)

# Predicting the output for the validation data set : mean and std
predMLE <- predict(fitMLE, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = TRUE)$mean
sdMLE <- predict(fitMLE, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = TRUE)$sd
var_nugget = fitMLE$varNoise

# Evaluation metrics : RMSE and Accuracy score
QMLE <- 1 - sum((dataVal$y - predMLE)^2)/sum((dataVal$y - mean(dataVal$y))^2)
print(paste("Accuracy Q2:", QMLE))

# Printing the quasi-Gaussian Percentile on training set
alpha = 0.10
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitMLE)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitMLE)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitMLE, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitMLE, nugget = TRUE)
P_upper - P_lower


# Do the same for alpha = 0.05
alpha = 0.05
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitMLE)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitMLE)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitMLE, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitMLE, nugget = TRUE)
P_upper - P_lower

# Do the same for alpha = 0.20
alpha = 0.01
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitMLE)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitMLE)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitMLE, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitMLE, nugget = TRUE)
P_upper - P_lower

# Z_tilde <- ( y - influence(fitMLE)$mean)/influence(fitMLE)$sd
# shapiro.test(Z_tilde)

################################### MSE Cross-Validation approach ###################################

# Initializing a GP model by Cross-Validation method
covFitCV = covRadial(inputs = inputs, d = d, k1Fun1 = k1Fun1Matern3_2, cov = "homo")
fitCV <- gp(formula = y ~ 1, data = dataFit, cov = covFitCV, estim = FALSE, varNoise = var_nugget)

source("Multistart_par.R")
param_quadr <- CVMSE_nugget_quadr(covObject = covFitCV, GP = fitCV, X = dataFit[ ,1:d], 
                                  y = dataFit$y, multistart = 8, sigma_eps = sqrt(var_nugget))

stopCluster(cl)
param_quadr

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



# Do the same for alpha = 0.20
alpha = 0.05
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitCV)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitCV)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitCV, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitCV, nugget = TRUE)
P_upper - P_lower

# Do the same for alpha = 0.20
alpha = 0.01
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, fitCV)
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, fitCV)
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, fitCV, nugget = TRUE)
P_lower = Coverage.Probability(alpha/2, dataVal, fitCV, nugget = TRUE)
P_upper - P_lower

################################### Full-Bayesian approach ###################################

# Define prior on both sigma2 and vtheta
library(mcmc)

covObject = CovModel
GP <- gp(formula = y ~ 1, dataFit, cov = covObject, estim = FALSE)
par =  fitMLE$covariance@par
sigma2_eps = fitMLE$varNoise


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
  
  logLik <-  -0.5 * (2 * sum(log(diag(L))) + sseStar )
  
  logLik <- as.numeric(logLik)
  
  return(logLik)
  
}

# run MCMC for the posterior distribution
out = metrop(logLikFun0, initial= fitMLE$covariance@par, nbatch = 1e5, scale = 2.5)
out$accept
out$time
out$final
plot(out$batch[ , 1])

# Build a sample of Y(x_new) given X,y for all MCMC realizations
burn.in = 4e4
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
  if ( dataVal$y[k] <= PI_boundsVal[k, 2] && PI_boundsVal[k, 1] <= dataVal$y[k] ){
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


# Deduce the quantile Y_alpha/2 and Y_1-alpha/2 
PI_boundsVal = t(apply(ypred_Val, 2, quantile, probs = c(0.025,0.975)))
PI_boundsTrain = t(apply(ypred_Train, 2, quantile, probs = c(0.025,0.975)))

CP_95 <- c()
n_Val <- length(dataVal$y)
for (k in 1:n_Val) {
  if ( dataVal$y[k] <= PI_boundsVal[k, 2] && PI_boundsVal[k, 1] <= dataVal$y[k] ){
    CP_95 <- c(CP_95, 1)
  }else {
    CP_95 <- c(CP_95, 0)
  }
}
CP_95 <- mean(CP_95)
print(CP_95)

Ptile_95 <- c()
n_Val <- length(dataFit$y)
for (k in 1:n_Val) {
  if ( dataFit$y[k] <= PI_boundsTrain[k, 2] && PI_boundsTrain[k, 1] <= dataFit$y[k] ){
    Ptile_95 <- c(Ptile_95, 1)
  }else {
    Ptile_95 <- c(Ptile_95, 0)
  }
}
Ptile_95 <- mean(Ptile_95)
print(Ptile_95)



# Deduce the quantile Y_alpha/2 and Y_1-alpha/2 
PI_boundsVal = t(apply(ypred_Val, 2, quantile, probs = c(0.005,0.995)))
PI_boundsTrain = t(apply(ypred_Train, 2, quantile, probs = c(0.005,0.995)))

CP_99 <- c()
n_Val <- length(dataVal$y)
for (k in 1:n_Val) {
  if ( dataVal$y[k] <= PI_boundsVal[k, 2] && PI_boundsVal[k, 1] <= dataVal$y[k] ){
    CP_99 <- c(CP_99, 1)
  }else {
    CP_99 <- c(CP_99, 0)
  }
}
CP_99 <- mean(CP_99)
print(CP_99)

Ptile_99 <- c()
n_Val <- length(dataFit$y)
for (k in 1:n_Val) {
  if ( dataFit$y[k] <= PI_boundsTrain[k, 2] && PI_boundsTrain[k, 1] <= dataFit$y[k] ){
    Ptile_99 <- c(Ptile_99, 1)
  }else {
    Ptile_99 <- c(Ptile_99, 0)
  }
}
Ptile_99 <- mean(Ptile_99)
print(Ptile_99)

