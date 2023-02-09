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


source("Wassersteinpar.R")

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
                           data = dataFit, CovModel = CovDist, tol = 1e-4)
  }else{
    sigma_min <- MinSearchLow(func = Sigma2_Percent, low = min, high = max, fit = GP, 
                              data = dataFit, CovModel = CovDist, tol = 1e-4)
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
fit = fitMLE
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



# Evaluating the performance of the MLE P5+P95 model 
predP95MLE <- predUpperMLE + qnorm(1-alpha/2)*sdUpperMLE
predP05MLE <- predLowerMLE + qnorm(alpha/2)*sdLowerMLE
MeanP50MLE <- (predP05MLE + predP95MLE)/2

Q <- data.frame(MLE = 1 - sum((dataVal$y - predMLE)^2)/sum((dataVal$y - mean(dataVal$y ))^2),
                GP_Mean_MLE = 1 - sum((dataVal$y - MeanP50MLE)^2)/sum((dataVal$y - mean(dataVal$y ))^2),
                row.names = "Accuracy")
print(Q)

# Printing the quasi-Gaussian Percentile on training set
alpha = 0.10
P_upper_CV = quasi.Gaussian.CV(1-alpha/2, dataFit, GP_quantile_MLE95)
print(paste("Percentile of ",1-alpha/2,"on training set:", P_upper_CV))
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, GP_quantile_MLE5)
print(paste("Percentile of ",alpha/2,"on training set:", P_lower_CV))
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, GP_quantile_MLE95, nugget = TRUE)
print(paste("Percentile of ",1-alpha/2,"on validation set:", P_upper))
P_lower = Coverage.Probability(alpha/2, dataVal, GP_quantile_MLE5, nugget = TRUE)
print(paste("Percentile of ",alpha/2,"on validation set:", P_lower))
P_upper - P_lower

# MPIW and SdPIW MLE model
mean(abs(qnorm(1-alpha/2)*sdUpperMLE - qnorm(alpha/2)*sdLowerMLE))
sd(abs(qnorm(1-alpha/2)*sdUpperMLE - qnorm(alpha/2)*sdLowerMLE))


##################################################################


## Shifting lambda for wasserstein distance
fit = fitCV
alpha = 0.95
result <- optimize(wasserstein_shift, interval = c(0.1,3))
lambda_opt_CV <- result$minimum
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
fit = fitCV
alpha = 0.05

coef(CovRef) <- fit$covariance@par

## Searching for the best shifting parameter
result5 <- optimize(wasserstein_shift, interval = c(0.1,2))
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
print(paste("Percentile of ",1-alpha/2,"on training set:", P_upper_CV))
P_lower_CV = quasi.Gaussian.CV(alpha/2, dataFit, GP_quantile_CV5)
print(paste("Percentile of ",alpha/2,"on training set:", P_lower_CV))
P_upper_CV - P_lower_CV

# Printing the quasi-Gaussian Percentile on validation set
P_upper = Coverage.Probability(1-alpha/2, dataVal, GP_quantile_CV95, nugget = TRUE)
print(paste("Percentile of ",1-alpha/2,"on validation set:", P_upper))
P_lower = Coverage.Probability(alpha/2, dataVal, GP_quantile_CV5, nugget = TRUE)
print(paste("Percentile of ",alpha/2,"on validation set:", P_lower))
P_upper - P_lower

# MPIW and SdPIW CV model
mean(abs(qnorm(1-alpha/2)*sdUpperCV - qnorm(alpha/2)*sdLowerCV))
sd(abs(qnorm(1-alpha/2)*sdUpperCV - qnorm(alpha/2)*sdLowerCV))
