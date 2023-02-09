h_plus <- function(x, delta=10e-2){
  if (x>delta) {
    h <- 1
  } else if ( 0<x & x<=delta){
    h <- x/delta
  } else {
    h <- 0
  }
  return(h)
}

h_minus <- function(x, delta=10e-2){
  if (x>=0) {
    h <- 1
  } else if ( -delta<=x & x<0){
    h <- 1+x/delta
  } else {
    h <- 0
  }
  return(h)
}

quasi.Gaussian.CV <- function(a, dataFit, GP){
  y_pred <- influence(GP)$mean
  sd_pred <- influence(GP)$sd
  quant_pred <- (dataFit$y - y_pred)/sd_pred
  n <- length(y_pred)
  if (a>1/2){
    psi_a <- 1/n*sum(sapply(qnorm(a)-quant_pred, h_plus))
  } else {
    psi_a <- 1/n*sum(sapply(qnorm(a)-quant_pred, h_minus))
  }
  return(psi_a)
}

Coverage.Probability <- function(alpha, dataVal, GP, nugget = FALSE){
  y_pred <- predict(GP, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = nugget)$mean
  sd_pred <- predict(GP, newdata = as.matrix(dataVal[ , inputs]), type = "UK", forceInterp = nugget)$sd
  quant_pred <- (dataVal$y - y_pred)/sd_pred
  CP <- 0
  n <- length(y_pred)
  for (k in 1:n) {
    if ( quant_pred[k] <= qnorm(alpha))
      CP <- CP + 1
    else
      CP <- CP + 0
  }
  CP <- CP/n
  return(CP)
}

score.fun <- function(IC) {
  score <- 0
  ktest <- nrow(IC)
  for (k in 1:ktest) {
    if ( (IC$y[k] > IC$UpperQ[k]) | ( IC$y[k]  < IC$LowerQ[k]) )  
      score <- score + 0
    else
      score <- score + 1 
  }
  score <- score/ktest
  return(score)
}

# Fitting GP variance noise by Cross-Validation
Sigma2_Percent <- function(sigma, alpha, fit, data, CovModel){
  kergp::coef(CovModel) <- fit$covariance@par
  kergp::coef(CovModel)[d+1] <- sigma
  
  GP <- gp(formula = y ~ 1, data = data, cov = CovModel, estim = FALSE, varNoise = fit$varNoise)
  percent <- quasi.Gaussian.CV(alpha, data, GP)
  percent_alpha <- percent
  return(percent_alpha)
}


# binary search of sigma_min (theta) (OK)
MinSearch <- function(func = Sigma2_Percent, low, high, fit, data, CovModel, tol = (CovModel$covariance@par[ncol(data)])/100 ) {
  n <- nrow(data)
  k0 <- floor(alpha*n)
  k1 <- floor(alpha*n)+1
  mid <- (low + high) / 2
  value <- func(mid, alpha, fit, data, CovModel)
  value0 <- func(high, alpha, fit, data, CovModel)
  if ( abs(high - low) < tol) {
    print(paste("Sigma.Min = ", high, "Percent.Value =", value0 ))
    return(high)
  } else {
    if ( value > k1/n )
      MinSearch(func, low, mid, fit, data, CovModel, tol)
    else if ( value < k0/n )
      MinSearch(func, mid, high, fit, data, CovModel, tol)
    else if (value == k0) 
      MinSearch(func, mid, high, fit, data, CovModel, tol)
    else 
      MinSearch(func, low, mid, fit, data, CovModel, tol)
  }
}

# binary search of sigma_min for lower quantiles (theta)
MinSearchLow <- function(func = Sigma2_Percent, low, high, fit, data, CovModel, tol = (CovModel$covariance@par[ncol(data)])/100 ) {
  n <- nrow(data)
  k0 <- floor(alpha*n)+1
  k1 <- floor(alpha*n)
  mid <- (low + high) / 2
  value <- func(mid, alpha, fit, data, CovModel)
  value0 <- func(high, alpha, fit, data, CovModel)
  if ( abs(high - low) < tol) {
    print(paste("Sigma.Min = ", low, "Percent.Value =", value0 ))
    return(high)
  } else {
    if ( value > k1/n )
      MinSearchLow(func, mid, high, fit, data, CovModel, tol)
    else if ( value < k0/n )
      MinSearchLow(func, low, mid, fit, data, CovModel, tol)
    else if (value == k0/n) 
      MinSearchLow(func, low, mid, fit, data, CovModel, tol)
    else 
      MinSearchLow(func, mid, high, fit, data, CovModel, tol)
  }
}
