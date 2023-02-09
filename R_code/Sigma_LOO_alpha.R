# Define the Heaviside functions h_delta^plus as described in the paper
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

# Define the Heaviside function h_delta^minus as described in the paper
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

# Compute the quasi Gaussian percentile psi_a^(delta) i.e. the number of Leave-One-Out predictions falling below the Gaussian quantile of level a
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

# Compute the empirical coverage probability on a validation set dataVal using a GP model
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

# Compute the quasi Gaussian percentile psi_a for a given model GP
Sigma2_Percent <- function(sigma, alpha, GP, data, CovModel){
  kergp::coef(CovModel) <- GP$covariance@par
  kergp::coef(CovModel)[d+1] <- sigma
  
  GPModel <- gp(formula = y ~ 1, data = data, cov = CovModel, estim = FALSE, varNoise = GP$varNoise)
  percent <- quasi.Gaussian.CV(alpha, data, GPModel)
  percent_alpha <- percent
  return(percent_alpha)
}


# Function for finding the minimal sigma2 that satisfies the coverage propability equal to a, for a quantile a>1/2 using binary search
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

# Function for finding the minimal sigma2 that satisfies the coverage propability equal to a, for a quantile a>1/2
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
