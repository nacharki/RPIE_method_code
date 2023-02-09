sqrtmatrix.tol <- function(mat, tol)
  {
    # Checking the symmetry of the matrix
    if (!isSymmetric(mat))
      stop("The matrix must be symmetric")
    
    ep<-eigen(mat,symmetric=TRUE);
    # Checking the positivity of the eigenvalues
    i.neg=which(ep$values < - tol)
    if (length(i.neg) > 0)
    {stop("All eigenvalues of the matrix must be greater or equal to zero")}
    
    # Transforming the small eigenvalues to zero
    i.small<-which(abs(ep$values) < tol)
    if (length(i.small)>0)
    {ep$values[i.small] = 0}
    
    valp = diag(sqrt(abs(ep$values)));
    vecp = ep$vectors;
    vecp %*% valp %*% t(vecp)
  }

wasserstein.tol <- function(mean1,var1,mean2,var2, check = FALSE, tol = 1e-8) {
    p <- length(mean1)
    d <- mean1-mean2
    vars <- var1+var2
    
    if (p == 1)
    {
      # Univariate distributions:
      if(check)
      {if(abs(var1) < .Machine$double.eps | abs(var2) < .Machine$double.eps)
      {stop("At least one variance is zero")
      }
      }
      return(sqrt( d^2 + var1 + var2 - 2*sqrt(var1*var2) ))
    } else
    {
      # Multivariate distributions:
      if(check)
      {
        if(abs(det(var1)) < .Machine$double.eps | abs(det(var2)) < .Machine$double.eps)
        {
          stop("One of the sample variances is degenerate")
        }
      }
      sqrtvar2 <- sqrtmatrix.tol(var2, tol)
      sqrtvars <- sqrtmatrix.tol(sqrtvar2%*%var1%*%sqrtvar2, tol)
      tracevar <- sum(diag(vars - 2*sqrtvars))
      
      return( sqrt( sum(d^2) + tracevar ) )
    }
  }

