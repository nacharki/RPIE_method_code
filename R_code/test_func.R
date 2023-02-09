# First industrial case : STEEL COLUMN or SULFUR MODEL or BareHole or WinWeight function

zhou98 <- function(xx) # Zhou function, d = 10
{
  d <- length(xx)
  
  xxa <- 10 * (xx-1/3)
  xxb <- 10 * (xx-2/3)
  
  norma <- sqrt(sum(xxa^2))
  normb <- sqrt(sum(xxb^2))
  
  phi1 <- (2*pi)^(-d/2) * exp(-0.5*(norma^2))
  phi2 <- (2*pi)^(-d/2) * exp(-0.5*(normb^2))
  
  y <- (10^d)/2 * (phi1 + phi2)
  return(y)
}

loepetal13 <- function(xx) # d = 10
{
  
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  x5 <- xx[5]
  x6 <- xx[6]
  x7 <- xx[7]
  
  term1 <- 6*x1 + 4*x2
  term2 <- 5.5*x3 + 3*x1*x2
  term3 <- 2.2*x1*x3 + 1.4*x2*x3
  term4 <- x4 + 0.5*x5
  term5 <- 0.2*x6 + 0.1*x7
  
  y <- term1 + term2 + term3 + term4 + term5
  return(y)
}


soblev99 <- function(xx, b, c0=0) # d = 10
{
  # xx = c(x1, x2, ..., xd)
  # b  = d-dimensional vector (optional), with default value
  #      c(2, 1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.65, 0.4228, 0.3077, 0.2169,
  #        0.1471, 0.0951, 0.0577, 0.0323, 0.0161, 0.0068, 0.0021, 0.0004, 0)
  #      (when d<=20)
  # c0 = constant term (optional), with default value 0  

  
  d <- length(xx)
  
  if (missing(b)) {
    if (d <= 20){
      b <- c(2, 1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.65, 0.4228, 0.3077, 0.2169, 0.1471, 0.0951, 0.0577, 0.0323, 0.0161, 0.0068, 0.0021, 0.0004, 0)
    }
    else {
      stop('Value of the d-dimensional vector b is required.')
    }
  }
  
  Id <- 1
  for (ii in 1:d) {
    bi  <- b[ii]
    new <- (exp(bi)-1) / bi
    Id  <- Id * new
  }
  
  sum <- 0
  for (ii in 1:d) {
    bi  <- b[ii]
    xi  <- xx[ii]
    sum <- sum + bi*xi
  }
  
  y <- exp(sum) - Id + c0
  return(y)
}

bratleyetal92 <- function(xx) # d = 10
{
  d <- length(xx)
  ii <- c(1:d)
  
  xxmat <- matrix(rep(xx,times=d), d, d, byrow=TRUE)
  xxmatlow <- xxmat
  xxmatlow[upper.tri(xxmatlow)] <- 1
  
  prod <- apply(xxmatlow, 1, prod)
  xxmatlow[upper.tri(xxmatlow)] <- 0
  sum <- sum(prod*(-1)^ii)
  
  y <- sum
  return(y)
}

roosarn63 <- function(xx) # d=10
{
  prod <- prod(abs(4*xx - 2))
  
  y <- prod
  return(y)
}

wingweight <- function(xx)
{

  Sw      <- xx[1]
  Wfw     <- xx[2]
  A       <- xx[3]
  LamCaps <- xx[4] * (pi/180)
  q       <- xx[5]
  lam     <- xx[6]
  tc      <- xx[7]
  Nz      <- xx[8]
  Wdg     <- xx[9]
  Wp      <- xx[10]
  
  fact1 <- 0.036 * Sw^0.758 * Wfw^0.0035
  fact2 <- (A / ((cos(LamCaps))^2))^0.6
  fact3 <- q^0.006 * lam^0.04
  fact4 <- (100*tc / cos(LamCaps))^(-0.3)
  fact5 <- (Nz*Wdg)^0.49
  
  term1 <- Sw * Wp
  
  y <- fact1*fact2*fact3*fact4*fact5 + term1
  return(y)
}

welchetal92 <- function(xx)
{
  
  x1  <- xx[1]
  x2  <- xx[2]
  x3  <- xx[3]
  x4  <- xx[4]
  x5  <- xx[5]
  x6  <- xx[6]
  x7  <- xx[7]
  x8  <- xx[8]
  x9  <- xx[9]
  x10 <- xx[10]
  x11 <- xx[11]
  x12 <- xx[12]
  x13 <- xx[13]
  x14 <- xx[14]
  x15 <- xx[15]
  x16 <- xx[16]
  x17 <- xx[17]
  x18 <- xx[18]
  x19 <- xx[19]
  x20 <- xx[20]
  
  term1 <- 5*x12 / (1+x1)
  term2 <- 5 * (x4-x20)^2
  term3 <- x5 + 40*x19^3 - 5*x19
  term4 <- 0.05*x2 + 0.08*x3 - 0.03*x6
  term5 <- 0.03*x7 - 0.09*x9 - 0.01*x10
  term6 <- -0.07*x11 + 0.25*x13^2 - 0.04*x14
  term7 <- 0.06*x15 - 0.01*x17 - 0.03*x18
  
  y <- term1 + term2 + term3 + term4 + term5 + term6 + term7
  return(y)
}

detpep108d <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  ii <- c(4:8)
  
  term1 <- 4 * (x1 - 2 + 8*x2 - 8*x2^2)^2
  term2 <- (3 - 4*x2)^2
  term3 <- 16 * sqrt(x3+1) * (2*x3-1)^2
  
  xxmat <- matrix(rep(xx[3:8],times=6), 6, 6, byrow=TRUE)
  xxmatlow <- xxmat
  xxmatlow[upper.tri(xxmatlow)] <- 0
  
  inner <- rowSums(xxmatlow)
  inner <- inner[2:6]
  outer <- sum(ii*log(1+inner))
  
  y <- term1 + term2 + term3 + outer
  return(y)
}


grlee09 <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  x5 <- xx[5]
  x6 <- xx[6]
  
  term1 <- exp(sin((0.9*(x1+0.48))^10))
  term2 <- x2 * x3
  term3 <- x4
  
  y <- term1 + term2 + term3
  return(y)
}

gfunc <- function(xx, a=(c(1:length(xx))-1)/2)
{
  ##########################################################################
  #
  # G-FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # a = c(a1, a2, ..., ad) (optional), with default values given by Crestaux
  #     et al. (2007)  
  #
  ##########################################################################
  
  new1 <- abs(4*xx-2) + a
  new2 <- 1 + a
  prod <- prod(new1/new2)
  
  y <- prod
  return(y)
}

morcaf95a <- function(xx)
{
  ##########################################################################
  #
  # MOROKOFF & CAFLISCH (1995) FUNCTION 1
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, ..., xd)
  #
  ##########################################################################
  
  d <- length(xx)
  fact1 <- (1 + 1/d)^d
  
  prod <- prod(xx^(1/d))
  
  y <- fact1 * prod
  return(y)
}

