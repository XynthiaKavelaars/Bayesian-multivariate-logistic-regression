#### 1. ComputeRho: Compute correlation ####
# Function to compute the correlation between two outcome variables
# Input: 
# phi: Vector of four joint response probabilities. Order of response combinations: (11),(10),(01),(00)

# Output:
# Correlation between two binary proportions. Outcome 1: (11) + (10), outcome 2: (11) + (01).

ComputeRho <- function(phi){
  (phi[1] * phi[4] - phi[2] * phi[3]) / 
    sqrt(sum(phi[c(1,2)]) * sum(phi[c(1,3)]) * sum(phi[c(3,4)]) * sum(phi[c(2,4)]))
}

#### 2. Theta2Phi: Transform theta to phi ####
# Function to transform a bivariate proportion to its multinomial proportions.
# Input: 
# theta: Vector of two proportions
# rho: Correlation between proportions

# Output: 
# phi: Vector of four multinomial proportions
Theta2Phi <- function(theta,rho){
  phi11 <- rho * sqrt(prod(theta)*prod(1-theta)) + prod(theta)
  phi <- c(phi11, theta[1]-phi11, theta[2]-phi11,1-theta[1]-theta[2]+phi11)
  if(!all(phi>= 0)) {
    warning("One or more negative probabilities. Rho adjusted.")
    while(!all(phi >= 0)){
      if(rho < 0) rho <- rho + 0.01
      if(rho > 0) rho <- rho - 0.01
      phi11 <- rho * sqrt(prod(theta)*prod(1-theta)) + prod(theta)
      phi <- c(phi11, theta[1]-phi11, theta[2]-phi11,1-theta[1]-theta[2]+phi11)
    }
  }
  return(phi)
}

#### 3. Phi2Theta: Transform phi to theta ####
# Function to transform a vector of four multinomial proportion to a bivariate proportion.
# Input: 
# phi: Vector of four multinomial proportions

# Output: 
# theta: Vector of two proportions
# rho: Correlation between proportions
Phi2Theta <- function(phi){
  c(sum(phi[c(1,2)]), sum(phi[c(1,3)]))
}

#### 4. Theta2DeltaW: Transform success probabilities to weighted treatment difference (E-C) ####
# Input:
# thetaE: Vector of two success probabilities of treatment E
# thetaC: Vector of two success probabilities of treatment C
# Weights: Vector of two weights of outcomes 1 and 2

# Output:
# deltaW: Weighted treatment difference (E-C)

Theta2DeltaW <- function(thetaE, thetaC, weights){
  deltaW <- (thetaE - thetaC) %*% weights
  return(deltaW)}

#### 5. Phi2DeltaW: Transform joint response probabilities to weighted treatment difference (E-C)####
# Input:
# phiE: Vector of four joint response probabilities of treatment E
# phiC: Vector of four joint response probabilities of treatment C
# Weights: Vector of two weights of outcomes 1 and 2

# Output:
# deltaW: Weighted treatment difference (E-C)

Phi2DeltaW <- function(phiE, phiC, weights){
  delta <- phiE - phiC
  deltaTheta <- cbind(rowSums(delta[,c(1,2)]), rowSums(delta[,c(1,3)]))
  deltaW <- deltaTheta %*% weights
  return(deltaW)
}

#### 6. FindTrueBeta: Transform phi to regression coefficients ####
# Function to compute regression coefficients from success probabilities for a model with treatment indicator, covariate, and interaction between treatment and covariate.
# Input:
# TrueThetaC_Lo: Vector of two success probabilities of treatment C for a population with low value of x
# TrueThetaE_Lo: Vector of two success probabilities of treatment E for a population with low value of x
# TrueThetaC_Hi: Vector of two success probabilities of treatment C for a population with high value of x
# TrueThetaE_Hi: Vector of two success probabilities of treatment E for a population with high value of x
# TrueRho_Lo: Scalar correlation between outcomes for a population with low value of x
# TrueRho_Hi: Scalar correlation between outcomes for a population with high value of x
# xLo: Scalar: low value of x
# xHi: Scalar: high value of x

# Output: 
# TrueBeta: 4 x 4 matrix with regression coefficients per joint response category (in columns)
FindTrueBeta <- function(TrueThetaC_Lo, TrueThetaE_Lo, TrueThetaC_Hi, TrueThetaE_Hi, TrueRho_Lo, TrueRho_Hi, xLo, xHi){
  PhiC_Lo <- Theta2Phi(TrueThetaC_Lo,TrueRho_Lo)
  PhiC_Hi <- Theta2Phi(TrueThetaC_Hi,TrueRho_Hi)
  PhiE_Lo <- Theta2Phi(TrueThetaE_Lo,TrueRho_Lo)
  PhiE_Hi <- Theta2Phi(TrueThetaE_Hi,TrueRho_Hi)
  
  PsiC_Lo <- log(PhiC_Lo) - log(PhiC_Lo[4])
  PsiC_Hi <- log(PhiC_Hi) - log(PhiC_Hi[4])
  PsiE_Lo <- log(PhiE_Lo) - log(PhiE_Lo[4])
  PsiE_Hi <- log(PhiE_Hi) - log(PhiE_Hi[4])
  
  Beta0 <- (xHi * PsiC_Lo - xLo * PsiC_Hi) / (xHi - xLo)
  Beta1 <- (xLo * (PsiC_Hi - PsiE_Hi) + xHi * (PsiE_Lo - PsiC_Lo)) / (xHi - xLo)
  Beta2 <- (PsiC_Hi - PsiC_Lo) / (xHi - xLo)
  Beta3 <- ((PsiE_Hi - PsiC_Hi) - (PsiE_Lo - PsiC_Lo)) / (xHi - xLo)
  
  TrueBeta <- rbind(Beta0, Beta1, Beta2, Beta3)
  return(TrueBeta)
}

#### 7. eTrueValues: True parameters over range of covariate ####
# Function to perform numerical integration with a single covariate
# Input: 
# TrueBeta: P x Q matrix with true regression coefficients 
# pXD: Scalar (0 =< pXD =< 1) with proportion of discrete x
# MuX: Scalar with mean of continuous x
# SigmaX: Scalar with standard deviation of continuous x
# RangeX: Vector with lower and upper bound of integration
# Continuous: Logical. TRUE if x is continuous; FALSE if x is discrete.
# Weights: Vector of two positively valued weights of weighted linear combination. Weights sum to 1.

# Output: 
# ThetaE: Vector of two success probabilities of treatment E
# ThetaC: Vector of two success probabilities of treatment C
# RhoE: Scalar correlation of treatment E
# RhoC: Scalar correlation of treatment C
# PhiE: Vector of four joint response probabilities of treatment E
# PhiC: Vector of four joint response probabilities of treatment C
# Delta: Vector of two treatment differences (E-C)
# DeltaW: Scalar weighted treatment difference (E-C)
eTrueValues <- function(TrueBeta, pXD = NULL, MuX = NULL, SigmaX = NULL, RangeX, Continuous, Weights){
  
  if(Continuous){
    integrand <- function(x, bs, MuX, SigmaX, RangeX, Trt, q){
      X <- cbind(1,Trt,x,Trt*x)
      Psi <- X %*% bs
      pX <- msm::dtnorm(x, mean = MuX, sd = SigmaX, lower = RangeX[1], upper = RangeX[2], log=TRUE)
      pXY <- Psi -  log(rowSums(exp(Psi))) + pX
      pXYq <- exp(pXY[,q])
      
      return(pXYq)
    }
    PhiE <- PhiC <- rep(NA, dim(TrueBeta)[2])
    for(q in 1:(dim(TrueBeta)[2])){
      PhiE[q] <- integrate(integrand, lower = RangeX[1], upper = RangeX[2], 
                           bs = TrueBeta, Trt = 1, q = q,
                           MuX = MuX, SigmaX = SigmaX, RangeX = RangeX)$value
      PhiC[q] <- integrate(integrand, lower = RangeX[1], upper = RangeX[2], 
                           bs = TrueBeta, Trt = 0, q = q,
                           MuX = MuX, SigmaX = SigmaX, RangeX = RangeX)$value
    }
  } else {
    pXY_E <- pXY_C <- matrix(NA, nrow = length(RangeX), ncol = dim(TrueBeta)[2])
    PhiE <- PhiC <- rep(NA, dim(TrueBeta)[2])

        for(x in 1:length(RangeX)){
      XE <- cbind(1, 1, RangeX[x], RangeX[x])
      XC <- cbind(1, 0, RangeX[x], 0)
      PsiE <- XE %*% TrueBeta
      PsiC <- XC %*% TrueBeta
      pX <- log(pXD[x]) - log(sum(pXD))
      pXY_E[x,] <- PsiE -  log(rowSums(exp(PsiE))) + pX
      pXY_C[x,] <- PsiC -  log(rowSums(exp(PsiC))) + pX
    }
    
    PhiE <- colSums(exp(pXY_E))
    PhiC <- colSums(exp(pXY_C))
    
  }
  
  ThetaE <- c(sum(PhiE[c(1,2)]), sum(PhiE[c(1,3)]))
  ThetaC <- c(sum(PhiC[c(1,2)]), sum(PhiC[c(1,3)]))
  RhoE <- ComputeRho(PhiE)
  RhoC <- ComputeRho(PhiC)
  
  Delta <- ThetaE - ThetaC
  DeltaW <- Delta %*% weights
  
  return(list(ThetaE = ThetaE, ThetaC = ThetaC, RhoE = RhoE, RhoC = RhoC, PhiE = PhiE, PhiC = PhiC, Delta = Delta, DeltaW = DeltaW))
}

#### 8. sTrueValues: True parameters for value of covariate ####  
# Function to compute true parameters for a value of a single covariate
# Input: 
# TrueBeta: P x Q matrix with true regression coefficients 
# ValueX: Scalar value of x
# Weights: Vector of two positively valued weights of weighted linear combination. Weights sum to 1.

# Output: 
# ThetaE: Vector of two success probabilities of treatment E
# ThetaC: Vector of two success probabilities of treatment C
# RhoE: Scalar correlation of treatment E
# RhoC: Scalar correlation of treatment C
# PhiE: Vector of four joint response probabilities of treatment E
# PhiC: Vector of four joint response probabilities of treatment C
# Delta: Vector of two treatment differences (E-C)
# DeltaW: Scalar weighted treatment difference (E-C)
sTrueValues <- function(TrueBeta, ValueX, Weights){
  XE <- cbind(1, 1, ValueX, ValueX)
  XC <- cbind(1, 0, ValueX, 0)
  
  PsiE <- XE %*% TrueBeta
  PsiC <- XC %*% TrueBeta
  
  PhiE <- exp(PsiE) / rowSums(exp(PsiE))
  PhiC <- exp(PsiC) / rowSums(exp(PsiC))
  
  ThetaE <- c(rowSums(PhiE[,c(1,2), drop=FALSE]), rowSums(PhiE[,c(1,3), drop=FALSE]))
  ThetaC <- c(rowSums(PhiC[,c(1,2), drop=FALSE]), rowSums(PhiC[,c(1,3), drop=FALSE]))
  
  RhoE <- ComputeRho(PhiE)
  RhoC <- ComputeRho(PhiC)
  
  Delta <- ThetaE - ThetaC
  DeltaW <- Delta %*% weights
  
  return(list(ThetaE = ThetaE, ThetaC = ThetaC, RhoE = RhoE, RhoC = RhoC, PhiE = PhiE, PhiC = PhiC, Delta = Delta, DeltaW = DeltaW))
}

#### 9. Multinomial2Multivariate: Transform multinomial response to bivariate binomial response ####
# Function to transform a matrix with multinomial responses to bivariate binomial responses
# Input: 
# yMult: n x 4 matrix with multinomial responses

# Output: 
# yMV: n x 2 matrix with bivariate binomial responses
Multinomial2Multivariate <- function(yMult){
  Answers <- rev(expand.grid(rev(list(c(0,1), c(0,1)))))
  Answers <- Answers[nrow(Answers):1,]
  yMV <- Answers[apply(yMult, 1, function(n) which(n == 1)),]
  return(yMV)
}

#### 10. Multivariate2Multinomial: Transform bivariate binomial response to multinomial response ####
# Input: 
# yMV: n x 2 matrix with bivariate binomial responses

# Output: 
# yMult: n x 4 matrix with multinomial responses
Multivariate2Multinomial <- function(yBV){
  Answers <- rev(expand.grid(rev(list(c(0,1), c(0,1)))))
  Answers <- Answers[nrow(Answers):1,]
  yMultVec <- apply(yBV, 1, function(x) which(apply(Responses, 1, function(y) all(y == x))))
  yMult <- matrix(0, length(yMultVec), 2^K)
  for(i in 1:length(yMultVec)){yMult[i,yMultVec[i]] <- 1}
  return(yMult) 
}

#### 11. RangeRho: Compute possible range of correlation given bivariate proportions ####
# Input: 
# theta: Vector of two proportions

# Output (list): 
# lb: Lower bound of correlation
# ub: Upper bound of correlation
RangeRho <- function(theta){
  lb <- round(max(-1 * sqrt(prod(theta) / prod(1-theta)), 
                  -1 * sqrt(prod(1-theta)/prod(theta))), 3)
  ub <- round(min((theta[1]*(1-theta[2]))/(theta[2]*(1-theta[1])), 
                  (theta[2]*(1-theta[1]))/(theta[1]*(1-theta[2]))), 3)
  return(list(lb = lb, ub = ub))
}