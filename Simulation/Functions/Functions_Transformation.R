#### 1. ComputeRhoMultivariate: Compute correlation ####
# Function to compute the correlation between two outcome variables
# Input: 
# Phi: Vector of four joint response probabilities. Order of response combinations: (11),(10),(01),(00)

# Output:
# Sigma: Correlation matrix

ComputeRhoMultivariate <- function(Phi){
  indicesTheta <- IndicesTheta(length(Phi))
  Theta <- Phi2Theta(Phi, indicesTheta)
  ThetaProducts <- t(Theta) %*% Theta
  ThetaProductsInv <- t(1-Theta) %*% (1-Theta)
  indicesPhi <- IndicesPhi(length(Phi))
  PairwisePhi <- Sigma <- diag(log2(length(Phi)))
  l <- 1
  for(i in 1:log2(length(Phi))){
    for(j in 1:log2(length(Phi))){
      if(i > j){
        PairwisePhi[i,j] <- sum(Phi[indicesPhi[[l]]])
        Sigma[i,j] <- Sigma[j,i] <- (PairwisePhi[i,j] - ThetaProducts[i,j]) / sqrt(ThetaProducts[i,j] * (ThetaProductsInv[i,j]))
        l <- l + 1
      }
    }
  }
  return(Sigma)
}


#### 2. Theta2Phi: Transform K-variate binomial proportion to multinomial proportions ####
# Wrapper function to transform a k-variate proportion to its multinomial proportions. Currently supported for K = 1,2,3
# Input: 
# Theta: Vector of K binomial proportions
# Rho: KxK matrix of pairwise correlations between binomial proportions. Only needed when K > 1.
# Phi111: Multinomial proportion of trivariate co-occurrence. Only needed when K=3.

# Output: 
# Phi: Vector of 2^K multinomial proportions

Theta2Phi <- function(Theta, Rho = NULL, Phi111 = NULL){
  if(length(Theta) == 1){
    Phi <- c(Theta, 1 - Theta)
  }
  
  if(length(Theta) == 2){
    Phi <- Theta2PhiBivariate(Theta, Rho)
  }
  
  if(length(Theta) == 3){
    Phi <- Theta2PhiTrivariate(Theta, Rho, Phi111)
  }
  return(Phi)
}
#### 2.1 Theta2PhiBivariate: Transform bivariate binomial proportion to multinomial proportions ####
# Function to transform a bivariate binomial proportion to its multinomial proportions.
# Input: 
# Theta: Vector of two binomial proportions
# Rho: 2x2 matrix of correlations between binomial proportions

# Output: 
# Phi: Vector of 2^2 multinomial proportions
Theta2PhiBivariate <- function(Theta,Rho){
  Phi11 <- Rho[2,1] * sqrt(prod(Theta)*prod(1-Theta)) + prod(Theta)
  Phi <- c(Phi11, Theta[1]-Phi11, Theta[2]-Phi11,1-Theta[1]-Theta[2]+Phi11)
  
  if(any(Phi<0)){
    stop("Error. At least one parameter < 0")
  }else{
    return(Phi)
  }
 }

#### 2.2 Theta2PhiTrivariate: Transform trivariate binomial proportion to multinomial proportions ####
# Function to transform a trivariate binomial proportion to its multinomial proportions.
# Input: 
# Theta: Vector of three binomial proportions
# Rho: 3x3 matrix of correlation between binomial proportions 
# Phi111: Multinomial proportion of trivariate co-occurrence 

# Output: 
# Phi: Vector of 2^3 multinomial proportions

Theta2PhiTrivariate <- function(Theta, Rho, Phi111){
  
  PhiIJ <- matrix(NA, length(Theta),length(Theta))
  PhiII <- rep(NA, length(Theta))
  for(i in 1:length(Theta)){
    for(j in 1:length(Theta)){
      if(j != i){
        PhiIJ[i,j] <- Rho[i,j] * sqrt(Theta[i] * (1-Theta[i]) * Theta[j] * (1-Theta[j])) + Theta[i] * Theta[j] - Phi111 
      }
    }
    PhiII[i] <- Theta[i] - Phi111 - sum(PhiIJ[i,], na.rm=TRUE)
  }
  PhiRes <- rbind(cbind(PhiIJ,PhiII), c(PhiII,NA))
  Phi <- c(Phi111, c(PhiRes[lower.tri(PhiRes)]), 1 - Phi111 - sum(c(PhiRes[lower.tri(PhiRes)])))
  
  if(any(Phi<0)){
    stop("Error. At least one parameter < 0")
  }else{
    return(Phi)
  }
}

#### 3. Phi2Theta: Transform multinomial proportions to K-variate binomial proportion ####
# Function to transform a Q-variate multinomial proportion to a K-variate binomial proportion.
# Input: 
# Phi: Vector of 2^K multinomial proportions
# Indices: List of length K with indices of the multinomial proportions that need to be summed. 
# K: Number of outcome variables

# Output: 
# Theta: Vector of K binomial proportions
Phi2Theta <- function(Phi, Indices, K = length(Indices)){
  if(!is.matrix(Phi)){Phi <- matrix(Phi, 1)}
  Theta <- matrix(NA, nrow = nrow(Phi), ncol = K)
  for(k in 1:K){
    Theta[,k] <- rowSums(Phi[,Indices[[k]],drop=FALSE])
  }
  return(Theta)
}

#### 4. IndicesTheta: Find multinomial proportions that sum to binomial proportions ####
# Function to find indices of Q multinomial proportions that form K binomial proportions. 
# Input:
# Q: Number of multinomial proportions

# Output:
# indices: List of length K = log2(Q) with indices per outcome variable
IndicesTheta <- function(Q){
  K <- log2(Q)
  indices <- start.indices <- end.indices <- vector("list", K)
  end.seqs <- seq(1, length.out = 2^(K-1), by = 2)
  for(k in 1:K){
    start.indices[[k]] <- seq(1,2^K-1,2^(K-k+1))
    end.indices[[k]] <- end.seqs[1:(2^(k-1))] * 1/2^k * Q
    indices[[k]] <- unlist(lapply(1:2^(k-1), function(l) seq(start.indices[[k]][l], end.indices[[k]][l])))
  }
  return(indices)
}

#### 5. IndicesPhi: Find co-occurring multinomial proportions (for computation of pairwise correlation) ####
# Function to find indices of the multinomial proportion that reflects co-occurence per set of two binomial proportions (e.g., 011 for K = 2 and K = 3 or ). 
# Used to compute pairwise correlations. 
# Input:
# Q: Number of multinomial proportions

# Output:
# indices: List of length K = log2(Q) with indices per outcome variable
IndicesPhi <- function(R){
  K <- log2(R)
  PhiCombs <- rev(expand.grid(rev(lapply(1:K, function(k) 1:0))))
  PhiComb <- as.list(as.data.frame(t(PhiCombs)))
  PhiMat <- which(upper.tri(diag(K)), arr.ind=TRUE)
  PhiMatOrdered <- PhiMat[do.call(order, as.data.frame(PhiMat)),,drop=FALSE]
  Phi <- vector("list", nrow(PhiMatOrdered))
  if(K==1){
    Phi[[1]] <- 1:0
  }else if(K > 1){
  for(i in 1:nrow(PhiMatOrdered)){
    Phi[[i]] <- which(unlist(lapply(PhiComb, function(x) all(x[PhiMatOrdered[i,]] == 1))))
  }
  }
  return(Phi)
}

#### 6. Phi2Beta: Transform multinomial proportions to regression coefficients ####
# Function to compute regression coefficients from success probabilities for a model with treatment indicator, covariate, and interaction between treatment and covariate.
# Input:
# PhiC_Lo: Vector of Q joint response probabilities of treatment C for a population with low value of x
# PhiE_Lo: Vector of Q joint response probabilities of treatment E for a population with low value of x
# PhiC_Hi: Vector of Q joint response probabilities of treatment C for a population with high value of x
# PhiE_Hi: Vector of Q joint response probabilities of treatment E for a population with high value of x
# xLo: Scalar: low value of x
# xHi: Scalar: high value of x

# Output: 
# TrueBeta: 4 x Q matrix with regression coefficients per joint response category 

Phi2Beta <- function(PhiE_Lo, PhiC_Lo, PhiE_Hi, PhiC_Hi, xLo = 0, xHi = 1){
  if(length(PhiE_Lo) == 1){
    PhiC_Lo <- c(TrueThetaC_Lo, 1-TrueThetaC_Lo)
    PhiC_Hi <- c(TrueThetaC_Hi, 1-TrueThetaC_Hi)
    PhiE_Lo <- c(TrueThetaE_Lo, 1-TrueThetaE_Lo)
    PhiE_Hi <- c(TrueThetaE_Hi, 1-TrueThetaE_Hi)
    }

  PsiC_Lo <- log(PhiC_Lo) - log(PhiC_Lo[length(PhiC_Lo)])
  PsiC_Hi <- log(PhiC_Hi) - log(PhiC_Hi[length(PhiC_Hi)])
  PsiE_Lo <- log(PhiE_Lo) - log(PhiE_Lo[length(PhiE_Lo)])
  PsiE_Hi <- log(PhiE_Hi) - log(PhiE_Hi[length(PhiE_Hi)])
  
  
  Beta0 <- (xHi * PsiC_Lo - xLo * PsiC_Hi) / (xHi - xLo)
  Beta1 <- (xLo * (PsiC_Hi - PsiE_Hi) + xHi * (PsiE_Lo - PsiC_Lo)) / (xHi - xLo)
  Beta2 <- (PsiC_Hi - PsiC_Lo) / (xHi - xLo)
  Beta3 <- ((PsiE_Hi - PsiC_Hi) - (PsiE_Lo - PsiC_Lo)) / (xHi - xLo)
  
  TrueBeta <- rbind(Beta0, Beta1, Beta2, Beta3)
  return(TrueBeta)
}

#### 7. Psi2Phi: Transform log-odds to multinomial proportions ####
# Function to transform log-odds to multinomial proportions via logistic function.
# Input: 
# xDat: n x P design matrix
# x: P x Q matrix with regression coefficients

# Output: 
# Phi: n x Q matrix of multinomial proportions

Psi2Phi <- function(xDat, x){
  psi <- exp(xDat %*% x)
  return(psi / rowSums(psi))
}

#### 8. eTrueValues: True parameters over range of covariate ####
# Function to perform numerical integration with a single covariate to find true values when subpopulation is defined by an interval.
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
  indicesTheta <- IndicesTheta(length(PhiE))
  ThetaE <- Phi2Theta(PhiE, indicesTheta, log2(length(PhiE)))
  ThetaC <- Phi2Theta(PhiC, indicesTheta, log2(length(PhiC)))
  RhoE <- ComputeRhoMultivariate(PhiE)
  RhoC <- ComputeRhoMultivariate(PhiC)
  
  Delta <- ThetaE - ThetaC
  DeltaW <- Delta %*% Weights
  
  Out <- list(ThetaE = ThetaE, ThetaC = ThetaC, RhoE = RhoE, RhoC = RhoC, PhiE = PhiE, PhiC = PhiC, Delta = Delta, DeltaW = DeltaW)
  
  return(Out)
}

#### 9. sTrueValues: True parameters for value of covariate ####  
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
  
  indicesTheta <- IndicesTheta(ncol(PhiE))
  ThetaE <- Phi2Theta(PhiE, Indices = indicesTheta, log2(ncol(PhiE)))
  ThetaC <- Phi2Theta(PhiC, Indices = indicesTheta, log2(ncol(PhiC)))
  
  RhoE <- ComputeRhoMultivariate(PhiE)
  RhoC <- ComputeRhoMultivariate(PhiC)
  
  Delta <- ThetaE - ThetaC
  DeltaW <- Delta %*% Weights

  Out <- list(ThetaE = ThetaE, ThetaC = ThetaC, RhoE = RhoE, RhoC = RhoC, PhiE = PhiE, PhiC = PhiC, Delta = Delta, DeltaW = DeltaW)

  return(Out)
  }
