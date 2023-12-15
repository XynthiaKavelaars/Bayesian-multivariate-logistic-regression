#### 1. Compensatory rule: Compute sample size for compensatory rule ####
# Function to compute sample size for a treatment difference (E-C) as a weighted linear combination of two (correlated) proportions
# Input:
# phiE: Vector with four joint response probabilities of treatment E
# phiC: Vector with four joint response probabilities of treatment C
# w1: Weight of outcome 1
# w2: Weight of outcome 2
# nMax: Maximum sample size when treatment difference (E-C) =< 0 
# alpha: Targeted Type I error rate
# beta: Targeted Type II error rate
# one.sided: TRUE for a right-sided test and FALSE for a two-sided test

# Output: 
# n: Sample size per treatment group
ComputeSampleSizeComp <- function(phiE, phiC, weights, nMax, alpha=0.05, beta=0.20, one.sided = TRUE){
   sigma2 <- function(theta, phi, weights){
    indicesPhi <- IndicesPhi(length(phi))
    WeightProducts <- outer(weights, weights, "*")
    Weights <- c(WeightProducts[lower.tri(WeightProducts)])
    ThetaProducts <- t(theta) %*% theta
    Var <- sum(weights^2 * c(theta * (1 - theta))) + 
      sum(2 * Weights * (sapply(indicesPhi, function(x) sum(phi[x])) - c(ThetaProducts[lower.tri(ThetaProducts)])))
    return(Var)
  }
  
   indicesTheta <- IndicesTheta(length(phiE))
   thetaE <- Phi2Theta(phiE, Indices = indicesTheta, K = length(indicesTheta))
   thetaC <- Phi2Theta(phiC, Indices = indicesTheta, K = length(indicesTheta))

  muE <- thetaE %*% weights
  muC <- thetaC %*% weights
  
  varE <- sigma2(thetaE, phiE, weights)
  varC <- sigma2(thetaC, phiC, weights)
  
  if(one.sided & muE - muC > 1e-2){
    n <- min((varE+varC)*((qnorm(1-alpha)+qnorm(1-beta))/(muE-muC))^2 , nMax)
  } else if (!one.sided & abs(muE - muC) > 1e-2){
    n <- min((varE+varC)*((qnorm(1-alpha/2)+qnorm(1-beta))/(muE-muC))^2, nMax)
  } else {n <- nMax}
  return(round(n))
}
#### 2. ComputeSampleSizeAll: Compute sample size for all rule ####
# Function to compute sample size for a positive treatment difference (E-C) on two (correlated) proportions
# Input:
# phiE: Vector with four joint response probabilities of treatment E
# phiC: Vector with four joint response probabilities of treatment C
# nMax: Maximum sample size when treatment difference (E-C) =< 0 
# alpha: Targeted Type I error rate
# beta: Targeted Type II error rate
# one.sided: TRUE for a right-sided test and FALSE for a two-sided test

# Output: 
# n: Sample size per treatment group

ComputeSampleSizeAll <- function(phiE, phiC, nMax, alpha=0.05, beta=0.20, one.sided = TRUE){
  indicesTheta <- IndicesTheta(length(phiE))
  thetaE <- Phi2Theta(phiE, Indices = indicesTheta, K = length(indicesTheta))
  thetaC <- Phi2Theta(phiC, Indices = indicesTheta, K = length(indicesTheta))
  
  varE <- thetaE * (1-thetaE)
  varC <- thetaC * (1-thetaC)
  
  PairwiseRhoE <- ComputeRhoMultivariate(Phi = phiE)
  PairwiseRhoC <- ComputeRhoMultivariate(Phi = phiC)
  
  Sigma <- diag(length(thetaE))
  for(i in 1:length(thetaE)){
    for(j in 1:length(thetaE)){
      if(i>j){
        Sigma[i,j] <- Sigma[j,i] <- 
          (PairwiseRhoE[i,j] * sqrt(varE[i] * varE[j]) + 
             PairwiseRhoC[i,j] * sqrt(varC[i] * varC[j])) / 
          sqrt((varE[i] + varC[i]) * (varE[j] + varC[j]))
      }
    }
  }
  
  
  pwr <- rep(NA, nMax)
  theta <- colMeans(rbind(thetaE, thetaC))
  delta <- thetaE - thetaC
  
  for(n in 1:nMax){
    se0 <- sqrt(2/n * theta * (1-theta))
    se <- sqrt((varE + varC) / n)
    
    
      if(one.sided){ 
        z.r <- c(delta / se - qnorm(1 - alpha))
 
    pwr[n] <- mvtnorm::pmvnorm(lower=-Inf, upper = z.r, 0, corr=Sigma) 
  } else if (!one.sided ){
     z.r <- c(delta / se - qnorm(1 - alpha/2))
    z.l <- c(delta / se - qnorm(alpha/2))
    
    pwr[n] <- mvtnorm::pmvnorm(lower=-Inf, upper = z.r, mean = 0, corr=Sigma) + mvtnorm::pmvnorm(lower = z.l, upper = Inf, mean = 0, corr = Sigma)
  }
    
  }
  n <- min(ifelse(any(pwr>=(1-beta)), which(pwr>=(1-beta))[1], nMax), nMax) 
   return(n)
}

#### 3. ComputeSampleSizeAny: Compute sample size for any rule ####
# Function to compute sample size for a positive treatment difference (E-C) on at least one of two (correlated) proportions
# Input:
# phiE: Vector with four joint response probabilities of treatment E
# phiC: Vector with four joint response probabilities of treatment C
# nMax: Maximum sample size when treatment difference (E-C) =< 0 
# alpha: Targeted Type I error rate. 
# beta: Targeted Type II error rate
# one.sided: TRUE for a right-sided test and FALSE for a two-sided test

# Output: 
# n: Sample size per treatment group

ComputeSampleSizeAny <- function(phiE, phiC, nMax, alpha=0.05, beta=0.20, one.sided = TRUE){
  indicesTheta <- IndicesTheta(length(phiE))
  thetaE <- Phi2Theta(phiE, Indices = indicesTheta, K = length(indicesTheta))
  thetaC <- Phi2Theta(phiC, Indices = indicesTheta, K = length(indicesTheta))

  varE <- thetaE * (1-thetaE)
  varC <- thetaC * (1-thetaC)
  
  PairwiseRhoE <- ComputeRhoMultivariate(Phi = phiE)
  PairwiseRhoC <- ComputeRhoMultivariate(Phi = phiC)
  
  Sigma <- diag(length(thetaE))
  for(i in 1:length(thetaE)){
    for(j in 1:length(thetaE)){
      if(i>j){
        Sigma[i,j] <- Sigma[j,i] <- 
          (PairwiseRhoE[i,j] * sqrt(varE[i] * varE[j]) + 
           PairwiseRhoC[i,j] * sqrt(varC[i] * varC[j])) / 
          sqrt((varE[i] + varC[i]) * (varE[j] + varC[j]))
      }
    }
  }
          
  
  pwr <- rep(NA, nMax)
    theta <- colMeans(rbind(thetaE, thetaC))
    delta <- thetaE - thetaC
    
    for(n in 1:nMax){
      se0 <- sqrt(2/n * theta * (1-theta))
      se <- sqrt((varE + varC) / n)
      
      if(one.sided){ 
      z.r <- c(qnorm(1 - alpha / length(theta)) - delta / se)
      pwr[n] <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=z.r, corr=Sigma)
      
    } else if (!one.sided ){
      z.l <- c(qnorm(alpha/(length(theta)*2)) - delta / se)
      z.r <- c(qnorm(1 - alpha/(length(theta)*2)) - delta / se)
      
      pwr[n] <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=z.r, corr=Sigma) + 1-mvtnorm::pmvnorm(lower=z.l, upper=Inf, corr=Sigma)
      
    }
  }
    
    n <- min(ifelse(any(pwr>=(1-beta)), which(pwr>=(1-beta))[1], nMax), nMax)
  return(n)
}

