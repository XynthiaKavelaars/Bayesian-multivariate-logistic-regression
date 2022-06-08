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
ComputeSampleSizeComp <- function(phiE, phiC, w1, w2, nMax, alpha=0.05, beta=0.20, one.sided = TRUE){
  
  theta1E <- sum(phiE[c(1,2)])
  theta2E <- sum(phiE[c(1,3)])
  theta1C <- sum(phiC[c(1,2)])
  theta2C <- sum(phiC[c(1,3)])
  
  sigma2 <- function(theta1, theta2, phi, w1, w2){
    w1^2 * theta1 * (1-theta1) + 
      w2^2 * theta2 * (1-theta2) +
      2 * w1 * w2 * (phi[1] - theta1 * theta2)
  }
  
  mu <- function(theta1, theta2, w1, w2){
    w1 * theta1 + w2 * theta2
  }
  
  muE <- mu(theta1E, theta2E, w1, w2)
  muC <- mu(theta1C, theta2C, w1, w2)
  
  varE <- sigma2(theta1E, theta2E, phiE, w1, w2)
  varC <- sigma2(theta1C, theta2C, phiC, w1, w2)
  
  if(one.sided & muE - muC > 1e-2){
    n <- min((varE+varC)*((qnorm(1-alpha)+qnorm(1-beta))/(muE-muC))^2 , nMax)
  } else if (!one.sided & abs(muE - muC) > 1e-2){
    n <- min((varE+varC)*((qnorm(1-alpha/2)+qnorm(1-beta))/(muE-muC))^2, nMax)
  } else {n <- nMax}
  return(ceiling(n))
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
  
  theta1E <- sum(phiE[c(1,2)])
  theta2E <- sum(phiE[c(1,3)])
  theta1C <- sum(phiC[c(1,2)])
  theta2C <- sum(phiC[c(1,3)])
  
  var1E <- theta1E * (1-theta1E)
  var1C <- theta1C * (1-theta1C)
  var2E <- theta2E * (1-theta2E)
  var2C <- theta2C * (1-theta2C)
  
  rho12E <- (phiE[1] - theta1E * theta2E)/
    sqrt(theta1E * (1-theta1E) * 
           theta2E * (1-theta2E))
  rho12C <- (phiC[1] - theta1C * theta2C)/
    sqrt(theta1C * (1-theta1C) * 
           theta2C * (1-theta2C))
  rho.nml <- (rho12E * sqrt(var1E * var2E) + 
                rho12C * sqrt(var1C * var2C)) / 
    sqrt((var1E + var1C) * (var2E + var2C))
  
  Sigma <- matrix(c(1,rho.nml,rho.nml,1), nrow=2, ncol=2, byrow=TRUE)
  
  pwr <- rep(NA, nMax)
  #mu <- matrix(NA, nrow=nMax, ncol=2)
  
  theta <- 1/2 * c(theta1E + theta1C, theta2E + theta2C)
  delta <- c(theta1E - theta1C, theta2E - theta2C)
  
  for(n in 1:nMax){
    
    se0 <- sqrt(2/n * theta * (1-theta))
    se <- sqrt(c(var1E + var1C, var2E + var2C) / n)
    
    
      if(one.sided){ 
        z.r <- delta / se - qnorm(1 - alpha)
 
    pwr[n] <- mvtnorm::pmvnorm(lower=-Inf, upper = z.r, 0, corr=Sigma) 
  } else if (!one.sided ){
   # z.r <- (delta - se0 * qnorm(1 - alpha/2)) / se
  #  z.l <- (delta - se0 * qnorm(alpha/2)) / se
    z.r <- delta / se - qnorm(1 - alpha/2)
    z.l <- delta / se - qnorm(alpha/2)
    
    pwr[n] <- mvtnorm::pmvnorm(lower=-Inf, upper = z.r, mean = 0, corr=Sigma) + mvtnorm::pmvnorm(lower = z.l, upper = Inf, mean = 0, corr = Sigma)
  }
    
  }
  n <- min(ifelse(any(pwr>=(1-beta)), which(pwr>=(1-beta))[1], nMax), nMax) 
   return(n)
}

#### 3. COmpueSampleSizeAny: Compute sample size for any rule ####
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
  theta1E <- sum(phiE[c(1,2)])
  theta2E <- sum(phiE[c(1,3)])
  theta1C <- sum(phiC[c(1,2)])
  theta2C <- sum(phiC[c(1,3)])
  
  var1E <- theta1E * (1-theta1E)
  var1C <- theta1C * (1-theta1C)
  var2E <- theta2E * (1-theta2E)
  var2C <- theta2C * (1-theta2C)
  
  rho12E <- (phiE[1] - theta1E * theta2E)/
    sqrt(theta1E * (1-theta1E) * 
           theta2E * (1-theta2E))
  rho12C <- (phiC[1] - theta1C * theta2C)/
    sqrt(theta1C * (1-theta1C) * 
           theta2C * (1-theta2C))
  rho.nml <- (rho12E * sqrt(var1E * var2E) + 
                rho12C * sqrt(var1C * var2C)) / 
    sqrt((var1E + var1C) * (var2E + var2C))
  
  Sigma <- matrix(c(1,rho.nml,rho.nml,1), nrow=2, ncol=2, byrow=TRUE)
  
  
  pwr <- rep(NA, nMax)
    theta <- 1/2 * c(theta1E + theta1C, theta2E + theta2C)
    delta <- c(theta1E - theta1C, theta2E - theta2C)
    
    for(n in 1:nMax){
      se0 <- sqrt(2/n * theta * (1-theta))
      se <- sqrt(c(var1E + var1C, var2E + var2C) / n)
      
      if(one.sided){ 
      #z.r <- (se0 * qnorm(1 - alpha / 2) - delta) / se
      z.r <- qnorm(1 - alpha / 2) - delta / se
      pwr[n] <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=z.r, corr=Sigma)
      
    } else if (!one.sided ){
      z.l <- qnorm(alpha/4) - delta / se
      z.r <- qnorm(1 - alpha/4) - delta / se
      
      pwr[n] <- 1-mvtnorm::pmvnorm(lower=-Inf, upper=z.r, corr=Sigma) + 1-mvtnorm::pmvnorm(lower=z.l, upper=Inf, corr=Sigma)
      
      #z.l <- (delta - se0 * qnorm(alpha/2)) / se
      #z.r <- (delta - se0 * qnorm(1 - alpha/2)) / se
   
  #    pwr[n] <- mvtnorm::pmvnorm(lower=-Inf, upper = z.r, mean = 0, corr=Sigma) + mvtnorm::pmvnorm(lower = z.l, upper = Inf, mean = 0, corr = Sigma)
      
    }
  }
    
    n <- min(ifelse(any(pwr>=(1-beta)), which(pwr>=(1-beta))[1], nMax), nMax)
  return(n)
}

