#### Functions ####
#### 1 SampleBeta PG: Estimate multinomial regression coefficients ####
# Function to estimate multinomial regression coefficients using Polya-Gamma gibbs sampling. 
# Input:
# - X: n x P design matrix
# - Y: n x Q matrix of multinomial response data
# - nBurn: Scalar. Number of burnin iterations
# - nIt: Scalar. Number of iterations
# - Start: Vector of two starting values, each used for a different chain
# - bMu0: Scalar. Prior mean.
# - bSigma0: Scalar. Prior variance of regression coefficients

# Output: 
# A list with:
# - ConvergencePG: Multivariate Gelman-Rubin statistic to asses convergence
# - bDrawPG: List of length nIt with P x Q matrices with estimated regression coefficients

SampleBetaPG <- function(X, Y, nBurn, nIt, Start, bMu0, bSigma0, filename, indicesE, indicesC, indicesK){
  
  P <- ncol(X);
  n <- nrow(Y);
  Q <- ncol(Y);
  
  #bDrawPG <- array(0, dim = c(P,Q,nIt));
  betaPG <- array(Start, dim = c(P,Q));
  betaPG[,Q] <- rep(0,P);
  thetaE <- array(NA, dim = c(length(indicesE), log2(Q), nIt))
  thetaC <- array(NA, dim = c(length(indicesE), log2(Q), nIt))
  kappa <- Y - 1/2 # Note: As opposed to the paper: kappa here is not divided by omega to improve programming efficiency.
  
  b0 <- rep(bMu0,P)
  B0 <- diag(bSigma0,P,P) # Prior precision matrix
  P0 <- B0 %*% b0;
  

  # Start Gibbs sampler
  for(i in 2:(nBurn + nIt)){
    for(q in 1:(Q-1)){
      # Compute linear predictor
      C <- log(rowSums(exp(X %*% betaPG[,-q])));
      eta <- X %*% betaPG[,q] - C;
      
      # Draw Polya-Gamma variable
      omegaDraw <- pgdraw::pgdraw(1,eta)
      
      # Draw regression coefficients
      bSigma <- chol2inv(chol(t(X) %*% (X * omegaDraw) + B0))
      bMu <- bSigma %*% (t(X) %*% (kappa[,q] + omegaDraw * C) + P0)
      betaPG[,q] <- bMu + t(chol(bSigma)) %*% rnorm(P) 
    } 
    if(i > nBurn){
      append <- ifelse(i == nBurn + 1, FALSE, TRUE)
      write.table(t(as.vector(betaPG)), file= filename, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, append=append)}
    
    phi <- Psi2Phi(xData, betaPG)
    thetaE[,,i-nBurn] <- Phi2Theta(Phi = phi[indicesE,], Indices = indicesK, K = length(indicesK))
    thetaC[,,i-nBurn] <- Phi2Theta(Phi = phi[indicesC,], Indices = indicesK, K = length(indicesK))
    
    
  }
  
  return(list(thetaE = thetaE, thetaC = thetaC))
}
