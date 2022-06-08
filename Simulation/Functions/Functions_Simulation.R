#### 1. GenerateData: Generate covariate and response data ####
# Function to generate covariate and response data 
# Input: 
# - TrueBeta: (P x Q) array with true parameters (regression coefficients)
# - n: Scalar sample size
# - pXD: Scalar probability that discrete covariate = 1
# - MuX: Scalar mean of continuous covariate
# - SigmaX: Scalar standard deviation of continuous covariate 
# - RangeX: Vector with values of interval of discrete variables or (2) vector with limits of low interval of continuous covariate 
# - Covariate: logical. TRUE if covariate is included
# - MeasurementLevel: "Continuous" or "Discrete"
# - GenerateY: Logical. If TRUE, response data are generated# - Interaction: logical. TRUE if interaction between covariate and treatment is included
# - GenerateMult: Logical. Generates multinomial response data if TRUE; Currently no other options supported.

# Output:
# - Data: 
# -- List with:
# -- X: n x P design matrix (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# -- yMult: n x Q (if multinomial) response data. Only if GenerateY = TRUE.
GenerateData <- function(TrueBeta = NULL, TrueSigma = NULL, n, pXD = pXD, MuX = MuX, SigmaX = SigmaX, Covariate = Covariate, MeasurementLevel = MeasurementLevel, Interaction = Interaction, RangeX = RangeX, GenerateY = TRUE, GenerateMult = GenerateMult){
  
  # Generate covariate data
  xTrt <- rep(0:1,n)
  X <- cbind(rep(1, 2*n), xTrt)
  if(Covariate){
    if(MeasurementLevel == "Continuous"){
      x1 <- msm::rtnorm(2*n, MuX, SigmaX, lower = min(RangeX), upper = max(RangeX))
    } else if (MeasurementLevel == "Discrete") {
      x1 <- apply(rmultinom(2*n, 1, pXD), 2, function(x) which(x == 1) - 1)}
    
    if(Interaction){
      X <- cbind(rep(1,2*n), xTrt, x1, xTrt * x1)
    } else {
      X <- cbind(rep(1,2*n), xTrt, x1)
    }
  } 
  
  # Generate response data
  if(GenerateY){
    if(GenerateMult){
    Phi <- exp(X %*% TrueBeta[1:ncol(X),]) / rowSums(exp(X %*% TrueBeta[1:ncol(X),]))
    yMult <- t(apply(Phi, 1, function(x) rmultinom(1, 1, x)))
    yMV <- Multinomial2Multivariate(yMult)
    }
    
    return(list(X = X, yMult = yMult))
  } else {
    return(X)
  }
}

#### 2. SamplePostBeta: Generate data and obtain regression coefficients ####
# Input:
# - nSim: Scalar number of simulations
# - TrueBeta: (P x Q) array with true parameters (regression coefficients)
# - n: Scalar sample size
# - pXD: Scalar probability that discrete covariate = 1
# - MuX: Scalar mean of continuous covariate
# - SigmaX: Scalar standard deviation of continuous covariate 
# - RangeX: Vector with values of interval of discrete variables or (2) vector with limits of low interval of continuous covariate 
# - Covariate: logical. TRUE if covariate is included
# - MeasurementLevel: "Continuous" or "Discrete"
# - Interaction: logical. TRUE if interaction between covariate and treatment is included
# - nBurn: number of burnin iterations
# - nIt: number of MCMC-iterations
# - Start: vector starting values of chains. Each chain can have one starting value.
# - b0: scalar. Prior mean.
# - B0: scalar. Prior variance of regression coefficients
# - Seed: Scalar.
# - Method: "PG" for Polya-Gamma sampling. Other options not supported yet.
# - GenerateMult: Logical. Generates multinomial response data if TRUE; generates bivariate binomial data if FALSE.
# - Parallel: Logical. Runs function in parallel if TRUE and sequential if FALSE.

# Output:
# A nested list with:
# - Data: 
# -- List with:
# -- X: n x P design matrix (Intercept, Treatment indicator, Covariate, Covariate x Treatment)
# -- yMult: n x Q (if multinomial) response data. Only if GenerateY = TRUE.
# - Parameters:
# -- List with: 
# -- 
SamplePostBeta <- function(nSim, TrueBeta, n, pXD = NULL, MuX = NULL, SigmaX = NULL, RangeX,
                           Covariate = TRUE, MeasurementLevel, Interaction = TRUE, nBurn, nIt, Start, b0, B0, Seed, Method = "PG", GenerateMult = TRUE, Parallel = TRUE){
 
  set.seed(Seed)
  if(Parallel){
  Res <- foreach(sim = 1:nSim,
                 .export = ls(globalenv()),
                 .packages = packages)%dopar%{ #('coda', 'pgdraw', 'TruncatedNormal', 'MCMCpack', 'abind', 'msm', 'nnet')) %dopar% {
                   Data <- GenerateData(TrueBeta = TrueBeta, TrueSigma = NULL, n = n, pXD = pXD, MuX = MuX, SigmaX = SigmaX, Covariate = Covariate, MeasurementLevel = MeasurementLevel, RangeX = RangeX, Interaction = Interaction, GenerateMult = GenerateMult)
                   if("PG" %in% Method){ParametersPG <- EstimateParametersPG(X = Data$X, Y = Data$yMult, nBurn = nBurn, nIt = nIt, Start = Start, bMu0 = b0, bSigma0 = B0)}
                     c(list(Data = Data), sapply(Method, function(x) get0(paste0("Parameters", x, collapse="")), simplify = FALSE, USE.NAMES = TRUE))
                                  }
} else {
  Res <- vector("list", nSim)
  for(sim in 1:nSim){  
     Data <- GenerateData(TrueBeta = TrueBeta, TrueSigma = NULL, n = n, pXD = pXD, MuX = MuX, SigmaX = SigmaX, Covariate = Covariate, MeasurementLevel = MeasurementLevel, Interaction = Interaction, RangeX = RangeX, GenerateMult = GenerateMult)
     if("PG" %in% Method){ParametersPG <- EstimateParametersPG(X = Data$X, Y = Data$yMult, nBurn = nBurn, nIt = nIt, Start = Start, bMu0 = b0, bSigma0 = B0)}
 
     Res[[sim]] <- c(list(Data = Data), sapply(Method, function(x) get0(paste0("Parameters", x, collapse="")), simplify = FALSE, USE.NAMES = TRUE))
        }
}
  return(Res)
}

#### 3. EstimateParametersPG: Wrapper for MCMC procedure to estimate regression coefficients - Polya-Gamma ####
# Function to estimate regression coefficients using Polya-Gamma gibbs sampling. 
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
EstimateParametersPG <- function(X, Y, nBurn, nIt, Start, bMu0, bSigma0){
tryCatch({
  P <- ncol(X)
  Q <- ncol(Y)
  
  Chain1PG <- SampleBetaPG(X=X, Y=Y, nBurn=nBurn, nIt=nIt, Start=Start[1], bMu0=bMu0, bSigma0=bSigma0)
  Chain2PG <- SampleBetaPG(X=X, Y=Y, nBurn=nBurn, nIt=nIt, Start=Start[2], bMu0=bMu0, bSigma0=bSigma0)
  
  ChainsPG <- mcmc.list(mcmc(do.call(rbind, lapply(Chain1PG[["bDrawPG"]], function(i) as.vector(i[,-seq(Q, Q * P, Q)])))),
                        mcmc(do.call(rbind, lapply(Chain2PG[["bDrawPG"]], function(i) as.vector(i[,-seq(Q, Q * P, Q)])))))
  
  ConvergencePG <- gelman.diag(ChainsPG)$mpsrf
  bDrawPG <- Chain1PG[["bDrawPG"]]  
  res <- list(ConvergencePG = ConvergencePG, bDrawPG = bDrawPG)
},
warning = function(warn) {
  print(paste("Warning: ", warn))
  res <- list(NULL, NULL)
},
error = function(err) {
  print(paste("Error: ", err))
  res <- list(NULL,NULL)
},
finally = function(f) {
  return(res)
})

}

#### 4. SampleBetaPG: Sample regression coefficients via Polya-Gamma logistic regression ####
# Function to sample regression coefficients via Polya-Gamma multinomial logistic regression
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
# - bDrawPG: List of length nIt with P x Q matrices with estimated regression coefficients

SampleBetaPG <- function(X, Y, nBurn, nIt, Start, bMu0, bSigma0){
  
  P <- ncol(X);
  n <- nrow(Y);
  Q <- ncol(Y);
  
  bDrawPG <- list("vector", nIt);
  betaPG <- array(Start, dim = c(P,Q));
  betaPG[,Q] <- rep(0,P);
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
      omegaDraw <- pgdraw(1,eta)
      
      # Draw regression coefficients
      bSigma <- chol2inv(chol(t(X) %*% (X * omegaDraw) + B0))
      bMu <- bSigma %*% (t(X) %*% (kappa[,q] + omegaDraw * C) + P0)
      betaPG[,q] <- bMu + t(chol(bSigma)) %*% rnorm(P) 
    } 
    if(i > nBurn){bDrawPG[[i-nBurn]] <- betaPG}
  }
  
  return(list(bDrawPG = bDrawPG))
}


###################################################################################
#            FUNCTIONS TRANSFORMATION: POSTERIOR PREDICTIVE PROBABILITIES         #
###################################################################################
#### 5. SamplePopulation: Function to (sub)sample covariate data ####
# Input: 
# - X: n x P design matrix with covariate data
# - Method: "Empirical" for empirical marginalization, "Value" for vector of fixed values, "MvB" for unconditional multivariate Bernoulli (reference approach)
# - Values: Scalar. Value of x representing subpopulation.
# - Range: Vector of lower and upper bound of covariate that represents subpopulation.

# Output: 
# - xE: nE x P design matrix with covariate data for treatment E
# - xC: nC x P design matrix with covariate data for treatment C
SamplePopulation <- function(X, Method, Values, Range){

  if("Value" %in% Method){
      xE <- cbind(1, 1, Values, Values)
      xC <- cbind(1, 0, Values, 0)
     } else if(any(c("MvB", "Empirical") %in% Method)){
        xE <- X[X[,2] == 1 & X[,3] >= min(Range) & X[,3] <= max(Range),]
        xC <- X[X[,2] == 0 & X[,3] >= min(Range) & X[,3] <= max(Range),]
         }
  
  return(list(xE = xE, xC = xC))
}

#### 6. EstimateThetaMvB: Function to estimate theta with Multivariate Bernoulli distribution ####
# Function to estimate success probailities with a Multivariate Bernoulli distribution (reference approach)
# Input: 
# - Y: n x Q matrix with multinomial responses.
# - PriorAlpha: Vector of length Q with prior hyperparameters.
# - nIt: Scalar. Number of draws from posterior distributions.

# Output:
# - mTheta: nIt x K matrix with bivariate Bernoulli probabilities. Currently supported for K=2 only.

EstimateThetaMvB <- function(Y, PriorAlpha, nIt){
  mPhi <- MCMCpack::rdirichlet(nIt, colSums(Y) + PriorAlpha)
  mTheta <- cbind(rowSums(mPhi[,c(1,2)]), rowSums(mPhi[,c(1,3)]))
  return(list(MVB = mTheta))
}

#### 7. EstimateThetaEmpirical: Function to estimate theta via empirical marginalization ####
# Function to estimate success probabilities via empirical marginalization
# Input:
# - BetaPG: P x Q x nIt array of posterior draws of regression coefficients
# - X: n x P matrix with covariate data

# Output:
# - mTheta: nIt x K matrix with bivariate probabilities. Currently supported for K=2 only.

EstimateThetaEmpirical <- function(BetaPG, X){
  nIt <- length(BetaPG)
 
    mPhi <- lapply(1:nIt, function(i) colMeans(exp(X %*% BetaPG[[i]]) / rowSums(exp(X %*% BetaPG[[i]]))))
    mThetaPG <- lapply(1:nIt, function(i) c(sum(mPhi[[i]][c(1,2)]), sum(mPhi[[i]][c(1,3)])))
 
  return(mTheta = mThetaPG)
}

#### 8. IntegrandPG: Function to integrate in numerical marginalization ####
# Integrand function for integration over a range of a covariate in numerical marginalization. 
# Input:
# - x: value of covariate x
# - Beta: P x Q matrix of regression coefficients
# - MuX: Mean of distribution of covariate
# - SigmaX: Standard deviation of distribution of covariate
# - RangeX: Vector of lower and upper bound of range to integrate over.
# - Trt: Scalar. Value of treatment indicator
# - q: Scalar. Response category in 1 to Q

# Output:
# - Phi[,q]: Scalar. Joint response of response category q.

IntegrandPG <- function(x, Beta, MuX, SigmaX, RangeX, Trt, q){
  xInt <- cbind(1,Trt,x,Trt*x)
  Psi <- xInt %*% Beta
  pX <- msm::dtnorm(x, mean = MuX, sd = SigmaX, lower = RangeX[1], upper = RangeX[2], log=TRUE)
  Phi <- exp(Psi - log(rowSums(exp(Psi))) + pX)
  return(Phi[,q])
}

#### 9. EstimateThetaAnalytical: Function to integrate over a range of x ####
# Function to integrate over a range of x
# Input: 
# - BetaPG: P x Q x nIt array of posterior regression coefficients
# - X: n x P matrix of covariate data, where the covariate of interest is in the third column
# - Trt: Scalar. Value of treatment indicator.
# - RangeX: Vector of lower and upper bound of range to integrate over.

# Output: 
# - mTheta: nIt x K matrix with bivariate probabilities. Currently supported for K=2 only.

EstimateThetaAnalytical <- function(BetaPG, X, Trt, RangeX){
  nIt <- length(BetaPG)
    
    mThetaPG <- vector("list", nIt)
    mPhi <- array(NA, dim = c(1, ncol(BetaPG[[1]])))
    for(i in 1:nIt){
      for(q in 1:(ncol(BetaPG[[1]]) - 1)){
        mPhi[,q] <- integrate(IntegrandPG, lower = RangeX[1], upper = RangeX[2],
                              Beta = BetaPG[[i]], MuX = mean(X[,3]), SigmaX = sd(X[,3]), RangeX = RangeX, Trt = Trt, q = q)$value
      }
      mThetaPG[[i]] <- c(rowSums(mPhi[,c(1,2),drop=FALSE]), rowSums(mPhi[,c(1,3),drop=FALSE]))
    }
    
  return(mTheta = mThetaPG)
}

#### 10. Transform2Theta: Function to compute theta from regression coefficients via various methods ####
# Function to compute theta from regression coefficients via various methods.
# Input: 
# - BetaDrawPG: P x Q x nIt array with posterior regression coefficients
# - X: n x P matrix with covariate data
# - Y: n x Q matrix with multinomial response data
# - Types: R x 4 matrix with R conditions, including columns with 
#	-- MeasurementLevels ("Discrete", "Continuous"),
#	-- Methods ("Value" for vector of fixed values, "Empirical" for empirical marginalization, "Analytical" for numerical marginalization, "MvB" for Multivariate Bernoulli (unconditional), 
#	-- Populations ("Trial" for study population, "Intra_Lo" for small subpopulation within trial), 
#	--Scopes ("Old" for existing covariate data, "New" for a new sample of covariate data to be drawn within the function.)
# - Ranges: List of S ranges that defines the population of interest by a vector of a lower and an upper bound. 
# - Values: List of S values that defines the population of interest by a scalar value.
# - PriorAlpha: Vector of length Q with prior hyperparameters. Currently supported for Q = 4 only.
# - MeasurementLevel: String. "Concrete" or "Discrete"

# Output:
# - List of:
# -- mTheta.E: nIt x K matrix with multivariate probabilities for experimental treatment (T=1). Currently supported for K=2 only.
# -- mTheta.C: nIt x K matrix with multivariate probabilities for control treatment (T=0). Currently supported for K=2 only.

Transform2Theta <- function(BetaDrawPG, X, Y = NULL, Types, Ranges, Values, PriorAlpha, MeasurementLevel){
  IndicesML <- which(Types[,"MeasurementLevels"] == MeasurementLevel)
 
  Empirical <- which(apply(Types, 1, function(x) any(c("Empirical", "Value") %in% x["Methods"])))
  Analytical <- which(apply(Types, 1, function(x) c("Analytical") %in% x["Methods"]))
  MvB <- which(apply(Types, 1, function(x) c("MvB") %in% x["Methods"]))
  
  xEmpirical <- lapply(1:length(IndicesML), function(x) vector("list", 2))
  mTheta.E <- mTheta.C <- vector("list", length(IndicesML))
  
  for(i in IndicesML){
    
    if(i %in% Empirical){
      xEmpirical[[which(IndicesML == i)]] <- SamplePopulation(X = X, Method = Types[i,"Methods"], 
                                          Values = Values[[Types[i,"Populations"]]], Range = Ranges[[Types[i,"Populations"]]])
      mTheta.E[[which(IndicesML == i)]] <- EstimateThetaEmpirical(BetaPG = BetaDrawPG, X = xEmpirical[[which(IndicesML == i)]][["xE"]])
      mTheta.C[[which(IndicesML == i)]] <- EstimateThetaEmpirical(BetaPG = BetaDrawPG, X = xEmpirical[[which(IndicesML == i)]][["xC"]])
      
    } else if(i %in% Analytical){
      mTheta.E[[which(IndicesML == i)]] <- EstimateThetaAnalytical(BetaPG = BetaDrawPG, X = X, Trt = 1, RangeX = Ranges[[Types[i,"Populations"]]])
      mTheta.C[[which(IndicesML == i)]] <- EstimateThetaAnalytical(BetaPG = BetaDrawPG, X = X, Trt = 0, RangeX = Ranges[[Types[i,"Populations"]]])
      
    } else if(i %in% MvB){
      if(Types[i,"MeasurementLevels"] == "Discrete"){
        Indices <- which(X[,3] %in% Ranges[[Types[i,"Populations"]]])
      }
      else if(Types[i,"MeasurementLevels"] == "Continuous"){
        Indices <- which(X[,3] > min(Ranges[[Types[i,"Populations"]]]) & X[,3] < max(Ranges[[Types[i,"Populations"]]]))
      }
      mTheta.E[[which(IndicesML == i)]] <- EstimateThetaMvB(Y[intersect(Indices, which(X[,2] == 1)),,drop=FALSE], PriorAlpha, nIt = nIt)
      mTheta.C[[which(IndicesML == i)]] <- EstimateThetaMvB(Y[intersect(Indices, which(X[,2] == 0)),,drop=FALSE], PriorAlpha, nIt = nIt)
      
    }
  }
  
  return(list(mTheta.E = mTheta.E, mTheta.C = mTheta.C))
}

###################################################################################
#            FUNCTIONS EVALUATION: BIAS AND POSTERIOR PROBABILITIES               #
###################################################################################
#### 11. EvaluateData: Function to evaluate data ####
# Function to compute bias and posterior probabilities 
# Input:
# - Data: nIt x K matrix with treatment differences (E - C). Currently supported for K=2 only. 
# - Weights: Vector of length K with weights for linear combination of treatment differences (Compensatory rule). Currently supported for K=2 only.
# - Alpha: Scalar. Desired Type I error rate. 
# - Alternative: String. Type of decision. "greater.than" for a right-sided test and "two.sided" for a two-sided test. 
# - Rule: String. Options: "All", "Any", "Compensatory" 
# - Truth: List of length 2 with true values: "Delta" is a vector of multivariate treatment differences and "DeltaW" is a scalar of a weighted treatment difference.
# - Types: R x 4 matrix with R conditions, including columns with 
#	-- MeasurementLevels ("Discrete", "Continuous"),
#	-- Methods ("Value" for vector of fixed values, "Empirical" for empirical marginalization, "Analytical" for numerical marginalization, "MvB" for Multivariate Bernoulli (unconditional), 
#	-- Populations ("Trial" for study population, "Intra_Lo" for small subpopulation within trial), 
#	-- Scopes ("Old" for existing covariate data. Currently no other option supported.)


# Output: 
# - Res: List of:
# -- Bias: List of length 2 with one vector of bias on individual treatment differences and one scalar of bias on the weighted linear combination.
# -- Pop: List of length 2 with one vector of posterior probabilities on individual treatment differences and one scalar of a posterior probability on the weighted linear combination. 
# -- Decision: Logical vector of S decisions.

EvaluateData <- function(Data, Weights, Alpha, Alternative = "greater.than", Rule, Truth, Types){
  DecisionTypes <- c("DecisionAny.RS", "DecisionAll.RS", "DecisionCompensatory.RS")
  BiasTypes <- c("BiasMultivariate", "BiasWeighted")
  PopTypes <- c("PopMultivariate", "PopWeighted")

 if(any(c("Any", "All") %in% Rule)){
    BiasMultivariate <- t(apply(Data, 1, function(x) x - Truth[["Delta"]]))
    PopMultivariate <- colMeans(Data > 0)
    
    if("Any" %in% Rule){
      if("greater.than" %in% Alternative){DecisionAny.RS <- max(PopMultivariate) > (1 - Alpha / length(PopMultivariate))} 
      if("two.sided" %in% Alternative){DecisionAny.TS <- max(PopMultivariate) > (1 - Alpha / (length(PopMultivariate) * 2)) | min(PopMultivariate) < (Alpha / (length(PopMultivariate) * 2))}
    }
    
    if("All" %in% Rule){
      if("greater.than" %in% Alternative){DecisionAll.RS <- min(PopMultivariate) > (1 - Alpha)} 
      if("two.sided" %in% Alternative){DecisionAll.TS <- min(PopMultivariate) > (1 - Alpha / 2) | max(PopMultivariate) < (Alpha / 2)}
    }
  }
  
  if("Compensatory" %in% Rule){
    DataWeighted <- Data %*% Weights
    BiasWeighted <- as.matrix(apply(DataWeighted, 1, function(x) x - Truth[["DeltaW"]]))
    PopWeighted <- colMeans(as.matrix(DataWeighted) > 0)
    
    if("greater.than" %in% Alternative){DecisionCompensatory.RS <- PopWeighted > (1 - Alpha)} 
    if("two.sided" %in% Alternative){DecisionCompensatory.TS <- PopWeighted > (1 - Alpha / 2)| PopWeighted < (Alpha / 2)}
  }
  
  Bias <- sapply(BiasTypes, function(x){ y <- get0(x); if(is.array(y)){colMeans(y)}}, USE.NAMES = TRUE)
  Pop <- sapply(PopTypes, function(x) get0(x), USE.NAMES = TRUE)
  Decision <- sapply(DecisionTypes, function(x) get0(x), USE.NAMES = TRUE)
  
  Res <- list(Bias = Bias[!sapply(Bias,is.null)], Pop = Pop[!sapply(Pop,is.null)], Decision = Decision[!sapply(Decision,is.null)])
  
  return(Res)
}
