####  DataGeneratingMechanisms ####
#### True success probabilities ####
TrueTheta <- list(
  list(ThetaE_Lo = c(0.35,0.45,0.35),
       ThetaC_Lo = c(0.35,0.45,0.35),
       ThetaE_Hi = c(0.65,0.55,0.65),
       ThetaC_Hi = c(0.65,0.55,0.65)),
  
  list(ThetaE_Lo = c(0.625,0.575,0.625),
       ThetaC_Lo = c(0.375,0.425,0.375),
       ThetaE_Hi = c(0.375,0.425,0.375),
       ThetaC_Hi = c(0.625,0.575,0.625)),
  
  list(ThetaE_Lo = c(0.65,0.60,0.65),
       ThetaC_Lo = c(0.35,0.40,0.35),
       ThetaE_Hi = c(0.45,0.40,0.45),
       ThetaC_Hi = c(0.55,0.60,0.55)),
  
  list(ThetaE_Lo = c(0.60,0.50,0.60),
       ThetaC_Lo = c(0.40,0.50,0.40),
       ThetaE_Hi = c(0.75,0.50,0.75),
       ThetaC_Hi <- c(0.25,0.50,0.25)))

#### True correlations ####
Rho <- c(-0.20,0.00,0.20)
names(Rho) <- c("Neg", "Zero", "Pos")

# As matrix
RhoXX <- lapply(1:length(Rho), function(r) {
  Sigma <- diag(1,max(K))
  Sigma[c(which(lower.tri(Sigma)),which(upper.tri(Sigma)))] <- Rho[r]
  Sigma}
)


#### Joint response probabilities ####
nDgm <- length(TrueTheta) * 2
TruePhi <- vector("list", length(K))
TrueBeta <- lapply(1:length(K), function(k) vector("list", nDgm))

for(k in K){
TruePhi[[k]] <- lapply(1:length(TrueTheta), function(es){
  sapply(1:length(TrueTheta[[1]]), function(subpop){
    if(k == 3){
    if(all(TrueTheta[[es]][[subpop]] <= 0.50)){
      if(max(TrueTheta[[es]][[subpop]]) - min(TrueTheta[[es]][[subpop]]) >= 0.11){
        Phi111Prop <- c(1/16,3/16,6/16)
      }else{
        Phi111Prop <- c(2/16,4/16,6/16)
      }
    }else{
       if(max(TrueTheta[[es]][[subpop]]) > 0.71){
        Phi111Prop <- c(6/16,9/16,12/16)
      }else{
        Phi111Prop <- c(3/16,7/16,9/16)
      }
    }
  }else if(k < 3){
    Phi111Prop <- NULL
  }
    lapply(1:length(RhoXX), function(r){
      results <- tryCatch({
        Theta2Phi(TrueTheta[[es]][[subpop]][1:k], RhoXX[[r]][1:k,1:k], Phi111=Phi111Prop[r] * min(TrueTheta[[es]][[subpop]][1:k]))
      },
      error = function(e){
        message <- paste0("Error in dgm", es, ", subpopulation ", subpop, ", rho ", r, ". Error message: ", e$message)
        return(message)
      })
    })
  }, simplify=FALSE, USE.NAMES = TRUE)
})


### True regression parameters ####
TrueBetaD <- lapply(1:length(TrueTheta), function(es){
  lapply(1:length(RhoXX), function(r){
    Phi2Beta(TruePhi[[k]][[es]][[1]][[r]], TruePhi[[k]][[es]][[2]][[r]], TruePhi[[k]][[es]][[3]][[r]], TruePhi[[k]][[es]][[4]][[r]], xLo = Values[["Discrete"]][["Intra_Lo"]], xHi = Values[["Discrete"]][["Intra_Hi"]])
  })
})


TrueBetaC <- lapply(1:length(TrueTheta), function(es){
  lapply(1:length(RhoXX), function(r){
    Phi2Beta(TruePhi[[k]][[es]][[1]][[r]], TruePhi[[k]][[es]][[2]][[r]], TruePhi[[k]][[es]][[3]][[r]], TruePhi[[k]][[es]][[4]][[r]], xLo = Values[["Continuous"]][["Intra_Lo"]], xHi = Values[["Continuous"]][["Intra_Hi"]])
  })
})

#lapply(TrueBetaD, function(x) sapply(x, function(y) any(!is.finite(y))))
#lapply(TrueBetaC, function(x) sapply(x, function(y) any(!is.finite(y))))

for(i in 1:(length(TrueBetaD)+length(TrueBetaC))){
    if(i %% 2 == 1){
      TrueBeta[[k]][[i]] <- TrueBetaD[[(i + 1) %/% 2]]
     }else{
      TrueBeta[[k]][[i]] <- TrueBetaC[[i %/% 2]]
      }
}
}

#### True pairwise correlation ####
PairwiseRho <- lapply(1:length(TrueTheta), function(es){
  sapply(1:length(TrueTheta[[1]]), function(subpop){
    lapply(1:length(RhoXX), function(r){
      ComputeRhoMultivariate(Phi = TruePhi[[3]][[es]][[subpop]][[r]])
    })
  }, simplify=FALSE, USE.NAMES = TRUE)
})




eTrueVal <- lapply(K, function(k) lapply(1:length(TrueBeta[[k]]), function(es) {lapply(1:length(Rho), function(r) {x <- vector("list", length(Populations))
names(x) <- Populations;
return(x)})}))
sTrueVal <- lapply(K, function(k) lapply(1:length(TrueBeta[[k]]), function(es) {lapply(1:length(Rho), function(r) {x <- vector("list", length(Populations)) 
names(x) <- Populations;
return(x)})}))

Cont <- seq(2,length(TrueBeta[[k]]),2)
Disc <- seq(1,length(TrueBeta[[k]]),2)
for(k in K){
for(es in 1:length(TrueBeta[[k]])){
  for(r in 1:length(Rho)){
  for(p in 1:length(Populations)){
    
    # Compute true values
    if(es %in% Cont){
      sTrueVal[[k]][[es]][[r]][[p]] <- sTrueValues(TrueBeta[[k]][[es]][[r]], ValueX = Values[["Continuous"]][[Populations[p]]], Weights = Weights[[k]])
      eTrueVal[[k]][[es]][[r]][[p]] <- eTrueValues(TrueBeta[[k]][[es]][[r]], MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[["Continuous"]][[Populations[p]]], Continuous = TRUE, Weights = Weights[[k]]) 
      
    } else if(es %in% Disc){ 
             sTrueVal[[k]][[es]][[r]][[p]] <- sTrueValues(TrueBeta[[k]][[es]][[r]], ValueX = Values[["Discrete"]][[Populations[p]]], Weights = Weights[[k]])
        pXd <- pX[[Populations[p]]]
        eTrueVal[[k]][[es]][[r]][[p]] <- eTrueValues(TrueBeta[[k]][[es]][[r]], pXD = pXd, RangeX = Ranges[["Discrete"]][[Populations[p]]], Continuous = FALSE, Weights = Weights[[k]])
      }
    
    names(eTrueVal[[k]][[es]]) <- names(sTrueVal[[k]][[es]]) <- Populations
  }
}
}
}

DgmNames <- paste0(rep(1:length(TrueTheta), each = 2), ".", rep(1:2, by = length(TrueTheta)))
DgmList <- list(1:nDgm, 1:nDgm, seq(1,nDgm,4)) # Select DGMs to be included per dimension of outcome variables
rm("eTrueValues", "sTrueValues")