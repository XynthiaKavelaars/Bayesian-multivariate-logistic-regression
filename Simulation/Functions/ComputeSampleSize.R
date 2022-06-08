#### Compute sample sizes ####
nNeededComp_OS <- 
 nNeededAll_OS <-
 nNeededAny_OS <-array(NA, dim = c(length(TrueBeta), length(Populations)), dimnames = list(paste(rep(1:(1/2*length(TrueBeta)), each = 2), 1:2, sep="."), Populations))
eTrueVal <- lapply(1:length(TrueBeta), function(s) vector("list", length(Ranges[["Continuous"]])))
sTrueVal <- lapply(1:length(TrueBeta), function(s) {x <- vector("list", length(Values[["Continuous"]])); 
names(x) <- names(Values[["Continuous"]]);
return(x)})

phiE <- phiC <- array(NA, dim = c(length(TrueBeta), 4, length(Populations)))
Cont <- seq(2,length(TrueBeta),2)
Disc <- seq(1,length(TrueBeta),2)


for(s in 1:length(TrueBeta)){
  for(p in 1:length(Populations)){
   
	# Compute true values
    if(s %in% Cont){
       sTrueVal[[s]][[Populations[p]]] <- sTrueValues(TrueBeta[[s]], ValueX = Values[[MeasurementLevel[s]]][[Populations[p]]], Weights = weights)
       eTrueVal[[s]][[p]] <- eTrueValues(TrueBeta[[s]], MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[[MeasurementLevel[s]]][[Populations[p]]], Continuous = TRUE, Weights = weights) 
     
      } else if(s %in% Disc){ 
      if(any(c("Extra_Lo", "Extra_Hi") %in% Populations[p])){next
       }else{
          sTrueVal[[s]][[Populations[p]]] <- sTrueValues(TrueBeta[[s]], ValueX = Values[["Discrete"]][[Populations[p]]], Weights = weights)
          pXd <- pX[[Populations[p]]]
         eTrueVal[[s]][[p]] <- eTrueValues(TrueBeta[[s]], pXD = pXd, RangeX = Ranges[["Discrete"]][[Populations[p]]], Continuous = FALSE, Weights = weights)
       }
      }
    
     names(eTrueVal[[s]]) <- Populations
     # Transform joint responses to success probabilities 
    phiE[s,,p] <- Theta2Phi(eTrueVal[[s]][[p]]$ThetaE, eTrueVal[[s]][[p]]$RhoE)
    phiC[s,,p] <- Theta2Phi(eTrueVal[[s]][[p]]$ThetaC, eTrueVal[[s]][[p]]$RhoC)
    
    # Compute sample sizes for different decision rules
    nNeededComp_OS[s,p] <- ComputeSampleSizeComp(phiE[s,,p], phiC[s,,p], weights[1], weights[2], nMax = nMax, alpha = 0.05, beta = 0.2, one.sided = TRUE)
    nNeededAll_OS[s,p] <- ComputeSampleSizeAll(phiE[s,,p], phiC[s,,p], nMax = nMax, alpha = 0.05, beta = 0.2, one.sided = TRUE)
    nNeededAny_OS[s,p] <- ComputeSampleSizeAny(phiE[s,,p], phiC[s,,p], nMax = nMax, alpha = 0.05, beta = 0.2, one.sided = TRUE)
    
  }
}

# Approximate remaining proportion of sample after stratification
propEffective <- rbind(c(1,rep(pXD[1], length(Populations)-1)), sapply(Populations, function(pop) abs(diff(pnorm(sort(Ranges[["Continuous"]][[pop]], decreasing = TRUE), lower.tail = TRUE)))))
pEffective <- propEffective[rep(1:nrow(propEffective), times = 1/2*length(TrueBeta)), ]
pEffective


# Required sample size given power analysis
nNeeded <- lapply(1:length(Populations), function(pop){
  ss <- cbind(nNeededComp_OS[,pop], nNeededAll_OS[,pop], nNeededAny_OS[,pop])
  colnames(ss) <- c(TestSide)
  return(ss)
})
names(nNeeded) <- Populations


# Approximate remaining sample size after stratification
nEffective <- lapply(1:ncol(pEffective), function(x) ceiling(nNeeded[[1]] * pEffective[,x]))
names(nEffective) <- Populations


# Total sampling conditions
nSample <- apply(nNeeded[[1]], 1, function(dgm) {ind <- which(!duplicated(dgm))
                 n <- dgm[ind]
                 n})


