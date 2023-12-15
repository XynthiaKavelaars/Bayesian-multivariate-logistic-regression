#### Compute sample sizes ####
nNeeded <- lapply(K, function(k) lapply(1:length(Populations), function(p) {lapply(1:length(TrueBeta[[1]]), function(dgm) {lapply(1:length(Rho), function(r) rep(NA, 4))})}))
for(k in K){
  for(p in 1:length(Populations)){
  for(es in 1:length(TrueBeta[[k]])){
    for(r in 1:length(Rho)){
      TruePhiE <- eTrueVal[[k]][[es]][[r]][[p]][["PhiE"]]
      TruePhiC <- eTrueVal[[k]][[es]][[r]][[p]][["PhiC"]]
      nNeeded[[k]][[p]][[es]][[r]][which(Rules == "Single")] <- ComputeSampleSizeComp(TruePhiE, TruePhiC, c(1, rep(0,k-1)), nMax)
    if(k > 1){
     nNeeded[[k]][[p]][[es]][[r]][which(Rules == "All")] <- ComputeSampleSizeAll(TruePhiE, TruePhiC, nMax)
     nNeeded[[k]][[p]][[es]][[r]][which(Rules == "Any")] <- ComputeSampleSizeAny(TruePhiE, TruePhiC, nMax)
}
     nNeeded[[k]][[p]][[es]][[r]][which(Rules == "Comp")] <- ComputeSampleSizeComp(TruePhiE, TruePhiC, Weights[[k]], nMax)
     names(nNeeded[[k]][[p]][[es]][[r]]) <- Rules
     if(max(nNeeded[[k]][[p]][[es]][[r]], na.rm=TRUE) - min(nNeeded[[k]][[p]][[es]][[r]], na.rm=TRUE) < .01){
       nNeeded[[k]][[p]][[es]][[r]] <- nNeeded[[k]][[p]][[es]][[r]]["Single"]
     }
  }
  }
  }
}

rm("ComputeSampleSizeComp", "ComputeSampleSizeAll", "ComputeSampleSizeAny")