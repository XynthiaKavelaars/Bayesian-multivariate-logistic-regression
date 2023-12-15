for(k in 2:max(K)){
  set.seed(202311)
  nIt_MvB <- ifelse(k==3, 20e3, nIt)
  for(dgm in DgmList[[k]]){
    for(r in 1:length(Rho)){
      if(k>1 & dgm > 4){rule_vec <- Rules
      }else if(k == 1 | (k > 1 & dgm <= 4)){rule_vec <- Rules[[1]]
      }
      for(rule in 1:length(rule_vec)){
        Res <-  foreach(sim = 1:nSim)%dopar%{
          indicesK <- IndicesTheta(2^k)
          load(file = paste(wd, "/Workspaces/Data/Data_K", k, "/xData_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "sim", sim, ".Rdata", sep=""))
            
          # Average treatment effect
          m.thetaE <- e.thetaE <- s.thetaE <- array(NA, dim = c(nIt_MvB, k))
          m.thetaC <- e.thetaC <- s.thetaC <- array(NA, dim = c(nIt_MvB, k))
          
          Phi <- Psi2Phi(xData, TrueBeta[[k]][[dgm]][[r]])
          Y <- t(apply(Phi, 1, function(x) rmultinom(1,1,x)))
          
          m.sE <- colSums(Y[xData[,2] == 1,])
          m.sC <- colSums(Y[xData[,2] == 0,])
          
          m.PhiE <- MCMCpack::rdirichlet(nIt_MvB, m.sE + rep(a0, 2^k))
          m.PhiC <- MCMCpack::rdirichlet(nIt_MvB, m.sC + rep(a0, 2^k))
          
          m.thetaE <- Phi2Theta(Phi = m.PhiE, Indices = indicesK, K = length(indicesK))
          m.thetaC <- Phi2Theta(Phi = m.PhiC, Indices = indicesK, K = length(indicesK))
          
          # Conditional average treatment effect: interval
          if(dgm%%2==0){
            e.sE <- colSums(Y[xData[,2] == 1 & xData[,3] >= min(Ranges$Continuous$Intra_Lo) & xData[,3] <= max(Ranges$Continuous$Intra_Lo),,drop=FALSE])
            e.sC <- colSums(Y[xData[,2] == 0 & xData[,3] >= min(Ranges$Continuous$Intra_Lo) & xData[,3] <= max(Ranges$Continuous$Intra_Lo),,drop=FALSE])
              
            }else if(dgm%%2 == 1){
              e.sE <- colSums(Y[xData[,2] == 1 & xData[,3] == Values$Discrete$Intra_Lo,,drop=FALSE])
              e.sC <- colSums(Y[xData[,2] == 0 & xData[,3] == Values$Discrete$Intra_Lo,,drop=FALSE])
            } 
        
           e.PhiE <- MCMCpack::rdirichlet(nIt_MvB, e.sE + rep(a0, 2^k))
           e.PhiC <- MCMCpack::rdirichlet(nIt_MvB, e.sC + rep(a0, 2^k))
          
           e.thetaE <- Phi2Theta(Phi = e.PhiE, Indices = indicesK, K = length(indicesK))
           e.thetaC <- Phi2Theta(Phi = e.PhiC, Indices = indicesK, K = length(indicesK))
            
            # Conditional average treatment effect: value
           # if(dgm%%2==0){
          #    s.sE <- colSums(Y[xData[,2] == 1 & xData[,3] >= min(Values$Continuous$Intra_Lo) & xData[,3] <= max(Values$Continuous$Intra_Lo),,drop=FALSE])
          #    s.sC <- colSums(Y[xData[,2] == 0 & xData[,3] >= min(Values$Continuous$Intra_Lo) & xData[,3] <= max(Values$Continuous$Intra_Lo),,drop=FALSE])
              
          #  }else if(dgm%%2 == 1){
          #    s.sE <- colSums(Y[xData[,2] == 1 & xData[,3] == Values$Discrete$Intra_Lo,,drop=FALSE])
          #    s.sC <- colSums(Y[xData[,2] == 0 & xData[,3] == Values$Discrete$Intra_Lo,,drop=FALSE])
          #  } 
            
          #  s.PhiE <- MCMCpack::rdirichlet(nIt_MvB, s.sE + rep(a0, 2^k))
          #  s.PhiC <- MCMCpack::rdirichlet(nIt_MvB, s.sC + rep(a0, 2^k))
            
          #  s.thetaE <- Phi2Theta(Phi = s.PhiE, Indices = indicesK, K = length(indicesK))
          #  s.thetaC <- Phi2Theta(Phi = s.PhiC, Indices = indicesK, K = length(indicesK))
        
          
          m.delta <- m.thetaE - m.thetaC
          m.deltaW <- apply(m.delta, 1, function(x) x %*% Weights[[k]])
          
          m.pSingle <- colMeans(m.delta > 0)
          m.pAll <- mean(apply(m.delta, 1, function(x) min(x) > 0))
          m.pAny <- mean(apply(m.delta, 1, function(x) max(x) > 0))
          m.pComp <- mean(m.deltaW > 0)
          
          Trial <- list(ThetaE = colMeans(m.thetaE), ThetaC = colMeans(m.thetaC), Delta = colMeans(m.delta), DeltaW = mean(m.deltaW),
                        pSingle = m.pSingle, pAll = m.pAll, pAny = m.pAny, pComp = m.pComp)
          
          e.delta <- e.thetaE - e.thetaC
          e.deltaW <- apply(e.delta, 1, function(x) x %*% Weights[[k]])
          
          e.pSingle <- colMeans(e.delta > 0)
          e.pAll <- mean(apply(e.delta, 1, function(x) min(x) > 0))
          e.pAny <- mean(apply(e.delta, 1, function(x) max(x) > 0))
          e.pComp <- mean(e.deltaW > 0)
          
          eIntra_Lo <- list(ThetaE = colMeans(e.thetaE), ThetaC = colMeans(e.thetaC), Delta = colMeans(e.delta), 
                            pSingle = e.pSingle, pAll = e.pAll, pAny = e.pAny, pComp = e.pComp)
          
          #s.delta <- s.thetaE - s.thetaC
          #s.deltaW <- apply(s.delta, 1, function(x) x %*% Weights[[k]])
          
          #s.pSingle <- colMeans(s.delta > 0)
          #s.pAll <- mean(apply(s.delta, 1, function(x) min(x) > 0))
          #s.pAny <- mean(apply(s.delta, 1, function(x) max(x) > 0))
          #s.pComp <- mean(s.deltaW > 0)
          
          #sIntra_Lo <- list(ThetaE = colMeans(s.thetaE), ThetaC = colMeans(s.thetaC), Delta = colMeans(s.delta), 
          #                  pSingle = s.pSingle, pAll = s.pAll, pAny = s.pAny, pComp = s.pComp)
          
          
          
          c(Trial = Trial, eIntra_Lo = eIntra_Lo)#, sIntra_Lo = sIntra_Lo)
        }
        saveRDS(Res, file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_MvB.rds", sep=""))  
      }
    }
  }}
