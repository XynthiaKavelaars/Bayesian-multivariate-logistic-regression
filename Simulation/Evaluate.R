for(k in 2:max(K)){
  for(dgm in DgmList[[k]]){
    for(r in 1:length(Rho)){
      if(k>1 & dgm > 4){rule_vec <- Rules
      }else if(k == 1 | (k > 1 & dgm <= 4)){rule_vec <- Rules[[1]]
      }
      for(rule in 1:length(rule_vec)){
        Res <-  foreach(sim = 1:nSim)%dopar%{
          Pars1 <- read.table(file = paste(wd, "/Workspaces/Pars/Pars_K", k, "_mLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_Sim", sim, "Chain1.txt", sep=""), sep=",")
          Pars2 <- read.table(file = paste(wd, "/Workspaces/Pars/Pars_K", k, "_mLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_Sim", sim, "Chain1.txt", sep=""), sep=",")
          Pars <- rbind(Pars1, Pars2)
          rm(Pars1, Pars2)
          load(file = paste(wd, "/Workspaces/Data/Data_K", k, "/xData_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "sim", sim, ".Rdata", sep=""))
          indicesK <- IndicesTheta(2^k)
          indicesE <- which(xData[,2] == 1)
          indicesC <- which(xData[,2] == 0)
          m.thetaE <- e.thetaE <- s.thetaE <- array(NA, dim = c(nrow(Pars), k))
          m.thetaC <- e.thetaC <- s.thetaC <- array(NA, dim = c(nrow(Pars), k))
          
          for(i in 1:nrow(Pars)){
            betaPG <- matrix(as.numeric(Pars[i,]), ncol = 2^k)
            
            m.phi <- Psi2Phi(xData, betaPG)
            m.thetaE[i,] <- colMeans(Phi2Theta(Phi = m.phi[indicesE,], Indices = indicesK, K = length(indicesK)))
            m.thetaC[i,] <- colMeans(Phi2Theta(Phi = m.phi[indicesC,], Indices = indicesK, K = length(indicesK)))
            
            
            if(dgm%%2==0){
              e.xE <- xData[xData[,2] == 1 & xData[,3] > Ranges$Continuous$Intra_Lo[1] & xData[,3] < Ranges$Continuous$Intra_Lo[2],]
              e.xC <- xData[xData[,2] == 0 & xData[,3] > Ranges$Continuous$Intra_Lo[1] & xData[,3] < Ranges$Continuous$Intra_Lo[2],]
              
            }else if(dgm%%2 == 1){
              e.xE <- xData[xData[,2] == 1 & xData[,3] == Values$Discrete$Intra_Lo,]
              e.xC <- xData[xData[,2] == 0 & xData[,3] == Values$Discrete$Intra_Lo,]
            } 
            
            e.phiE <- Psi2Phi(e.xE, betaPG)
            e.phiC <- Psi2Phi(e.xC, betaPG)
            
            e.thetaE[i,] <- colMeans(Phi2Theta(Phi = e.phiE, Indices = indicesK, K = length(indicesK)))
            e.thetaC[i,] <- colMeans(Phi2Theta(Phi = e.phiC, Indices = indicesK, K = length(indicesK)))
            
            
            if(dgm%%2==0){
              s.xE <- matrix(c(1,1,Values$Continuous$Intra_Lo, Values$Continuous$Intra_Lo * 1), nrow=1)
              s.xC <- matrix(c(1,0,Values$Continuous$Intra_Lo, Values$Continuous$Intra_Lo * 0), nrow=1)
              
            }else if(dgm%%2 == 1){
              s.xE <- matrix(c(1,1,Values$Discrete$Intra_Lo, Values$Discrete$Intra_Lo * 1), nrow=1)
              s.xC <- matrix(c(1,0,Values$Discrete$Intra_Lo, Values$Discrete$Intra_Lo * 0), nrow=1)
            } 
            s.phiE <- Psi2Phi(s.xE, betaPG)
            s.phiC <- Psi2Phi(s.xC, betaPG)
            
            s.thetaE[i,] <- colMeans(Phi2Theta(Phi = s.phiE, Indices = indicesK, K = length(indicesK)))
            s.thetaC[i,] <- colMeans(Phi2Theta(Phi = s.phiC, Indices = indicesK, K = length(indicesK)))
          }
          
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
          
          s.delta <- s.thetaE - s.thetaC
          s.deltaW <- apply(s.delta, 1, function(x) x %*% Weights[[k]])
          
          s.pSingle <- colMeans(s.delta > 0)
          s.pAll <- mean(apply(s.delta, 1, function(x) min(x) > 0))
          s.pAny <- mean(apply(s.delta, 1, function(x) max(x) > 0))
          s.pComp <- mean(s.deltaW > 0)
          
          sIntra_Lo <- list(ThetaE = colMeans(s.thetaE), ThetaC = colMeans(s.thetaC), Delta = colMeans(s.delta), 
                            pSingle = s.pSingle, pAll = s.pAll, pAny = s.pAny, pComp = s.pComp)
          
          
          
          c(Trial = Trial, eIntra_Lo = eIntra_Lo, sIntra_Lo = sIntra_Lo)
        }
        saveRDS(Res, file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], ".rds", sep=""))  
      }
    }
  }}
