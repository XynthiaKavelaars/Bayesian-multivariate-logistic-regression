
for(k in 2){
    for(dgm in 1:nDgm){
    for(r in 1:length(Rho)){
       if(k>1 & dgm > 4){rule_vec <- Rules
      }else if(k == 1 | (k > 1 & dgm <= 4)){rule_vec <- Rules[[1]]
      }
      for(rule in 1:length(rule_vec)){
     Qk <- Q[k]
    n <- nNeeded[[k]][[1]][[dgm]][[r]][rule]
    xT <- c(rep(0,n), rep(1,n))
    set.seed(Seeds[dgm,k,r,rule])

            Res <- foreach(sim = 1:nSim,
                  .errorhandling="pass", .verbose=TRUE)%dopar%{
                     if(dgm %in% Disc){
                     x <- rbinom(2*n,1,pXD)
                     }else if(dgm %in% Cont){
                      x <- rnorm(2*n)} 
                     xData <- cbind(rep(1, 2*n),xT, x, xT * x)
                     save(xData, file = paste(wd, "/Workspaces/Data/Data_K", k, "_uLR/xData_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "sim", sim, ".Rdata", sep=""))
                     
                     TruePhi <- Psi2Phi(xData, TrueBeta[[k]][[dgm]][[r]][,1:Qk])
                     
                     yData <- t(apply(TruePhi, 1, function(x) rmultinom(1,1,x)))
                     #Data <- list(X = xData, yMult = yData)
                     indicesK <- IndicesTheta(Qk)
                     indicesE <- which(xData[,2] == 1)
                     indicesC <- which(xData[,2] == 0)
                     
                    #### 3. Estimate parameters ####
                    filename.c1_1 <- paste(wd, "/Workspaces/Pars/Pars_K", k, "_uLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_Sim", sim, "Chain1_1.txt", sep="")
                    filename.c2_1 <- paste(wd, "/Workspaces/Pars/Pars_K", k, "_uLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule],  "_Sim", sim, "Chain2_1.txt", sep="")
                   
     			  filename.c1_2 <- paste(wd, "/Workspaces/Pars/Pars_K", k, "_uLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule],  "_Sim", sim, "Chain1_2.txt", sep="")
                    filename.c2_2 <- paste(wd, "/Workspaces/Pars/Pars_K", k, "_uLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule],  "_Sim", sim, "Chain2_2.txt", sep="")

                    theta.c1_1 <- SampleBetaPG(X=xData, Y=cbind(rowSums(yData[,indicesK[[1]]]), 1 - rowSums(yData[,indicesK[[1]]])), nBurn=nBurn, nIt=nIt, Start=Start[1], bMu0=b0, bSigma0=B0, filename = filename.c1_1, indicesE = indicesE, indicesC = indicesC, indicesK = list(1))
                    theta.c2_1 <- SampleBetaPG(X=xData, Y=cbind(rowSums(yData[,indicesK[[1]]]), 1 - rowSums(yData[,indicesK[[1]]])), nBurn=nBurn, nIt=nIt, Start=Start[2], bMu0=b0, bSigma0=B0, filename = filename.c2_1, indicesE = indicesE, indicesC = indicesC, indicesK = list(1))
                    
   			  theta.c1_2 <- SampleBetaPG(X=xData, Y=cbind(rowSums(yData[,indicesK[[2]]]), 1 - rowSums(yData[,indicesK[[2]]])), nBurn=nBurn, nIt=nIt, Start=Start[1], bMu0=b0, bSigma0=B0, filename = filename.c1_2, indicesE = indicesE, indicesC = indicesC, indicesK = list(1))
                    theta.c2_2 <- SampleBetaPG(X=xData, Y=cbind(rowSums(yData[,indicesK[[2]]]), 1 - rowSums(yData[,indicesK[[2]]])), nBurn=nBurn, nIt=nIt, Start=Start[2], bMu0=b0, bSigma0=B0, filename = filename.c2_2, indicesE = indicesE, indicesC = indicesC, indicesK = list(1))
               	
		        ThetaEstMean.E <- colMeans(do.call(rbind, lapply(1:nIt, function(i) rbind(colMeans(cbind(theta.c1_1[["thetaE"]][,,i], theta.c1_2[["thetaE"]][,,i])), colMeans(cbind(theta.c2_1[["thetaE"]][,,i], theta.c2_2[["thetaE"]][,,i]))))))
                    ThetaEstMean.C <- colMeans(do.call(rbind, lapply(1:nIt, function(i) rbind(colMeans(cbind(theta.c1_1[["thetaC"]][,,i], theta.c1_2[["thetaC"]][,,i])), colMeans(cbind(theta.c2_1[["thetaC"]][,,i], theta.c2_2[["thetaC"]][,,i]))))))
                    delta.c1 <- do.call(rbind, lapply(1:nIt, function(i) colMeans(cbind(theta.c1_1[["thetaE"]][,,i], theta.c1_2[["thetaE"]][,,i]) - cbind(theta.c1_1[["thetaC"]][,,i], theta.c1_2[["thetaC"]][,,i]))))
                    delta.c2 <- do.call(rbind, lapply(1:nIt, function(i) colMeans(cbind(theta.c2_1[["thetaE"]][,,i], theta.c2_2[["thetaE"]][,,i]) - cbind(theta.c2_1[["thetaC"]][,,i], theta.c2_2[["thetaC"]][,,i]))))

                    rm("filename.c1", "filename.c2", "theta.c1", "theta.c2")
                  
                    #### 5. Evaluate parameters ####
                  
                    DeltaEstMean <- colMeans(rbind(delta.c1, delta.c2))
                  
                            #### 6. Make decision ####
                     pSingle <- colMeans(rbind(delta.c1, delta.c2) > 0)
                     pAll <- mean(apply(rbind(delta.c1, delta.c2), 1, function(x) min(x) > 0))
                     pAny <- mean(apply(rbind(delta.c1, delta.c2), 1, function(x) max(x) > 0))
                     pComp <- mean(apply(rbind(delta.c1, delta.c2), 1, function(x) Weights[[k]] %*% x > 0))
                  
                     list(ThetaE = ThetaEstMean.E, ThetaC = ThetaEstMean.C, Delta = DeltaEstMean, pSingle = pSingle, pAll = pAll, pAny = pAny, pComp = pComp)
}

     saveRDS(Res, file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_uLR.rds", sep=""))
  gc()
}
}
}
  }

