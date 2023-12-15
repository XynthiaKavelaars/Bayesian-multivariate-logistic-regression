
for(k in K){
    for(dgm in DgmList[[k]]){
    for(r in 1:length(Rho)){
         for(rule in 1:length(nNeeded[[k]][[1]][[dgm]][[r]])){
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
                     save(xData, file = paste(wd, "/Workspaces/Data/Data_K", k, "/xData_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "sim", sim, ".Rdata", sep=""))
                     
                     TruePhi <- Psi2Phi(xData, TrueBeta[[k]][[dgm]][[r]][,1:Qk])
                     
                     yData <- t(apply(TruePhi, 1, function(x) rmultinom(1,1,x)))
                     indicesK <- IndicesTheta(Qk)
                     indicesE <- which(xData[,2] == 1)
                     indicesC <- which(xData[,2] == 0)
                     
                    #### 3. Estimate parameters ####
                    filename.c1 <- paste(wd, "/Workspaces/Pars/Pars_K", k, "_mLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_Sim", sim, "Chain1.txt", sep="")
                    filename.c2 <- paste(wd, "/Workspaces/Pars/Pars_K", k, "_mLR/Pars_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_Sim", sim, "Chain2.txt", sep="")
                   
                    theta.c1 <- SampleBetaPG(X=xData, Y=yData, nBurn=nBurn, nIt=nIt, Start=Start[1], bMu0=b0, bSigma0=B0, filename = filename.c1, indicesE = indicesE, indicesC = indicesC, indicesK = indicesK)
                    theta.c2 <- SampleBetaPG(X=xData, Y=yData, nBurn=nBurn, nIt=nIt, Start=Start[2], bMu0=b0, bSigma0=B0, filename = filename.c2, indicesE = indicesE, indicesC = indicesC, indicesK = indicesK)
                    
                    ThetaEstMean.E <- colMeans(do.call(rbind, lapply(1:nIt, function(i) rbind(colMeans(theta.c1[["thetaE"]][,,i]), colMeans(theta.c2[["thetaE"]][,,i])))))
                    ThetaEstMean.C <- colMeans(do.call(rbind, lapply(1:nIt, function(i) rbind(colMeans(theta.c1[["thetaC"]][,,i]), colMeans(theta.c2[["thetaC"]][,,i])))))
                     delta.c1 <- do.call(rbind, lapply(1:nIt, function(i) colMeans(theta.c1[["thetaE"]][,,i] - theta.c1[["thetaC"]][,,i])))
                    delta.c2 <- do.call(rbind, lapply(1:nIt, function(i) colMeans(theta.c2[["thetaE"]][,,i] - theta.c2[["thetaC"]][,,i])))
                    
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

     saveRDS(Res, file = paste(wd, "/Workspaces/Results/Results_K2/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], ".rds", sep=""))
  gc()
}
}
}
}
