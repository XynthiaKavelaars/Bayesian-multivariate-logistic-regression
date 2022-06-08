#### 1. Decisions and bias of treatment differences ####
IndicesML <- eTruths <- sTruths <- lapply(1:nDgm, function(dgm) c())
BiasRC <- Decision <- Convergence <- Converged <- lapply(1:nDgm, function(dgm) lapply(1:length(nSample[[dgm]]), function(ss) c() ))
Convergence <- lapply(1:nDgm, function(dgm) vector("list"))


for(dgm in 1:nDgm){
  IndicesML[[dgm]] <- which(Types[,"MeasurementLevels"] == MeasurementLevel[dgm])
  eTruths <- sapply(Populations, function(pop) eTrueVal[[dgm]][[pop]][c("Delta", "DeltaW")], simplify = FALSE, USE.NAMES = TRUE)
  sTruths <- sapply(Populations, function(pop) sTrueVal[[dgm]][[pop]][c("Delta", "DeltaW")], simplify = FALSE, USE.NAMES = TRUE)
  
  for(ss in 1:length(nSample[[dgm]])){
    load(paste("Workspaces/Pars_DGM", dgm, "_n", names(nSample[[dgm]])[ss], ".RData", sep=""))
  
    Convergence[[dgm]][[ss]] <- sapply(lapply(Pars, "[[", "PG"), "[[", "ConvergencePG")
    Converged[[dgm]][[ss]] <- which(Convergence[[dgm]][[ss]] < GR.Cut)[1:nSim_Eval]

    BiasRC[[dgm]][[ss]] <- Reduce("+", lapply(Converged[[dgm]][[ss]], function(sim){
      Reduce("+", lapply(1:nIt, function(i) Pars[[sim]][["PG"]][["bDrawPG"]][[i]] - TrueBeta[[dgm]])) / nIt
      })) / nSim_Eval
 

    Decisions <- foreach(sim = Converged[[dgm]][[ss]],
                         .packages = packages, .verbose=TRUE)%dopar%{
                          
                           Theta <- Transform2Theta(BetaDrawPG = if(is.list(Pars[[sim]][["PG"]])){Pars[[sim]][["PG"]][["bDrawPG"]]} else{ NULL}, 
                                                    X = Pars[[sim]]$Data$X, 
                                                    Y = Pars[[sim]]$Data$yMult, 
                                                    Types = Types, 
                                                    Ranges = Ranges[[MeasurementLevel[dgm]]], Values = Values[[MeasurementLevel[dgm]]], 
                                                    PriorAlpha = A0,
                                                    MeasurementLevel = MeasurementLevel[dgm])
                           
                           lapply(IndicesML[[dgm]], function(i){
                             if(c("Value") %in% Types[i,"Methods"]){
                               Truths <- sTruths[[Types[i,"Populations"]]]
                               } else if(any(c("MvB", "Empirical", "Analytical") %in% Types[i,"Methods"])){
                                 Truths <- eTruths[[Types[i,"Populations"]]]
                                 }
                             x <- EvaluateData(Data = do.call(rbind, Map("-", Theta[["mTheta.E"]][[which(IndicesML[[dgm]] == i)]],  
                                                          Theta[["mTheta.C"]][[which(IndicesML[[dgm]] == i)]])), 
                                               Weights = weights, Alpha = Alpha, Alternative = c("greater.than", "two.sided"), 
                                               Rule = c("Any", "All", "Compensatory"), Truth = Truths, Types = Types[i,])
                             return(x)
                             
})
}                         
    Decision[[dgm]][[ss]] <- Decisions
    save(Decisions, file = paste("Workspaces/Decisions_DGM", dgm, "_n", names(nSample[[dgm]])[ss], ".RData", sep="")) 
    rm(Decisions)
 }
}

#### 2. Median, bias and convergence of regression coefficients ####
MedianRC <- lapply(1:nDgm, function(dgm) vector("list", length(nSample[[dgm]])))
for(dgm in 1:nDgm){
  for(ss in 1:length(nSample[[dgm]])){
    
    load(paste("Workspaces/Pars_DGM", dgm, "_n", names(nSample[[dgm]])[ss], ".RData", sep=""))
    Convergence[[dgm]][[ss]] <- sapply(lapply(Pars, "[[", "PG"), "[[", "ConvergencePG")
    Converged[[dgm]][[ss]] <- which(Convergence[[dgm]][[ss]] < GR.Cut)[1:nSim_Eval]

    BiasRC[[dgm]][[ss]] <- Reduce("+", lapply(Converged[[dgm]][[ss]], function(sim){
      cbind(do.call(rbind, lapply(1:P, function(p) sapply(1:(Q-1), function(q)
        median(sapply(1:nIt, function(i) Pars[[sim]][["PG"]][["bDrawPG"]][[i]][p,q] - TrueBeta[[dgm]][p,q]))
        ))), 0)})) / nSim_Eval
    
    MedianRC[[dgm]][[ss]] <- lapply(Converged[[dgm]][[ss]], function(sim){
      cbind(do.call(rbind, lapply(1:P, function(p) sapply(1:(Q-1), function(q)
        median(sapply(1:nIt, function(i) Pars[[sim]][["PG"]][["bDrawPG"]][[i]][p,q]))
        ))), 0)
    })
    }}

save(Decision, file = "Workspaces/Decisions.RData")
save(Convergence, file = "Workspaces/Convergence.RData")
save(MedianRC, file = "Workspaces/MedianRC.RData")
save(BiasRC, file = "Workspaces/BiasRC.RData")



