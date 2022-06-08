load("Workspaces/Convergence.RData")

ConvergenceRange <- lapply(1:nDgm, function(dgm) lapply(1:length(nSample[[dgm]]), function(ss) range(Convergence[[dgm]][[ss]])))
NonConvergence <- lapply(1:nDgm, function(dgm) names(nSample[[dgm]])[which(max(unlist(ConvergenceRange[[dgm]])) > GR.Cut)])
names(NonConvergence) <- paste0("DGM ", 1:nDgm)

load("Workspaces/BiasRC.RData")
BiasRange <-  lapply(1:nDgm, function(dgm) lapply(1:length(nSample[[dgm]]), function(ss) range(BiasRC[[dgm]][[ss]])))
nLarge <- lapply(nSample, function(x) which(x >= 2000))
MaxBias <- lapply(1:nDgm, function(dgm) if(length(BiasRange[[dgm]][nLarge[[dgm]]]) > 0){max(abs(unlist(BiasRange[[dgm]][nLarge[[dgm]]])))})

BiasDelta <- lapply(1:2, function(type){x <- as.matrix(rbind(BiasTrial_4.1[[type]], BiasIntraLo_4.1[[type]], BiasTrial_4.2[[type]], BiasIntraLo_4.2[[type]]))
return(range(x[which(!is.na(x))]))
})

Out <- list(
  paste0("Maximum bias of regression coefficients among large sample size (n leq 1000) = ", round(max(unlist(MaxBias[!is.null(MaxBias)])), 3), collapse = ""),
  paste0("Maximum bias of treatment differences among DGM 4.1 & 4.2 = ", round(max(unlist(BiasDelta)), 3), collapse = ""),
  paste0("Maximum Gelman-Rubin statistic  = ", round(max(unlist(ConvergenceRange)), 3), collapse = ""),
  paste0("Gelman-Rubin was above ", GR.Cut, " in DGM ", DgmNames, " in conditions ", NonConvergence)
)
Out

