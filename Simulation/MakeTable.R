#### Table 1. Effect sizes ####
Delta.Trial <- DeltaW.Trial <- Rho.Trial <- 
  Delta.Intra_Lo <- DeltaW.Intra_Lo <- Rho.Intra_Lo <- 
  lapply(1:length(Rho), function(r) vector("list", nDgm))
DeltaTab <- lapply(K, function(k) vector("list", length(Rho)))
for(k in 2:max(K)){
  for(r in 1:length(Rho)){
    Delta.Trial[[r]] <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]))
    DeltaW.Trial[[r]] <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[k]][[dgm]][[r]][["Trial"]][["DeltaW"]])) 
    Rho.Trial[[r]] <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[k]][[dgm]][[r]][["Trial"]][["RhoE"]][2,1]))
    
    Delta.Intra_Lo[[r]] <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]))
    DeltaW.Intra_Lo[[r]] <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["DeltaW"]])) 
    Rho.Intra_Lo[[r]] <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["RhoE"]][2,1]))
    
    DeltaTab[[k]][[r]] <- cbind(Delta.Trial[[r]], DeltaW.Trial[[r]], ifelse(abs(Rho.Trial[[r]]) < 1e-3, abs(Rho.Trial[[r]]), Rho.Trial[[r]]),
                           Delta.Intra_Lo[[r]], DeltaW.Intra_Lo[[r]], ifelse(abs(Rho.Intra_Lo[[r]]) < 1e-3, abs(Rho.Intra_Lo[[r]]), Rho.Intra_Lo[[r]]))
  }
}
#### 1.1 K=2 ####
Ind_Pairs.Delta.k2 <- matrix(1:ncol(DeltaTab[[2]][[1]]), nrow = 4)
Ind_Truths.k2 <- t(apply(DeltaTab[[2]][[1]][,c(1:3, 5:7)], 1, function(x) 
  c(x[1] > 1e-2, all(x[c(1,2)] > 1e-2), any(x[c(1,2)] > 1e-2), x[3] > 1e-2,
    x[4] > 1e-2, all(x[c(4,5)] > 1e-2), any(x[c(4,5)] > 1e-2), x[6] > 1e-2)))
colnames(Ind_Truths.k2) <- paste0(rep(c(RuleNames), times = 2), rep(c("Trial", "Sub"), each = length(Rules) ))
rownames(Ind_Truths.k2) <- DgmNames

Delta.Print.k2 <- do.call(rbind, lapply(1:nDgm, function(dgm){
  do.call(rbind, lapply(1:length(Rho),function(r){
  c(apply(Ind_Pairs.Delta.k2[1:2,], 2, function(z){
     if(r == 1){
       ES <- sapply(DeltaTab[[2]][[r]][dgm,z], function(x) ifelse(abs(x) < 1e-3, abs(x), x))
       c(NA, paste0("(", paste(sprintf("%5.3f", ES[1]),
                            sprintf("%7.3f", ES[2]),
                            sep=",", collapse = ""), ")"), 
               sprintf("%6.3f", DeltaTab[[2]][[r]][dgm,max(z)+1:2]))
    }else{
      c(rep(NA, 3), sprintf("%6.3f", DeltaTab[[2]][[r]][dgm,max(z)+2]))
  }
   }))
}))
}))

tabDelta.k2 <- cbind.data.frame(c(rbind(DgmNames, rep(" ", nDgm), rep(" ", nDgm))), 
                             rep(c("D", rep("", length(Rho) - 1),  "C", rep("", length(Rho) - 1)), times = nDgm/2), 
                             Delta.Print.k2)

AddToRow.Delta.k2 <- list()
AddToRow.Delta.k2$pos <- list(0, 0, nrow(tabDelta.k2))
AddToRow.Delta.k2$command <- c(paste0("
& ", paste0(" & & \\multicolumn{", ncol(Delta.Print.k2)/2 - 1, "}{c}{", c("ATE", "CATE"), "}", collapse = ""), " \\\\ 
                                      \\cmidrule(lr){3-", 3 + ncol(Delta.Print.k2)/2 - 1,"} 
                                      \\cmidrule(l){", 3 + ncol(Delta.Print.k2)/2, "-", ncol(tabDelta.k2),"} \n"),
                            paste0("ES & ", paste0(" & ", rep(paste0(" & \\multicolumn{1}{l}{", c("$(\\delta_{1}, \\delta_{2})$", "$\\delta (\\bm{\\mathsf{w}})$", "$\\rho_{\\theta^{k},\\theta^{l}}$"), "}", collapse = ""), 2), collapse = ""), " \\\\\n ", collapse = ""),
                            paste0("\\midrule \n \\multicolumn{", ncol(tabDelta.k2), "}{l}{Es = Effect size, D = Discrete covariate, C = Continuous covariate} \\\\\n"))

Align.Delta.k2 <- paste0("lll", paste(rep(paste0("p{0.02cm}", paste0(rep("r", ncol(Delta.Print.k2)/2 - 1), collapse = ""), collapse = ""), 2), collapse = ""))

print(xtable(tabDelta.k2, 
             caption="Parameters of average treatment effects (ATEs) in the trial and conditional average treatment effects (CATEs) in a subpopulation for two outcome variables. ",
             label="tab:Delta.k2", align=Align.Delta.k2,
             digits=c(1,1,1,rep(0,ncol(Delta.Print.k2)))),
      add.to.row=AddToRow.Delta.k2, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp",
      NA.string = " ", booktabs=TRUE)


#### 1.2 K=3 ####
Ind_Pairs.Delta.k3 <- matrix(1:ncol(DeltaTab[[3]][[1]]), nrow = 5)
Ind_Truths.k3 <- t(apply(DeltaTab[[3]][[1]][,c(1:4, 6:9)], 1, function(x) 
  c(x[1] > 1e-2, all(x[1:3] > 1e-2), any(x[1:3] > 1e-2), x[4] > 1e-2,
    x[5] > 1e-2, all(x[5:7] > 1e-2), any(x[5:7] > 1e-2), x[8] > 1e-2)))
colnames(Ind_Truths.k3) <- paste0(rep(c(Rules), times = 2), rep(c("Trial", "Sub"), each = length(Rules) ))
rownames(Ind_Truths.k3) <- DgmNames

Delta.Print.k3 <- do.call(rbind, lapply(1:nDgm, function(dgm){
  do.call(rbind, lapply(1:length(Rho),function(r){
    c(apply(Ind_Pairs.Delta.k3[1:max(K),], 2, function(z){
      if(r == 1){
        ES <- sapply(DeltaTab[[3]][[r]][dgm,z], function(x) ifelse(abs(x) < 1e-3, abs(x), x))
        c(NA, paste0("(", paste(sprintf("%5.3f", ES[1]),
                                sprintf("%7.3f", ES[2]),
                                sprintf("%7.3f", ES[3]), sep=",", collapse = ""), ")"), 
          sprintf("%6.3f", DeltaTab[[3]][[r]][dgm,max(z)+1:2]))
      }else{
        c(rep(NA, 3), sprintf("%6.3f", DeltaTab[[3]][[r]][dgm,max(z)+2]))
      }
    }))
  }))
}))

tabDelta.k3 <- cbind.data.frame(c(rbind(DgmNames, rep(" ", nDgm), rep(" ", nDgm))), 
                             rep(c("D", rep("", length(Rho) - 1),  "C", rep("", length(Rho) - 1)), times = nDgm/2), 
                             Delta.Print.k3)[c(1:3,13:15),]

AddToRow.Delta.k3 <- list()
AddToRow.Delta.k3$pos <- list(0, 0, nrow(tabDelta.k3))
AddToRow.Delta.k3$command <- c(paste0("
& ", paste0(" & & \\multicolumn{", ncol(Delta.Print.k3)/2 - 1, "}{c}{", c("ATE", "CATE"), "}", collapse = ""), " \\\\ 
                                      \\cmidrule(lr){3-", 3 + ncol(Delta.Print.k3)/2 - 1,"} 
                                      \\cmidrule(l){", 3 + ncol(Delta.Print.k3)/2, "-", ncol(tabDelta.k3),"} \n"),
                               paste0("ES & ", paste0(" & ", rep(paste0(" & \\multicolumn{1}{l}{", c("$(\\delta_{1}, \\delta_{2}, \\delta_{3})$", "$\\delta (\\bm{\\mathsf{w}})$", "$\\rho_{\\theta^{k},\\theta^{l}}$"), "}", collapse = ""), 2), collapse = ""), " \\\\\n ", collapse = ""),
                               paste0("\\midrule \n \\multicolumn{", ncol(tabDelta.k3), "}{l}{ES = Effect size, D = Discrete covariate, C = Continuous covariate} \\\\\n"))

Align.Delta.k3 <- paste0("lll", paste(rep(paste0("p{0.02cm}", paste0(rep("r", ncol(Delta.Print.k3)/2 - 1), collapse = ""), collapse = ""), 2), collapse = ""))

print(xtable(tabDelta.k3, 
             caption="Parameters of average treatment effects (ATEs) in the trial and conditional average treatment effects (CATEs) in a subpopulation for tree outcome variables. ",
             label="tab:Delta.k3", align=Align.Delta.k3,
             digits=c(1,1,1,rep(0,ncol(Delta.Print.k3)))),
      add.to.row=AddToRow.Delta.k3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp",
      NA.string = " ", booktabs=TRUE)

#### Table 2. Sample sizes ####
#### 2.2 K=2 ####
nNeeded.Trial <- do.call(rbind, lapply(1:nDgm, function(dgm){
  if(length(nNeeded[[2]][[1]][[dgm]][[1]]) > 1){
  do.call(rbind, nNeeded[[2]][[1]][[dgm]][1:length(Rho)])
  }else if(length(nNeeded[[2]][[1]][[dgm]][[1]]) == 1){
     matrix(nNeeded[[2]][[1]][[dgm]][1:length(Rho)], ncol = length(Rules), nrow = length(Rho)) 
    }
}))

nNeeded.Intra_Lo <- do.call(rbind, lapply(1:nDgm, function(dgm){
  if(length(nNeeded[[2]][[2]][[dgm]][[1]]) > 1){
    do.call(rbind, nNeeded[[2]][[2]][[dgm]][1:length(Rho)])
  }else if(length(nNeeded[[2]][[2]][[dgm]][[1]]) == 1){
    matrix(nNeeded[[2]][[2]][[dgm]][1:length(Rho)], ncol = length(Rules), nrow = length(Rho)) 
  }
}))

# Approximate remaining proportion of sample after stratification
propEffective <- rbind(c(1,rep(pXD[1], length(Populations)-1)), sapply(Populations, function(pop) abs(diff(pnorm(sort(Ranges[["Continuous"]][[pop]], decreasing = TRUE), lower.tail = TRUE)))))
pEffective <- propEffective[rep(rep(1:nrow(propEffective), each = length(Rho)), times = 1/2*nDgm), ]

pEffective

nEffective.Intra_Lo <- ceiling(do.call(cbind, lapply(1:length(Rules), function(rule) unlist(nNeeded.Trial[,rule]) * pEffective[,2])))

nDiff.Intra_Lo <- nEffective.Intra_Lo - unlist(nNeeded.Intra_Lo) 

SampleSizes <- do.call(cbind, lapply(1:length(Rules), function(rule){
  cbind.data.frame(unlist(nNeeded.Trial[,rule]), unlist(nNeeded.Intra_Lo[,rule]), unlist(nEffective.Intra_Lo[,rule]), NA)}))

SampleSizes.Print <- do.call(rbind.data.frame, lapply(1:(nDgm*length(Rho)), function(dgm){
   y <- x <- SampleSizes[dgm,-c(1:4)]
  y[which(x == nMax)] <- "-"
  ind <- which(unlist(x[seq(3,length(x),4)]) - unlist(x[seq(2,length(x),4)]) < 0 & unlist(y[seq(2,length(x),4)]) != "-")
  if(!(length(ind) == 0)){y[seq(3,length(y),4)[ind]] <- paste0("$\\textbf{", x[seq(3,length(x),4)[ind]],"}$")}
  y
}))


tabSampleSize <- cbind.data.frame(c(rbind(DgmNames, rep(" ", nDgm), rep(" ", nDgm))), 
                                  rep(paste0(c("$<0$", "$\\approx 0$", "$>0$")), times = nDgm), NA,
                                      SampleSizes.Print[,-ncol(SampleSizes.Print)])
rownames(tabSampleSize) <- NULL
AddToRow.SampleSizes <- list()
AddToRow.SampleSizes$pos <- list(0, 0, nrow(tabSampleSize))
AddToRow.SampleSizes$command <- c(paste0(" & ", paste0(" & & \\multicolumn{3}{c}{", RuleNames[-1], "}",collapse = ""), " \\\\\n 
                                         \\cmidrule(lr){4-6}
                                         \\cmidrule(l){8-10}
                                         \\cmidrule(l){12-14}"),
                                  paste0("ES & $\\rho_{k^{\\theta},l^{\\theta}}$ ", paste0(" &", rep(paste0(" & ", c("ATE", "CATE", "Sub"), collapse = ""), length(Rules)-1), collapse = ""), " \\\\\n ", collapse = ""),
                                  paste0("\\midrule \n
                                  \\multicolumn{", ncol(tabSampleSize), "}{l}{Sub = expected size of subsample} \\\\  
                                         \\multicolumn{", ncol(tabSampleSize), "}{l}{Bold-faced subsamples are smaller than required for estimation of the CATE} \\\\"))

Align.SampleSizes <- paste0("lll", paste0(rep("p{0.02cm}", length(Rules)-1), rep("rrr", length(Rules)-1), collapse = ""))

print(xtable(tabSampleSize, 
             caption="Required sample sizes to evaluate the average treatment effect (ATE) and conditional treatment effect (CATE) for two outcome variables.",
             label="tab:SampleSizes.k2", align=Align.SampleSizes),
             #digits=c(1,1,rep(0,ncol(SampleSizes.Print)))),
      add.to.row=AddToRow.SampleSizes, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE, sanitize.text.function = identity,
      table.placement="htbp",
      NA.string = " ", booktabs=TRUE)
#### 2.3 K=3 ####
nNeeded.Trial.k3 <- do.call(rbind, lapply(1:nDgm, function(dgm){
  if(length(nNeeded[[3]][[1]][[dgm]][[1]]) > 1){
    do.call(rbind, nNeeded[[3]][[1]][[dgm]][1:length(Rho)])
  }else if(length(nNeeded[[3]][[1]][[dgm]][[1]]) == 1){
    matrix(nNeeded[[3]][[1]][[dgm]][1:length(Rho)], ncol = length(Rules), nrow = length(Rho)) 
  }
}))

nNeeded.Intra_Lo.k3 <- do.call(rbind, lapply(1:nDgm, function(dgm){
  if(length(nNeeded[[3]][[2]][[dgm]][[1]]) > 1){
    do.call(rbind, nNeeded[[3]][[2]][[dgm]][1:length(Rho)])
  }else if(length(nNeeded[[3]][[2]][[dgm]][[1]]) == 1){
    matrix(nNeeded[[3]][[2]][[dgm]][1:length(Rho)], ncol = length(Rules), nrow = length(Rho)) 
  }
}))

# Approximate remaining proportion of sample after stratification
propEffective <- rbind(c(1,rep(pXD[1], length(Populations)-1)), sapply(Populations, function(pop) abs(diff(pnorm(sort(Ranges[["Continuous"]][[pop]], decreasing = TRUE), lower.tail = TRUE)))))
pEffective <- propEffective[rep(rep(1:nrow(propEffective), each = length(Rho)), times = 1/2*nDgm), ]

pEffective

nEffective.Intra_Lo.k3 <- ceiling(do.call(cbind, lapply(1:length(Rules), function(rule) unlist(nNeeded.Trial.k3[,rule]) * pEffective[,2])))

nDiff.Intra_Lo.k3 <- nEffective.Intra_Lo.k3 - unlist(nNeeded.Intra_Lo.k3) 

SampleSizes.k3 <- do.call(cbind, lapply(1:length(Rules), function(rule){
  cbind(unlist(nNeeded.Trial.k3[,rule]), unlist(nNeeded.Intra_Lo.k3[,rule]), unlist(nEffective.Intra_Lo.k3[,rule]), NA)}))
colnames(SampleSizes.k3) <- NULL
SampleSizes.k3

SampleSizes.Print.k3 <- do.call(rbind.data.frame, lapply(1:(nDgm*length(Rho)), function(dgm){
  y <- x <- SampleSizes.k3[dgm,-c(1:4)]
  y[which(x == nMax)] <- "-"
  ind <- which(unlist(x[seq(3,length(x),4)]) - unlist(x[seq(2,length(x),4)]) < 0 & unlist(y[seq(2,length(x),4)]) != "-")
  if(!(length(ind) == 0)){y[seq(3,length(y),4)[ind]] <- paste0("$\\bm{\\mathsf{", x[seq(3,length(x),4)[ind]],"}}$")}
  y
}))


tabSampleSize.k3 <- cbind.data.frame(c(rbind(DgmNames, rep(" ", nDgm), rep(" ", nDgm))), 
                                  rep(paste0(c("$<0$", "$\\approx 0$", "$>0$")), times = nDgm), NA,
                                  SampleSizes.Print.k3[,-ncol(SampleSizes.Print.k3)])[c(1:3,13:15),]
rownames(tabSampleSize.k3) <- NULL
AddToRow.SampleSizes.k3 <- list()
AddToRow.SampleSizes.k3$pos <- list(0, 0, nrow(tabSampleSize.k3))
AddToRow.SampleSizes.k3$command <- c(paste0(" & ", paste0(" & & \\multicolumn{3}{c}{", RuleNames[-1], "}",collapse = ""), " \\\\\n 
                                         \\cmidrule(lr){4-6}
                                         \\cmidrule(l){8-10}
                                         \\cmidrule(l){12-14}"),
                                     paste0("ES & $\\rho_{k^{\\theta},l^{\\theta}}$ ", paste0(" &", rep(paste0(" & ", c("ATE", "CATE", "Sub"), collapse = ""), length(Rules)-1), collapse = ""), " \\\\\n ", collapse = ""),
                                     paste0("\\midrule \n
                                  \\multicolumn{", ncol(tabSampleSize), "}{l}{Sub = expected size of subsample} \\\\ \n"))

Align.SampleSizes.k3 <- paste0("lll", paste0(rep("p{0.02cm}", length(Rules)-1), rep("rrr", length(Rules)-1), collapse = ""))

print(xtable(tabSampleSize.k3, 
             caption="Required sample sizes to evaluate the average treatment effect (ATE) and conditional treatment effect (CATE) for three outcome variables.",
             label="tab:SampleSizes.k3", align=Align.SampleSizes.k3),
             #digits=c(1,1,rep(0,ncol(SampleSizes.Print.k3)))),
      add.to.row=AddToRow.SampleSizes.k3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE, sanitize.text.function = identity,
      table.placement="htbp",
      NA.string = " ", booktabs=TRUE)

#### Table 3. pReject ####
pReject.ATE <- pReject.eCATE <- pReject.sCATE <- pReject.CATE <- 
  vector("list", max(K))

for(k in 2:max(K)){
  pSingle <- pAny <- pAll <- pComp <- 
    pSingle_uLR <- pAny_uLR <- pAll_uLR <- pComp_uLR <-
    pSingle_MvB <- pAny_MvB <- pAll_MvB <- pComp_MvB <- 
   lapply(1:nDgm, function(dgm) vector("list", length(Rho)))
  
  for(dgm in DgmList[[k]]){
    for(r in 1:length(Rho)){
      if(k > 1 & dgm > 4){rule_vec <- Rules
      }else if(k == 1 | (k > 1 & dgm <= 4)){rule_vec <- Rules[[1]]
      }
      for(rule in 1:length(rule_vec)){
        
Res <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], ".rds", sep=""))
Res_MvB <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_MvB.rds", sep=""))
if(k == 2){
   Res_uLR <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_uLR.rds", sep=""))
}else{
Res_uLR <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[rule], "_MvB.rds", sep=""))
}

if(length(rule_vec) == 1){
  if(k == 2){
  pRejectAll_uLR <- lapply(1:length(Res_uLR), function(sim) {
    if(all(is.finite(unlist(Res_uLR[[sim]])))){
      apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
      apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))
      }else{NA}})
  }else{
    pRejectAll_uLR <- lapply(1:length(Res_uLR), function(sim) {apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
        apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))})
    }
    pAll_uLR[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAll_uLR[which(sapply(pRejectAll_uLR, function(x) all(is.finite(x))))[1:nSim_eval]]))
  
  pRejectAll_MvB <- lapply(1:length(Res_MvB), function(sim) {apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
      apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))})
  pAll_MvB[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAll_MvB[which(sapply(pRejectAll_MvB, function(x) all(is.finite(x))))[1:nSim_eval]]))
  
  pRejectAll <- lapply(1:length(Res), function(sim) {apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
      apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))})
  pAll[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAll[which(sapply(pRejectAll, function(x) all(is.finite(x))))[1:nSim_eval]]))
  
  if(k == 2){
    pRejectAny_uLR <- lapply(1:length(Res_uLR), function(sim) {
      if(all(is.finite(unlist(Res_uLR[[sim]])))){
        apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
        apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))
        }else{NA}})
  }else{
    pRejectAny_uLR <- lapply(1:length(Res_uLR), function(sim) {apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
        apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))})
  }
  pAny_uLR[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAny_uLR[which(sapply(pRejectAny_uLR, function(x) all(is.finite(x))))[1:nSim_eval]]))
  
  pRejectAny_MvB <- lapply(1:length(Res_MvB), function(sim) {apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
      apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))})
  pAny_MvB[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAny_MvB[which(sapply(pRejectAny_MvB, function(x) all(is.finite(x))))[1:nSim_eval]]))
  
  pRejectAny <- lapply(1:length(Res), function(sim) {apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
      apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))})
  pAny[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAny[which(sapply(pRejectAny, function(x) all(is.finite(x))))[1:nSim_eval]]))
  
  if(k == 2){
    pRejectComp_uLR <- lapply(1:length(Res_uLR), function(sim) {
      if(all(is.finite(unlist(Res_uLR[[sim]])))){
        do.call(rbind, Res_uLR[[sim]][c("pComp")]) > max(pCut[[k]][4,]) | 
        do.call(rbind, Res_uLR[[sim]][c("pComp")]) < min(pCut[[k]][4,])
        }else{NA}})
  }else{
    pRejectComp_uLR <- lapply(1:length(Res_uLR), function(sim) {do.call(rbind, Res_uLR[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) > max(pCut[[k]][4,]) | 
        do.call(rbind, Res_uLR[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) < min(pCut[[k]][4,])})
  }
  pComp_uLR[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectComp_uLR[which(sapply(pRejectComp_uLR, function(x) all(is.finite(x))))[1:nSim_eval]]))

pRejectComp_MvB <- lapply(1:length(Res_MvB), function(sim) {do.call(rbind, Res_MvB[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) > max(pCut[[k]][4,]) | 
    do.call(rbind, Res_MvB[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) < min(pCut[[k]][4,])})
pComp_MvB[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectComp_MvB[which(sapply(pRejectComp_MvB, function(x) all(is.finite(x))))[1:nSim_eval]]))

pRejectComp <- lapply(1:length(Res), function(sim) {do.call(rbind, Res[[sim]][c("Trial.pComp", "eIntra_Lo.pComp", "sIntra_Lo.pComp")]) > max(pCut[[k]][4,]) | 
    do.call(rbind, Res[[sim]][c("Trial.pComp", "eIntra_Lo.pComp", "sIntra_Lo.pComp")]) < min(pCut[[k]][4,])})
pComp[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectComp[which(sapply(pRejectComp, function(x) all(is.finite(x))))[1:nSim_eval]]))

}else if(length(rule_vec) > 1){
 if(rule == 2){
   if(k == 2){
     pRejectAll_uLR <- lapply(1:length(Res_uLR), function(sim) {
       if(all(is.finite(unlist(Res_uLR[[sim]])))){
         apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
           apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))
       }else{NA}})
   }else{
     pRejectAll_uLR <- lapply(1:length(Res_uLR), function(sim) {apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
         apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))})
   }
   pAll_uLR[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAll_uLR[which(sapply(pRejectAll_uLR, function(x) all(is.finite(x))))[1:nSim_eval]]))
   
   pRejectAll_MvB <- lapply(1:length(Res_MvB), function(sim) {apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
       apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))})
   pAll_MvB[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAll_MvB[which(sapply(pRejectAll_MvB, function(x) all(is.finite(x))))[1:nSim_eval]]))
   
  pRejectAll <- lapply(1:length(Res), function(sim) {apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) min(x) > max(pCut[[k]][2,])) | 
      apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) max(x) < min(pCut[[k]][2,]))})
  pAll[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAll[which(sapply(pRejectAll, function(x) all(is.finite(x))))[1:nSim_eval]]))

     }else if(rule == 3){
       if(k == 2){
         pRejectAny_uLR <- lapply(1:length(Res_uLR), function(sim) {
           if(all(is.finite(unlist(Res_uLR[[sim]])))){
             apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
               apply(do.call(rbind, Res_uLR[[sim]][c("pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))
           }else{NA}})
       }else{
         pRejectAny_uLR <- lapply(1:length(Res_uLR), function(sim) {apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
             apply(do.call(rbind, Res_uLR[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))})
       }
       pAny_uLR[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAny_uLR[which(sapply(pRejectAny_uLR, function(x) all(is.finite(x))))[1:nSim_eval]]))
       
       pRejectAny_MvB <- lapply(1:length(Res_MvB), function(sim) {apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
           apply(do.call(rbind, Res_MvB[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))})
       pAny_MvB[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAny_MvB[which(sapply(pRejectAny_MvB, function(x) all(is.finite(x))))[1:nSim_eval]]))
       
       pRejectAny <- lapply(1:length(Res), function(sim) {apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) max(x) > max(pCut[[k]][3,])) | 
           apply(do.call(rbind, Res[[sim]][c("Trial.pSingle", "eIntra_Lo.pSingle", "sIntra_Lo.pSingle")]), 1, function(x) min(x) < min(pCut[[k]][3,]))})
       pAny[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectAny[which(sapply(pRejectAny, function(x) all(is.finite(x))))[1:nSim_eval]]))
     
       }else if(rule == 4){
   
         if(k == 2){
           pRejectComp_uLR <- lapply(1:length(Res_uLR), function(sim) {
             if(all(is.finite(unlist(Res_uLR[[sim]])))){
               do.call(rbind, Res_uLR[[sim]][c("pComp")]) > max(pCut[[k]][4,]) | 
                 do.call(rbind, Res_uLR[[sim]][c("pComp")]) < min(pCut[[k]][4,])
             }else{NA}})
         }else{
           pRejectComp_uLR <- lapply(1:length(Res_uLR), function(sim) {do.call(rbind, Res_uLR[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) > max(pCut[[k]][4,]) | 
               do.call(rbind, Res_uLR[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) < min(pCut[[k]][4,])})
         }
         pComp_uLR[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectComp_uLR[which(sapply(pRejectComp_uLR, function(x) all(is.finite(x))))[1:nSim_eval]]))
   
   pRejectComp_MvB <- lapply(1:length(Res_MvB), function(sim) {do.call(rbind, Res_MvB[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) > max(pCut[[k]][4,]) | 
       do.call(rbind, Res_MvB[[sim]][c("Trial.pComp", "eIntra_Lo.pComp")]) < min(pCut[[k]][4,])})
   pComp_MvB[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectComp_MvB[which(sapply(pRejectComp_MvB, function(x) all(is.finite(x))))[1:nSim_eval]]))
   
   pRejectComp <- lapply(1:length(Res), function(sim) {do.call(rbind, Res[[sim]][c("Trial.pComp", "eIntra_Lo.pComp", "sIntra_Lo.pComp")]) > max(pCut[[k]][4,]) | 
       do.call(rbind, Res[[sim]][c("Trial.pComp", "eIntra_Lo.pComp", "sIntra_Lo.pComp")]) < min(pCut[[k]][4,])})
   pComp[[dgm]][[r]] <- rowMeans(do.call(cbind, pRejectComp[which(sapply(pRejectComp, function(x) all(is.finite(x))))[1:nSim_eval]]))
 }
}
  
}
      }
    }
  
  pReject.ATE[[k]] <- do.call(rbind, 
    lapply(1:nDgm, function(dgm) {
      do.call(rbind, lapply(1:length(Rho), function(r) {
      c(pAll_uLR[[dgm]][[r]][1], pAll_MvB[[dgm]][[r]][1], pAll[[dgm]][[r]][1], 
        pAny_uLR[[dgm]][[r]][1], pAny_MvB[[dgm]][[r]][1], pAny[[dgm]][[r]][1], 
        pComp_uLR[[dgm]][[r]][1], pComp_MvB[[dgm]][[r]][1], pComp[[dgm]][[r]][1])
        
    }))
  }))
  pReject.eCATE[[k]] <- do.call(rbind, 
                              lapply(1:nDgm, function(dgm) {
                                do.call(rbind, lapply(1:length(Rho), function(r) {
                                  c(pAll_MvB[[dgm]][[r]][2], pAll[[dgm]][[r]][2],
                                    pAny_MvB[[dgm]][[r]][2], pAny[[dgm]][[r]][2],
                                    pComp_MvB[[dgm]][[r]][2], pComp[[dgm]][[r]][2])
                                  
                                }))
                              }))
  
  pReject.sCATE[[k]] <- do.call(rbind, 
                                lapply(1:nDgm, function(dgm) {
                                  do.call(rbind, lapply(1:length(Rho), function(r) {
                                    c(pAll[[dgm]][[r]][3], 
                                      pAny[[dgm]][[r]][3], 
                                      pComp[[dgm]][[r]][3])
                                    
                                  }))
                                }))
}

Effects <- c("ATEs with two outcome variables", "CATEs with two outcome variables", "three outcome variables")
Captions <- paste0("Proportions of superiority decisions for ", Effects, 
                   " by data-generating mechanism, correlation, and decision rule.")
names(Captions) <- c("k2.ATE", "k2.CATE", "k3")

#### 3.1 K=2 - ATE ####
pReject.ATE.k2.list <- lapply(1:ncol(pReject.ATE[[2]]), function(dgm) {t(matrix(pReject.ATE[[2]][,dgm], nrow = length(Rho)))})
pReject.ATE.k2 <- rbind(matrix(do.call(rbind, pReject.ATE.k2.list[1:3]), nrow = nrow(pReject.ATE.k2.list[[1]])),
                     matrix(do.call(rbind, pReject.ATE.k2.list[4:6]), nrow = nrow(pReject.ATE.k2.list[[4]])),
                     matrix(do.call(rbind, pReject.ATE.k2.list[7:9]), nrow = nrow(pReject.ATE.k2.list[[7]])))
label.ATE.k2 <- "tab:pReject.ATE.k2" 

pReject.ATE.Print <- do.call(rbind.data.frame, lapply(1:nrow(pReject.ATE.k2), function(dgm) {
  if(dgm %in% 1:nDgm){Truths <- Ind_Truths.k2[,2]
  }else if(dgm %in% c(nDgm+1:nDgm)){Truths <- Ind_Truths.k2[,3]
  }else if(dgm %in% c(2*nDgm+1:nDgm)){Truths <- Ind_Truths.k2[,4]}
  
  if(Truths[ifelse(dgm %% length(Truths) == 0, length(Truths), dgm %% length(Truths))]){
    x <- paste0("$\\textbf{", sprintf("%.3f", pReject.ATE.k2[dgm,]), "}$")
  }else{
    x <- sprintf("%.3f", pReject.ATE.k2[dgm,])
   }
  unlist(rbind.data.frame(NA, matrix(x, ncol = 3)))
}))
colnames(pReject.ATE.Print) <- NULL
pReject.ATE.Print

tabPRejectATE <- cbind.data.frame(rep(DgmNames, times = length(Rho)), pReject.ATE.Print)

AddToRow.pReject.ATE <- list()
AddToRow.pReject.ATE$pos <- as.list(c(0,0,seq(0,nrow(tabPRejectATE), nDgm)))
AddToRow.pReject.ATE$command <- c("\\toprule \n  & & \\multicolumn{3}{c}{$\\rho < 0$} & & \\multicolumn{3}{c}{$\\rho = 0$} & & \\multicolumn{3}{c}{$\\rho > 0$} \\\\
                                  \\cmidrule(l){3-5} 
                                  \\cmidrule(l){7-9} 
                                  \\cmidrule(l){11-13} \n",
                                   paste0("ES ", paste0(rep(" & & uLR & mB & mLR", 3), collapse = ""), " \\\\\n", collapse = ""),
                                   paste0(rep(paste0("\\cmidrule(l){3-", ncol(tabPRejectATE), "} \n"), length(Rules)-1), 
                                          "\\multicolumn{", ncol(tabPRejectATE), "}{c}{Rule = ", RuleNames[-1], "} \\\\\n",
                                          rep(paste0("\\cmidrule(l){3-", ncol(tabPRejectATE), "} \n"), length(Rules)-1) 
                                   ),
                                   paste0("\\midrule \n
                                           \\multicolumn{", ncol(tabPRejectATE),"}{l}{uLR = uLRvariate logistic regression} \\\\\n" , 
                                          "\\multicolumn{", ncol(tabPRejectATE),"}{l}{mB = Multivariate Bernoulli} \\\\\n", 
                                          "\\multicolumn{", ncol(tabPRejectATE),"}{l}{mLR = Multivariate logistic regression} \\\\\n",
                                          "\\multicolumn{", ncol(tabPRejectATE),"}{l}{Bold-faced entries have effect sizes that should lead to a superiority conclusion} \\\\\n
                                          \\bottomrule \n"))

Align.pReject.ATE <- paste0("ll", paste0(rep("p{0.02cm}rrr", 3), collapse = ""), collapse = "")

print(xtable(tabPRejectATE, 
             caption= Captions["k2.ATE"],
             label=label.ATE.k2, align=Align.pReject.ATE),
      add.to.row=AddToRow.pReject.ATE, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE,
      hline.after=NULL)


#### 3.2 K=2 - CATE ####
pReject.eCATE.list <- lapply(1:ncol(pReject.eCATE[[2]]), function(dgm) {t(matrix(pReject.eCATE[[2]][,dgm], nrow = length(Rho)))})
pReject.sCATE.list <- lapply(1:ncol(pReject.sCATE[[2]]), function(dgm) {t(matrix(pReject.sCATE[[2]][,dgm], nrow = length(Rho)))})
pReject.CATE.k2 <- rbind(matrix(do.call(rbind, c(pReject.eCATE.list[1:2], pReject.sCATE.list[1])), nrow = nrow(pReject.eCATE.list[[1]])),
                      matrix(do.call(rbind, c(pReject.eCATE.list[3:4], pReject.sCATE.list[2])), nrow = nrow(pReject.eCATE.list[[1]])),
                      matrix(do.call(rbind, c(pReject.eCATE.list[5:6], pReject.sCATE.list[3])), nrow = nrow(pReject.eCATE.list[[1]])))

label.CATE <- "tab:pReject.CATE.k2" 

pReject.CATE.Print <- do.call(rbind.data.frame, lapply(1:nrow(pReject.CATE.k2), function(dgm) {
  if(dgm %in% 1:nDgm){Truths <- Ind_Truths.k2[,6]
  }else if(dgm %in% c(nDgm+1:nDgm)){Truths <- Ind_Truths.k2[,7]
  }else if(dgm %in% c(2*nDgm+1:nDgm)){Truths <- Ind_Truths.k2[,8]}
  
  if(Truths[ifelse(dgm %% length(Truths) == 0, length(Truths), dgm %% length(Truths))]){
    x <- paste0("$\\textbf{", sprintf("%.3f", pReject.CATE.k2[dgm,]), "}$")
  }else{
    x <- sprintf("%.3f", pReject.CATE.k2[dgm,])
  }
  if(dgm %% 2 == 1){
    x[seq(2,length(x),length(Rho))] <- "-"
  }
  unlist(rbind.data.frame(NA, matrix(as.character(x), ncol = 3)))
  }))

colnames(pReject.CATE.Print) <- NULL
pReject.CATE.Print

tabPRejectCATE <- cbind.data.frame(rep(DgmNames, times = length(Rho)), pReject.CATE.Print)

AddToRow.pReject.CATE <- list()
AddToRow.pReject.CATE$pos <- as.list(c(0,0,seq(0,nrow(tabPRejectCATE), nDgm)))
AddToRow.pReject.CATE$command <- c("\\toprule \n  & & \\multicolumn{3}{c}{$\\rho < 0$} & & \\multicolumn{3}{c}{$\\rho = 0$} & & \\multicolumn{3}{c}{$\\rho > 0$} \\\\
                                  \\cmidrule(l){3-5} 
                                  \\cmidrule(l){7-9} 
                                  \\cmidrule(l){11-13} \n",
                                   paste0("ES ", paste0(rep(" & & mB & mLR-S & mLR-V", 3), collapse = ""), " \\\\\n", collapse = ""),
                                   paste0(rep(paste0("\\cmidrule(l){3-", ncol(tabPRejectCATE), "} \n"), length(Rules)-1), 
                                          "\\multicolumn{", ncol(tabPRejectCATE), "}{c}{Rule = ", RuleNames[-1], "} \\\\\n",
                                          rep(paste0("\\cmidrule(l){3-", ncol(tabPRejectCATE), "} \n"), length(Rules)-1) 
                                   ),
                                    paste0("\\midrule \n
                                           \\multicolumn{", ncol(tabPRejectCATE),"}{l}{mB = Multivariate Bernoulli} \\\\\n" , 
                                           "\\multicolumn{", ncol(tabPRejectCATE),"}{l}{mLR-S = Multivariate logistic regression - sample} \\\\\n",
                                           "\\multicolumn{", ncol(tabPRejectCATE),"}{l}{mLR-V = Multivariate logistic regression - value} \\\\\n",
                                           "\\multicolumn{", ncol(tabPRejectATE),"}{l}{Bold-faced entries have effect sizes that should lead to a superiority conclusion} \\\\\n
                                           \\bottomrule \n"))

Align.pReject.CATE <- paste0("ll", paste0(rep("p{0.02cm}rrr", 3), collapse = ""), collapse = "")

print(xtable(tabPRejectCATE, 
             caption= Captions["k2.CATE"],
             label=label.CATE, align=Align.pReject.CATE),
        add.to.row=AddToRow.pReject.CATE, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE, hline.after = NULL)

#### 3.3 K=3 ####
pReject.ATE.k3.list <- lapply(sort(c(seq(2,ncol(pReject.ATE[[3]]),3), seq(3,ncol(pReject.ATE[[3]]),3)), decreasing = FALSE), function(dgm) {
  t(matrix(pReject.ATE[[3]][,dgm], nrow = length(Rho)))})
pReject.eCATE.k3.list <- lapply(1:ncol(pReject.eCATE[[3]]), function(dgm) {t(matrix(pReject.eCATE[[3]][,dgm], nrow = length(Rho)))})
pReject.k3 <- rbind(matrix(do.call(rbind, pReject.ATE.k3.list[1:2]), nrow = nrow(pReject.ATE.k3.list[[1]])),
                    matrix(do.call(rbind, pReject.eCATE.k3.list[1:2]), nrow = nrow(pReject.eCATE.k3.list[[1]])),
                    matrix(do.call(rbind, pReject.ATE.k3.list[3:4]), nrow = nrow(pReject.ATE.k3.list[[1]])),
                    matrix(do.call(rbind, pReject.eCATE.k3.list[3:4]), nrow = nrow(pReject.eCATE.k3.list[[1]])),
                    matrix(do.call(rbind, pReject.ATE.k3.list[5:6]), nrow = nrow(pReject.ATE.k3.list[[1]])),
                    matrix(do.call(rbind, pReject.eCATE.k3.list[5:6]), nrow = nrow(pReject.eCATE.k3.list[[1]])))

label.k3 <- "tab:pReject.k3.RS" 

pReject.k3.Print <- do.call(rbind.data.frame, lapply(1:nrow(pReject.k3), function(dgm) {
  if(dgm %in% 1:length(DgmList[[3]])){Truths <- Ind_Truths.k3[DgmList[[k]],2]
  }else if(dgm %in% c(length(DgmList[[3]])+1:length(DgmList[[3]]))){Truths <- Ind_Truths.k3[DgmList[[3]],6]
  }else if(dgm %in% c(2*length(DgmList[[3]])+1:length(DgmList[[3]]))){Truths <- Ind_Truths.k3[DgmList[[3]],3]
  }else if(dgm %in% c(3*length(DgmList[[3]])+1:length(DgmList[[3]]))){Truths <- Ind_Truths.k3[DgmList[[3]],7]
  }else if(dgm %in% c(4*length(DgmList[[3]])+1:length(DgmList[[3]]))){Truths <- Ind_Truths.k3[DgmList[[3]],4]
  }else if(dgm %in% c(5*length(DgmList[[3]])+1:length(DgmList[[3]]))){Truths <- Ind_Truths.k3[DgmList[[3]],8]
  }
  
  if(Truths[ifelse(dgm %% length(Truths) == 0, length(Truths), dgm %% length(Truths))]){
    x <- paste0("$\\textbf{", sprintf("%.3f", pReject.k3[dgm,]), "}$")
  }else{
    x <- sprintf("%.3f", pReject.k3[dgm,])
  }
  unlist(rbind.data.frame(NA, matrix(x, ncol = 3)))
}))
colnames(pReject.k3.Print) <- NULL
pReject.k3.Print

tabPReject.k3 <- cbind.data.frame(rep(rep(c("1.1", "3.1"), times = 2), times = 3), rep(rep(c(" ATE", " CATE"), each = 2), times = 3), pReject.k3.Print)

AddToRow.pReject.k3 <- list()
AddToRow.pReject.k3$pos <- as.list(c(0,0, seq(0,nrow(tabPReject.k3), 4)))
AddToRow.pReject.k3$command <- c("\\toprule \n & & & \\multicolumn{2}{c}{$\\rho < 0$} & & \\multicolumn{2}{c}{$\\rho = 0$} & & \\multicolumn{2}{c}{$\\rho > 0$} \\\\
                                  \\cmidrule(l){4-5} 
                                  \\cmidrule(l){7-8} 
                                  \\cmidrule(l){10-11} \n",
                                  paste0("ES & Type", paste0(rep(" & & mB & mLR", 3), collapse = ""), " \\\\\n", collapse = ""),
                                  paste0(rep(paste0("\\cmidrule(l){4-", ncol(tabPReject.k3), "} \n"), length(Rules)-1), 
                                         "\\multicolumn{", ncol(tabPReject.k3), "}{c}{Rule = ", RuleNames[-1], "} \\\\\n",
                                         rep(paste0("\\cmidrule(l){4-", ncol(tabPReject.k3), "} \n"), length(Rules)-1) 
                                  ),
                                  paste0("\\midrule \n
                                        \\multicolumn{", ncol(tabPReject.k3),"}{l}{mB = Multivariate Bernoulli} \\\\\n", 
                                         "\\multicolumn{", ncol(tabPReject.k3),"}{l}{mLR = Multivariate logistic regression} \\\\\n",
                                         "\\multicolumn{", ncol(tabPReject.k3),"}{l}{Bold-faced entries should lead to a superiority conclusion} \\\\\n
                                          \\bottomrule \n"))

Align.pReject.k3 <- paste0("lll", paste0(rep("p{0.02cm}rr", 3), collapse = ""), collapse = "")

print(xtable(tabPReject.k3, 
             caption= Captions["k3"],
             label=label.k3, align=Align.pReject.k3),
       add.to.row=AddToRow.pReject.k3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE, hline.after = NULL)


#### Table 4. Bias - treatment difference delta ####
bias.ATE <- bias.eCATE <- bias.sCATE <-
  vector("list", max(K))


for(k in 2){
  biasSingle <- biasAny <- biasAll <- biasComp <- 
    biasSingle_uLR <- biasAny_uLR <- biasAll_uLR <- biasComp_uLR <-
    biasSingle_MvB <- biasAny_MvB <- biasAll_MvB <- biasComp_MvB <- lapply(1:nDgm, function(dgm) {lapply(1:length(Rho), function(r) vector("list", 3))})
  
  for(dgm in DgmList[[k]]){
  for(r in 1:length(Rho)){
      for(rule in 1:length(Rules)){
        if(k > 1 & dgm <= 4){
          Rule_dat <- 1
        }else{Rule_dat <- rule}
        
        tryCatch({
          Res <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[Rule_dat], ".rds", sep=""))
        Res_MvB <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[Rule_dat], "_MvB.rds", sep=""))
        if(k == 2){
          Res_uLR <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[Rule_dat], "_uLR.rds", sep=""))
        }else{
        Res_uLR <- readRDS(file = paste(wd, "/Workspaces/Results/Results_K", k, "/Results_DGM", dgm, "K", k, "Rho", names(Rho)[r], Rules[Rule_dat], "_MvB.rds", sep=""))
        }
          if(rule == 2){
            if(k == 2){
          bias1 <- lapply(1:length(Res_uLR), function(sim) {Res_uLR[[sim]][["Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]})
          #bias4 <- lapply(1:length(Res_uLR), function(sim) {Res_uLR[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
            }else{
              bias1 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]})
           #   bias4 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
             
            }
         biasAll_uLR[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias1[which(sapply(bias1, function(x) all(is.finite(x))))[1:nSim_eval]]))
        #   biasAll_uLR[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias4[which(sapply(bias4, function(x) all(is.finite(x))))[1:nSim_eval]]))
         
           bias2 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]})
          biasAll_MvB[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias2[which(sapply(bias2, function(x) all(is.finite(x))))[1:nSim_eval]]))
           bias5 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
          biasAll_MvB[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias5[which(sapply(bias5, function(x) all(is.finite(x))))[1:nSim_eval]]))
          bias3 <- lapply(1:length(Res), function(sim) {Res[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]})
          biasAll[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias3[which(sapply(bias3, function(x) all(is.finite(x))))[1:nSim_eval]]))
          bias6 <- lapply(1:length(Res), function(sim) {Res[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
          biasAll[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias6[which(sapply(bias6, function(x) all(is.finite(x))))[1:nSim_eval]]))
    
           }else if(rule == 3){
             if(k == 2){
               bias1 <- lapply(1:length(Res_uLR), function(sim) {Res_uLR[[sim]][["Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]}) 
            #   bias4 <- lapply(1:length(Res_uLR), function(sim) {Res_uLR[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
             }else{
               bias1 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]}) 
            #   bias4 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
             }
             
             biasAny_uLR[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias1[which(sapply(bias1, function(x) all(is.finite(x))))[1:nSim_eval]]))
             #biasAny_uLR[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias4[which(sapply(bias4, function(x) all(is.finite(x))))[1:nSim_eval]]))
        
              bias2 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]})
             biasAny_MvB[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias2[which(sapply(bias2, function(x) all(is.finite(x))))[1:nSim_eval]]))
              bias5 <- lapply(1:length(Res_MvB), function(sim) {Res_MvB[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
             biasAny_MvB[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias5[which(sapply(bias5, function(x) all(is.finite(x))))[1:nSim_eval]]))
             bias3 <- lapply(1:length(Res), function(sim) {Res[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["Delta"]]})
             biasAny[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias3[which(sapply(bias3, function(x) all(is.finite(x))))[1:nSim_eval]]))
             bias6 <- lapply(1:length(Res), function(sim) {Res[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["Delta"]]})
             biasAny[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias6[which(sapply(bias6, function(x) all(is.finite(x))))[1:nSim_eval]]))
     
           }else if(rule == 4){
             if(k == 2){
             bias1 <- lapply(1:length(Res_uLR), function(sim) {
               if(!is.null(Res_uLR[[sim]][["Delta"]])){Weights[[k]] %*% Res_uLR[[sim]][["Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["DeltaW"]]
               }else{NULL}})
            # bias4 <- lapply(1:length(Res_uLR), function(sim) {Weights[[k]] %*% Res_uLR[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["DeltaW"]]})
            }else{ 
               bias1 <- lapply(1:length(Res_MvB), function(sim) {Weights[[k]] %*% Res_MvB[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["DeltaW"]]})
            #   bias4 <- lapply(1:length(Res_MvB), function(sim) {Weights[[k]] %*% Res_MvB[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["DeltaW"]]})
               
             }
             biasComp_uLR[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias1[which(sapply(bias1, function(x) all(is.finite(x))))[1:nSim_eval]]))
             #biasComp_uLR[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias4[which(sapply(bias4, function(x) all(is.finite(x))))[1:nSim_eval]]))
             
            bias2 <- lapply(1:length(Res_MvB), function(sim) {Weights[[k]] %*% Res_MvB[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["DeltaW"]]})
            biasComp_MvB[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias2[which(sapply(bias2, function(x) all(is.finite(x))))[1:nSim_eval]]))
            bias5 <- lapply(1:length(Res_MvB), function(sim) {Weights[[k]] %*% Res_MvB[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["DeltaW"]]})
            biasComp_MvB[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias5[which(sapply(bias5, function(x) all(is.finite(x))))[1:nSim_eval]]))
            bias3 <- lapply(1:length(Res), function(sim) {Weights[[k]] %*% Res[[sim]][["Trial.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Trial"]][["DeltaW"]]})
            biasComp[[dgm]][[r]][[1]] <- colMeans(do.call(rbind, bias3[which(sapply(bias3, function(x) all(is.finite(x))))[1:nSim_eval]]))
            bias6 <- lapply(1:length(Res), function(sim) {Weights[[k]] %*% Res[[sim]][["eIntra_Lo.Delta"]] - eTrueVal[[k]][[dgm]][[r]][["Intra_Lo"]][["DeltaW"]]})
            biasComp[[dgm]][[r]][[2]] <- colMeans(do.call(rbind, bias6[which(sapply(bias6, function(x) all(is.finite(x))))[1:nSim_eval]]))
           }},
        warning = function(w) {print(paste0("Warning in dgm ", dgm, "correlation ", r, "rule ", rule, ": ", w))}
        )
           }
      }
  } 
        bias.ATE[[k]] <- c(lapply(1:nDgm, function(dgm) {list(
          do.call(rbind, lapply(1:length(Rho), function(r) {c(cbind(biasAll_uLR[[dgm]][[r]][[1]], biasAll_MvB[[dgm]][[r]][[1]], biasAll[[dgm]][[r]][[1]]))})),
          do.call(rbind, lapply(1:length(Rho), function(r) {c(cbind(biasAny_uLR[[dgm]][[r]][[1]], biasAny_MvB[[dgm]][[r]][[1]], biasAny[[dgm]][[r]][[1]]))})),
          do.call(rbind, lapply(1:length(Rho), function(r) {c(cbind(biasComp_uLR[[dgm]][[r]][[1]], biasComp_MvB[[dgm]][[r]][[1]], biasComp[[dgm]][[r]][[1]]))}))
          )}))
        
        #bias.eCATE[[k]] <- c(lapply(1:nDgm, function(dgm) {list(
        #  do.call(rbind, lapply(1:length(Rho), function(r) {c(cbind(biasAll_uLR[[dgm]][[r]][[2]], biasAll_MvB[[dgm]][[r]][[2]], biasAll[[dgm]][[r]][[2]]))})),
        #  do.call(rbind, lapply(1:length(Rho), function(r) {c(cbind(biasAny_uLR[[dgm]][[r]][[2]], biasAny_MvB[[dgm]][[r]][[2]], biasAny[[dgm]][[r]][[2]]))})),
        #  do.call(rbind, lapply(1:length(Rho), function(r) {c(cbind(biasComp_uLR[[dgm]][[r]][[2]], biasComp_MvB[[dgm]][[r]][[2]], biasComp[[dgm]][[r]][[2]]))}))
        #)}))
  }
    
#### 4.1 K=2 ####
bias.ATE.k2 <- bias.ATE[[2]][7:8]
label.bias.k2 <- "tab:Bias.k2" 

tabBias.ATE.k2 <- do.call(rbind.data.frame, 
  lapply(1:3, function(rule) {
  do.call(rbind, lapply(1:2, function(dgm) {
  do.call(rbind, lapply(1:length(Rho), function(r) {
    if(rule %in% c(1:2)){
      c(paste0("(", paste0(sprintf("%6.3f", bias.ATE.k2[[dgm]][[rule]][r,1:2]), collapse = ", "), ")", collapse = ""),
  paste0("(", paste0(sprintf("%6.3f", bias.ATE.k2[[dgm]][[rule]][r,3:4]), collapse = ", "), ")", collapse = ""),
  paste0("(", paste0(sprintf("%6.3f", bias.ATE.k2[[dgm]][[rule]][r,5:6]), collapse = ", "), ")", collapse = ""))
      }else if(rule == 3){
       paste0("\\multicolumn{1}{c}{", sprintf("%6.3f", bias.ATE.k2[[dgm]][[rule]][r,]), "}")
          }
  }))
  }))})
)

Effects <- c("Average", "Conditional")
tabBiasDelta.k2 <- cbind.data.frame(rep(c("4.1", "", "", "4.2", "", ""), 3), 
                                    rep(c("$<0$", "$\\approx 0$", "$>0$"), 6),
                                    tabBias.ATE.k2)

AddToRow.BiasDelta.k2 <- list()
AddToRow.BiasDelta.k2$pos <- as.list(c(0,seq(0,nrow(tabBiasDelta.k2),6)))
AddToRow.BiasDelta.k2$command <- c("\\toprule \n ES & $\\rho_{\\theta^{k},\\theta^{l}}$ & \\multicolumn{1}{c}{uLR} & 
                                          \\multicolumn{1}{c}{mB} & \\multicolumn{1}{c}{mLR} \\\\\n",
                                   paste0(c(paste0("\\cmidrule(l){1-", ncol(tabBiasDelta.k2), "}"), rep(paste0("\\cmidrule(l){2-", ncol(tabBiasDelta.k2), "}"), length(RuleNames)-2)),
                                          "\\multicolumn{", ncol(tabBiasDelta.k2), "}{c}{", RuleNames[-1], " rule} \\\\ \\midrule \n", 
                                          rep(paste0("\\cmidrule(l){2-", ncol(tabBiasDelta.k2), "} \n"), length(RuleNames)-1)),
                                          paste0("\\midrule 
                                   \\multicolumn{", ncol(tabBiasDelta.k2), "}{l}{uLR = uLRvariate logistic regression} \\\\
                                   \\multicolumn{", ncol(tabBiasDelta.k2), "}{l}{mB = multivariate Bernoulli} \\\\
                                   \\multicolumn{", ncol(tabBiasDelta.k2), "}{l}{mLR = multivariate logistic regression} \\\\
                                    \\bottomrule \n"))

AddToRow.BiasDelta.k2.pars_mv <- rep(paste0("ES & $\\rho_{\\theta^{k},\\theta^{l}}$", paste0(rep(c(" & \\multicolumn{1}{c}{($\\delta^{1},\\delta^{2}$)}"), length(Rho)), collapse = ""), "\\\\\n", collapse = ""), 2)
AddToRow.BiasDelta.k2.pars_w <- paste0("ES & $\\rho_{\\theta^{k},\\theta^{l}}$", paste0(rep(c(" & \\multicolumn{1}{c}{$\\delta (\\textbf{w})$}"), length(Rho)), collapse = ""), "\\\\\n", collapse = "")


AddToRow.BiasDelta.k2 <- list()
AddToRow.BiasDelta.k2$pos <- as.list(c(0,seq(0,nrow(tabBiasDelta.k2),6)))
AddToRow.BiasDelta.k2$command <- c("\\toprule ",
                                   paste0(c("", rep(paste0("\\midrule \n"), 2)),
                                          "\\multicolumn{", ncol(tabBiasDelta.k2), "}{c}{", RuleNames[-1], " rule} \\\\ \n",
                                          paste0("\\cmidrule(l){3-", ncol(tabBiasDelta.k2), "}"),
                                          c(" & & \\multicolumn{1}{c}{uLR} & 
                                            \\multicolumn{1}{c}{mB} & \\multicolumn{1}{c}{mLR} \\\\\n"),
                                          c("\\cmidrule(lr){3} \n \\cmidrule(lr){4} \n \\cmidrule(lr){5} \n"),
                                          c(AddToRow.BiasDelta.k2.pars_mv, AddToRow.BiasDelta.k2.pars_w), 
                                          rep(paste0("\\cmidrule(l){2-", ncol(tabBiasDelta.k2), "} \n"), length(RuleNames)-1)),
                                   paste0("\\midrule 
                                   \\multicolumn{", ncol(tabBiasDelta.k2), "}{l}{uLR = uLRvariate logistic regression} \\\\
                                   \\multicolumn{", ncol(tabBiasDelta.k2), "}{l}{mB = multivariate Bernoulli} \\\\
                                   \\multicolumn{", ncol(tabBiasDelta.k2), "}{l}{mLR = multivariate logistic regression} \\\\
                                    \\bottomrule \n"))

Align.BiasDelta.k2 <- paste0("ll", paste0(rep("c",ncol(tabBiasDelta.k2) - 1), collapse = ""), collapse = "")

print(xtable(tabBiasDelta.k2, 
             caption="Bias in average treatment differences of effect size (ES) $4.1$ and $4.2$ by decision rule.",
             label=label.bias.k2, align=Align.BiasDelta.k2),
     add.to.row=AddToRow.BiasDelta.k2, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " - ", booktabs=TRUE, hline.after = NULL)

