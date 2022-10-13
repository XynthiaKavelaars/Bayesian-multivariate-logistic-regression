#### 1. Simulation setup ####
#### 1.1 Table: Sample size ####
nAny_rs <- cbind(nNeededAny_OS[,c("Trial","Intra_Lo")], nEffective[["Intra_Lo"]][,"Any_rs"])
nAll_rs <- cbind(nNeededAll_OS[,c("Trial","Intra_Lo")], nEffective[["Intra_Lo"]][,"All_rs"])
nComp_rs <- cbind(nNeededComp_OS[,c("Trial","Intra_Lo")], nEffective[["Intra_Lo"]][,"Comp_rs"])

SampleSizes.Print <- cbind(nAny_rs, NA, nAll_rs, NA, nComp_rs)
SampleSizes.Print1 <- do.call(rbind.data.frame, lapply(1:nDgm, function(dgm){
  y <- x <- SampleSizes.Print[dgm,]
  y[which(x == nMax)] <- "-"
  ind <- which(x[seq(3,length(x),4)] - x[seq(2,length(x),4)] < 0 & y[seq(2,length(x),4)] != "-")
if(!(length(ind) == 0)){y[seq(3,length(y),4)[ind]] <- paste0("$\\bm{", x[seq(3,length(x),4)[ind]],"}$")}
  y
}))


tabSampleSize <- cbind.data.frame(DgmNames, SampleSizes.Print1)
AddToRow.SampleSizes <- list()
AddToRow.SampleSizes$pos <- list(0, 0, nrow(tabSampleSize))
AddToRow.SampleSizes$command <- c(paste0(paste0(c("", rep(" & ", length(Rules) - 1)), " &  ", c(paste0("\\multicolumn{3}{l}{Any}"),
                                                 paste0("\\multicolumn{3}{l}{All}"),
                                                 paste0("\\multicolumn{3}{l}{Compensatory}")), collapse = ""), " \\\\\n "),
                                  paste0("DGM ", paste0(c(""," & "," & "), rep(paste0(" & ", c("ATE", "CTE", "Sub"), collapse = ""), length(Rules)), collapse = ""), " \\\\\n ", collapse = ""),
                                  paste0("\\midrule \n
                                  \\multicolumn{", ncol(tabSampleSize), "}{l}{Sub = expected size of subsample} \\\\  
                                         \\multicolumn{", ncol(tabSampleSize), "}{l}{Bold-faced subsamples are smaller than required for estimation of the CTE} \\\\"))
  
Align.SampleSizes <- paste0("ll", paste0(c("", rep("p{0.02cm}", length(Rules)-1)), rep("rrr", length(Rules)), collapse = ""))

print(xtable(tabSampleSize, 
             caption="Required sample sizes to evaluate the average treatment effect (ATE) and conditional treatment effect (CTE).",
             label="tab:SampleSizes", align=Align.SampleSizes,
             digits=c(1,1,rep(0,ncol(SampleSizes.Print)))),
      add.to.row=AddToRow.SampleSizes, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE, sanitize.text.function = identity,
      table.placement="htbp",
      NA.string = " ", booktabs=TRUE)

#### 1.2 Table: Effect size ####
Delta.Trial <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[dgm]][["Trial"]][["Delta"]]))
DeltaW.Trial <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[dgm]][["Trial"]][["DeltaW"]])) 
Rho.Trial <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[dgm]][["Trial"]][["RhoE"]]))

Delta.Intra_Lo <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[dgm]][["Intra_Lo"]][["Delta"]]))
DeltaW.Intra_Lo <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[dgm]][["Intra_Lo"]][["DeltaW"]])) 
Rho.Intra_Lo <- do.call(rbind, lapply(1:nDgm, function(dgm) eTrueVal[[dgm]][["Intra_Lo"]][["RhoE"]]))

DeltaTab <- cbind(Delta.Trial, DeltaW.Trial, Rho.Trial, Delta.Intra_Lo, DeltaW.Intra_Lo, Rho.Intra_Lo)
Ind_Pairs.Delta <- matrix(1:ncol(DeltaTab), nrow = 4)
Ind_Truths <- t(apply(DeltaTab[,c(1:3, 5:7)], 1, function(x) 
  c(any(x[c(1,2)] > 1e-2), all(x[c(1,2)] > 1e-2), x[3] > 1e-2,
    any(x[c(4,5)] > 1e-2), all(x[c(4,5)] > 1e-2), x[6] > 1e-2)))
colnames(Ind_Truths) <- paste0(rep(c(Rules), times = 2), rep(c("Trial", "Sub"), each = length(Rules) ))
rownames(Ind_Truths) <- DgmNames

Delta.Print <- do.call(rbind, lapply(1:nrow(DeltaTab), function(x) {
  c(apply(Ind_Pairs.Delta[1:2,], 2, function(z) {
    c(NA, paste0("(", paste(sprintf("%5.3f", DeltaTab[x,z[1]]),
                     sprintf("%7.3f", DeltaTab[x,z[2]]), sep=",", collapse = ""), ")"), sprintf("%6.3f", DeltaTab[x,max(z)+1:2]))
  }))
}))

tabDelta <- cbind.data.frame(DgmNames, MeasurementLevel, Delta.Print)

AddToRow.Delta <- list()
AddToRow.Delta$pos <- list(0, 0, nrow(tabDelta))
AddToRow.Delta$command <- c(paste0(" & ", paste0(" & & \\multicolumn{", ncol(Delta.Print)/2 - 1, "}{l}{", c("Average treatment effect", "Conditional treatment effect"), "}", collapse = ""), " \\\\\n \\midrule \n"),
                            paste0("DGM & Covariate", paste0(" & ", rep(paste0(" & \\multicolumn{1}{l}{", c("$(\\delta_{1}, \\delta_{2})$", "$\\delta (\\bm{w})$", "$\\rho_{\\theta^{k},\\theta^{l}}$"), "}", collapse = ""), 2), collapse = ""), " \\\\\n ", collapse = ""),
                            paste0("\\midrule \n \\multicolumn{", ncol(tabDelta), "}{l}{DGM = Data-generating mechanism} \\\\"))

Align.Delta <- paste0("lll", paste(rep(paste0("p{0.02cm}", paste0(rep("r", ncol(Delta.Print)/2 - 1), collapse = ""), collapse = ""), 2), collapse = ""))

print(xtable(tabDelta, 
             caption="Parameters of average treatment effects in the trial and conditional treatment effects in a subpopulation. ",
             label="tab:Delta", align=Align.Delta,
             digits=c(1,1,1,rep(0,ncol(Delta.Print)))),
      add.to.row=AddToRow.Delta, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp",
      NA.string = " ", booktabs=TRUE)


#### 1.3 Compute prior means (Appendix B) ####
theta1_Lo <- c(0.60,0.70)
theta1_Hi <- c(0.40,0.30)
theta0_Lo <- c(0.40,0.30)
theta0_Hi <- c(0.60,0.70)

rho_Lo <- -0.30
rho_Hi <- -0.30

BetaExample <- FindTrueBeta(TrueThetaC_Lo = theta0_Lo, TrueThetaC_Hi = theta0_Hi,
                      TrueThetaE_Lo = theta1_Lo, TrueThetaE_Hi = theta1_Hi, 
                      TrueRho_Lo = rho_Lo, TrueRho_Hi = rho_Hi,
                      xLo = -1, xHi = 1)

RC.Names <- paste0("$\\beta^{q}_{", 0:(P-1), "}$", collapse = "")
tabRCs <- cbind(RC.Names, apply(BetaExample, 2, function(x) sprintf("%1.3f", x)))

Align.PriorMeans <- paste0("l", paste0(rep("r", ncol(tabRCs)), collapse = ""), collapse = "")
AddToRow.PriorMeans <- list()
AddToRow.PriorMeans$pos <- as.list(c(0))
AddToRow.PriorMeans$command <- paste0(paste0("& $q = ",1:Q, "$", collapse = ""), "\\\\ \n", collapse = "")


print(xtable(tabRCs, 
             caption = "Example of means of the prior distribution of regression coefficients.",
             label = "tab:appPriorMeans",align =Align.PriorMeans,
             digits=c(1,1,rep(3,ncol(tabRCs)-1))),
      
      add.to.row=AddToRow.PriorMeans, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)

#### 2. pReject ####
load("Workspaces/Decisions.RData")
Effects <- c("ATEs", "CTEs")
Captions <- paste0("Proportions of superiority decisions for ", Effects, 
" by data-generating mechanism, estimation method, and decision rule.")
names(Captions) <- c("ATE", "CTE")

#### 2. Initialization ####
DecisionTypes <- c("DecisionAny.RS", "DecisionAny.LS", "DecisionAny.TS", 
                   "DecisionAll.RS", "DecisionAll.LS", "DecisionAll.TS", 
                   "DecisionCompensatory.RS", "DecisionCompensatory.LS", "DecisionCompensatory.TS")
BiasTypes <- c("BiasMultivariate", "BiasWeighted")
PopTypes <- c("PopMultivariate", "PopWeighted")
DecisionSelect_RS <- c("DecisionAny.RS", "DecisionAll.RS", "DecisionCompensatory.RS")
Bias <- sapply(BiasTypes, function(x){ y <- get0(x); if(is.array(y)){colMeans(y)}}, USE.NAMES = TRUE)
Pop <- sapply(PopTypes, function(x) get0(x), USE.NAMES = TRUE)

SupSpaces <- c("Any", "All", "Comp")
Rules <- c("Any", "All", "Compensatory")

Ind_Comp_rs <- do.call(rbind, lapply(1:nDgm, function(dgm) ifelse(any(grepl("Comp_rs", names(nSample[[dgm]]))), grep("Comp_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))
Ind_All_rs <- do.call(rbind, lapply(1:nDgm, function(dgm) ifelse(any(grepl("All_rs", names(nSample[[dgm]]))), grep("All_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))
Ind_Any_rs <- do.call(rbind, lapply(1:nDgm, function(dgm)ifelse(any(grepl("Any_rs", names(nSample[[dgm]]))), grep("Any_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))

IndicesML <- eTruths <- sTruths <- lapply(1:nDgm, function(dgm) c())

for(dgm in 1:nDgm){
  IndicesML[[dgm]] <- which(Types[,"MeasurementLevels"] == MeasurementLevel[dgm])
  eTruths <- sapply(Populations, function(pop) eTrueVal[[dgm]][[pop]][c("Delta", "DeltaW")], simplify = FALSE, USE.NAMES = TRUE)
  sTruths <- sapply(Populations, function(pop) sTrueVal[[dgm]][[pop]][c("Delta", "DeltaW")], simplify = FALSE, USE.NAMES = TRUE)
}

nSample_RS <- lapply(nSample, function(x) c(x[1], x[grep("rs", names(x))]))

#### 2.1 Average treatment effect (trial)  ####
# Compute proportion of rejections among samples without signs of non-convergence
pRejectTrial_RS <- lapply(1:nDgm, function(dgm){ 
  lapply(1:length(nSample[[dgm]]), function(ss){
    ar <- aperm(do.call(abind, c(sapply(seq_along(DecisionSelect_RS), function(s){
     
      # Proportion of rejected samples
       p <- rowMeans(do.call(cbind, lapply(1:nSim_Eval, function(sim) {
         tryCatch({ 
       
       do.call(rbind, sapply(rownames(Types[intersect(IndicesML[[dgm]], which(Types[,"Populations"] == "Trial")),]), function(i){
          Decision[[dgm]][[ss]][[sim]][[i]][["Decision"]][[DecisionSelect_RS[s]]]
        }, simplify = FALSE, USE.NAMES = TRUE))
        },
        error = function(e){
          message(paste0(e, " in sim ", sim, collapse = ""))
          })
      })))
      
      # Standard error
      se <- sqrt(p * (1-p) / nSim_Eval)
      return(cbind(p,se))
          }, simplify = FALSE, USE.NAMES = TRUE), along = 3)), c(1,3,2))
    dimnames(ar) <- list(rownames(Types[intersect(IndicesML[[dgm]], which(Types[,"Populations"] == "Trial")),]),
                         DecisionSelect_RS,
                         c("p", "se"))
    return(ar)
  })
})
pRejectTrial_RS

# Extract proportions per decision rule
pRejectTrial.Any_RS <- lapply(1:nrow(Ind_Any_rs), function(dgm){
  pRejectTrial_RS[[dgm]][[Ind_Any_rs[dgm]]][,"DecisionAny.RS",]
})

pRejectTrial.All_RS <- lapply(1:nrow(Ind_All_rs), function(dgm){
  pRejectTrial_RS[[dgm]][[Ind_All_rs[dgm]]][,"DecisionAll.RS",]
})

pRejectTrial.Comp_RS <- lapply(1:nrow(Ind_Comp_rs), function(dgm){
  pRejectTrial_RS[[dgm]][[Ind_Comp_rs[dgm]]][,"DecisionCompensatory.RS",]
})

# Conditions with discrete covariate
pRejectTrial.Disc_RS <- do.call(rbind, sapply(seq(1,nDgm,2), function(dgm) 
  unlist(cbind.data.frame(t(pRejectTrial.Any_RS[[dgm]]), 
                          t(pRejectTrial.All_RS[[dgm]]), 
                          t(pRejectTrial.Comp_RS[[dgm]]))) , simplify = FALSE, USE.NAMES = FALSE))

# Conditions with continuous covariate
pRejectTrial.Cont_RS <- do.call(rbind, sapply(seq(2,nDgm,2), function(dgm) 
  unlist(cbind.data.frame(t(pRejectTrial.Any_RS[[dgm]]), 
                          t(pRejectTrial.All_RS[[dgm]]), 
                          t(pRejectTrial.Comp_RS[[dgm]]))) , simplify = FALSE, USE.NAMES = FALSE))

# Attribute row/col names. Split, because continuous covariate has additional conditions (empirical, numerical marginalization)
rownames(pRejectTrial.Disc_RS) <- paste0(1:(nDgm/2), ".", rep(1,nDgm/2), sep = "")
rownames(pRejectTrial.Cont_RS) <- paste0(1:(nDgm/2), ".", rep(2,nDgm/2), sep = "")
names.Trial.Disc_RS <- Types[intersect(IndicesML[[1]], which(Types[,"Populations"] == "Trial")),"Methods"]
names.Trial.Cont_RS <- Types[intersect(IndicesML[[2]], which(Types[,"Populations"] == "Trial")),"Methods"]
colnames(pRejectTrial.Disc_RS) <- colnames(pRejectTrial.Cont_RS) <- paste0(rep(rep(names.Trial.Cont_RS, each = 2), times = length(SupSpaces)), "-", rep(SupSpaces, each = 2 * length(names.Trial.Cont_RS)))

# Rbind and reorder rows
pReject.Trial_RS <- rbind(pRejectTrial.Disc_RS, pRejectTrial.Cont_RS)
pReject.Trial_RS <- pReject.Trial_RS[c(matrix(1:nrow(pReject.Trial_RS), nrow = 2, byrow = TRUE)),]

# Reorder columns
Ind_Val <- grep("Value", colnames(pReject.Trial_RS))
Ind_Ana <- grep("Analytical", colnames(pReject.Trial_RS))
Ind_Emp <- grep("Empirical", colnames(pReject.Trial_RS))
Ind_Cols <- replace(1:ncol(pReject.Trial_RS), c(Ind_Val,Ind_Emp), c(Ind_Emp,Ind_Val))[-Ind_Ana]
#Ind_Cols <- replace(1:ncol(pReject.Trial_RS), c(Ind_Val,Ind_Emp,Ind_Ana), c(Ind_Emp,Ind_Ana,Ind_Val))
Names.Cols <- c("Reference", "Empirical",# "Analytical", 
                "Value")
pReject.Trial_RS <- pReject.Trial_RS[,Ind_Cols]


#### 2.1.1 Average treatment effect (trial) - Latex table ####
pReject.Trial <- pReject.Trial_RS
label.ATE <- "tab:pReject.ATE.RS" 

pReject.Trial.Print <- do.call(rbind, lapply(seq_along(Rules), function(rule) {
  pReject.Trial[,((rule-1) * ncol(pReject.Trial) / length(Rules) + 1):(rule * ncol(pReject.Trial) / length(Rules))]}))

pReject.Trial.Print <- t(sapply(1:(nDgm * length(Rules)), function(dgm) {
  x <- matrix(pReject.Trial.Print[dgm,], ncol = 2, byrow = TRUE)
  y <- matrix(NA, nrow = nrow(x), ncol = 3, byrow = TRUE);
  if((dgm %% 2 == 0)){
    if(c(Ind_Truths[1:nDgm,grep("Trial",colnames(Ind_Truths))])[dgm]){
      y[,2] <- paste0("$\\bm{", sprintf("%.3f", x[,1]), "}$")
    }else{
      y[,2] <- sprintf("%.3f", x[,1])}
    y[,3] <- paste0("(", sprintf("%.3f", x[,2]), ")")
  }else if((dgm %% 2) == 1){
    if(c(Ind_Truths[1:nDgm,grep("Trial",colnames(Ind_Truths))])[dgm]){ 
      y[,2] <- paste0("$\\bm{", sprintf("%.3f", x[,1]), "}$")
    }else{
      y[,2] <- sprintf("%.3f", x[,1])}
    y[,3] <- paste0("(", sprintf("%.3f", x[,2]), ")")}
   return(c(t(y)))}, simplify = TRUE, USE.NAMES = TRUE))

tabPRejectTrial <- cbind.data.frame(DgmNames, pReject.Trial.Print)

AddToRow.pReject.Trial <- list()
AddToRow.pReject.Trial$pos <- as.list(c(0,seq(0,nDgm * length(Rules)-1, nDgm)))
AddToRow.pReject.Trial$command <- c(paste0(paste0(paste0(" & & \\multicolumn{2}{l}{", Names.Cols, "}", collapse = ""), " \\\\\n "),
                                    paste0("DGM", 
                                           paste0(
                                                  " & ", rep(c(" & \\multicolumn{1}{l}{p} & \\multicolumn{1}{l}{SE}"), 
                                                             times = max(sapply(1:nDgm, function(dgm) 
                                                               length(intersect(IndicesML[[2]], which(Types[,"Populations"] == "Trial" & Types[,"Methods"] != "Analytical")))))), collapse = ""), " \\\\\n ", collapse = ""), collapse = ""),
                                    paste0(rep("\\midrule \n ", length(Rules)), 
                                           "\\multicolumn{", ncol(pReject.Trial.Print) + 
                                             1, "}{c}{Rule = ", Rules, "} \\\\\n",
                                             c("", rep("\\midrule \n", length(Rules)-1)) 
                                     )
                                    #paste0("\\midrule \n \\multicolumn{", ncol(pReject.Trial.Print) + 1, "}{l}{p = proportion of samples that concluded superiority} \\\\ \n 
                                    #       \\multicolumn{", ncol(pReject.Trial.Print) + 1, "}{l}{SE = standard error} \\\\ \n 
                                    #       \\multicolumn{", ncol(pReject.Trial.Print) + 1, "}{l}{Bold-faced proportions represent correct rejections (i.e. power).} \\\\ \n ")
                                    )
Align.pReject.Trial <- paste0("ll", paste(rep("p{0.02cm}rr", length(intersect(IndicesML[[2]], which(Types[,"Populations"] == "Trial" & Types[,"Methods"] != "Analytical")))), collapse = ""))

print(xtable(tabPRejectTrial, 
             caption= Captions["ATE"],
             label=label.ATE, align=Align.pReject.Trial,
             digits=c(1,1,rep(3,ncol(pReject.Trial.Print)))),
      add.to.row=AddToRow.pReject.Trial, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)
 
#### 2.2 Conditional treatment effect (low) ####
# Compute proportion of rejections among samples without signs of non-convergence
pRejectIntra_Lo_RS <- lapply(1:nDgm, function(dgm){ 
  lapply(1:length(nSample[[dgm]]), function(ss){
    ar <- aperm(do.call(abind, c(sapply(seq_along(DecisionSelect_RS), function(s){
      
      # Proportion of rejected samples
      p <- rowMeans(do.call(cbind, lapply(1:nSim_Eval, function(sim) {
          do.call(rbind, sapply(
            rownames(Types[intersect(IndicesML[[dgm]], which(Types[,"Populations"] == "Intra_Lo")),]), function(i){
            Decision[[dgm]][[ss]][[sim]][[i]][["Decision"]][[DecisionSelect_RS[s]]]
          }, simplify = FALSE, USE.NAMES = TRUE))
        }))) 
      
      # Standard error
      se <- sqrt(p * (1-p) / nSim_Eval)
      return(cbind(p,se))
    }, simplify = FALSE, USE.NAMES = TRUE), along = 3)), c(1,3,2))
    dimnames(ar) <- list(rownames(Types[intersect(IndicesML[[dgm]], which(Types[,"Populations"] == "Intra_Lo")),]),
                         DecisionSelect_RS,
                         c("p", "se"))
    return(ar)
  })
}) 
pRejectIntra_Lo_RS

# Extract proportions per decision rule
pRejectIntra_Lo.Any_RS <- lapply(1:nrow(Ind_Any_rs), function(dgm){
  pRejectIntra_Lo_RS[[dgm]][[Ind_Any_rs[dgm]]][,"DecisionAny.RS",]
})

pRejectIntra_Lo.All_RS <- lapply(1:nrow(Ind_All_rs), function(dgm){
  pRejectIntra_Lo_RS[[dgm]][[Ind_All_rs[dgm]]][,"DecisionAll.RS",]
})

pRejectIntra_Lo.Comp_RS <- lapply(1:nrow(Ind_Comp_rs), function(dgm){
  pRejectIntra_Lo_RS[[dgm]][[Ind_Comp_rs[dgm]]][,"DecisionCompensatory.RS",]
})

# Conditions with discrete covariate
pRejectIntra_Lo.Disc_RS <-  do.call(rbind, sapply(seq(1,nDgm,2), function(dgm) unlist(cbind.data.frame(
  cbind(t(pRejectIntra_Lo.Any_RS[[dgm]]), rep(NA, nrow(t(pRejectIntra_Lo.Any_RS[[dgm]]))), rep(NA, nrow(t(pRejectIntra_Lo.Any_RS[[dgm]])))), 
  cbind(t(pRejectIntra_Lo.All_RS[[dgm]]), rep(NA, nrow(t(pRejectIntra_Lo.All_RS[[dgm]]))), rep(NA, nrow(t(pRejectIntra_Lo.All_RS[[dgm]])))),
  cbind(t(pRejectIntra_Lo.Comp_RS[[dgm]]), rep(NA, nrow(t(pRejectIntra_Lo.Comp_RS[[dgm]]))), rep(NA, nrow(t(pRejectIntra_Lo.Comp_RS[[dgm]])))))),
  simplify = FALSE, USE.NAMES = FALSE))

# Conditions with continuous covariate
pRejectIntra_Lo.Cont_RS <- do.call(rbind, sapply(seq(2,nDgm,2), function(dgm) unlist(cbind.data.frame(
  t(pRejectIntra_Lo.Any_RS[[dgm]]), 
  t(pRejectIntra_Lo.All_RS[[dgm]]), 
  t(pRejectIntra_Lo.Comp_RS[[dgm]]))) , simplify = FALSE, USE.NAMES = FALSE))

# Attribute row/col names
rownames(pRejectIntra_Lo.Cont_RS) <- paste0(1:(nDgm/2), ".", rep(2,(nDgm/2)), sep = "")
rownames(pRejectIntra_Lo.Disc_RS) <- paste0(1:(nDgm/2), ".", rep(1,(nDgm/2)), sep = "")
names.Intra_Lo.Disc_RS <- Types[intersect(IndicesML[[1]], which(Types[,"Populations"] == "Intra_Lo")),"Methods"]
names.Intra_Lo.Cont_RS <- Types[intersect(IndicesML[[2]], which(Types[,"Populations"] == "Intra_Lo")),"Methods"]
colnames(pRejectIntra_Lo.Disc_RS) <- colnames(pRejectIntra_Lo.Cont_RS) <- paste0(rep(rep(names.Intra_Lo.Cont_RS, each = 2), times = length(SupSpaces)), "-", rep(SupSpaces, each = 2 * length(names.Intra_Lo.Cont_RS)))

# Rbind and reorder rows
pReject.Intra_Lo_RS <- rbind(pRejectIntra_Lo.Disc_RS, pRejectIntra_Lo.Cont_RS)
pReject.Intra_Lo_RS <- pReject.Intra_Lo_RS[c(matrix(1:nrow(pReject.Intra_Lo_RS), nrow = 2, byrow = TRUE)),]

# Reorder columns
Ind_Val <- grep("Value", colnames(pReject.Intra_Lo_RS))
Ind_Ana <- grep("Analytical", colnames(pReject.Intra_Lo_RS))
Ind_Emp <- grep("Empirical", colnames(pReject.Intra_Lo_RS))
Ind_Cols <- replace(1:ncol(pReject.Intra_Lo_RS), c(Ind_Val,Ind_Emp), c(Ind_Emp,Ind_Val))[-Ind_Ana]
Names.Cols <- c("Reference", "Empirical",# "Analytical", 
                "Value")

pReject.Intra_Lo_RS <- pReject.Intra_Lo_RS[,Ind_Cols]

#### 2.2.1 Conditional treatment effect (low) - Latex table ####
pReject.Intra_Lo <- pReject.Intra_Lo_RS
label.CTE <- "tab:pReject.CTE.RS"

Names.Cols <- c( "Reference", "Empirical", # "Analytical", 
                 "Value")
pReject.Intra_Lo.Print <- do.call(rbind, lapply(seq_along(Rules), function(rule) {
  pReject.Intra_Lo[,((rule-1) * ncol(pReject.Intra_Lo) / length(Rules) + 1):(rule * ncol(pReject.Intra_Lo) / length(Rules))]}))

pReject.Intra_Lo.Print <- t(sapply(1:(nDgm * length(Rules)), function(dgm) {
    x <- matrix(pReject.Intra_Lo.Print[dgm,], ncol = 2, byrow = TRUE)
    y <- matrix(NA, nrow = nrow(x), ncol = 3, byrow = TRUE);
    if((dgm %% 2 == 0)){
      if(c(Ind_Truths[1:nDgm,grep("Sub",colnames(Ind_Truths))])[dgm]){
        y[,2] <- paste0("$\\bm{", sprintf("%.3f", x[,1]), "}$")
      }else{
          y[,2] <- sprintf("%.3f", x[,1])}
      y[,3] <- paste0("(", sprintf("%.3f", x[,2]), ")")
      }else if((dgm %% 2) == 1){
        if(c(Ind_Truths[1:nDgm,grep("Sub",colnames(Ind_Truths))])[dgm]){ 
          y[-which(is.na(x[,1])),2] <- paste0("$\\bm{", sprintf("%.3f", x[-which(is.na(x[,1])),1]), "}$")
          }else{
    y[-which(is.na(x[,1])),2] <- sprintf("%.3f", x[-which(is.na(x[,1])),1])}
        y[-which(is.na(x[,1])),3] <- paste0("(", sprintf("%.3f", x[-which(is.na(x[,1])),2]), ")")}
    return(c(t(y)))}, simplify = TRUE, USE.NAMES = TRUE))

tabPRejectIntra_Lo <- cbind.data.frame(DgmNames, pReject.Intra_Lo.Print)


AddToRow.pReject.Intra_Lo <- list()
AddToRow.pReject.Intra_Lo$pos <- as.list(c(0,seq(0,nDgm * length(Rules)-1, nDgm)))
AddToRow.pReject.Intra_Lo$command <- c(paste0(paste0(paste0(" & & \\multicolumn{2}{l}{", Names.Cols, "}", collapse = ""), " \\\\\n "),
                                           paste0("DGM", 
                                                  paste0(
                                                    " & ", rep(c(" & \\multicolumn{1}{l}{p} & \\multicolumn{1}{l}{SE}"), 
                                                               times = max(sapply(1:nDgm, function(dgm) 
                                                                 length(intersect(IndicesML[[2]], which(Types[,"Populations"] == "Intra_Lo" & Types[,"Methods"] != "Analytical")))))), collapse = ""), " \\\\\n ", collapse = ""), collapse = ""),
                                    paste0(rep("\\midrule \n ", length(Rules)), 
                                           "\\multicolumn{", ncol(pReject.Intra_Lo.Print) + 
                                             1, "}{c}{Rule = ", Rules, "} \\\\\n", 
                                             c("", rep(" \\midrule \n", length(Rules)-1)) 
                                    )
                                    #paste0("\\midrule \n \\multicolumn{", ncol(pReject.Intra_Lo.Print) + 1, "}{l}{p = proportion of samples that concluded superiority} \\\\ \n 
                                    #       \\multicolumn{", ncol(pReject.Intra_Lo.Print) + 1, "}{l}{SE = standard error} \\\\ \n 
                                    #       \\multicolumn{", ncol(pReject.Intra_Lo.Print) + 1, "}{l}{Bold-faced proportions represent correct rejections (i.e. power).} \\\\ \n ")
)
Align.pReject.Intra_Lo <- paste0("ll", paste(rep("p{0.02cm}rr", length(intersect(IndicesML[[2]], which(Types[,"Populations"] == "Intra_Lo" & Types[,"Methods"] != "Analytical")))), collapse = ""))


print(xtable(tabPRejectIntra_Lo, 
             caption=Captions["CTE"],
             label=label.CTE, align=Align.pReject.Intra_Lo,
             digits=c(1,1,rep(3,ncol(pReject.Intra_Lo.Print)))),
      add.to.row=AddToRow.pReject.Intra_Lo, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)

#### 3. Bias - treatment difference delta ####
#### 3.0 Initialization ####
DecisionTypes <- c("DecisionAny.RS", "DecisionAny.LS", "DecisionAny.TS", 
                   "DecisionAll.RS", "DecisionAll.LS", "DecisionAll.TS", 
                   "DecisionCompensatory.RS", "DecisionCompensatory.LS", "DecisionCompensatory.TS")
BiasTypes <- c("BiasMultivariate", "BiasWeighted")
PopTypes <- c("PopMultivariate", "PopWeighted")
DecisionSelect_RS <- c("DecisionAny.RS", "DecisionAll.RS", "DecisionCompensatory.RS")
Bias <- sapply(BiasTypes, function(x){ y <- get0(x); if(is.array(y)){colMeans(y)}}, USE.NAMES = TRUE)
Pop <- sapply(PopTypes, function(x) get0(x), USE.NAMES = TRUE)

SupSpaces <- c("Any", "All", "Comp")
Rules <- c("Any", "All", "Compensatory")
Names.Cols <- c("Reference", "Empirical",# "Analytical", 
                "Value")
Ind_Comp_rs <- do.call(rbind, lapply(1:nDgm, function(dgm) ifelse(any(grepl("Comp_rs", names(nSample[[dgm]]))), grep("Comp_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))
Ind_All_rs <- do.call(rbind, lapply(1:nDgm, function(dgm) ifelse(any(grepl("All_rs", names(nSample[[dgm]]))), grep("All_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))
Ind_Any_rs <- do.call(rbind, lapply(1:nDgm, function(dgm)ifelse(any(grepl("Any_rs", names(nSample[[dgm]]))), grep("Any_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))

#### 3.1 Average treatment effect (trial) ####
names.Trial.Disc <- Types[intersect(IndicesML[[1]], which(Types[,"Populations"] == "Trial")),"Methods"]
names.Trial.Cont <- Types[intersect(IndicesML[[2]], which(Types[,"Populations"] == "Trial")),"Methods"]


BiasTrial <- lapply(1:nDgm, function(dgm){
  lapply(1:length(nSample[[dgm]]), function(ss){
    sapply(seq_along(BiasTypes), function(type){
      l <- lapply(1:nSim_Eval, function(sim) {
        tryCatch({ 
        do.call(rbind, sapply(rownames(Types[intersect(IndicesML[[dgm]], which(Types[,"Populations"] == "Trial")),]), function(i){
          Decision[[dgm]][[ss]][[sim]][[i]][["Bias"]][[type]]
        }, simplify = FALSE, USE.NAMES = TRUE))
        },
        error = function(e){
          message(paste0(e, " in sim ", sim, collapse = ""))
        })
        })
       Reduce("+", l[!sapply(l, is.null)]) / length(l)
    }, simplify = FALSE, USE.NAMES = TRUE)
   })
})

# Extract bias per decision rule
BiasTrial.Any_RS <- lapply(1:nrow(Ind_Any_rs), function(dgm){
  BiasTrial[[dgm]][[Ind_Any_rs[dgm]]]
})

BiasTrial.All_RS <- lapply(1:nrow(Ind_All_rs), function(dgm){
  BiasTrial[[dgm]][[Ind_All_rs[dgm]]]
})

BiasTrial.Comp_RS <- lapply(1:nrow(Ind_Comp_rs), function(dgm){
  BiasTrial[[dgm]][[Ind_Comp_rs[dgm]]]
})

# Conditions with discrete covariate
BiasTrial.Disc_RS <- sapply(seq_along(BiasTypes), function(type){ 
  do.call(rbind, sapply(seq(1,nDgm,2), function(dgm) 
    c(t(cbind(BiasTrial.Any_RS[[dgm]][[type]], BiasTrial.All_RS[[dgm]][[type]], BiasTrial.Comp_RS[[dgm]][[type]]))), simplify = FALSE, USE.NAMES = FALSE))})

BiasTrial.Disc_RS <- lapply(BiasTrial.Disc_RS, function(type){rownames(type) <- paste0(1:4, ".", rep(1,4), sep = ""); 
colnames(type) <- paste0(
  rep(paste0(rep(Rules, times = length(names.Trial.Cont)), "-", rep(names.Trial.Cont, each = length(Rules))),
      each = ncol(type) / (length(Rules) * length(names.Trial.Cont))),
  rep(1:(ncol(type) / (length(Rules) * length(names.Trial.Cont))), times = (length(Rules) * length(names.Trial.Cont))));
type})

# Conditions with continuous covariate
BiasTrial.Cont_RS <- sapply(seq_along(BiasTypes), function(type){
  do.call(rbind, sapply(seq(2,nDgm,2), function(dgm) 
    c(t(cbind(BiasTrial.Any_RS[[dgm]][[type]], BiasTrial.All_RS[[dgm]][[type]], BiasTrial.Comp_RS[[dgm]][[type]]))), simplify = FALSE, USE.NAMES = FALSE))})

BiasTrial.Cont_RS <- lapply(BiasTrial.Cont_RS, function(type){rownames(type) <- paste0(1:4, ".", rep(2,4), sep = ""); 
colnames(type) <- paste0(
  rep(paste0(rep(Rules, times = length(names.Trial.Cont)), "-", rep(names.Trial.Cont, each = length(Rules))),
      each = ncol(type) / (length(Rules) * length(names.Trial.Cont))),
  rep(1:(ncol(type) / (length(Rules) * length(names.Trial.Cont))), times = (length(Rules) * length(names.Trial.Cont))));
type})

# Rbind and reorder rows
Bias.Trial_RS <- mapply(rbind, BiasTrial.Disc_RS, BiasTrial.Cont_RS)
Bias.Trial_RS <- lapply(Bias.Trial_RS, function(x) x[c(matrix(1:nrow(x), nrow = 2, byrow = TRUE)),])

names(Bias.Trial_RS) <- BiasTypes

# Reorder columns
Ind_Val <- lapply(Bias.Trial_RS, function(x) grep("Value", colnames(x)))
Ind_Ana <- lapply(Bias.Trial_RS, function(x) grep("Analytical", colnames(x)))
Ind_Emp <- lapply(Bias.Trial_RS, function(x) grep("Empirical", colnames(x)))
Ind_Cols <- lapply(1:length(Bias.Trial_RS), function(x) replace(1:ncol(Bias.Trial_RS[[x]]),
                                                                   c(Ind_Val[[x]],Ind_Emp[[x]],Ind_Ana[[x]]
                                                                     ), 
                                                                   c(Ind_Emp[[x]],Ind_Ana[[x]],
                                                                     Ind_Val[[x]]))[-Ind_Ana[[x]]])
Bias.Trial_RS <- lapply(1:length(Bias.Trial_RS), function(x) Bias.Trial_RS[[x]][,Ind_Cols[[x]]])

#### 3.2 Conditional treatment effect (low)  ####
names.Intra_Lo.Disc <- Types[intersect(IndicesML[[1]], which(Types[,"Populations"] == "Intra_Lo")),"Methods"]
names.Intra_Lo.Cont <- Types[intersect(IndicesML[[2]], which(Types[,"Populations"] == "Intra_Lo")),"Methods"]

# Compute bias among samples without signs of non-convergence
BiasIntra_Lo <- lapply(1:nDgm, function(dgm){
  lapply(1:length(nSample[[dgm]]), function(ss){
    sapply(seq_along(BiasTypes), function(type){
      l <- lapply(1:nSim_Eval, function(sim) {
        tryCatch({ 
          do.call(rbind, sapply(rownames(
Types[intersect(IndicesML[[dgm]], which(Types[,"Populations"] == "Intra_Lo")),]), function(i){
            Decision[[dgm]][[ss]][[sim]][[i]][["Bias"]][[type]]
          }, simplify = FALSE, USE.NAMES = TRUE))
        },
        error = function(e){
          message(paste0(e, " in sim ", sim, collapse = ""))
        })
      })
      Reduce("+", l[!sapply(l, is.null)]) / length(l)
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
})

# Extract bias per decision rule
BiasIntra_Lo.Any_RS <- lapply(1:nrow(Ind_Any_rs), function(dgm){
  BiasIntra_Lo[[dgm]][[Ind_Any_rs[dgm]]]
})

BiasIntra_Lo.All_RS <- lapply(1:nrow(Ind_All_rs), function(dgm){
  BiasIntra_Lo[[dgm]][[Ind_All_rs[dgm]]]
})

BiasIntra_Lo.Comp_RS <- lapply(1:nrow(Ind_Comp_rs), function(dgm){
  BiasIntra_Lo[[dgm]][[Ind_Comp_rs[dgm]]]
})

# Conditions with discrete covariate
BiasIntra_Lo.Disc_RS <- sapply(seq_along(BiasTypes), function(type){
  do.call(rbind, sapply(seq(1,nDgm,2), function(dgm) c(t(cbind(
    rbind(BiasIntra_Lo.Any_RS[[dgm]][[type]], 
          t(matrix(NA, ncol = nrow(BiasIntra_Lo.Any_RS[[dgm + 1]][[type]]) - nrow(BiasIntra_Lo.Any_RS[[dgm]][[type]]), 
                 nrow = nrow(t(BiasIntra_Lo.Any_RS[[dgm]][[type]]))))), 
    rbind(BiasIntra_Lo.All_RS[[dgm]][[type]], 
          t(matrix(NA, ncol = nrow(BiasIntra_Lo.All_RS[[dgm + 1]][[type]]) - nrow(BiasIntra_Lo.All_RS[[dgm]][[type]]), 
                                                      nrow = nrow(t(BiasIntra_Lo.All_RS[[dgm]][[type]]))))),
    rbind(BiasIntra_Lo.Comp_RS[[dgm]][[type]], 
          t(matrix(NA, ncol = nrow(BiasIntra_Lo.Comp_RS[[dgm + 1]][[type]]) - nrow(BiasIntra_Lo.Comp_RS[[dgm]][[type]]), 
                                                       nrow = nrow(t(BiasIntra_Lo.Comp_RS[[dgm]][[type]])))))))),
    simplify = FALSE, USE.NAMES = FALSE))
})

# Conditions with continuous covariate
BiasIntra_Lo.Cont_RS <- sapply(seq_along(BiasTypes), function(type){
  do.call(rbind, sapply(seq(2,nDgm,2), function(dgm){
    c(t(cbind(BiasIntra_Lo.Any_RS[[dgm]][[type]], 
              BiasIntra_Lo.All_RS[[dgm]][[type]], 
              BiasIntra_Lo.Comp_RS[[dgm]][[type]])))
  }, simplify = FALSE, USE.NAMES = FALSE))
})

BiasIntra_Lo.Disc_RS <- lapply(BiasIntra_Lo.Disc_RS, function(type){rownames(type) <- paste0(1:4, ".", rep(1,4), sep = ""); 
colnames(type) <- paste0(
  rep(paste0(rep(Rules, times = length(names.Intra_Lo.Cont)), "-", rep(names.Intra_Lo.Cont, each = length(Rules))),
      each = ncol(type) / (length(Rules) * length(names.Intra_Lo.Cont))),
  rep(1:(ncol(type) / (length(Rules) * length(names.Intra_Lo.Cont))), times = (length(Rules) * length(names.Intra_Lo.Cont))));
type})

BiasIntra_Lo.Cont_RS <- lapply(BiasIntra_Lo.Cont_RS, function(type){rownames(type) <- paste0(1:4, ".", rep(2,4), sep = ""); 
colnames(type) <- paste0(
  rep(paste0(rep(Rules, times = length(names.Intra_Lo.Cont)), "-", rep(names.Intra_Lo.Cont, each = length(Rules))),
      each = ncol(type) / (length(Rules) * length(names.Intra_Lo.Cont))),
  rep(1:(ncol(type) / (length(Rules) * length(names.Intra_Lo.Cont))), times = (length(Rules) * length(names.Intra_Lo.Cont))));
type})

# Rbind and reorder rows
Bias.Intra_Lo_RS <- mapply(rbind,BiasIntra_Lo.Disc_RS, BiasIntra_Lo.Cont_RS)
Bias.Intra_Lo_RS <- lapply(Bias.Intra_Lo_RS, function(x) x[c(matrix(1:nrow(x), nrow = 2, byrow = TRUE)),])

names(Bias.Intra_Lo_RS) <- BiasTypes

# Reorder columns
Ind_Val <- lapply(Bias.Intra_Lo_RS, function(x) grep("Value", colnames(x)))
Ind_Ana <- lapply(Bias.Intra_Lo_RS, function(x) grep("Analytical", colnames(x)))
Ind_Emp <- lapply(Bias.Intra_Lo_RS, function(x) grep("Empirical", colnames(x)))
Ind_Cols <- lapply(1:length(Bias.Intra_Lo_RS), function(x) replace(1:ncol(Bias.Intra_Lo_RS[[x]]),
                                                                c(Ind_Val[[x]],Ind_Emp[[x]]),#,Ind_Ana[[x]]), 
                                                                c(Ind_Emp[[x]],#Ind_Ana[[x]],
                                                                  Ind_Val[[x]]))[-Ind_Ana[[x]]])
Bias_Intra_Lo_RS <- lapply(1:length(Bias.Intra_Lo_RS), function(x) Bias.Intra_Lo_RS[[x]][,Ind_Cols[[x]]])

#### 3.3 Select DGMs 4.1 & 4.2 ####
  Bias.Trial <- Bias.Trial_RS
  Bias_Intra_Lo <- Bias_Intra_Lo_RS
label.bias <- "tab:Bias4.2.RS" 


BiasTrial_4.1 <- lapply(seq_along(BiasTypes), function(type) {
  x <- as.data.frame(matrix(Bias.Trial[[type]][seq(7,nrow(Bias.Trial[[type]]),nDgm),], nrow = ncol(Bias.Trial[[2]]), byrow = TRUE))
  rownames(x) <- colnames(Bias.Trial[[2]])
  return(x)})

BiasIntraLo_4.1 <- lapply(seq_along(BiasTypes), function(type) {
  x <- as.data.frame(matrix(Bias_Intra_Lo[[type]][seq(7,nrow(Bias_Intra_Lo[[type]]),nDgm),], nrow = ncol(Bias_Intra_Lo[[2]]), byrow = TRUE))
  rownames(x) <- colnames(Bias_Intra_Lo[[2]])
  return(x)})

BiasTrial_4.2 <- lapply(seq_along(BiasTypes), function(type) {
  x <- as.data.frame(matrix(Bias.Trial[[type]][seq(nDgm,nrow(Bias.Trial[[type]]),nDgm),], nrow = ncol(Bias.Trial[[2]]), byrow = TRUE))
  rownames(x) <- colnames(Bias.Trial[[2]])
  return(x)})

BiasIntraLo_4.2 <- lapply(seq_along(BiasTypes), function(type) {
  x <- as.data.frame(matrix(Bias_Intra_Lo[[type]][seq(nDgm,nrow(Bias_Intra_Lo[[type]]),nDgm),], nrow = ncol(Bias_Intra_Lo[[2]]), byrow = TRUE))
  rownames(x) <- colnames(Bias_Intra_Lo[[2]])
  return(x)})



tabBiasTrial_4.1 <- cbind.data.frame(apply(BiasTrial_4.1[[1]], 1, function(x) paste0("(", paste0(sprintf("%6.3f", x), collapse = ", "), ")", collapse = "")), 
                                     sprintf("%6.3f", unlist(BiasTrial_4.1[[2]])))
tabBiasIntraLo_4.1 <- cbind.data.frame(apply(BiasIntraLo_4.1[[1]], 1, function(x) 
  if(all(!is.na(x))){paste0("(", paste0(sprintf("%6.3f", x), collapse = ", "), ")", collapse = "")
  }else{NA}),
  apply(BiasIntraLo_4.1[[2]], 1, function(x) 
    if(all(!is.na(x))){sprintf("%6.3f", x)
    }else{NA}))
 

tabBiasTrial_4.2 <- cbind.data.frame(apply(BiasTrial_4.2[[1]], 1, function(x) paste0("(", paste0(sprintf("%6.3f", x), collapse = ", "), ")", collapse = "")),
                                     sprintf("%6.3f", unlist(BiasTrial_4.2[[2]])))
tabBiasIntraLo_4.2 <- cbind.data.frame(apply(BiasIntraLo_4.2[[1]], 1, function(x) 
  if(all(!is.na(x))){paste0("(", paste0(sprintf("%6.3f", x), collapse = ", "), ")", collapse = "")
  }else{NA}),
  apply(BiasIntraLo_4.2[[2]], 1, function(x) 
    if(all(!is.na(x))){sprintf("%6.3f", x)
    }else{NA}))

# Discrete covariate
tabBiasTrial_4.1a <- cbind(tabBiasTrial_4.1[seq(1,nrow(tabBiasTrial_4.1),3),1],
                           tabBiasTrial_4.1[seq(2,nrow(tabBiasTrial_4.1),3),1],
                           tabBiasTrial_4.1[seq(3,nrow(tabBiasTrial_4.1),3),2])
tabBiasIntraLo_4.1a <- cbind(tabBiasIntraLo_4.1[seq(1,nrow(tabBiasIntraLo_4.1),3),1],
                           tabBiasIntraLo_4.1[seq(2,nrow(tabBiasIntraLo_4.1),3),1],
                           tabBiasIntraLo_4.1[seq(3,nrow(tabBiasIntraLo_4.1),3),2])

# Continuous covariate
tabBiasTrial_4.2a <- cbind(tabBiasTrial_4.2[seq(1,nrow(tabBiasTrial_4.2),3),1],
                           tabBiasTrial_4.2[seq(2,nrow(tabBiasTrial_4.2),3),1],
                           tabBiasTrial_4.2[seq(3,nrow(tabBiasTrial_4.2),3),2])
tabBiasIntraLo_4.2a <- cbind(tabBiasIntraLo_4.2[seq(1,nrow(tabBiasIntraLo_4.2),3),1],
                             tabBiasIntraLo_4.2[seq(2,nrow(tabBiasIntraLo_4.2),3),1],
                             tabBiasIntraLo_4.2[seq(3,nrow(tabBiasIntraLo_4.2),3),2])


colnames(tabBiasTrial_4.1a) <- colnames(tabBiasIntraLo_4.1a) <- 
  colnames(tabBiasTrial_4.2a) <- colnames(tabBiasIntraLo_4.2a) <- 
  SupSpaces 
rownames(tabBiasTrial_4.1a) <- rownames(tabBiasIntraLo_4.1a) <- 
  rownames(tabBiasTrial_4.2a) <- rownames(tabBiasIntraLo_4.2a) <- 
  Names.Cols

tabBias4.1 <- rbind.data.frame(tabBiasTrial_4.1a, tabBiasIntraLo_4.1a)
tabBias4.2 <- rbind.data.frame(tabBiasTrial_4.2a, tabBiasIntraLo_4.2a)

#### 3.3.1 Latex table DGMs 4.1 & 4.2 ####
Effects <- c("Average", "Conditional")
tabBiasDelta <- cbind.data.frame(rep(Names.Cols, times = 4), rbind(tabBias4.1, tabBias4.2))

ATR.DGM <- c(rbind(paste0(c("", rep("\n \\midrule \n", 3)), "\\multicolumn{", ncol(tabBiasDelta), "}{c}{DGM", rep(c("4.1 Discrete", "4.2 Continuous"), each = 2), " covariate - ", rep(Effects, times = 2), " Treatment Effect} \\\\ \n", c("", rep("\\midrule \n", 3)))))
ATR.Rules <- paste0("Method", paste0(" & ", "\\multicolumn{1}{l}{", SupSpaces, "}", collapse = ""), " \\\\ \n",
                               paste0(c(rep(" & \\multicolumn{1}{l}{$\\bm{\\delta}(\\bm{x})$}", 2), "& \\multicolumn{1}{l}{$\\delta (\\bm{w}, \\bm{x})$} \\\\ \n \\midrule \n  "), collapse = "")
                        )
   
AddToRow.BiasDelta <- list()
AddToRow.BiasDelta$pos <- as.list(c(0, 0, seq(from = length(Names.Cols), to = nrow(tabBiasDelta)-1, by = length(Names.Cols))))
AddToRow.BiasDelta$command <- c(ATR.Rules,
                                ATR.DGM)
                               
  
Align.BiasDelta <- paste0("ll", paste0(rep("r",ncol(tabBiasDelta) - 1), collapse = ""), collapse = "")

print(xtable(tabBiasDelta, 
             caption="Bias in average and conditional treatment differences of DGM $4.2$ by decision rule.",
             label="tab:BiasDelta", align=Align.BiasDelta,
             digits=c(1,1,rep(3,ncol(tabBiasDelta) - 1))),
      add.to.row=AddToRow.BiasDelta, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " - ", booktabs=TRUE)



#### 4.0 Figure median bias regression coefficients ####
Ind_Comp_rs <- do.call(rbind, lapply(1:nDgm, function(dgm) ifelse(any(grepl("Comp_rs", names(nSample[[dgm]]))), grep("Comp_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))
Ind_All_rs <- do.call(rbind, lapply(1:nDgm, function(dgm) ifelse(any(grepl("All_rs", names(nSample[[dgm]]))), grep("All_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))
Ind_Any_rs <- do.call(rbind, lapply(1:nDgm, function(dgm)ifelse(any(grepl("Any_rs", names(nSample[[dgm]]))), grep("Any_rs", names(nSample[[dgm]])), grep(nMax, nSample[[dgm]]))))

load("Workspaces/MedianRC.RData")

# Extract median regression coefficients for each of the decision rules per simulated dataset.
MedianRC.Comp_RS <- lapply(1:nDgm, function(dgm){
  MedianRC[[dgm]][[Ind_Comp_rs[dgm]]]
})

# Select dgm 4.2 and sample size of compensatory decision rule
  dgm <- 8
  ss <- nSample[[dgm]][which(nSample[[dgm]] == "Comp_rs")]
  
  # Draw figure
  postscript(file=paste0("Plots/Median_hist_RC_", dgm, ".eps"), horizontal = FALSE, width = 4.75, height = 4.75, font = "Times")

  layout(matrix(c(1:((P+1)*Q+1), rep((P+1)*Q+2, Q-1)), nrow= P+2, byrow = TRUE), 
         heights = c(2,rep(P,6),2), widths = c(1, rep(4, Q-1)))
  
  xlim <- matrix(NA, Q-1, 2)
  ylim <- matrix(NA, P, 2)
  
  par(mar = rep(0.01, 4), family = "Times")
  plot.new()
  
  for(q in 1:(Q-1)){
    par(mar = rep(0.01, 4))
    plot.new()
    text(x=0.5,y=0.5,labels=paste0("q = ", q), pos = 3)
    xlim[q,] <-c(floor(min(sapply(1:nSim_Eval, function(sim) MedianRC.Comp_RS[[dgm]][[sim]][,q]))), 
                 ceiling(max(sapply(1:nSim_Eval, function(sim) MedianRC.Comp_RS[[dgm]][[sim]][,q]))))
 }
  
  for(p in 1:P){
    par(mar = rep(0.01,4))
    plot.new()
    text(x = 0.05, y = 0.5, labels = bquote(beta[.(p)]), pos = 4)
    ylim[p,] <- c(0,2)
     
    for(q in 1:(Q-1)){
      par(mar = c(0.5,1,0.5,1))
      
      plot(NULL, xlim = xlim[q,], ylim=ylim[p,], xaxs = "i", yaxs = "i", xaxt="n", yaxt="n")
     
      if(q == 1){axis(2, at = c(ylim[p,1],0,ylim[p,2]), las = 1)}
      if(p == P){axis(1, at = c(min(xlim[q,]), 0, max(xlim[q,])),las = 1)}
      if(p == 1){axis(3, at = c(min(xlim[q,]), 0, max(xlim[q,])),las = 1)}
      
      hist(sapply(1:nSim_Eval, function(sim) MedianRC.Comp_RS[[dgm]][[sim]][p,q]), freq = FALSE, col = NULL, border = "gray30", add = TRUE, breaks = seq(min(xlim[q,]), max(xlim[q,]), 0.25))
      abline(v = TrueBeta[[dgm]][p,q], lty = 2)
      box()
    }
  }

  par(mar = rep(0.01,4))
  plot.new()
  plot.new()
  legend(x = 0.5, y = 0.5, xjust = 0.5, yjust = 1, legend = "True value", lty = 2, col = "black")
  
  dev.off()

