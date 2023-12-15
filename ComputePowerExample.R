#### 0. Initialization ####
RS <- TRUE
TS <- FALSE
RhoNames <- c("Neg", "Zero", "Pos")
nMax <- 1e3

#### 1. Two outcomes ####
Weights.Bi <- c(0.75,0.25)

Rho.Bi <- list(c(-0.20, 0.00, 0.20),
                c(-0.40, 0.00, 0.40))


Delta.Bi <- list(c(0.20, 0.10), 
                  c(0.30, 0.20))
ThetaE.Bi <- lapply(1:length(Delta.Bi), function(es) {0.50 + 1/2 * Delta.Bi[[es]]})
ThetaC.Bi <- lapply(1:length(Delta.Bi), function(es) {0.50 - 1/2 * Delta.Bi[[es]]})

PhiE.Bi <- lapply(1:length(Rho.Bi), function(es) {lapply(1:length(Rho), function(r) {
  tryCatch({Theta2Phi(ThetaE.Bi[[es]], Rho = diag(3)+Rho.Bi[[es]][r]-diag(Rho.Bi[[es]][r],3))
  }, error = function(e) {paste0("Error in effect size ", es, ", correlation structure ", RhoNames[r], ". \n
                                 Error message ", e)})
})
})

PhiC.Bi <- lapply(1:length(Rho.Bi), function(es) {lapply(1:length(Rho), function(r) {
  tryCatch({Theta2Phi(ThetaC.Bi[[es]], Rho = diag(3)+Rho.Bi[[es]][r]-diag(Rho.Bi[[es]][r],3))
  }, error = function(e) {paste0("Error in effect size ", es, ", correlation structure ", RhoNames[r], ". \n
                                 Error message ", e)})
})
})


RhoE.Bi <- lapply(1:length(Rho.Bi), function(es) {lapply(1:length(Rho), function(r) {ComputeRhoMultivariate(PhiE.Bi[[es]][[r]])})})
RhoC.Bi <- lapply(1:length(Rho.Bi), function(es) {lapply(1:length(Rho), function(r) {ComputeRhoMultivariate(PhiC.Bi[[es]][[r]])})})


#### 1.1 One sided ####
nAllBi.RS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAll(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nMax = nMax, one.sided = RS)})})
pwrAllBi.RS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAllBi.RS[[es]][[r]], one.sided = RS)
})})
pwrAllBi.RS_Ind <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAllBi.RS[[es]][[2]], one.sided = RS)
})})

nAnyBi.RS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAny(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nMax = nMax, one.sided = RS)})})
pwrAnyBi.RS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAnyBi.RS[[es]][[r]], one.sided = RS)
})})
pwrAnyBi.RS_Ind <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAnyBi.RS[[es]][[2]], one.sided = RS)
})})

nCompBi.RS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeComp(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], weights = Weights.Bi, nMax = nMax, one.sided = RS)})})
pwrCompBi.RS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], weights = Weights.Bi, nCompBi.RS[[es]][[r]], one.sided = RS)
})})
pwrCompBi.RS_Ind <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], weights = Weights.Bi, nCompBi.RS[[es]][[2]], one.sided = RS)
})})
#### 1.2 Two-sided ####
nAllBi.TS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAll(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nMax = nMax, one.sided = TS)})})
pwrAllBi.TS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAllBi.TS[[es]][[r]], one.sided = TS)
})})
pwrAllBi.TS_Ind <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAllBi.TS[[es]][[2]], one.sided = TS)
})})

nAnyBi.TS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAny(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nMax = nMax, one.sided = TS)})})
pwrAnyBi.TS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAnyBi.TS[[es]][[r]], one.sided = TS)
})})
pwrAnyBi.TS_Ind <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], nAnyBi.TS[[es]][[2]], one.sided = TS)
})})

nCompBi.TS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeComp(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], weights = Weights.Bi, nMax = nMax, one.sided = TS)})})
pwrCompBi.TS <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], weights = Weights.Bi, nCompBi.TS[[es]][[r]], one.sided = TS)
})})
pwrCompBi.TS_Ind <- lapply(1:length(Rho.Bi), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Bi[[es]][[r]], phiC = PhiC.Bi[[es]][[r]], weights = Weights.Bi, nCompBi.TS[[es]][[2]], one.sided = TS)
})})

nRS.Bi <- list(do.call(rbind, nAllBi.RS), do.call(rbind, nAnyBi.RS), do.call(rbind, nCompBi.RS))
nTS.Bi <- list(do.call(rbind, nAllBi.TS), do.call(rbind, nAnyBi.TS), do.call(rbind, nCompBi.TS))

pwrRS.Bi <- list(do.call(rbind, pwrAllBi.RS), do.call(rbind, pwrAnyBi.RS), do.call(rbind, pwrCompBi.RS))
pwrTS.Bi <- list(do.call(rbind, pwrAllBi.TS), do.call(rbind, pwrAnyBi.TS), do.call(rbind, pwrCompBi.TS))

pwrRS.Bi_Ind <- list(do.call(rbind, pwrAllBi.RS_Ind), do.call(rbind, pwrAnyBi.RS_Ind), do.call(rbind, pwrCompBi.RS_Ind))
pwrTS.Bi_Ind <- list(do.call(rbind, pwrAllBi.TS_Ind), do.call(rbind, pwrAnyBi.TS_Ind), do.call(rbind, pwrCompBi.TS_Ind))
#### 2. Three outcomes ####
Weights.Tri <- c(0.50,0.25, 0.25)
Rho.Tri <- list(c(-0.20, 0.00, 0.20),
                c(-0.40, 0.00, 0.40))
Phi111.E <- list(c(0.05, 0.20, 0.30),
                 c(0.01, 0.20, 0.40))
Phi111.C <- list(c(0.00, 0.10, 0.15),
                 c(0.00, 0.05, 0.10))

Delta.Tri <- list(c(0.20, 0.10, 0.20), 
                  c(0.30, 0.20, 0.30))
ThetaE.Tri <- lapply(1:length(Delta.Tri), function(es) {0.50 + 1/2 * Delta.Tri[[es]]})
ThetaC.Tri <- lapply(1:length(Delta.Tri), function(es) {0.50 - 1/2 * Delta.Tri[[es]]})

PhiE.Tri <- lapply(1:length(Rho.Tri), function(es) {lapply(1:length(Rho), function(r) {
  tryCatch({Theta2Phi(ThetaE.Tri[[es]], Rho = diag(3)+Rho.Tri[[es]][r]-diag(Rho.Tri[[es]][r],3), Phi111 = Phi111.E[[es]][r])
  }, error = function(e) {paste0("Error in effect size ", es, ", correlation structure ", RhoNames[r], ". \n
                                 Error message ", e)})
})
  })

PhiC.Tri <- lapply(1:length(Rho.Tri), function(es) {lapply(1:length(Rho), function(r) {
  tryCatch({Theta2Phi(ThetaC.Tri[[es]], Rho = diag(3)+Rho.Tri[[es]][r]-diag(Rho.Tri[[es]][r],3), Phi111 = Phi111.C[[es]][r])
  }, error = function(e) {paste0("Error in effect size ", es, ", correlation structure ", RhoNames[r], ". \n
                                 Error message ", e)})
})
})


RhoE.Tri <- lapply(1:length(Rho.Tri), function(es) {lapply(1:length(Rho), function(r) {ComputeRhoMultivariate(PhiE.Tri[[es]][[r]])})})
RhoC.Tri <- lapply(1:length(Rho.Tri), function(es) {lapply(1:length(Rho), function(r) {ComputeRhoMultivariate(PhiC.Tri[[es]][[r]])})})


#### 2.1 One sided ####
nAllTri.RS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAll(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nMax = nMax, one.sided = RS)})})
pwrAllTri.RS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAllTri.RS[[es]][[r]], one.sided = RS)
})})
pwrAllTri.RS_Ind <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAllTri.RS[[es]][[2]], one.sided = RS)
})})

nAnyTri.RS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAny(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nMax = nMax, one.sided = RS)})})
pwrAnyTri.RS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAnyTri.RS[[es]][[r]], one.sided = RS)
})})
pwrAnyTri.RS_Ind <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAnyTri.RS[[es]][[2]], one.sided = RS)
})})

nCompTri.RS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeComp(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], weights = Weights.Tri, nMax = nMax, one.sided = RS)})})
pwrCompTri.RS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], weights = Weights.Tri, nCompTri.RS[[es]][[r]], one.sided = RS)
})})
pwrCompTri.RS_Ind <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], weights = Weights.Tri, nCompTri.RS[[es]][[2]], one.sided = RS)
})})

#### 2.2 Two-sided ####
nAllTri.TS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAll(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nMax = nMax, one.sided = TS)})})
pwrAllTri.TS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAllTri.TS[[es]][[r]], one.sided = TS)
})})
pwrAllTri.TS_Ind <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAll(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAllTri.TS[[es]][[2]], one.sided = TS)
})})

nAnyTri.TS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeAny(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nMax = nMax, one.sided = TS)})})
pwrAnyTri.TS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAnyTri.TS[[es]][[r]], one.sided = TS)
})})
pwrAnyTri.TS_Ind <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerAny(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], nAnyTri.TS[[es]][[2]], one.sided = TS)
})})

nCompTri.TS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputeSampleSizeComp(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], weights = Weights.Tri, nMax = nMax, one.sided = TS)})})
pwrCompTri.TS <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], weights = Weights.Tri, nCompTri.TS[[es]][[r]], one.sided = TS)
})})
pwrCompTri.TS_Ind <- lapply(1:length(Rho.Tri), function(es) {sapply(1:length(Rho), function(r) {
  ComputePowerComp(phiE = PhiE.Tri[[es]][[r]], phiC = PhiC.Tri[[es]][[r]], weights = Weights.Tri, nCompTri.TS[[es]][[2]], one.sided = TS)
})})

nRS.Tri <- list(do.call(rbind, nAllTri.RS), do.call(rbind, nAnyTri.RS), do.call(rbind, nCompTri.RS))
nTS.Tri <- list(do.call(rbind, nAllTri.TS), do.call(rbind, nAnyTri.TS), do.call(rbind, nCompTri.TS))

pwrRS.Tri <- list(do.call(rbind, pwrAllTri.RS), do.call(rbind, pwrAnyTri.RS), do.call(rbind, pwrCompTri.RS))
pwrTS.Tri <- list(do.call(rbind, pwrAllTri.TS), do.call(rbind, pwrAnyTri.TS), do.call(rbind, pwrCompTri.TS))

pwrRS.Tri_Ind <- list(do.call(rbind, pwrAllTri.RS_Ind), do.call(rbind, pwrAnyTri.RS_Ind), do.call(rbind, pwrCompTri.RS_Ind))
pwrTS.Tri_Ind <- list(do.call(rbind, pwrAllTri.TS_Ind), do.call(rbind, pwrAnyTri.TS_Ind), do.call(rbind, pwrCompTri.TS_Ind))

#### 3. Create x-table ####
tabSampleSizesRS.Bi <- lapply(1:length(nRS.Bi), function(rule) {matrix(rbind(nRS.Bi[[rule]], pwrRS.Bi[[rule]], pwrRS.Bi_Ind[[rule]]), nrow = nrow(nRS.Bi[[rule]]))})
tabSampleSizesRS.Tri <- lapply(1:length(nRS.Tri), function(rule) {matrix(rbind(nRS.Tri[[rule]], pwrRS.Tri[[rule]], pwrRS.Tri_Ind[[rule]]), nrow = nrow(nRS.Tri[[rule]]))})
tabSampleSizesTS.Bi <- lapply(1:length(nTS.Bi), function(rule) {matrix(rbind(nTS.Bi[[rule]], pwrTS.Bi[[rule]], pwrTS.Bi_Ind[[rule]]), nrow = nrow(nTS.Bi[[rule]]))})
tabSampleSizesTS.Tri <- lapply(1:length(nTS.Tri), function(rule) {matrix(rbind(nTS.Tri[[rule]], pwrTS.Tri[[rule]], pwrTS.Tri_Ind[[rule]]), nrow = nrow(nTS.Tri[[rule]]))})

listSampleSizes <- Reduce(function(...) Map(rbind, ...), list(tabSampleSizesRS.Bi, tabSampleSizesRS.Tri))
tabSampleSizes <- do.call(rbind.data.frame, lapply(listSampleSizes, function(x) cbind.data.frame(
  paste0(rep(c("S", "L"), times = length(Delta.Bi)), rep(c(2,3), each = length(Delta.Bi))),
                                          do.call(rbind, lapply(1:nrow(x), function(y) {
                                          x[y,sort(c(seq(2,ncol(x),3), seq(3,ncol(x),3)), decreasing = FALSE)] <- sprintf("%.3f", x[y,sort(c(seq(2,ncol(x),3), seq(3,ncol(x),3)), decreasing = FALSE)])
                                          unlist(rbind.data.frame(NA, matrix(x[y,],nrow = 3)))
                                          })))))

AddToRow.SampleSizeEx <- list()
AddToRow.SampleSizeEx$pos <- as.list(c(0,seq(0,nrow(tabSampleSizes)-1, 4)))
AddToRow.SampleSizeEx$command <- c(" & & \\multicolumn{4}{l}{True $\\rho < 0$} & \\multicolumn{4}{l}{True $\\rho = 0$} & \\multicolumn{3}{l}{True $\\rho > 0$} \\\\\n",
                                   paste0("\\midrule  
                                          \\multicolumn{", ncol(tabSampleSizes), "}{c}{", RuleNames[-1], " rule} \\\\\n", 
                                          paste0(" DGM ", paste0(rep(c("& & n & $p_{\\rho_{T}}$ & $p_{\\rho_{U}}$"), length(Rho)), collapse = ""), "\\\\\n"),
                                          c("", rep("\\midrule \n", length(RuleNames)-2))))

Align.SampleSizeEx <- paste0("ll", paste0(rep("p{0.02cm}rrr", length(Rho)), collapse = ""), collapse = "")
 
print(xtable(tabSampleSizes, 
             caption= "Example of required sample sizes (n) for analysis with correlated data and anticipated probabilities to conclude superiority when sample size computations use the true correlation ($p_{\\rho_{T}}$) vs. assume uncorrelated dependent variables ($p_{\\rho_{U}}$) under four different data-generating mechanisms (DGM)",
             label="tab:SampleSize_example", align=Align.SampleSizeEx),
      add.to.row=AddToRow.SampleSizeEx, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)                                  
                                   
                                   
