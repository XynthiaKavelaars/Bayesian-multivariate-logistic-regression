load("Workspaces/Application3.RData")

pReject.3 <- cbind.data.frame(Types[,c("Methods", "Populations")], 
                              do.call(rbind, sapply(Decisions.3, function(x) unlist(x[["Bias"]]), simplify = FALSE, USE.NAMES = TRUE)),
                              do.call(rbind, sapply(Decisions.3, function(x) unlist(x[["Pop"]]), simplify = FALSE, USE.NAMES = TRUE)),
                              do.call(rbind, sapply(Decisions.3, function(x) c(x[["Decision"]][grep(".LS", names(x[["Decision"]]))]), simplify = FALSE, USE.NAMES = TRUE)))

pReject.3[pReject.3 == "MvB"] <- "Reference"

#### 1 Make Latex table ####

tab.pReject.3 <- data.frame(matrix(NA, nrow = nrow(Types), ncol = 10))

# Layout table, with brackets etc.
for(i in 1:nrow(pReject.3)){
  tab.pReject.3[i,] <-
    unlist(c(pReject.3[i,c("Methods")], NA,
             paste0("(", paste0(sprintf("%6.3f", pReject.3[i,grep("BiasM", colnames(pReject.3))]), collapse = ", "), ")", collapse = ""),
             paste0("(", paste0(sprintf("%5.3f", pReject.3[i,grep("PopM", colnames(pReject.3))]), collapse = ", "), ")", collapse = ""),
             if(max(pReject.3[i,grep("PopM", colnames(pReject.3))]) > (1 - Alpha / 4) 
                & min(pReject.3[i,grep("PopM", names(pReject.3))]) < (Alpha / 4)){"$\\bm{<} \\& \\bm>}$"
             } else if(max(pReject.3[i,grep("PopM", colnames(pReject.3))]) > (1 - Alpha / 4)){"$\\bm{<}$"
             } else if(min(pReject.3[i,grep("PopM", names(pReject.3))]) < (Alpha / 4)){"$\\bm{>}$"
             }  else {"-"},
             if(min(pReject.3[i,grep("PopM", colnames(pReject.3))]) > (1 - Alpha / 2)){"$\\bm{<}$"
             } else if(max(pReject.3[i,grep("PopM", names(pReject.3))]) < (Alpha / 2)){"$\\bm{>}$"
             } else {"-"},
             NA,
             sprintf("%6.3f", pReject.3[i,grep("BiasW", colnames(pReject.3))]),
             sprintf("%5.3f", pReject.3[i,grep("PopW", colnames(pReject.3))]),
             if(pReject.3[i,grep("PopW", colnames(pReject.3))] > (1 - Alpha / 2)){"$\\bm{<}$"
             } else if(pReject.3[i,grep("PopW", colnames(pReject.3))] < (Alpha / 2)){"$\\bm{>}$"} else {"-"}))
}

#### 1.1 Range of blood pressure ####
tab.pReject.Range.3 <- tab.pReject.3[tab.pReject.3[,"X1"] %in% c("Reference", "Empirical", "Analytical"),]
AddToRow.pReject.Range.3 <- list()
AddToRow.pReject.Range.3[["pos"]] <- as.list(cumsum(c(0, sapply(Populations, function(pop) sum(Types.Range[,"Populations"] == pop)))))
AddToRow.pReject.Range.3[["command"]] <- c(paste0("\\midrule \n \\multicolumn{4}{l}{\\textbf{", names.Populations, "} (",
                                                  #"} & \\multicolumn{3}{l}{", 
                                                  range.Populations[[3]], ") 
                                                  } & \\multicolumn{", ncol(tab.pReject.Range.3) - 4, "}{l}{$n_{H+A} = ", nRange[["RSBP_C"]][,"nE"], "$, $n_{A} = ", nRange[["RSBP_C"]][,"nC"],"$} \\\\\n  Method  & & \\multicolumn{1}{l}{$\\bm{\\delta}(Bp)$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[1:2], collapse = ""), 
                                                  " & & \\multicolumn{1}{l}{$\\delta(\\bm{w},Bp)$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[3], collapse = ""), " \\\\\n \\midrule \n"),
                                           paste0("\\midrule \n
                                           \\multicolumn{", ncol(tab.pReject.Range.3), "}{l}{pp = posterior probability, Comp = Compensatory rule} \\\\ \n
                                           \\multicolumn{", ncol(tab.pReject.Range.3), "}{l}{$\\bm{>}$ = superiority concluded, $\\bm{<}$ = inferiority concluded} \\\\ \n
                                          ", collapse = ""))

Align.pReject.Range.3 <- paste0("ll", paste0(rep("r", ncol(tab.pReject.Range.3) - 1), collapse = ""))

print(xtable(tab.pReject.Range.3, 
             caption="Average and conditional treatment differences in the IST-data, by range of blood pressure.",
             label="tab:Application.Range_rbsp", align=Align.pReject.Range.3),
      add.to.row=AddToRow.pReject.Range.3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE)                                   


#### 1.2 Values of blood pressure ####
names.Values.3 <- paste0(rep(c("-", "+"), each = length(Values[["RSBP_C"]]) / 2), c((length(Values[["RSBP_C"]]) / 2):1, 1:(length(Values[["RSBP_C"]]) / 2)), " SD")# (" , 

tab.pReject.Value.3 <- cbind.data.frame(names.Values.3, subset(tab.pReject.3[tab.pReject.3[,"X1"] %in% c("Value"),-1]))

AddToRow.pReject.Value.3 <- list()
AddToRow.pReject.Value.3[["pos"]] <- as.list(c(0,nrow(tab.pReject.Value.3)))# sapply(Populations[-length(Populations)], function(pop) sum(Types.Value[,"Populations"] == pop)))))
AddToRow.pReject.Value.3[["command"]] <- c(paste0("Value  & & \\multicolumn{1}{l}{$\\bm{\\delta} (\\bm{Bp})$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[1:2], collapse = ""), 
                                                  " & & \\multicolumn{1}{l}{$\\delta(\\bm{w}, Bp)$} & \\multicolumn{1}{l}{pp} ", paste0(" & ", names.Rules[3], collapse = ""), " \\\\\n \\midrule \n"),
                                           paste0("\\midrule \n
                                                   \\multicolumn{", ncol(tab.pReject.Value.3), "}{l}{pp = posterior probability, Comp = Compensatory} \\\\ \n
                                                  \\multicolumn{", ncol(tab.pReject.Value.3), "}{l}{$\\bm{>}$ = superiority concluded, $\\bm{<}$ = inferiority concluded} \\\\ ", collapse = ""))

Align.pReject.Value.3 <- paste0("ll", paste0(rep("r", ncol(tab.pReject.Value.3) - 1), collapse = ""))

print(xtable(tab.pReject.Value.3, 
             caption="Conditional treatment differences in the IST-data, by value of blood pressure.",
             label="tab:Application.Value_rbsp", align=Align.pReject.Value.3),
       add.to.row=AddToRow.pReject.Value.3, caption.placement="top", 
      include.colnames=FALSE, include.rownames=FALSE,
      table.placement="htbp", sanitize.text.function = identity,
      NA.string = " ", booktabs=TRUE) 




#### 2 Make Figure ####
StandardizedData <- cbind(Data[,c("TRT")], (Data[,"RSBP"] - mean(Data[,"RSBP"])) / sd(Data[,"RSBP"]), Y12)
colnames(StandardizedData) <- c("TRT", "RSBP", colnames(Y12))


xPoints <- seq(-3,3,0.25)
yPoints.E <- do.call(rbind, lapply(1:(length(xPoints)-1), function(x) {
  ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 1)
  phi <- rdirichlet(1e4, colSums(Y12[ind,]) + Alpha0)
  theta <- colMeans(cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)])))
  return(theta)
}))

yPoints.C <- do.call(rbind, lapply(1:(length(xPoints)-1), function(x) {
  ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 0)
  phi <- rdirichlet(1e4, colSums(Y12[ind,]) + Alpha0)
  theta <- colMeans(cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)])))
  return(theta)
}))


thetaEmp.E <- Reduce("+", lapply(1:nIt, function(i){
  do.call(rbind, lapply(1:(length(xPoints)-1), function(x) {
    ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 1)
    phi <- exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]) / rowSums(exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]))
    theta <- colMeans(cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)])))
    return(theta)
  }))
})) / nIt

thetaEmp.C <- Reduce("+", lapply(1:nIt, function(i){
  do.call(rbind, lapply(1:(length(xPoints)-1), function(x) {
    ind <- which(StandardizedData[,"RSBP"] > xPoints[x] & StandardizedData[,"RSBP"] < xPoints[x+1] & StandardizedData[,"TRT"] == 0)
    phi <- exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]) / rowSums(exp(X3[ind,] %*% Pars.3[["bDrawPG"]][[i]]))
    theta <- colMeans(cbind(rowSums(phi[,c(1,2)]), rowSums(phi[,c(1,3)])))
    return(theta)
  }))
})) / nIt


thetaAna.E <- do.call(rbind, lapply(1:(length(xPoints)-1), function(x) {
  theta <- EstimateThetaAnalytical(BetaPG = Pars.3[["bDrawPG"]], X = X3, Trt = 1, RangeX = xPoints[x+0:1] * sdBp)
  return(colMeans(do.call(rbind, theta)))
}))

thetaAna.C <- do.call(rbind, lapply(1:(length(xPoints)-1), function(x) {
  theta <- EstimateThetaAnalytical(BetaPG = Pars.3[["bDrawPG"]], X = X3, Trt = 0, RangeX = xPoints[x+0:1] * sdBp)
  return(colMeans(do.call(rbind, theta)))
}))

save(thetaEmp.E, thetaEmp.C, thetaAna.E, thetaAna.C, file = "Workspaces/Application_figure.RData")

load("Workspaces/Application_figure.RData")

postscript("Plots/Application.eps", width = 6.5, height = 4, horizontal = FALSE, font = "Times")
layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE), heights = c(5,2))
par(mar = c(4,4,2,1), family = "Times")
plot(x = xPoints[-1] * sdBp + meanBp, y = yPoints.E[,1] - yPoints.C[,1], xlim = c(RoundChoose(min(xPoints[-1] * sdBp + meanBp),10,0), RoundChoose(max(xPoints[-1] * sdBp + meanBp),10,1)), ylim = c(-0.5,0.5), 
     type = "l", main = "Recurring stroke", xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
lines(x = xPoints[-1] * sdBp + meanBp, y = thetaEmp.E[,1] - thetaEmp.C[,1], lty = 2)
lines(x = xPoints[-1] * sdBp + meanBp, y = thetaAna.E[,1] - thetaAna.C[,1], lty = 3)

plot(x = xPoints[-1] * sdBp + meanBp, y = yPoints.E[,2] - yPoints.C[,2], xlim = c(RoundChoose(min(xPoints[-1] * sdBp + meanBp),10,0), RoundChoose(max(xPoints[-1] * sdBp + meanBp),10,1)), ylim = c(-0.5,0.5), 
     type = "l", main = "Dependency",  xlab = "Blood pressure", ylab = bquote(theta[H+A] - theta[A]), las = 1)
lines(x = xPoints[-1] * sdBp + meanBp, y = thetaEmp.E[,2] - thetaEmp.C[,2], lty = 2)
lines(x = xPoints[-1] * sdBp + meanBp, y = thetaAna.E[,2] - thetaAna.C[,2], lty = 3)

par(mar = rep(0.01, 4))
plot.new()
legend(x = "center", legend = c("Reference", "Empirical", "Analytical"), lty = 1:3, ncol = 3) #col = c("black", "red", "darkgreen"))
dev.off()

