#### 1. Sample posterior regression coefficients ####
for(dgm in 1:nDgm){
  for(ss in 1:length(nSample[[dgm]])){
    Pars <- SamplePostBeta(nSim = nSim, TrueBeta = TrueBeta[[dgm]], n = nSample[[dgm]][ss], 
                           pXD = pXD, 
                           MuX = MuX, SigmaX = SigmaX, RangeX = Ranges[[MeasurementLevel[dgm]]][["Trial"]],
                           Covariate = TRUE, MeasurementLevel = MeasurementLevel[dgm], Interaction = Interaction, 
                           nBurn = nBurn, nIt = nIt, Start = Start, b0 = b0, B0 = B0,
                           Seed = Seeds[dgm,ss], Method = MethodFit, GenerateMult = TRUE, Parallel = TRUE)
    
    save(Pars, file = paste("Workspaces/Pars_DGM", dgm, "_n", names(nSample[[dgm]])[ss], ".RData", sep=""))
    
  }}