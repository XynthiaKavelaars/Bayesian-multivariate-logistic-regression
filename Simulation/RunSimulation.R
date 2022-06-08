#### Generate data ####
rm(list=ls())

#### 0. Initialization ####
#### 0.1 Load packages ####
packages <- c('MCMCpack', 'coda', 'pgdraw', 'TruncatedNormal', 'pracma', 'msm', 'abind', 'foreach', 'doParallel', 'nnet', 'xtable')

package.check <- lapply(packages, function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }})

#### 0.2 Initialize parallellization ####
nCores <- detectCores()
cl <- makeCluster(nCores - 1)
registerDoParallel(cl)


#### 0.3 Set working directory ####
setwd('')

#### 0.4 Load functions ####
source("Functions/Functions_Transformation.R")
source("Functions/Functions_Simulation.R")
source("Functions/Functions_SampleSize.R")
source("Functions/VariableDefinitions.R")
source("Functions/DataGeneratingMechanisms.R")
source("Functions/ComputeSampleSize.R")


#### 1. Run simulation ####
source("SampleParameters.R")
source("Evaluate.R")
source("MakeTableFigure.R")
source("TextualOutput.R")

