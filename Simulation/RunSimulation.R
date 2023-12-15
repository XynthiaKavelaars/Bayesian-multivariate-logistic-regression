#### Run simulation ####
rm(list = ls())
#setwd('')
wd <- getwd()

#### 0. Initialization ####
# Load packages after installing if needed
RequiredPackages <- c('MCMCpack', 'coda', 'pgdraw', 'TruncatedNormal', 'pracma', 'msm', 'abind', 'foreach', 'doParallel', 'nnet', 'xtable')
package.check <- lapply(RequiredPackages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }})

# Initiate parallel computing
nCores <- detectCores() - 1
registerDoParallel(cores = nCores)

# Load functions 
source("Functions/Functions_Simulation.R")
source("Functions/Functions_SampleSize.R")
source("Functions/Functions_Transformation.R")

# Load variable definitions and data generating mechanisms
source("Functions/VariableDefinitions.R")
source("Functions/DefineDataGeneratingMechanisms.R")

# Compute required sample sizes 
source("Functions/ComputeSampleSizes.R")

#### 1. Run simulation ####
source("SimulateMvLR.R")
source("SimulateUvLR.R")
source("SimulateMvB.R")

#### 2. Evaluate results ####
source("Evaluate.R")

#### 3. Make tables and figures ####
source("MakeTable.R")

#### 4. Compute power example ####
source("Functions/Functions_SampleSize.R")
source("Functions/Functions_Power.R")
source("ComputePowerExample.R")
