This folder contains all necessary files to replicate the simulation study for the manuscript ”Bayesian multivariate logistic regression for superiority and inferiority decision-making under treatment heterogeneity." by Xynthia Kavelaars, Joris Mulder, and Maurits Kaptein.


########################################################################################################
#					Simulation 						       #
########################################################################################################
Description of files and folders:
\Simulation:				Folder containing subfolders with scripts, functions, simulation conditions.
Simulation.zip:				Zipped folder with all files and (empty) folders in folder \Simulation. This folder was added since (empty) (sub)folders to store workspaces are not uploaded to GitHub.
RunSimulation.R				R-script to source all functions, run simulation, and create output in tables and figures
SampleParameters_mLR.R			R-script to sample multivariate logistic regression parameters via Gibbs sampler
SampleParameters_uLR.R			R-script to sample univariate logistic parameters via Gibbs sampler
Simulate_MvB.R				R-script to run simulationm of multivariate Bernoulli model
Evaluate.R				R-script to extract regression parameters from output and 
MakeTableFigure.R			R-script to create Latex tables and plots presented in manuscript
TextualOutput.R				R-script to create in-text results
ComputePowerExample.R			R-script to compute illustrative example of power, sample size, and effect size

\\Functions:				Folder containing functions used by scripts in folder \Simulation:
ComputeSampleSizes.R			Script to compute sample sizes for all conditions using "FunctionsSampleSize.R"
DefineDataGeneratingMechanisms.R	Script to find parameters of defined data generating mechanisms
FunctionsPower.R			Functions to compute sample sizes and power, used for illustrative example
FunctionsSampleSize.R			Functions to compute sample sizes for various decison rules
FunctionsSimulation.R			Various functions to estimate parameters 
FunctionsTransformation.R		Various functions to perform transformations of data and parameters
VariableDefinitions.R			Various variable definitions

\\Workspaces:				Empty folder, where workspaces will be stored in subfolder after running "RunSimulation.R"

\\\Pars:				Empty folder, where regression parameters will be stored in subfolder

\\\\Pars_K2_mLR:			Empty folder, where multivariate logistic regression parameters will be stored for two-dimensional outcome variable
\\\\Pars_K2_uLR:			Empty folder, where univariate logistic regression parameters will be stored for two-dimensional outcome variable
\\\\Pars_K3:				Empty folder, where multivariate regression parameters will be stored for three-dimensional outcome variable

\\\Data:				Empty folder, where covariate data will be stored in subfolder

\\\\Data_K2:				Empty folder, where covariate data will be stored for multivariate logistic regression with two-dimensional outcome variable
\\\\Data_K2:				Empty folder, where covariate data will be stored for univariate logistic regression with two-dimensional outcome variable
\\\\Data_K3:				Empty folder, where covariate data will be stored for multivariate logistic regression with three-dimensional outcome variable

\\\Results:				Empty folder, where simulation results will be stored in subfolder

\\\\Results_K2:				Empty folder, where simulation results will be stored for two-dimensional outcome variable
\\\\Results_K3:				Empty folder, where simulation results will be stored for three-dimensional outcome variable

\\Plots:				Empty folder, where figures will be stored after running "RunSimulation.R"

Instructions:
1. Set path to working directory to folder "\Simulation" in call "setwd('')" in "RunSimulation.R"
2. Execute "RunSimulation.R"
########################################################################################################
#					Application 						       #
########################################################################################################
Description of files and folders:
\Application:				Folder containing subfolders with scripts, functions, and dataset of the IST.
Application.zip:				Zipped folder with all files and (empty) folders in folder \Application. This folder was added since (empty) (sub)folders to store workspaces are not uploaded to GitHub.
RunApplication.R			R-script to run application
SampleParameters.R			R-script to sample parameters via Gibbs sampler
MakeTableFigure.R			R-script to create Latex tables and plots presented in manuscript
TextualOutput.R				R-script to create in-text results

\\Functions:				Folder containing functions used by scripts in folder \Application:
FunctionsApplication.R			Functions to run application, used by "RunApplication.R"
VariableDefinitions.R			Various variable definitions

\\Workspaces:				Empty folder, where workspaces will be stored after running "RunApplication.R"

\\Plots:				Empty folder, where figures will be stored after running "RunApplication.R"

Instructions:
1. Set path to working directory to folder "\Application" in call "setwd('')" in "RunApplication.R"
2. Download data from International Stroke Trial from https://doi.org/10.1186/1745-6215-12-101 to folder "\Application"
3. Name the dataset "IST_Data.csv".
4. Execute "RunApplication.R"
########################################################################################################
#					Further information					       #
########################################################################################################
Warning: This simulation takes a long time to run (a few months over multiple cores) and requires a lot of memory!

These scripts are based on parallel computation. The default number of cores is set at the number of cores of the PC minus one. 
To change the number of cores, change parameter(s) "nCores" in "\Simulation\SampleParameters.R" and "\Application\SampleParameters.R".

For any help with the files in this archive, please contact Xynthia Kavelaars (xynthia.kavelaars@ou.nl). 
