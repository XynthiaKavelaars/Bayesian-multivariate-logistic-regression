This folder contains all necessary files to replicate the simulation study for the manuscript ”Bayesian multivariate logistic regression for superiority and inferiority decision-making under treatment heterogeneity." by Xynthia Kavelaars, Joris Mulder, and Maurits Kaptein.


########################################################################################################
#					Simulation 						       #
########################################################################################################
Description of files and folders:
\Simulation:				Folder containing subfolders with scripts, functions, simulation conditions.
RunSimulation.R				R-script to source all functions, run simulation, and create output in tables and figures
SampleParameters.R			R-script to sample parameters via Gibbs sampler
Evaluate.R				R-script to extract regression parameters from output and 
MakeTableFigure.R			R-script to create Latex tables and plots presented in manuscript
TextualOutput.R				R-script to create in-text results

\\Functions:				Folder containing functions used by scripts in folder \Simulation:
ComputeSampleSize.R			Script to compute sample sizes for all conditions using "FunctionsSampleSize.R"
DefineScenarios.R			Script to find parameters of defined data generating mechansms
FunctionsComputeSampleSize.R		Functions to compute sample sizes for various decison rules under a fixed design
FunctionsSimulation.R			Various functions to estimate parameters 
FunctionsTransformation.R		Various functions to perform transformations of data and parameters
VariableDefinitions.R			Various variable definitions

\\Workspaces:				Empty folder, where workspaces will be stored after running "RunSimulation.R"

\\Plots:				Empty folder, where figures will be stored after running "RunSimulation.R"

Instructions:
1. Set path to working directory to folder "\Simulation" in call "setwd('')" in "RunSimulation.R"
2. Execute "RunSimulation.R"
########################################################################################################
#					Application 						       #
########################################################################################################
Description of files and folders:
\Application:				Folder containing subfolders with scripts, functions, and dataset of the IST.
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
These scripts are based on parallel computation. The default number of cores is set at the number of cores of the PC minus one. 
To change the number of cores, change parameter(s) "nCores" in "\Simulation\SampleParameters.R" and "\Application\SampleParameters.R".

For any help with the files in this archive, please contact Xynthia Kavelaars (x.m.kavelaars@tilburguniversity.edu). 
