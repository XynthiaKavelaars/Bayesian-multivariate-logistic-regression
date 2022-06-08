#### Variable definitions ####

#### 1. Simulation parameters ####
nSim <- 1e3 + 1e2		# Number of samples to run
nSim_Eval <- 1e3		# Number of samples to evaluate (lower than nSim, so non-converging samples can be excluded)
nMax <- 2e3			# Maximum sample size per treatment group

#### 2. Covariate data ####
Interaction <- TRUE		# Include interaction between x and treatment
P <- 4				# Number of covariates
nDgm <- 8			# Number of data generating mechanisms
set.seed(2021)

MeasurementLevel <- rep(c("Discrete", "Continuous"), 1/2 * nDgm)
Q <- 4;				# Number of joint response categories
K <- log2(Q);			# Number of outcome variables

#### 2.1 Covariate data - Discrete ####
pXD <- c(0.5,0.5)		# Probability of x

#### 2.2 Covariate data - Continuous ####
MuX <- 0			# Mean
SigmaX <- 1			# Standard deviation

#### 2.3 Define subpopulations ####
# Range defining subpopulations
Ranges <- list()
Ranges[["Continuous"]] <- list(Trial = c(-Inf, Inf),
                               Intra_Lo = c(-1, 0),
                               Intra_Hi = c(0, 1))
Ranges[["Discrete"]] <- list(Trial = c(0, 1),
                             Intra_Lo = c(0),
                             Intra_Hi = c(1))

# Values defining subpopulations
Values <- list()
Values[["Continuous"]] <- list(Trial = c(0),
				Intra_Lo = c(-1),
				Intra_Hi = c( 1))
Values[["Discrete"]] <- list(Trial = c(0.5),
				     Intra_Lo = c(0),
				     Intra_Hi = c(1))

# Probability of discrete covariate
pX <- list()
pX[["Trial"]] <- pXD
pX[["Intra_Lo"]] <- pXD[1]
pX[["Intra_Hi"]] <- pXD[2]

#### 3. Parameters Gibbs sampler ####
nBurn <- 1e3			# No. of burnin iterations
nIt <- 1e3#1e4			# No. of iterations
Start <- c(0,0.5)		# Starting values for chains 1 and 2
GR.Cut <- 1.10			# Cutoff value for (non-)convergence

# Prior parameters 
b0 <- 0;			# Prior mean regression coefficient
B0 <- 1e-2;			# Prior precision of regression coefficients
A0 <- rep(0,Q)			# Prior parameters multivariate Bernoulli distribution

#### 3. Parameter superiority decision ####
weights <- c(0.50,0.50)		# Weights compensatory rule
Alpha <- 0.05			# Type I error rate

#### 4. Types of populations ####
# Specify conditions per data generating mechanism as function of the measurement level,  method to estimate parameters, population of interest, and scope of data
MeasurementLevels <- c("Discrete", "Continuous")
Methods <- c("MvB", "Value", "Empirical", "Analytical")
Populations <- c("Trial", "Intra_Lo")
Scopes <- c("Old")
MethodFit <- c("PG")
TypesSpace <- as.data.frame(t(expand.grid(MeasurementLevels, Methods, Populations, Scopes, stringsAsFactors = FALSE)))
Exclude <- lapply(TypesSpace, function(x) 
  any(c("Intra_Hi", "Extra_Lo", "Extra_Hi") %in% x) | 
    any(c("New") %in% x) |
    all(c("Discrete", "Analytical", "Intra_Lo") %in% x) | 
    all(c("Discrete", "Empirical", "Intra_Lo") %in% x))
 

Types <- t(TypesSpace[,which(do.call(rbind, Exclude) == FALSE)])
colnames(Types) <- rownames(TypesSpace) <- c("MeasurementLevels", "Methods", "Populations", "Scopes")

Seeds <- matrix(sample(1:1e4, 9 * nDgm, replace = FALSE),
nrow = nDgm)

Rules <- c("Any", "All", "Compensatory")		# Name rules
TestSide <- c("Comp_rs", "All_rs", "Any_rs")		# Names of performed tests

