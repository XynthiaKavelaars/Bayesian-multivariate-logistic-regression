# Variable Definitions

# Simulation parameters 
nSim <- 1e3+1e2 # Number of simulations to run
nSim_eval <- 1e3  # Number of simulations to evaluate
nBurn <- 1e3			# No. of burnin iterations
nIt <- 1e4			# No. of iterations
Start <- c(0,0.5)		# Starting values for chains 1 and 2
GR.Cut <- 1.10			# Cutoff value for (non-)convergence

# Prior hyperparameters
b0 <- 0       # Prior mean of regression coefficient (univariate and multivariate logistic regression)
B0 <- 1e-1    # Prior precision of regression coefficient (univariate and multivariate logistic regression)
a0 <- 1e-2    # Prior observation per response category (multivariate Bernoulli)

# Covariate data
P <- 4		  	# Number of covariates
nDgm <- 8 # Number of data generating mechanisms
Rho <- c(-0.20,0,0.20)
names(Rho) <- c("Neg", "Zero", "Pos")
set.seed(2021)
Discrete <- TRUE

K <- 1:3;			  # Number of outcome variables
Q <- 2^(K);	# Number of joint response categories

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
                               Intra_Hi = c(1))
Values[["Discrete"]] <- list(Trial = c(0.5),
                             Intra_Lo = c(0), 
                             Intra_Hi = c(1))

pX <- list()
pX[["Trial"]] <- pXD
pX[["Intra_Lo"]] <- pXD[1]
pX[["Intra_Hi"]] <- pXD[2]

#### 3. Decision-making ####
# Weights for compensatory rule
Weights <- list(1, 
                c(0.75, 0.25),
                c(0.50,0.25,0.25))

nMax <- 1e3 # Maximum sample size per treatment group
alpha <- 0.05 # Type I error rate
pwr <- 0.80 # Statistical power
Rules <- c("Single", "All", "Any", "Comp")
RuleNames <- c("Single", "All", "Any", "Compensatory")

# Cutoff for superiority
pCut <- lapply(K, function(k){
  x <- matrix(NA, nrow = length(Rules), ncol = 2)
  x[which(Rules == "Single"),] <- c(0, 1 - alpha);
  x[which(Rules == "All"),] <- c(0, 1 - alpha);
  x[which(Rules == "Any"),] <- c(0, 1 - (alpha/k));
  x[which(Rules == "Comp"),] <- c(0, 1 - alpha);
  return(x)})

Populations <- c("Trial", "Intra_Lo")

Seeds <- array(sample(1:1e4, nDgm * max(K) * length(Rho) * length(Rules), replace = FALSE),
                dim = c(nDgm, max(K), length(Rho), length(Rules)))
