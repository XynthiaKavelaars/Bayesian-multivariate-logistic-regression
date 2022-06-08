#### Data Generating Mechanisms ####
Rho <- Rho_Lo <- Rho_Hi <- -0.20

#### 1. Delta = 0 - homogeneous delta ####
Theta1E_Lo <- c(0.35,0.45)
Theta1C_Lo <- c(0.35,0.45)
Theta1E_Hi <- c(0.65,0.55)
Theta1C_Hi <- c(0.65,0.55)

TrueBeta1.1 <- FindTrueBeta(Theta1C_Lo, Theta1E_Lo, Theta1C_Hi, Theta1E_Hi, Rho, Rho, xLo = Values[["Discrete"]][["Intra_Lo"]], xHi = Values[["Discrete"]][["Intra_Hi"]])
TrueBeta1.2 <- FindTrueBeta(Theta1C_Lo, Theta1E_Lo, Theta1C_Hi, Theta1E_Hi, Rho, Rho, xLo = Values[["Continuous"]][["Intra_Lo"]], xHi = Values[["Continuous"]][["Intra_Hi"]])

#### 2. Delta = 0 - heterogeneous delta ####
Theta2E_Lo <- c(0.625,0.575)
Theta2C_Lo <- c(0.375,0.425)
Theta2E_Hi <- c(0.375,0.425)
Theta2C_Hi <- c(0.625,0.575)

TrueBeta2.1 <- FindTrueBeta(Theta2C_Lo, Theta2E_Lo, Theta2C_Hi, Theta2E_Hi, Rho, Rho, xLo = Values[["Discrete"]][["Intra_Lo"]], xHi = Values[["Discrete"]][["Intra_Hi"]])
TrueBeta2.2 <- FindTrueBeta(Theta2C_Lo, Theta2E_Lo, Theta2C_Hi, Theta2E_Hi, Rho, Rho, xLo = Values[["Continuous"]][["Intra_Lo"]], xHi = Values[["Continuous"]][["Intra_Hi"]])

#### 3. Delta > 0 - heterogeneous delta ####
Theta3E_Lo <- c(0.70,0.65)
Theta3C_Lo <- c(0.30,0.35)
Theta3E_Hi <- c(0.45,0.40)
Theta3C_Hi <- c(0.55,0.60)

TrueBeta3.1 <- FindTrueBeta(Theta3C_Lo, Theta3E_Lo, Theta3C_Hi, Theta3E_Hi, Rho, Rho, xLo = Values[["Discrete"]][["Intra_Lo"]], xHi = Values[["Discrete"]][["Intra_Hi"]])
TrueBeta3.2 <- FindTrueBeta(Theta3C_Lo, Theta3E_Lo, Theta3C_Hi, Theta3E_Hi, Rho, Rho, xLo = Values[["Continuous"]][["Intra_Lo"]], xHi = Values[["Continuous"]][["Intra_Hi"]])

#### 4. Delta > 0 - heterogeneous delta ####
Theta4E_Lo <- c(0.60,0.50)
Theta4C_Lo <- c(0.40,0.50)
Theta4E_Hi <- c(0.80,0.50)
Theta4C_Hi <- c(0.20,0.50)

TrueBeta4.1 <- FindTrueBeta(Theta4C_Lo, Theta4E_Lo, Theta4C_Hi, Theta4E_Hi, Rho, Rho, xLo = Values[["Discrete"]][["Intra_Lo"]], xHi = Values[["Discrete"]][["Intra_Hi"]])
TrueBeta4.2 <- FindTrueBeta(Theta4C_Lo, Theta4E_Lo, Theta4C_Hi, Theta4E_Hi, Rho, Rho, xLo = Values[["Continuous"]][["Intra_Lo"]], xHi = Values[["Continuous"]][["Intra_Hi"]])

#### TrueBeta: list of true regression coefficients for all DGMs ####
TrueBeta <- list(TrueBeta1.1 = TrueBeta1.1, TrueBeta1.2 = TrueBeta1.2,
                 TrueBeta2.1 = TrueBeta2.1, TrueBeta2.2 = TrueBeta2.2,
                 TrueBeta3.1 = TrueBeta3.1, TrueBeta3.2 = TrueBeta3.2,
                 TrueBeta4.1 = TrueBeta4.1, TrueBeta4.2 = TrueBeta4.2)

# Names of data generating mechanisms
DgmNames <- paste(rep(1:(nDgm/2), each = 2), rep(c(".1", ".2"), times = length(TrueBeta)/2), sep = "")[1:nDgm]