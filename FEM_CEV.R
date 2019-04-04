## Option pricing using Finite-Element Method (FEM)

## Packages

library(Matrix)
library(gaussquad)
library(pracma)
library(ggplot2)
library(latex2exp)
library(RQuantLib)

## Required functions

source("spTridiags.R")
source("Stiff.R")
source("BSFormula.R")

## Numerical approximation of CEV Model

n <- 2^10 - 1      # Number of spacial nodes
r <- 0.06          # Interest rate
Mat <- 1           # Maturity
sig <- 0.6         # Volatility
K <- 1             # Strike price
theta <- 0.5       # Theta scheme parameter

Rho <- c(0.35, 0.60, 0.85, 1)
Mu <- c(0, -0.2, -0.4, 0)
#Rho <- 0.6; Mu <- -0.2
Names <- c("R35", "R60", "R85", "R100")
IV <- list()

for (j in 1:length(Rho)){
  rho <- Rho[j]
  mu <- Mu[j]
  
  ## Check if mu is valid:
  if (rho >= 0 & rho < 0.5)
    stopifnot(mu == 0)
  if (rho >= 0.5 & rho < 1)
    stopifnot(mu > -0.5 & mu < 0.5 - rho)
  
  payoff <- function(x){ # Put option
    pmax(x - K, 0)
  }
  
  ## Define alpha, beta, gamma
  alpha <- function(x){
    sig^2 / 2 * x^(2 * (rho + mu))
  }
  beta <- function(x){
    sig^2 * (rho + mu) * x^(2 * (rho + mu) - 1) - 
      r * x^(2 * mu + 1)
  }
  gamma <- function(x){
    r * x^(2 * mu)
  }
  
  ## Bounded domain
  R <- 60
  s <- seq(0, R, length = n + 2) # Elements
  
  ## Time/space-stepping variables
  h <- R / (n + 1)    # space
  k <- h              # time
  M <- floor(Mat / k) # number of time nodes
  k <- Mat / M        # adjust time-stepping
  dof <- 1:(n + 1)    # degrees of freedom (non-boundary)
  
  u <- numeric(n + 2)      # Initialize vector
  u[dof] <- payoff(s[dof]) # Initial data
  
  ## Mass and Stiffness matrix
  Am <- Stiff(s, function(x) 0, function(x) 0, 
              function(x) x^(2 * mu))
  A <- Stiff(s, alpha, beta, gamma)
  
  ## Initialize
  B <- (Am + k * theta * A)
  C <- (Am - k * (1 - theta) * A)
  for (i in 1:M){
    u[dof] <- solve(B[dof, dof], C[dof, dof] %*% u[dof])
  }
  I <- (s > 0.1) & (s < 1.5)  # Domain of interest
  pos <- which(I)
  
  ## Implied volatility
  # IV[[j]] <- sapply(pos, function(i){
  #   EuropeanOptionImpliedVolatility(
  #     "call",
  #     value=u[i],
  #     underlying=s[i],
  #     strike=K,
  #     dividendYield=0,
  #     riskFreeRate=r,
  #     maturity=Mat,
  #     volatility=sig
  # )[1]})
  
  ## Fill data.fram
  if (j == 1)
    data <- data.frame(u[I])
  else{
    data <- data.frame(data, u[I])
  }
}
data$Payoff <- payoff(s[I])
data$x <- s[I]
colnames(data) <- c(Names, "Payoff", "x")

## Prices

p <- ggplot(data) +
  geom_line(aes(x = x, y = R35, color=Names[1])) +
  geom_line(aes(x = x, y = R60, color=Names[2])) +
  geom_line(aes(x = x, y = R85, color=Names[3])) +
  geom_line(aes(x = x, y = R100, color=Names[4])) +
  geom_line(aes(x = x, y = Payoff, color="Payoff"))

p <- p +
  scale_color_manual(name = "Method",
                     values=c(R35="blue",
                              R60="red",
                              R85="green",
                              R100="purple",
                              Payoff="black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Underlying price") +
  ylab("Option/payoff value")

## Implied volatility

# data_IV <- data.frame(x = s[I],
#                    IV35 =  IV[[1]],
#                    IV60 =  IV[[2]],
#                    IV85 =  IV[[3]],
#                    IV100 = IV[[4]])
# 
# IVplot <- ggplot(data_IV) +
#   geom_line(aes(x = x, y = IV35, color="IV35")) +
#   geom_line(aes(x = x, y = IV60, color="IV60")) +
#   geom_line(aes(x = x, y = IV85, color="IV85")) +
#   geom_line(aes(x = x, y = IV100, color="IV100"))
# 
# IVplot <- IVplot +
#   scale_color_manual(name = "Implied Vol",
#                      values=c(IV35="blue",
#                               IV60="red",
#                               IV85="green",
#                               IV100="purple")) +
#   xlab("Underlying price") +
#   ylab("Implied volatility")

## Plot price + IV

# grid.arrange(gg, IVplot, ncol=2)



