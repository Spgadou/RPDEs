## Option pricing using Finite-Element Method (FEM)

## Packages

library(Matrix)
library(gaussquad)
library(pracma)
library(ggplot2)
library(latex2exp)

## Required functions

source("spTridiags.R")
source("Stiff.R")
source("BSFormula.R")

## Numerical approximation of Black-Scholes Model

n <- 2^10 - 1       # Number of spacial nodes
r <- 0.05          # Interest rate
Mat <- 5           # Maturity
sig <- 0.6         # Volatility
K <- 2             # Strike price
theta <- 0.5       # Theta scheme parameter
DoI <- c(0.1, 6)   # Domain of interest

payoff <- function(x){ # Payoff function
  pmax(exp(x) - K, 0)
}
solution <- function(s, t){
  BSFormula(s, r, t, K, sig, type = "call")
}

## Define alpha, beta, gamma
alpha <- function(x){
  sig^2 / 2
}
beta <- function(x){
  (sig^2 / 2 - r)
}
gamma <- function(x){
  r
}

## Bounded domain
a <- -10
b <- 10
x <- seq(a, b, length = n + 2) # Elements
s <- exp(x)

## Time/space-stepping variables
h <- x[2] - x[1]    # space
k <- h              # time
M <- floor(Mat / k) # number of time nodes
k <- Mat / M        # adjust time-stepping
dof <- 2:(n + 1)    # degrees of freedom (non-boundary)

# Problem:
#
# d(u, v)/dt + a(u, v) = 0,     on JxG
# u(0, x) = g(x),               on G
# u(t, x) = 0,                  on JxdG

u <- matrix(0, ncol = M+1, nrow = n+2)
u[dof,1] <- payoff(x[dof]) # Initial data

#u <- numeric(n + 2)      # Initialize vector
#u[dof] <- payoff(x[dof]) # Initial data

## Mass and Stiffness matrix
Am <- (h / 6) * spTridiags(rep(1, n+1),
                           rep(4, n+2),
                           rep(1, n+1))
A <- Stiff(x, alpha, beta, gamma)

## Initialize
B <- (Am + k * theta * A)
C <- (Am - k * (1 - theta) * A)
for (i in 1:M){
  u[dof,i+1] <- as.numeric(solve(B[dof, dof], C[dof, dof] %*% u[dof,i]))
}
time <- seq(0, Mat, k)
u <- data.frame(u)
colnames(u) <- time

I <- (s > DoI[1]) & (s < DoI[2])                          # Domain of interest
sol <- solution(s[I], Mat)    # True sol.
error <- sum((sol - u[I, M+1])^2 * h^2)                # L2 error

## Plot price

FEM <- u[I, M+1]
BS <- sol
POF <- payoff(x[I])

data_FEM <- data.frame(x = s[I], y = FEM)
data_BS <-  data.frame(x = s[I], y = BS)
data_POF <- data.frame(x = s[I], y = POF)

theme_set(theme_bw())
p <- ggplot(data_FEM, aes(x = x, y = y)) + 
  geom_line(aes(color="FEM")) + 
  geom_line(data=data_BS, aes(color="BS")) +
  geom_line(data=data_POF, aes(color="Payoff")) + 
  scale_color_manual(name = "Method",
                     values=c(FEM="red", BS="blue", Payoff="black")) + 
  ggtitle(TeX(paste("$e_{L^2}(N) \\approx ",
                    round(error, 4), sep = ""))) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Underlying price") + 
  ylab("Option/payoff value")
p

