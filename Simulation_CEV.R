## Option pricing via Simulations:
##   => CEV Model

m <- 50000         # Number of paths
dt <- 0.00001       # Timestep
r <- 0.06          # Interest rate
Mat <- 1           # Maturity
sig <- 0.6         # Volatility
K <- 1             # Strike price
theta <- 0.5       # Theta scheme parameter
rho <- 0.6         # CEV parameter
S0 <- 2     # Initial value of stock

n <- Mat / dt + 1
paths <- matrix(0, ncol = m, nrow = n)
paths[1,] <- S0 # Initial value

for (i in 2:n){
  print(i)
  Z <- rnorm(m, 0, sqrt(dt))
  paths[i,] <- paths[i-1,] + 
    r * paths[i-1,] * dt + 
    sig * paths[i-1,]^rho * Z
}

OptionValue <- exp(-r * Mat) * mean(pmax(K - paths[n,!is.nan(paths[n,])], 0))
OptionValue




