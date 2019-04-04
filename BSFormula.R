BSFormula <- Vectorize(function(s, r, Mat, K, sig, type = "call"){
  d1 <- (log(s / K) + (r + 0.5 * sig^2) * Mat) / (sig * sqrt(Mat))
  d2 <- d1 - sig * sqrt(Mat)
  call <- s * pnorm(d1) - K * exp(-r * Mat) * pnorm(d2)
  
  if (type == "call")
    call
  else
    call - (s - K * exp(-r * Mat))
}, 's')
