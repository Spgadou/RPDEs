spTridiags <- function(lower, main, upper){
  ## lower: lower diagonal vector
  ## main: middle diagonal vector
  ## upper: upper diagonal vector
  ##
  ## return: sparse, tridiagonal matrix
  
  diags <- list(lower, main, upper)
  bandSparse(length(main), k = c(-1,0,1), diag = diags)
}
