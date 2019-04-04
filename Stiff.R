Stiff <- function(vertices, alpha, beta, gamma){
  ## Compute the stiffness matrix
  
  n <- length(vertices)
  SDiag <- BDiag <- MDiag <- matrix(0, ncol = 3, nrow = n)
  Sloc <- Bloc <- Mloc <- matrix(0, 2, 2)
  h <- vertices[2] - vertices[1]
  
  basis <- function(x, xi){
      abs(x - xi) / h * (abs(x - xi) <= h)
  }
  d_basis <- function(x, xi){
      -sign(x - xi) / h * (abs(x - xi) <= h)
  }
  
  order <- 2
  # rule <- gaussquad::legendre.quadrature.rules(order)
  # integ <- function(f, l, u){
  #   gaussquad::legendre.quadrature(f, rule[[2]], l, u)
  # }
  
  vidx <- c(1, 2)
  for (i in 1:(n - 1)){
    
    rule <- pracma::gaussLegendre(order, vertices[max(1,i)], vertices[i+1])
    integ <- function(f, l, u){
      sum(f(rule$x) * rule$w)
    }
    
    ## Fill Sloc:
    Sloc[1,1] <- integ(function(x){
      alpha(x) * d_basis(x, vertices[i]) * 
        d_basis(x, vertices[i])
    }, vertices[max(1,i)], vertices[i+1])
    
    Sloc[1,2] <- integ(function(x){
      alpha(x) * d_basis(x, vertices[i+1]) * 
        d_basis(x, vertices[i])
    }, vertices[max(1,i)], vertices[i+1])
    
    Sloc[2,1] <- integ(function(x){
      alpha(x) * d_basis(x, vertices[i]) * 
        d_basis(x, vertices[i+1])
    }, vertices[max(1,i)], vertices[i+1])
    
    Sloc[2,2] <- integ(function(x){
      alpha(x) * d_basis(x, vertices[i+1]) * 
        d_basis(x, vertices[i+1])
    }, vertices[max(1,i)], vertices[i+1])
    
    ## Fill Bloc
    Bloc[1,1] <- integ(function(x){
      beta(x) * d_basis(x, vertices[i]) * 
        basis(x, vertices[i])
    }, vertices[max(1,i)], vertices[i+1])
    
    Bloc[1,2] <- integ(function(x){
      beta(x) * d_basis(x, vertices[i+1]) * 
        basis(x, vertices[i])
    }, vertices[max(1,i)], vertices[i+1])
    
    Bloc[2,1] <- integ(function(x){
      beta(x) * d_basis(x, vertices[i]) * 
        basis(x, vertices[i+1])
    }, vertices[max(1,i)], vertices[i+1])
    
    Bloc[2,2] <- integ(function(x){
      beta(x) * d_basis(x, vertices[i+1]) * 
        basis(x, vertices[i+1])
    }, vertices[max(1,i)], vertices[i+1])
    
    ## Fill Mloc
    Mloc[1,1] <- integ(function(x){
      gamma(x) * basis(x, vertices[i]) * 
        basis(x, vertices[i])
    }, vertices[max(1,i)], vertices[i+1])
    
    Mloc[1,2] <- integ(function(x){
      gamma(x) * basis(x, vertices[i+1]) * 
        basis(x, vertices[i])
    }, vertices[max(1,i)], vertices[i+1])
    
    Mloc[2,1] <- integ(function(x){
      gamma(x) * basis(x, vertices[i]) * 
        basis(x, vertices[i+1])
    }, vertices[max(1,i)], vertices[i+1])
    
    Mloc[2,2] <- integ(function(x){
      gamma(x) * basis(x, vertices[i+1]) * 
        basis(x, vertices[i+1])
    }, vertices[max(1,i)], vertices[i+1])
    
    ## Fill Diag matrices
    SDiag[vidx[1], 2] <- SDiag[vidx[1], 2] + Sloc[1,1]
    SDiag[vidx[2], 2] <- SDiag[vidx[2], 2] + Sloc[2,2]
    SDiag[vidx[2], 3] <- SDiag[vidx[2], 3] + Sloc[1,2]
    SDiag[vidx[1], 1] <- SDiag[vidx[1], 1] + Sloc[2,1]
    
    BDiag[vidx[1], 2] <- BDiag[vidx[1], 2] + Bloc[1,1]
    BDiag[vidx[2], 2] <- BDiag[vidx[2], 2] + Bloc[2,2]
    BDiag[vidx[2], 3] <- BDiag[vidx[2], 3] + Bloc[1,2]
    BDiag[vidx[1], 1] <- BDiag[vidx[1], 1] + Bloc[2,1]
    
    MDiag[vidx[1], 2] <- MDiag[vidx[1], 2] + Mloc[1,1]
    MDiag[vidx[2], 2] <- MDiag[vidx[2], 2] + Mloc[2,2]
    MDiag[vidx[2], 3] <- MDiag[vidx[2], 3] + Mloc[1,2]
    MDiag[vidx[1], 1] <- MDiag[vidx[1], 1] + Mloc[2,1]
    
    # Update vidx
    vidx <- vidx + 1
  }
  
  # Construct tridiagonal matrices
  S <- spTridiags(SDiag[1:(n-1),1],
                  SDiag[1:n,2],
                  SDiag[2:n,3])
  B <- spTridiags(BDiag[1:(n-1),1],
                  BDiag[1:n,2],
                  BDiag[2:n,3])
  M <- spTridiags(MDiag[1:(n-1),1],
                  MDiag[1:n,2],
                  MDiag[2:n,3])
  
  # Return stiffness matrix
  S + B + M
}