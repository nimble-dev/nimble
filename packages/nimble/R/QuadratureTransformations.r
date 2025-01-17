## Linear transformations to / from the standard and observed scale.

## This is the key transformation methods needed for transforming the quadrature grid the way that I want to.
quadGridTransformMethods = nimbleFunction(
  setup = function(){},
  run = function(){},
  methods = list(
    ## Spectral transform for a vector of z to theta.
    z_to_theta_spectral = function( z = double(1), 
                          mode = double(1), 
                          eigenVectors = double(2), 
                          eigenValues = double(1) ) {
      returnType(double(1))
      d <- dim(z)[1]
      theta <- numeric(value = 0, length = d)
      for( i in 1:d ){
        theta[i] <- mode[i] + sum(eigenVectors[,i] * z) / sqrt(eigenValues[i])
      }
      return(theta)
    },
    ## Cholesky transform for a vector of z to theta.
    z_to_theta_cholesky = function( z = double(1), 
                          mode = double(1), 
                          A = double(2) ) {
      ## Cholesky transformation.
      theta <- mode + backsolve(A, z)
      returnType(double(1))
      return(theta)
    },
    ## Spectral transform for a vector of theta to z.
    theta_to_z_spectral = function(theta = double(1), 
                          mode = double(1), 
                          eigenVectors = double(2), 
                          eigenValues = double(1) ){
      d <- dim(theta)[1]
      z <- numeric(value = 0, length = d)
      theta_mean <- theta - mode
      sqrtlambda <- sqrt(eigenValues)
      for( i in 1:d ){
        z[i] <- sum(eigenVectors[i,]*sqrtlambda*theta_mean)
      }
      returnType(double(1))
      return(z)
    },
    ## Cholesky transform for a vector of theta to z.    
    theta_to_z_cholesky = function(theta = double(1), mode = double(1), A = double(2)){
      z <- (A %*% (theta - mode))[,1]
      returnType(double(1))
      return(z)
    },
    ## Spectral transform for a matrix of z to theta.   
    z_to_theta_spectral_vec = function(z = double(2), 
                          mode = double(1), 
                          eigenVectors = double(2), 
                          eigenValues = double(1) ) {
      d <- dim(z)[2]
      nQ <- dim(z)[1]
      theta <- matrix(0, nrow = nQ, ncol = d)

      sqrtlambda <- sqrt(eigenValues)
      for( i in 1:nQ ) {
        for( j in 1:d ) {
          theta[i,j] <- mode[j] + sum(eigenVectors[,j] * z[i,]) / sqrtlambda[j]
        }
      }
      returnType(double(2))
      return(theta)
    },
    ## Cholesky transform for a matrix of z to theta.   
    z_to_theta_cholesky_vec = function(z = double(2), 
                          mode = double(1), 
                          A = double(2) ) {
      d <- dim(z)[2]
      nQ <- dim(z)[1]
      theta <- matrix(0, nrow = nQ, ncol = d)

      ## Cholesky transformation.
      for( i in 1:nQ ) {
        theta[i,] <- mode + backsolve(A, z[i,])
      }
      returnType(double(2))      
      return(theta)
    },
    ## Spectral transform for a matrix of theta to z.   
    theta_to_z_spectral_vec = function(theta = double(2), 
                          mode = double(1), 
                          eigenVectors = double(2), 
                          eigenValues = double(1) ) {
      d <- dim(theta)[2]
      nQ <- dim(theta)[1]
      z <- matrix(0, nrow = nQ, ncol = d)

      sqrtlambda <- sqrt(eigenValues)
      for( i in 1:nQ ){
        theta_mean <- theta[i,] - mode
        for( j in 1:d ){
          z[i,j] <- sum(eigenVectors[j,]*sqrtlambda*theta_mean)
        }
      }
      returnType(double(2))
      return(z)
    },
    ## Choleksy transform for a matrix of theta to z.      
    theta_to_z_cholesky_vec = function(theta = double(2), 
                          mode = double(1), 
                          A = double(2) ){
      d <- dim(theta)[2]
      nQ <- dim(theta)[1]
      z <- matrix(0, nrow = nQ, ncol = d)

      for( i in 1:nQ ) {
        z[i,] <- (A %*% (theta[i,] - mode))[,1]
      }
      returnType(double(2))
      return(z)
    },        
    ## Skew grid points:
    skew_z = function( z = double(1), skewSD = double(2) ){
      d <- dim(z)[1]
      znew <- z
      for( j in 1:d ) {
        if(znew[i,j] <= 0) 
          znew[i,j] <-  znew[i, j]*skewSD[j, 1]	# Negative Skew
        else
          znew[i, j] <- znew[i, j]*skewSD[j,2] # Positive Skew
      }
      returnType(double(2))
      return(znew)
    },
    ## Skew grid points:
    skew_z_vec = function( z = double(2), skewSD = double(2) ){
      d <- dim(z)[2]
      nQ <- dim(z)[1]
      znew <- z
      for( i in 1:nQ) {
        for( j in 1:d ) {
          if(znew[i,j] <= 0) 
            znew[i,j] <-  znew[i, j]*skewSD[j, 1]	# Negative Skew
          else
            znew[i, j] <- znew[i, j]*skewSD[j,2] # Positive Skew
        }
      }
      returnType(double(2))
      return(znew)
    }
  )
)
