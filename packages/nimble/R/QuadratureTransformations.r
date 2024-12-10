## Linear transformations to / from the standard and observed scale.

## This is the key transformation methods needed for transforming the quadrature grid the way that I want to.
quadGridTransformMethods = nimbleFunction(
  setup = function(){},
  run = function(){},
  methods = list(
    z_to_theta = function(z = double(1), 
                          mode = double(1), 
                          A = double(2), 
                          lambda = double(1),
                          transformation = character(0, default = "spectral") ) {
      returnType(double(1))
      if( transformation == "spectral" ) {
        d <- dim(z)[1]
        theta <- numeric(value = 0, length = d)
        for( i in 1:d ){
          theta[i] <- mode[i] + sum(A[,i] * z) / sqrt(lambda[i])
        }
      }else{
        ## Cholesky transformation.
        theta <- mode + backsolve(A, z)
      }
      return(theta)
    },
    theta_to_z = function(theta = double(1), 
                          mode = double(1), 
                          A = double(2), 
                          lambda = double(1),
                          transformation = character(0, default = "spectral") ){
      returnType(double(1))
      if( transformation == "spectral" ) {
        d <- dim(theta)[1]
        z <- numeric(value = 0, length = d)
        theta_mean <- theta - mode
        sqrtlambda <- sqrt(lambda)
        for( i in 1:d ){
          z[i] <- sum(A[i,]*sqrtlambda*theta_mean)
        }
      }else{
        ## Cholesky transformation.
        z <- (A %*% (theta - mode))[,1]
      }
      return(z)
    },
    z_to_theta_vec = function(z = double(2), 
                          mode = double(1), 
                          A = double(2), 
                          lambda = double(1),
                          transformation = character(0, default = "spectral") ) {
      returnType(double(2))
      d <- dim(z)[2]
      nQ <- dim(z)[1]
      theta <- matrix(0, nrow = nQ, ncol = d)

      if( transformation == "spectral" ) {
        sqrtlambda <- sqrt(lambda)
        for( i in 1:nQ ) {
          for( j in 1:d ) {
            theta[i,j] <- mode[j] + sum(A[,j] * z[i,]) / sqrtlambda[j]
          }
        }
      }else{
        ## Cholesky transformation.
        for( i in 1:nQ ) {
          theta[i,] <- mode + backsolve(A, z[i,])
        }
      }
      return(theta)
    },
    theta_to_z_vec = function(theta = double(2), 
                          mode = double(1), 
                          A = double(2), 
                          lambda = double(1),
                          transformation = character(0, default = "spectral") ){
      returnType(double(2))
      
      d <- dim(theta)[2]
      nQ <- dim(theta)[1]
      z <- matrix(0, nrow = nQ, ncol = d)

      if( transformation == "spectral" ) {
        sqrtlambda <- sqrt(lambda)
        for( i in 1:nQ ){
          theta_mean <- theta[i,] - mode
          for( j in 1:d ){
            z[i,j] <- sum(A[j,]*sqrtlambda*theta_mean)
          }
        }
      }else{
        ## Cholesky transformation.
        for( i in 1:nQ ) {
          z[i,] <- (A %*% (theta[i,] - mode))[,1]
        }
      }
      return(z)
    },    
    ## Skew grid points:
    skew_z = function( z = double(2), skewSD = double(2) ){
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