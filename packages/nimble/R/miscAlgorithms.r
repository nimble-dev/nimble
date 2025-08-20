#' Matrix Exponential times a vector
#'
#'   Compute the combined term expm(A) %*% v
#'   to avoid a full matrix exponentiation.
#' 
#' @name expmAv
#' 
#' @param A Infinitesimal Generator Matrix.
#' @param v vector to multiply by the matrix exponential exp(A) %*% v.
#' @param tol level of accuracy required (default = 1e-8).
#' @param rescaleFreq How frequently should the terms be scaled to avoid underflow/overflow (default = 10).
#' @param Nmax maximum number of iterations to compute (default = 10000).
#' @param sparse (logical) specify if the matrix may be sparse and to do sparse computation (default = FALSE).
#' @author Paul van Dam-Bates
#' @details For large matrix exponentials it is much more efficient to compute exp(A) %*% v, than to actually compute the entire matrix exponential. A, the generator matrix must have
#' negative value diagonals.
#'
#' This function follows the function `expAv` from the R package RTMB (Kristensen, 2025), and theory outlined in Sherlock (2021). It is developed for working with continuous times
#' Markov chains. If using the matrix exponential to create a transition probability matrix in a HMM context just once, 
#' this function may be slower than the one time call to compute the full matrix exponentiation. If a full matrix exponentiation is required, refer to `nimbleRcall` and 
#' functions such as `pracma::expm` or `expm::expm` to compute. Choosing sparse = TRUE will then check which matrix A values are non-zero and do sparse linear 
#' algebra on the non-zero terms only. Note that for computation efficiency matrix uniformization is always done by A* = A + rho I, where rho = max(abs(diag(A))), see Algorithm 2#' in Sherlock (2021).
#'
#' @return \code{expmAv} gives a vector that is ans = exp(A) %*% v.
#' @references 
#' Sherlock, C. (2021). Direct statistical inference for finite Markov jump processes via the matrix exponential. Computational Statistics, 36(4), 2863-2887.
#' 
#' Kristensen K (2025). _RTMB: 'R' Bindings for 'TMB'_. R package version 1.7, commit 6bd7a16403ccb4d3fc13ff7526827540bf27b352, 
#' <https://github.com/kaskr/RTMB>.
#' @examples
#' A <- rbind(c(-1, 0.25, 0.75), c(0, -2, 2), c(0.25, 0.25, -0.5))
#' v <- c(0.35, 0.25, 0.1)
#' expmAv(A, v)
NULL

#' @rdname expmAv
#' @export
expmAv <- nimbleFunction(
  run = function(A = double(2), v = double(1),
                 tol = double(0, default = 1e-8), 
                 rescaleFreq = double(0, default = 10),
                 Nmax = integer(0, default = 10000),
                 sparse = logical(0, default = FALSE)){
    if(any(diag(A) > 0))
      stop("The expAv algorithm is designed for generator matrices only, which must have negative or zero diagonals.")
    C <- max(abs(diag(A)))
    diag(A) <- diag(A) + C
    ans <- v
    term <- v
    log_scale <- 0
    m <- 1L
 
    ## Check and build sparse matrix if sparsity exists.
    if(sparse){
      q <- dim(A)
      m <- prod(q)
      indexi <- nimNumeric(value = 0, length = m)
      indexj <- nimNumeric(value = 0, length = m)
      values <- nimNumeric(value = 0, length = m)
      i <- 1L; j <- 1L; k <- 1L
      for( i in 1:q[1] ){
        for( j in 1:q[2] ){
          if(A[i,j] != 0){
            indexi[k] <- i
            indexj[k] <- j
            values[k] <- A[i,j]
            k <- k+1
          }
        }
      }
      m <- k-1
      indexi <- indexi[1:m]
      indexj <- indexj[1:m]
      values <- values[1:m]
    }
    maxterm <- which(term == max(term))[1]
    n <- 1L
    check <- rescaleFreq
    while(n < Nmax & term[maxterm] > tol) {
      if(sparse){
        tmp <- nimNumeric(0, length = q[2])
        for( k in 1:m ){
          tmp[indexi[k]] <- tmp[indexi[k]] + term[indexj[k]]*values[k]
        }
        term <- tmp / n
      }else{
        term <- (A %*% asCol(term/n))[,1]
      }
      ans <- ans + term
      if( check == n ) {
        s <- sum(abs(ans))
        term <- term / s
        ans <- ans / s
        log_scale <- log_scale + log(s)
        check <- check + rescaleFreq
      }
      n <- n+1
    }
    if(n == Nmax)
      cat("  Warning: Nmax in `expAv` is less than N required for the specified tolerance. Consider increasing Nmax.\n")

    ans <- exp(-C + log_scale) * ans
    returnType(double(1))
    return(ans)
  },
  buildDerivs = TRUE
)

#' Matrix Exponential
#'
#'   Compute the the matrix exponential expm(A)
#'   by scaling and squaring.
#' 
#' @name expmA
#' 
#' @param A Infinitesimal Generator Matrix.
#' @param tol level of accuracy required (default = 1e-8).
#' @author Paul van Dam-Bates
#' @details A, the generator matrix must have negative value diagonals.
#'
#' This function follows the scaling and squaring algorithm from Sherlock (2021), except that we compute the full
#' matrix exponential. It differs from the standard Taylor scaling and squaring algorithm reviewed by Ruiz et al (2016) and found in common texts. 
#' If using the matrix exponential to create a transition probability matrix in a HMM context just once, 
#' this function is good. If dimension is large, we recommend avoiding the matrix exponential and using `expmAv` instead. 
#' Note that for computation efficiency matrix uniformization is always done by A* = A + rho I, where rho = max(abs(diag(A))).
#'
#' @return \code{expmAv} gives a vector that is ans = exp(A) %*% v.
#' @references 
#' Sherlock, C. (2021). Direct statistical inference for finite Markov jump processes via the matrix exponential. Computational Statistics, 36(4), 2863-2887.
#' 
#' Ruiz, P., Sastre, J., Ibáñez, J., & Defez, E. (2016). High performance computing of the matrix exponential. 
#' Journal of computational and applied mathematics, 291, 370-379.
#'
#' @examples
#' A <- rbind(c(-1, 0.25, 0.75), c(0, -2, 2), c(0.25, 0.25, -0.5))
#' Lambda <- diag(c(0.25, 0.1, 0))
#' expmA((A-Lambda)*2.5)
NULL


#' @rdname expmA
#' @export
expmA <- nimbleFunction(
  run = function(A = double(2), tol = double(0, default = 1e-8)){
    returnType(double(2))
    n <- dim(A)[1]
    rho <- max(abs(diag(A)))
    log2 <- log(2)
    shat <- ceiling((log(rho)+log(log2))/log2)
    smin <- max(c(0, shat))
    smax <- smin+6;
    rhosmall <- rho/exp(smin*log(2))
    opt <- Inf
    for( k in smin:smax ){
      rhosmall <- rhosmall/2
      test <- qpois(tol, rhosmall, lower.tail = FALSE) + k 
      if(test < opt){
        opt <- test
        s <- k
      }
    }
    twos <- 2^s
    if(uniformization){
      Msmall <- (A + rho*diag(n))/twos
    }else{
      Msmall <- A/twos
    }
    m <- qpois(tol, rho/twos, lower.tail = FALSE)

    expMsmall <- diag(n) + Msmall
    expMpro <- Msmall 
    for( i in 2:m ) {
      expMpro <- expMpro %*% Msmall / i
      expMsmall <- expMsmall + expMpro
    }

    if(uniformization) 
      expMsmall <- expMsmall * exp(-rho/twos)

    expM <- expMsmall
    for( i in 1:s ) expM <- expM %*% expM
    return(expM)
  }
)
