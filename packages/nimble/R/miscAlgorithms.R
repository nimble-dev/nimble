
## Not exported function to get approximate number of iterations.
## Based on Sherlock 2021 Appendix A and code from:
## https://github.com/ChrisGSherlock/expQ/blob/master/rexpQ.cpp
## get_m(). Switching back to qpois for now. Will turn use when 
## buildDerivs = TRUE in next iteration.
# getNPrec = nimbleFunction(
  # run = function(rho = double(), prec = double(0, default= 1e-8)) {
    
    # logprec <- log(prec)
    # pi <- 3.14159265
    # logroot2pi <- 0.5*log(2*pi)
    
    # mhi <- rho + (-logprec)/3.0*(1.0+sqrt(1.0+18.0*rho/(-logprec)))
    # A <- 2*rho*(1.0-mhi/rho+mhi/rho*log(mhi/rho))
    # mlo <- rho+sqrt(2*rho)*sqrt(-logroot2pi-logprec-1.5*log(A)+log(A-1))
    # returnType(integer())
    # return(ceiling((mlo+mhi)/2))
  # } #, buildDerivs = TRUE
# )

## Not exported function to bound the spectrum of a square matrix:
## Based on eigenDisc from RTMB, which allows for the 
## generalization of the Sherlock 2021 algorithm for more than just 
## generator matrices.
spectrumBound <- nimbleFunction(
  run = function(A = double(2)){
    q <- dim(A)
    Min <- Inf
    Max <- -Inf
    for( i in 1:q[1] ){
      center <- A[i,i]
      radius <- sum(abs(A[i,])) - abs(center)
      minC <- center - radius
      maxC <- center + radius
      if(minC < Min) Min <- minC
      if(maxC > Max) Max <- maxC
    }
    ans <- nimNumeric(value = 0, length = 2)
    returnType(double(1))
    ans[1] <- (Max + Min)/2
    ans[2] <- (Max - Min)/2
    return(ans)
  }#, buildDerivs = TRUE
)

#' Matrix Exponential times a vector
#'
#' Compute the combined term \code{expm(A) \%*\% v} to avoid a full matrix exponentiation.
#' 
#' @name expAv
#' 
#' @param A Square matrix.
#' @param v vector to multiply by the matrix exponential \code{exp(A) \%*\% v}.
#' @param tol level of accuracy required (default = 1e-8).
#' @param rescaleFreq How frequently should the terms be scaled to avoid underflow/overflow (default = 10).
#' @param Nmax Maximum number of iterations to compute (default = 10000).
#' @param sparse (logical) specify if the matrix may be sparse and to do sparse computation (default = TRUE).
#' @author Paul van Dam-Bates
#' @details
#' For large matrix exponentials it is much more efficient to compute \code{exp(A) \%*\% v}, than to actually compute the entire matrix exponential.
#'
#' This function follows the function \code{expAv} from the R package \pkg{RTMB} (Kristensen, 2025), and theory outlined in Sherlock (2021). It is developed for working with continuous times
#' Markov chains. If using the matrix exponential to create a transition probability matrix in a HMM context just once, 
#' this function may be slower than the one time call to compute the full matrix exponentiation. If a full matrix exponentiation is required, refer to 
#' \code{expm} to compute. Choosing \code{sparse = TRUE} will check which values of \code{A} are non-zero and do sparse linear algebra.
#' Note that for computation efficiency matrix uniformization is done by \code{A* = A + rho I}, where \code{rho = max(abs(diag(A)))}; see Algorithm 2' in Sherlock (2021).
#' When the row sums of the matrix are not zero, then uniformization is not done, and the number of iterations to reach tolerance are approximated based on the
#' a bound of the spectrum, similar to \pkg{RTMB} (Kristensen, 2025).
#'
#' @return the result as a vector.
#'
#' @references 
#' Sherlock, C. (2021). Direct statistical inference for finite Markov jump processes via the matrix exponential. Computational Statistics, 36(4), 2863-2887.
#' 
#' Kristensen K (2025). _RTMB: 'R' Bindings for 'TMB'_. R package version 1.7, commit 6bd7a16403ccb4d3fc13ff7526827540bf27b352, 
#' <https://github.com/kaskr/RTMB>.
#'
#' @examples
#' A <- rbind(c(-1, 0.25, 0.75), c(0, -2, 2), c(0.25, 0.25, -0.5))
#' v <- c(0.35, 0.25, 0.1)
#' expAv(A, v)
NULL

#' @rdname expAv
#' @export
expAv <- nimbleFunction(
  run = function(A = double(2), v = double(1),
                 tol = double(0, default = 1e-8), 
                 rescaleFreq = double(0, default = 10),
                 Nmax = integer(0, default = 10000),
                 sparse = logical(0, default = TRUE)){
    q <- dim(A)
    if(q[2] != q[1])      stop("A must be a square matrix.")
    if(length(v) != q[1]) stop("Length of v must be equal to the number of columns of A.")

    ## Only center, perform uniformization, if no diag > 0.
    CR <- spectrumBound(A)
    diag(A) <- diag(A) - CR[1]
    ans <- v
    term <- v
    log_scale <- 0
    m <- 1L
    Niter <- 1L
    # R <- ADbreak(CR[2])
    # Niter <- getNPrec(CR[2], tol)
    Niter <- qpois(tol, CR[2], lower.tail = FALSE)
    if(Niter > Nmax) Niter <- Nmax
    
    ## Check and build sparse matrix if sparsity exists.
    ## Compressed sparse col format
    if(sparse){
      m <- prod(q)
      colptr <- nimNumeric(value = 0, length = q[2]+1)
      rowindex <- nimNumeric(value = 0, length = m)
      values <- nimNumeric(value = 0, length = m)
      i <- 1L; j <- 1L; k <- 1L
      colptr[1] <- 1
      for( j in 1:q[2] ){
        for( i in 1:q[1] ){
          if(A[i,j] != 0){
            rowindex[k] <- i
            values[k] <- A[i,j]
            k <- k+1
          }
        }
        colptr[j+1] <- k
      }
      m <- k-1
      rowindex <- rowindex[1:m]
      values <- values[1:m]
    }
    check <- rescaleFreq
    for(n in 1:Niter) {
      if(sparse){
        tmp <- nimNumeric(0, length = q[1]) ## Row Vector
        ## Sparse column matrix vector multiplication: term <- A %*% term/n
        ## for each column, add a[rowindices,j] * term[j]
        for( j in 1:q[2] ){
          if(colptr[j+1] > colptr[j])
            tmp[rowindex[colptr[j]:(colptr[j+1]-1)]] <- tmp[rowindex[colptr[j]:(colptr[j+1]-1)]] + term[j]*values[colptr[j]:(colptr[j+1]-1)]
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
    }
    if(Niter == Nmax)
      cat("  Warning: Nmax in `expAv` is less than N required for the specified tolerance. Consider increasing Nmax.\n")

    ans <- exp(CR[1] + log_scale) * ans
    returnType(double(1))
    return(ans)
  },# buildDerivs = TRUE
)

#' Matrix Exponential
#'
#' Compute the the matrix exponential \code{expm(A)} by scaling and squaring.
#' 
#' @name expm
#' 
#' @param A Square matrix.
#' @param tol level of accuracy required (default = 1e-8).
#' @author Paul van Dam-Bates
#'
#' @details
#' This function follows the scaling and squaring algorithm from Sherlock (2021), except that we compute the full
#' matrix exponential. It differs from the standard Taylor scaling and squaring algorithm reviewed by Ruiz et al (2016) and found in common texts,
#' by doing uniformization if the matrix is a generator matrix from a continuous time Markov chain. If using the matrix exponential to create a 
#' transition probability matrix in a HMM context just once, this function may be efficient. If dimension is large, we recommend avoiding the 
#' matrix exponential and using \code{expAv} instead. Note that for computation efficiency, when the columns are non-positive, matrix uniformization 
#' is done by \code{A* = A + rho I}, where \code{rho = max(abs(diag(A)))}. When the row sums of the matrix are not zero, then uniformization is not done, 
#' and the number of iterations to reach tolerance are approximated based on the a bound of the spectrum, similar to \pkg{RTMB} (Kristensen, 2025).
#'
#' @return a matrix that is ans = exp(A).
#'
#' @references 
#' Sherlock, C. (2021). Direct statistical inference for finite Markov jump processes via the matrix exponential. Computational Statistics, 36(4), 2863-2887.
#' 
#' Ruiz, P., Sastre, J., Ibáñez, J., & Defez, E. (2016). High performance computing of the matrix exponential. 
#' Journal of computational and applied mathematics, 291, 370-379.
#'
#' Kristensen K (2025). _RTMB: 'R' Bindings for 'TMB'_. R package version 1.7, commit 6bd7a16403ccb4d3fc13ff7526827540bf27b352, 
#' <https://github.com/kaskr/RTMB>.
#'
#' @examples
#' A <- rbind(c(-1, 0.25, 0.75), c(0, -2, 2), c(0.25, 0.25, -0.5))
#' Lambda <- diag(c(0.25, 0.1, 0))
#' expm((A-Lambda)*2.5)
NULL


#' @rdname expm
#' @export
expm <- nimbleFunction(
  run = function(A = double(2), tol = double(0, default = 1e-8)){
    returnType(double(2))
    
    if(dim(A)[1] != dim(A)[2]) stop("`A' must be a square matrix")
    # if(dim(tol)[1] > 1)  stop("`tol' must be a scalar.")

    ## Only center, perform uniformization, if no diag > 0.
    CR <- spectrumBound(A)
    Niter <- 1L
    C <- CR[1]
    R <- CR[2]
    # Niter <- getNPrec(R, tol)

    n <- dim(A)[1]
    log2 <- log(2)
    shat <- ceiling((log(R)+log(log2))/log2)
    smin <- max(c(0, shat))
    smax <- smin+6;
    rhosmall <- R/exp(smin*log(2))
    opt <- Inf
    for( k in smin:smax ){
      rhosmall <- rhosmall/2
      # test <- getNPrec(rhosmall, tol) + k
      test <- qpois(tol, rhosmall, lower.tail = FALSE) + k
      if(test < opt){
        opt <- test
        s <- k
      }
    }
    twos <- 2^s
    Msmall <- (A - C*diag(n))/twos
    # m <- getNPrec(R/twos, tol)
    m <- qpois(tol, R/twos, lower.tail = FALSE)    
    m <- max(2, m)
    
    expMsmall <- diag(n) + Msmall
    expMpro <- Msmall 
    for( i in 2:m ) {
      expMpro <- expMpro %*% Msmall / i
      expMsmall <- expMsmall + expMpro
    }

    expMsmall <- expMsmall * exp(C/twos)

    for( i in 1:s ) expMsmall <- expMsmall %*% expMsmall
    return(expMsmall)
  } #, buildDerivs = TRUE
)
