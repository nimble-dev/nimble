
#----------------------------------------------#
# Chinesse restaurant process distributions:
#----------------------------------------------#

# random number generator: rCRP
#   input:
#      n:  size of the sample
#      conc: concentration parameter of the DP (related to sampling a new value)
#   output:
#      x: vector of size n
#      size: 
#
#   a new value is sampled with probability conc/(i-1+conc), an existing
#   value is sampled otherwise.

# density evaluation: dCRP
#   input: 
#      x: a vector
#      conc: concentration parameter of the DP
#      log: density returned in log scale or not
#   output:
#      x: a value

#' The Chinesse Restaurant Process Distribution
#'
#'   Density and random generation for the Chinesse
#'   Restaurant Process distribution with concentration 
#'   parameter \code{conc}.
#' 
#' @name ChinesseRestaurantProcess 
#' 
#' @param x vector of values.
#' @param n number of observations (only n = 1 is handled currently).
#' @param conc scalar with concentration parameter.
#' @param size integer of length of vector x.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Claudia Wehrhahn
#' @details The chineese restaurant process distribution is a distribution
#' on the space of partitions of the positive integers. 
#' The chineese restaurant process
#' distribution with concentration parameter \eqn{=\alpha}{= conc} has 
#' probability function 
#' \deqn{f(x_i \mid x_1, \ldots, x_{i-1})=\frac{1}{i-1+\alpha}\sum_{j=1}^{i-1}\delta_{x_j}+
#' \frac{\alpha}{i-1+\alpha}\delta_{x^{new}},}
#' where \eqn{x^{new}} is a new integer not in \eqn{x_1, \ldots, x_{i-1}}.
#' 
#' If \code{conc} is not specified, it assumes the default value of 1. 
#' @return \code{dCRP} gives the density, and \code{rCRP} gives random generation.
#' @references Blackwell, D., & MacQueen, J. B. (1973) Ferguson distributions via 
#' Pólya urn schemes \emph{The annals of statistics}, 353-355.
#' 
#' Aldous, D. J. (1985) Exchangeability and related topics \emph{In École d'Été de 
#' Probabilités de Saint-Flour XIII} 1983 (pp. 1-198) Springer, Berlin, Heidelberg.
#' 
#' Pitman, Jim. Some developments of the Blackwell-MacQueen urn scheme 
#' \emph{Statistics, probability and game theory}, 245--267, Institute of
#'  Mathematical Statistics, Hayward, CA, 1996.
#' 
#' @examples
#' x <- rCRP(n=1, conc = 1, size=10)
#' dCRP(x, conc = 1, size=10)
NULL

#' @rdname ChinesseRestaurantProcess
#' @export
rCRP=nimbleFunction(
  run=function(n = integer(0), 
               conc = double(0, default=1),
               size = integer(0))
  {
    returnType(double(1))
    
    if(n>1){
      print(paste('Warning message: In rCRP(', n, ',', conc, ',', size, ') : rCRP only handles n = 1 at the moment'))
    }
    
    if( conc <= 0 ){
      concCond1 <- 1
      print('conc parameter has to be larger than zero')
    }else{
      concCond1 <- 0
    }
    
    if(concCond1 == 1){
      return(c(NaN))
    }else{
      x <- numeric(size)# returned vector:
      x[1] <- 1
      if(size>1){
        for(i in 2:size){
          if( runif(1)<=conc/(conc+i-1) ){# a new value
            x[i] <- max(x[1:(i-1)])+1 
          }else{# an old value
            index <- rcat(n=1, rep(1, i-1)) 
            x[i] <- x[index]
          }
        }
      }
      return(x)
    }
  }
)


#' @rdname ChinesseRestaurantProcess
#' @export
dCRP=nimbleFunction(
  run=function(x = double(1), 
               conc = double(0, default=1),
               size = integer(0),
               log = integer(0, default=0))
  {
    returnType(double(0))
    n <- length(x) 
    tmpden <- numeric(n) 
    
    if( conc <= 0 ){
      concCond1 <- 1
      print('conc parameter has to be larger than zero')
    }else{
      concCond1 <- 0
    }
    
    if(concCond1 == 1){
      return(NaN)
    }else{
      tmpden[1] <- 1
      if(n>1){
        for(i in 2:n){
          #counts=0 # replaces sum(x[i]==x[1:(i-1)]) (works in nimble?)
          #for(j in 1:(i-1)){
          #  if(x[i]==x[j]){
          #    counts=counts+1
          #  }
          #}
          counts <- sum(x[i]==x[1:(i-1)])
          if( counts>0 ){
            tmpden[i] <- 1/(i-1+conc)
          }else{
            tmpden[i] <- conc/(i-1+conc)
          }
        }
      }
      logProb <- sum(log(tmpden))
      if(log) return(logProb)
      else return(exp(logProb)) 
    }
  }
)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



#----------------------------------------------#
# stick_breaking nimbleFunction
#----------------------------------------------#

#' The Stick Breaking function
#'
#'   Based on the stick breaking construction, weights are computed 
#'   with the argument vector.
#' 
#' @name StickBreakingFunction 
#' 
#' @param z vector argument.
#' @param log logical; if TRUE, weights are returned on the log scale.
#' @author Claudia Wehrhahn
#' @details values in \code{z} have to be between \deqn{[0,1)]}. If one of the components 
#' is equal to 1, then the returned vector has length smaller than \code{z}. If one of the
#' components is smaller than 0 or greater than 1, NaNs are returned.
#' @references Sethuraman, J. (1994) A constructive definition of Dirichlet 
#' priors. \emph{Statistica sinica}, 639-650.
#' @examples
#' z <- rbeta(5, 1,1)
#' stick_breaking(z)
#' 
#' cstick_breaking <- compileNimble(stick_breaking)
#' cstick_breaking(z)
NULL

#' @rdname StickBreakingFunction
#' @export
stick_breaking=nimbleFunction(
  run=function(z=double(1),
               log=integer(0, default=0)){ # z values must be different of 1, otherwise the process is truncated to a smaller value tha N; never use z[N]!
    returnType(double(1))
    
    cond <- sum(z < 0)
    if(cond > 0){
      print('values in vector z have to be in [0,1)')
      cond1 <- 1
    }else{
      cond1 <- 0
    }
    cond <- sum(z > 1)
    if(cond > 0){
      print('values in vector z have to be in [0,1)')
      cond2 <- 1
    }else{
      cond2 <- 0
    }
    
    cond <- sum(z == 1)
    if(cond > 0){ # el vector de probabilidades es mas chico
      print('length of returned vector is less than length(z)')
      N <- 1
      while(z[N] != 1){
        N <- N + 1
      }
    }else{
      N<-length(z)
    }
    
    if(cond1 + cond2 > 0){
      return(c(NaN))
    }else{
      x<-numeric(N) # returned vector
      tmp<-0#1
      
      x[1]<-log(z[1]) #z[1]
      for(i in 2:(N-1)){
        tmp=tmp+log(1-z[i-1]) #tmp*(1-z[i-1])
        x[i]<-log(z[i])+tmp #z[i]*tmp
      }
      x[N]<-tmp+log(1-z[N-1])#tmp*(1-z[N-1])
      if(log) return(x)
      else return(exp(x))
    }
  }
)


registerDistributions(list(
    dCRP = list(
        BUGSdist = 'dCRP(conc, size)',
        discrete = TRUE,
        range = c(1, Inf),
        types = c('value = double(1)')
    )
), verbose = FALSE)
