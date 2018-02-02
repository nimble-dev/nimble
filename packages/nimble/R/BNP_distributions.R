
#----------------------------------------------#
# Chinesse restaurant process distributions:
#----------------------------------------------#

# random number generator: rCRP
#   input:
#      n: sample size
#      conc: concentration parameter of the DP (related to sampling a new value)
#   output:
#      x: vector of size n
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
#' @param n number of observations.
#' @param conc scalar with concentration parameter.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Claudia Wehrhahn
#' @details The chineese restaurant process distribution is a distribution
#' on the space of partitions of the positive integers. 
#' The chineese restaurant process
#' distribution with concentration parameter \eqn{=\alpha}{= conc} has 
#' probability function 
#' \deqn{f(x_i)=\frac{\alpha}{i-1+\alpha}}
#' 
#' If \code{conc} is not sepcified, it assumes the default value of 1. 
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
#' x <- rCRP(10, conc = 1)
#' dCRP(x, conc = 1)
NULL

#' @rdname ChinesseRestaurantProcess
#' @export
rCRP=nimbleFunction(
  run=function(n=integer(0), 
               conc=double(0, default=1))
  {
    returnType(double(1))
    
    n <- floor(n)
    if( length(n) > 1){
      n <- floor(n[1])
      print('length(n) > 1 only the first element will be used')
    }
    
    if( conc[1] <= 0 ){
      concCond1 <- 0
      print('conc parameter has to be larger than zero')
    }else{
      concCond1 <- 1
    }
    
    if( length(conc) > 1 ){
      concCond2 <- 0
      print('conc parameter has to be a scalar')
    }else{
      concCond2 <- 1
    }
    
    if(concCond1 == 0 ){
      return(NaN)
    }
    if( concCond2 == 0 ){
      return(NaN)
    }
    if( concCond1+concCond2 == 2){
      x=numeric(n)# returned vector:
      x[1]=1
      if(n>1){
        for(i in 2:n){
          if(runif(1)<=conc/(conc+i-1)){# a new value
            x[i]=max(x[1:(i-1)])+1 
          }else{# an old value
            index=rcat(n=1, rep(1, i-1)) 
            x[i]=x[index]
          }
        }
      }
      return(x)
    }
  }
)



dCRP=nimbleFunction(
  run=function(x=double(1), 
               conc=double(0, default=1), 
               log=integer(0, default=0))
  {
    returnType(double(0))
    n=length(x) 
    tmpden=numeric(n) 
    
    if( conc[1] <= 0 ){
      concCond1 <- 0
      print('conc parameter has to be larger than zero')
    }else{
      concCond1 <- 1
    }
    
    if( length(conc) > 1 ){
      concCond2 <- 0
      print('conc parameter has to be a scalar')
    }else{
      concCond2 <- 1
    }
    
    if(concCond1 == 0 ){
      return(NaN)
    }
    if( concCond2 == 0 ){
      return(NaN)
    }
    if( concCond1+concCond2 == 2){
      tmpden[1]=1
      if(n>1){
        for(i in 2:n){
          #counts=0 # replaces sum(x[i]==x[1:(i-1)]) (works in nimble?)
          #for(j in 1:(i-1)){
          #  if(x[i]==x[j]){
          #    counts=counts+1
          #  }
          #}
          counts=sum(x[i]==x[1:(i-1)])
          if(counts>0){
            tmpden[i]=1/(i-1+conc)
          }else{
            tmpden[i]=conc/(i-1+conc)
          }
        }
      }
      logProb <- sum(log(tmpden))
      if(log) return(logProb)
      else return(exp(logProb)) 
    }
  }
)


