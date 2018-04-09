
#----------------------------------------------#
# Chinese restaurant process distributions:
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

#' The Chinese Restaurant Process Distribution
#'
#'   Density and random generation for the Chinese
#'   Restaurant Process distribution
#' 
#' @name ChineseRestaurantProcess 
#' 
#' @param x vector of values.
#' @param n number of observations (only n = 1 is handled currently).
#' @param conc scalar concentration parameter.
#' @param size integer-valued length of \code{x}.
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Claudia Wehrhahn
#' @details The Chinese restaurant process distribution is a distribution
#' on the space of partitions of the positive integers. 
#' The distribution with concentration parameter \eqn{=\alpha}{= conc} has 
#' probability function 
#' \deqn{f(x_i \mid x_1, \ldots, x_{i-1})=\frac{1}{i-1+\alpha}\sum_{j=1}^{i-1}\delta_{x_j}+
#' \frac{\alpha}{i-1+\alpha}\delta_{x^{new}},}
#' where \eqn{x^{new}} is a new integer not in \eqn{x_1, \ldots, x_{i-1}}.
#' 
#' If \code{conc} is not specified, it assumes the default value of 1. 
#' @return \code{dCRP} gives the density, and \code{rCRP} gives random generation.
#' @references Blackwell, D., and MacQueen, J. B. (1973). Ferguson distributions via 
#' Pólya urn schemes. \emph{The Annals of Statistics}, 1: 353-355.
#' 
#' Aldous, D. J. (1985). Exchangeability and related topics. \emph{In École d'Été de 
#' Probabilités de Saint-Flour XIII -- 1983}, Pages 1-198. Springer.
#' 
#' Pitman, Jim. Some developments of the Blackwell-MacQueen urn scheme (1986). 
#' \emph{Statistics, Probability and Game Theory: Papers in Honor of David Blackwell},
#' 30: 245--267.
#' 
#' @examples
#' x <- rCRP(n=1, conc = 1, size=10)
#' dCRP(x, conc = 1, size=10)
NULL

#' @rdname ChineseRestaurantProcess
#' @export
rCRP <- nimbleFunction(
    run = function(n = integer(0), 
                   conc = double(0, default=1),
                   size = integer(0))
    {
        returnType(double(1))
        
        if(n != 1) {
            stop("rCRP only handles n = 1 at the moment")
        }
        
        if( conc <= 0 ) {
            nimCat("rCRP: value of concentration parameter is not positive. NaNs produced.\n")
            return(rep(NaN, size))
        }

        x <- nimNumeric(size) 
        x[1] <- 1
        if(size > 1) {
            numComponents <- 1
            ones <- rep(1, size)
            for(i in 2:size) {
                if(runif(1) <= conc / (conc+i-1) ) {  # a new value
                    numComponents <- numComponents + 1
                    x[i] <- numComponents
                } else {                             # an old value
                    index <- rcat(n=1, ones[1:(i-1)])
                    x[i] <- x[index]
                }
            }
        }
        
        return(x)
    }
)


#' @rdname ChineseRestaurantProcess
#' @export
dCRP=nimbleFunction(
    run=function(x = double(1), 
               conc = double(0, default=1),
               size = integer(0),
               log = integer(0, default=0))
    {
        returnType(double(0))
        n <- length(x)  
        dens <- nimNumeric(n) 
    
        if(n != size) {
            stop("length of x has to be equal to size")
        }
    
        if( conc <= 0 ) {
            nimCat("value of concentration parameter has to be larger than zero")
            return(NaN)
        }
        
        dens[1] <- 1
        if(n > 1) {
            for(i in 2:n) {
                counts <- sum(x[i] == x[1:(i-1)])
                if( counts > 0 ) {
                    ## Claudia: shouldn't this be 'counts / (i-1+conc)' ???
                    ## Answer: yes
                    dens[i] <- counts / (i-1+conc)
                } else {
                dens[i] <- conc / (i-1+conc)
                }
            }
        }
        
        logProb <- sum(log(dens))
        if(log) return(logProb)
        else return(exp(logProb)) 
    }
)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



#----------------------------------------------#
# stick_breaking nimbleFunction
#----------------------------------------------#

#' The Stick Breaking Function
#'
#' Computes probabilities based on stick breaking construction
#' 
#' @name StickBreakingFunction 
#' 
#' @param z vector argument.
#' @param log logical; if TRUE, weights are returned on the log scale.
#' @author Claudia Wehrhahn
#' @details
#' The stick breaking construction produces a vector of probabilities that sums to one,
#' based on a series of individual probabilities in \code{z} that have to be between \deqn{[0,1)]}.
#' The construction starts by breaking a piece of stick off of a stick of length one, based on
#' \code{z[1]}. The first element the output is the length of the piece that was broken off.
#' The construction then proceeds by breaking a piece of stick off of the remaining stick, based on
#' \code{z[2]} and so forth. If \code{z[k]}
#' is equal to 1, then the returned vector has length smaller than \code{z}. If one of the
#' components is smaller than 0 or greater than 1, NaNs are returned.
#' @references Sethuraman, J. (1994). A constructive definition of Dirichlet 
#' priors. \emph{Statistica Sinica}, 4: 639-650.
#' @examples
#' z <- rbeta(5, 1, 1)
#' stick_breaking(z)
#' 
#' cstick_breaking <- compileNimble(stick_breaking)
#' cstick_breaking(z)
NULL

#' @rdname StickBreakingFunction
#' @export
stick_breaking <- nimbleFunction(
    run = function(z = double(1),
                 log = integer(0, default=0)) 
    {
        returnType(double(1))
    
        N <- length(z)   
        cond <- sum(z < 0) | sum(z > 1)
        if(cond) {
            nimCat("values in 'z' have to be in (0,1)")
            return(rep(NaN, N+1))
        }
    
        x <- nimNumeric(N+1) 
        remainingLogProb <- 0 
    
        x[1] <- log(z[1]) 
        for(i in 2:N) {
            remainingLogProb <- remainingLogProb + log(1-z[i-1]) 
            x[i] <- log(z[i]) + remainingLogProb 
        }
        
        x[N+1] <- remainingLogProb + log(1-z[N]) 
        if(log) return(x)
        else return(exp(x))
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
