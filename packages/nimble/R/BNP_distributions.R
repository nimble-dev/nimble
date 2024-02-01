#' The Chinese Restaurant Process Distribution
#'
#'   Density and random generation for the Chinese
#'   Restaurant Process distribution.
#' 
#' @name ChineseRestaurantProcess 
#' 
#' @param x vector of values.
#' @param n number of observations (only n = 1 is handled currently).
#' @param conc scalar concentration parameter.
#' @param size integer-valued length of \code{x} (required).
#' @param log logical; if TRUE, probability density is returned on the log scale.
#' @author Claudia Wehrhahn
#' @details The Chinese restaurant process distribution is a distribution
#' on the space of partitions of the positive integers. 
#' The distribution with concentration parameter \eqn{\alpha} equal to \code{conc} 
#' has probability function 
#' \deqn{
#' f(x_i \mid x_1, \ldots, x_{i-1})=\frac{1}{i-1+\alpha}\sum_{j=1}^{i-1}\delta_{x_j}+
#' \frac{\alpha}{i-1+\alpha}\delta_{x^{new}},}
#' where \eqn{x^{new}} is a new integer not in \eqn{x_1, \ldots, x_{i-1}}.
#' 
#' If \code{conc} is not specified, it assumes the default value of 1. The \code{conc} 
#' parameter has to be larger than zero. Otherwise, \code{NaN} are returned.
#' @return \code{dCRP} gives the density, and \code{rCRP} gives random generation.
#' @references Blackwell, D., and MacQueen, J. B. (1973). Ferguson distributions via 
#' \enc{Pólya}{Polya} urn schemes. \emph{The Annals of Statistics}, 1: 353-355.
#' 
#' Aldous, D. J. (1985). Exchangeability and related topics. In \emph{\enc{École}{Ecole} d'\enc{Été}{Ete} 
#' de \enc{Probabilités}{Probabilites} de Saint-Flour XIII - 1983} (pp. 1-198). Springer, Berlin, 
#' Heidelberg.
#' 
#' Pitman, J. (1996). Some developments of the Blackwell-MacQueen urn scheme. \emph{IMS Lecture
#' Notes-Monograph Series}, 30: 245-267.
#' 
#' @examples
#' x <- rCRP(n=1, conc = 1, size=10)
#' dCRP(x, conc = 1, size=10)
NULL

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
    
        if(n != size) {
            stop("dCRP: length of 'x' has to be equal to 'size'.\n")
        }
    
        if( conc <= 0 | is.na(conc) ) {
        #    nimCat("dCRP: value of concentration parameter has to be larger than zero.\n")
            return(NaN)
        }
        if(any_na(x)) return(NaN)
        
        ldens <- 0 # log scale
        if(n > 1) {
            for(i in 2:n) {
                counts <- sum(x[i] == x[1:(i-1)])
                if( counts > 0) {
                    ldens <- ldens + log(counts / (i-1+conc))
                } else {
                    ldens <- ldens + log(conc / (i-1+conc))
                }
            }
        }
        
        if(log) return(ldens)
        else return(exp(ldens)) 
    }
)

#' @rdname ChineseRestaurantProcess
#' @export
rCRP <- nimbleFunction(
    run = function(n = double(0), 
                   conc = double(0, default=1),
                   size = integer(0))
    {
        returnType(double(1))
        
        if(n != 1) {
            stop("rCRP only handles n = 1 at the moment.\n")
        }
        
        if( conc <= 0 | is.na(conc) ) {
        #    nimCat("rCRP: value of concentration parameter is not positive. NaNs produced.\n")
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

#' The Stick Breaking Function
#'
#' Computes probabilities based on stick breaking construction.
#' 
#' @name StickBreakingFunction
#'
#' @aliases stickbreaking
#' 
#' @param z vector argument.
#' @param log logical; if \code{TRUE}, weights are returned on the log scale.
#' @author Claudia Wehrhahn
#' @details
#' The stick breaking function produces a vector of probabilities that add up to one,
#' based on a series of individual probabilities in \code{z}, which define the breaking
#' points relative to the remaining stick length. The first element of \code{z} determines
#' the first probability based on breaking a proportion \code{z[1]} from a stick of length one.
#' The second element of \code{z} determines the second probability based on breaking a
#' proportion \code{z[2]} from the remaining stick (of length \code{1-z[1]}), and so forth.
#' Each element of \code{z} should be in 
#' \eqn{(0,1)}.
#' The returned vector has length equal to the length of \code{z} plus 1. 
#' If \code{z[k]} is equal to 1 for any \code{k}, then the returned vector has length smaller than \code{z}. 
#' If one of the components is smaller than 0 or greater than 1, \code{NaN}s are returned.
#' @references Sethuraman, J. (1994). A constructive definition of Dirichlet priors.
#'  \emph{Statistica Sinica}, 639-650.
#' @examples
#' z <- rbeta(5, 1, 1)
#' stick_breaking(z)
#'
#' \dontrun{
#' cstick_breaking <- compileNimble(stick_breaking)
#' cstick_breaking(z)
#' }
NULL

#' @rdname StickBreakingFunction
#' @export
stick_breaking <- nimbleFunction(
    run = function(z = double(1),
                 log = integer(0, default=0)) 
    {
        returnType(double(1))
    
        N <- length(z)   
        cond <- any(z < 0 | z > 1)
        if(cond) {
            nimCat("  [Warning] stick_breaking: values in 'z' have to be in (0,1).\n")
            return(rep(NaN, N+1))
        }
    
        x <- nimNumeric(N+1) 
        remainingLogProb <- 0 
    
        x[1] <- log(z[1])
        if(N > 1) {
            for(i in 2:N) {
                remainingLogProb <- remainingLogProb + log(1-z[i-1]) 
                x[i] <- log(z[i]) + remainingLogProb 
            }
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
