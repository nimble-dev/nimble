double <- function(ndim, dims) {}

labelFunctionCreator <- function(lead, start = 1) {
  nextIndex <- start
  force(lead)
  labelGenerator <- function(reset = FALSE, count = 1, envName = "") {
    if(reset) {
      nextIndex <<- 1
      return(invisible(NULL))
    }
    lead <- paste(lead, envName , sep = '_')
    ans <- paste0(lead, nextIndex - 1 + (1:count))
    nextIndex <<- nextIndex + count
    ans
  }
  labelGenerator
}

nimbleUniqueID <- labelFunctionCreator("UID")

dimOrLength <- function(obj, scalarize = FALSE) {
    if(scalarize) if(length(obj) == 1) return(numeric(0))
    if(is.null(dim(obj))) return(length(obj))
    return(dim(obj))
}

#' return sizes of an object whether it is a vector, matrix or array
#'
#' R's regular \code{dim} function returns NULL for a vector.  It is useful to have this function that treats a vector similarly to a matrix or array.  Works in R and NIMBLE.  In NIMBLE \code{dim} is identical to \code{nimDim}, not to R's \code{dim}
#'
#' @param obj  objects for which the sizes are requested
#'
#' @return a vector of sizes in each dimension
#'
#' @author NIMBLE development team
#'
#' @examples
#' x <- rnorm(4)
#' dim(x)
#' nimDim(x)
#' y <- matrix(x, nrow = 2)
#' dim(y)
#' nimDim(y)
nimDim <- function(obj) {
    if(is.null(dim(obj))) return(length(obj))
    return(dim(obj))
}

getLoadingNamespace <- function() {
    if(!is.null(nimbleOptions()$notUsingPackage)) if(nimbleOptions()$notUsingPackage) return(globalenv())
    if(system.file(package = "nimble") == "")
         return(globalenv())
    
    nimbleenv <- getNamespace('nimble')
    if(environmentIsLocked(nimbleenv)) globalenv() else nimbleenv
}

getNimbleFunctionEnvironment <- function() {
    ## This is how we determine which environment will contain the reference class definition underlying a NIMBLE function.
    ## This implementation allows for:
    ## (1) nimbleFunctions to be defined in the namespace of another package (in a R/ source file), and
    ## (2) nimbleFunctions to be defined inside other package member functions
    env <- topenv(parent.frame(2))  # 2 gets out of nimble namespace and into frame where topenv gives the namespace in which the nf is being defined
    if(!environmentIsLocked(env)) return(env)
    return(globalenv())
}


#' Resizes a modelValues object
#'
#' Adds or removes rows to a modelValues object. If rows are added to a modelValues object, the default values are \code{NA}. Works in both R and NIMBLE.
#'
#' @param container		modelValues object
#' @param k				number of rows that modelValues is set to
#' @author				Clifford Anderson-Bergman
#' @details
#' See the User Manual or \code{help(modelValuesBaseClass)} for infomation about modelValues objects
#'
#' @examples
#' mvConf <- modelValuesConf(vars = c('a', 'b'),
#'              types = c('double', 'double'),
#'              sizes = list(a = 1, b = c(2,2) ) ) 
#' mv <- modelValues(mvConf)
#' as.matrix(mv)
#' resize(mv, 3)
#' as.matrix(mv)
resize <- function(container, k) {
    container$resize( as.integer(k) ) 
}

#' Returns number of rows of modelValues
#' 
#' Returns the number of rows of NIMBLE modelValues object. Works in R and NIMBLE. 
#' 
#' @param container		modelValues object
#' @author Clifford Anderson-Bergman
#' @export
#' @details
#' See the User Manual or \code{help(modelValuesBaseClass)} for information about modelValues objects 
#'
#' @examples
#'	mvConf <- modelValuesConf(vars = 'a', types = 'double', sizes = list(a = 1) )
#'	mv <- modelValues(mvConf)
#'  resize(mv, 10)
#'	getsize(mv)
getsize <- function(container) {
    container$getSize()
}

# simply adds width.cutoff = 500 as the default to deal with creation of long variable names from expressions
deparse <- function(...) {
    if("width.cutoff" %in% names(list(...))) {
        base::deparse(...)
    } else {
          base::deparse(..., width.cutoff = 500L)
      }
}

