double <- function(ndim, dims) {}

## Sequential label generation system:
## labelFunctionMetaCreator returns a function that returns a function.
## labelFunctionMetaCreator is only called once, immediately below, to create labelFunctionCreator
## The outer layer allows allLabelFunctionCreators to be in the closure of every function returned
## by labelFunctionCreator.  Each of those functions is registered as an element of allLableFunctionCreators.
## 
## This scheme allows the function resetLabelFunctionCreators below to work simply,
## resetting the count to 1 for all of the label generators.
##
## The motivation for resetLabelFunctionCreators is for testing: If we want to check
## that two pathways to code generation (one existing, one experimental) create identical
## code, it is helpful to have identical generated labels.  Resetting all label generators
## supports this goal.
labelFunctionMetaCreator <- function() {
    allLabelFunctionCreators <- list()

    creatorFun <- function(lead, start = 1) {
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
        allLabelFunctionCreators[[ length(allLabelFunctionCreators) + 1 ]] <<- labelGenerator
        labelGenerator
    }
    creatorFun
}

labelFunctionCreator <- labelFunctionMetaCreator()

resetLabelFunctionCreators <- function() {
    allLabelFunctionCreators <- environment(labelFunctionCreator)$allLabelFunctionCreators
    for(i in allLabelFunctionCreators) {
        i(reset = TRUE)
    }
}

nimbleUniqueID <- labelFunctionCreator("UID")
nimbleModelID  <- labelFunctionCreator("MID")
ADinfoLabel <- labelFunctionCreator("nimbleCppADinfo_")
StridedLabelMaker <- labelFunctionCreator("Strided")

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
#' @aliases dim
#'
#' @examples
#' x <- rnorm(4)
#' dim(x)
#' nimDim(x)
#' y <- matrix(x, nrow = 2)
#' dim(y)
#' nimDim(y)
#'
#' @export
nimDim <- function(obj) {
    if(is.null(dim(obj))) return(length(obj))
    return(dim(obj))
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
#' 
#' @export
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


## creates objects in the parent.frame(), named by names(lst), values are eval(lst[[i]])
## this is used for creating the conjugate sampler nimble function generators, among other things
createNamedObjectsFromList <- function(lst, writeToFile = NULL, envir = parent.frame()) {
    for(i in seq_along(lst)) {
        objName <- names(lst)[i]
        obj <- eval(lst[[i]])
        assign(objName, obj, envir = envir)
    }
    if(!is.null(writeToFile)) {
        write('', file = writeToFile)
        for(i in seq_along(lst)) {
            expr <- substitute(VAR <- VALUE, list(VAR = as.name(names(lst)[i]), VALUE = lst[[i]]))
            deparseExpr <- deparse(expr, control=c())
            deparseExpr <- gsub('\"', '\'', deparseExpr)
            write(deparseExpr, file = writeToFile, append = TRUE)
            write('\n\n\n', file = writeToFile, append = TRUE)
        }
    }
}

#' Print error messages after failed compilation
#' 
#' Retrieves the error file from R's tempdir and prints to the screen.
#'
#' @param excludeWarnings logical indicating whether compiler warnings should be printed; generally such warnings can be ignored.
#'
#' @author Christopher Paciorek
#' 
#' @export
printErrors <- function(excludeWarnings = TRUE) {
    if(exists('errorFile', nimbleUserNamespace)) {
        errors <- readLines(nimbleUserNamespace$errorFile)
        if(excludeWarnings)
            errors <- grep("^[Ww]arning", errors, value = TRUE, invert = TRUE)
        cat(errors, sep = "\n")
    } else cat("No error file found.")
}
