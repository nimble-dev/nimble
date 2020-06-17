#' set the size of a numeric variable in NIMBLE
#'
#' set the size of a numeric variable in NIMBLE.  This works in R and NIMBLE, but in R it usually has no effect.
#'
#' @param numObj    This is the object to be resized
#' @param ...       sizes, provided as scalars, in order, or as a single vector
#' @param copy      logical indicating whether values should be preserved (in column-major order)
#' @param fillZeros logical indicating whether newly allocated space should be initialized with zeros (in compiled code)
#' 
#' @author NIMBLE development team
#' @export
#' @details
#' Note that assigning the result of \code{numeric}, \code{integer}, \code{logical}, \code{matrix}, or \code{array} is often as good or better than using \code{setSize}.  For example, `x <- matrix(nrow = 5, ncol = 5)` is equivalent to `setSize(x, 5, 5)` but the former allows more control over initialization.
#' 
#' This function is part of the NIMBLE language.  Its purpose is to explicitly resize a multivariate object (vector, matrix or array), currently up to 4 dimensions.  Explicit resizing is not needed when an entire object is assigned to.  For example, in \code{Y <- A \%*\% B}, where A and B are matrices, \code{Y} will be resized automatically.  Explicit resizing is necessary when assignment will be by indexed elements or blocks, if the object is not already an appropriate size for the assignment.  E.g. prior to \code{Y[5:10] <- A \%*\% B}, one can use setSize to ensure that \code{Y} has a size (length) of at least 10.
#'
#' This does work in uncompiled (R) and well as compiled execution, but in some cases it is only necessary for compiled execution. During uncompiled execution, it may not catch bugs due to resizing because some R objects will be dynamically resized during assignments anyway.
#'
#' If preserving values in the resized object and/or initializing new values with 0 is not necessary, then setting these arguments to FALSE will yield slightly more efficient compiled code.
#' 
setSize <- function(numObj, ..., copy = TRUE, fillZeros = TRUE){ ## fillValues isn't used here but is included for consistency with the DSL
    thisCall <- as.list(match.call()[-1])
    if(length(thisCall) < 2) stop("No information provided to setSize")
    newDimsList <- list(...)
    if(is.numeric(numObj) | is.logical(numObj)) {
        targetVar <- deparse(thisCall[[1]])
        if(!exists(targetVar, envir = parent.frame()))
            stop(paste0("Variable ", targetVar, " does not exist."))
        targetIsLocal <- exists(targetVar, envir = parent.frame(), inherits = FALSE)
        oldClass <- class(numObj)
        oldDims <- nimDim(numObj)
        if(length(oldDims) != length(newDimsList)) {
            ## newDims could have been provided as a vector
            if(length(newDimsList) > 1)
                stop("Number of dimensions provided does not match object to change in setSize")
            newDims <- newDimsList[[1]]
        } else {
            newDims <- unlist(newDimsList)
        }
        if(any(is.na(newDims))) warning("Not sure what to do with NA dims in setSize")
        if(length(newDims) != length(oldDims))
            if(length(newDims) < length(oldDims))
                stop("Number of dimensions provided does not match object to change in setSize")
            else
                warning("Number of dimensions provided to setSize does not match object to change in setSize.  Sizes will be truncated.")
        newDims <- newDims[1:length(oldDims)]

        if(length(oldDims) == 1) {
            if(oldDims[1] < newDims[1]) {
                newObj <- as(rep(0, newDims[1]), oldClass)
                if(copy) newObj[seq_len(oldDims[1])] <- numObj
            } else {
                if(copy) newObj <- numObj[seq_len(newDims[1])]
                else newObj <- as(rep(0, newDims[1]), oldClass)
            }
        } else {
            newObj <- array(0, newDims)
            if(copy) {
                if(length(numObj) < length(newObj))
                    newObj[seq_len(length(numObj))] <- numObj
                else
                    newObj[seq_len(length(newObj))] <- numObj[seq_len(length(newObj))]
            }
        }
        assign("_SETSIZE_TEMP_VAL", newObj, parent.frame())
        if(targetIsLocal)
            assignCall <- substitute(A <- B, list(A = thisCall[[1]], B = as.name("_SETSIZE_TEMP_VAL")))
        else
            assignCall <- substitute(A <<- B, list(A = thisCall[[1]], B = as.name("_SETSIZE_TEMP_VAL")))
        eval(assignCall, envir = parent.frame())
        rm("_SETSIZE_TEMP_VAL", envir = parent.frame())
        return(invisible(NULL))
    }
    ## what follows is partly defunct -- to be examined
    dims = newDims
    dims <- dims[!is.na(dims)]
    if(inherits(numObj, "RNumericList") ) {
        oldDims = nimDim(numObj[[row]])
        if(length(oldDims) != length(dims) )
            stop("Wrong number of dimensions")
        numObj$Values[[row]] <- array(data = 0, dim = dims)
    }
}
	
