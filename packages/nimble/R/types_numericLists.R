NumericListBase <- setRefClass(Class = "NumericListBase")

RNumericList <- setRefClass(Class = "RNumericList",
                            contains = "NumericListBase",
                            fields = list(listType = "character",
                                nDim = "numeric", 
                                Values = "list",
                                nRow = "numeric"),
                            methods = list(initialize = function (listType = "double", nDim = 1, nRow = 1){
                                listType <<- listType
                                nDim <<- nDim
                                nRow <<- nRow
                                Values <<- list()
                                if(nRow >0)
                                    for(i in 1:nRow)
                                        Values[[i]] <<- array(data = 0, dim = rep(1, nDim) )
                            },
                                
                                show = function(){
                                    writeLines("R-Based Numeric List")
                                    cat("Number of Dimensions = ", nDim, "Data Type = ", listType, "\nNumber of Rows = ", nRow, "\n" )
                                },
                                resize = function(rows){
                                    if(nRow == rows)
                                        break
                                    else if(nRow > rows){
                                        Values <<- numList$Values[1:rows]
                                        nRow <<- rows
                                    }
                                    else if(nRow < rows){
                                        ndims = nDim
                                        for(i in (nRow+1):rows)
                                            Values[[i]] <<- array(0, dim = rep(1, ndims))
                                        nRow <<- rows
                                    }
                                },
                                getSize = function(){
                                    return(nRow)
                                }
                                )
                            
                            )

setMethod('[[', 'RNumericList',
          function(x, i)
          return(x$Values[[i]])
          )

setMethod('[[<-', 'RNumericList',
          function(x, i, value) 
          {
              if(i > x$nRow)
                  stop("Index is higher than max row of numericList")
              rowDims <- nimble:::dimOrLength(x[[i]])
              inputDims <- nimble:::dimOrLength(value)
              if(length(rowDims) != length(inputDims) )
                  stop("Incorrect number of dimensions for numericList")
              if(any(rowDims != inputDims) ) 
                  stop("Dimensions of row incorrect. Try resizing row first")
              x$Values[[i]] <- value
              return(x)
          })
				
numericList <- function(type = double(2), length = 1, buildType = 'R', extPtr = NA){
    typeSub <- substitute(type)
    if(as.character(typeSub[[1]]) == "double")
        type = "double"
    else if(as.character(typeSub[[1]]) == "integer")
        type = "int"
    else 
        stop("Specified type not currently supported")
    nDim = typeSub[[2]]
    if(buildType == 'R')
        return(RNumericList(listType = type, nDim = nDim, nRow = length))
    if(inherits(extPtr, 'externalptr') ) 
        return(CNumericList(extPtr = extPtr))
    return(CNumericList( extPtr = .Call('makeNumericList',  as.integer(nDim), type, as.integer(length) ) ) )
}

#' set the size of a numeric variable in NIMBLE
#'
#' set the size of a numeric variable in NIMBLE.  This works in R and NIMBLE, but in R it usually has no effect.
#'
#' @param numObj   This is the object to be resized
#' @param ...      sizes, provided as scalars, in order, with exactly as many as needed for the object
#' @param row      Optional argument that is not currently used
#'
#' @author NIMBLE development team
#' @export
#' @details
#' This function is part of the NIMBLE language.  Its purpose is to explicitly resize a multivariate object (vector, matrix or array), currently up to 4 dimensions.  Explicit resizing is not needed when an entire object is assigned to.  For example, in \code{Y <- A \%*\% B}, where A and B are matrices, \code{Y} will be resized automatically.  Explicit resizing is necessary when assignment will be by indexed elements or blocks, if the object is not already an appropriate size for the assignment.  E.g. prior to \code{Y[5:10] <- A \%*\% B}, one can use setSize to ensure that \code{Y} has a size (length) of at least 10.
#'
#' This does work in uncompiled (R) and well as compiled execution, but in some cases it is only necessary for compiled execution. During uncompiled execution, it may not catch bugs due to resizing because some R objects will be dynamically resized during assignments anyway.
setSize <- function(numObj, ..., row){
    thisCall <- as.list(match.call()[-1])
    if(length(thisCall) < 2) stop("No information provided to setSize")
    newDims <- unlist(list(...))
    if(any(is.na(newDims))) warning("Not sure what to do with NA dims in setSize")
    if(is.numeric(numObj)) {
        oldDims <- nimDim(numObj)
        if(length(oldDims) != length(newDims)) stop("Number of dimensions provided does not match object to change in setSize") 
        if(length(oldDims) == 1) {
            if(oldDims[1] < newDims[1]) {
                newObj <- rep(0, newDims[1])
                newObj[1:oldDims[1]] <- numObj
            } else {
                newObj <- numObj[1:newDims[1]]
            }
        } else {
            newObj <- array(0, newDims)
            if(length(numObj) < length(newObj))
                newObj[1:length(numObj)] <- numObj
            else
                newObj[1:length(newObj)] <- numObj[1:length(newObj)]
        }
        assign("_SETSIZE_TEMP_VAL", newObj, parent.frame())
        assignCall <- substitute(A <- B, list(A = thisCall[[1]], B = as.name("_SETSIZE_TEMP_VAL")))
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
    if(inherits(numObj, "CNumericList") ) 
        jnk <- .Call("resizeNumListRow", numObj$valuesPtr, as.integer(row), as.integer(dims) )
}

CNumericList <- setRefClass(Class = "CNumericList",
                            contains = "NumericListBase",
                            fields = list(listType = "character",
                                valuesPtr = "ANY", 
                                nRow = function(x)
                                {
                                    if(missing(x))
                                        getCRows(valuesPtr)
                                    else
                                        writeLines("Need to set nRow using setRows function")
                                }
                                ),
                            methods = list(initialize = function (extPtr){
                                valuesPtr <<- extPtr
                                
                            },
                                
                                show = function(){
                                    writeLines("C-Based Numeric List")
                                        #				cat("Number of Dimensions = ", nDim, "Data Type = ", listType, "\nNumber of Rows = ", getCRows(valuesPtr), "\n" )
                                }
                                
                                
                                ,
                                resize = function(rows){
                                    jnk <- .Call("setNumListRows", valuesPtr, as.integer(rows), FALSE ) 
                                },
                                getSize = function()
                                getCRows(valuesPtr)
                                )
                            )

setMethod('[[', 'CNumericList', 
          function(x, i){
              if(i < 0 | i > x$nRow)
                  stop("Invalid index for this numericList")
              getCNumericListValue(x$valuesPtr, i)
          }
          )
		
setMethod('[[<-', 'CNumericList',
          function(x, i, value){
              if(i < 0 | i > x$nRow)
                  stop("Invaid index for this numericList")
              setCNumericListValue(x$valuesPtr, value, i)
              return(x)
          }
          )

	
makeCNumericListPtr <- function(nDims = 1, type = "double", nRows = 1)
	.Call("makeNumericList", as.integer(nDims), as.character(type), as.integer(nRows))
	
setCNumericListValue <- function(listPtr, values, row)
	jnk <- .Call("setMVElement", listPtr, as.integer(row) , values ) 
	
getCNumericListValue <- function(list, row)
	.Call("getMVElement", list, as.integer(row))	
	
getCRows = function(list)
	.Call("getNRow", list)
	
