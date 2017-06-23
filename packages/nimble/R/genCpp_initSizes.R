
###############################################
## Section for initializing size information ##
###############################################
##
## format for sizeList entries is a vector whose length gives nDim and whose entries are NA for non-fixed size or a number for fixed size

## exprClasses_initSizes
## This function uses the symbol table (or alternatively the sizeList) to initialize typeEnv
## with sizes expressions of known object
##
## exprClasses_initSizes takes as input:
## code is an exprClass object
## symTab is a symbol table of known objects
## sizeList is an alternative way of providing known objects: a list of lists of types and size vectors
## The format of the size vectors is the same as for symTab$sizes:
## size = 0 indicates a scalar, size entries of NA indicate unknown size, size entries of numbers indicate fixed sizes
## The purpose of sizeList is to provide a quicker way to test the function than building a symbolTable for each test case.
##
## returns an environment (typeEnv) that is populated with exprTypeInfoClass objects.  These are lightweight type information, and we'll see if we should use symbolBase (and derived) objects instead.

## The purpose of this environment is to keep track of types *dynamically*.  That is, to the extent that known sizes change in known ways, this can be tracked as we later process each line of code
## Size information for each exprClass nested in the code will be generated later, based on the context in the order of code processing.
##
## It is possible that instead of typeEnv we should just use a symbolTable with copies of entries from the input symbolTable.
## There still needs to be processing because only NEEDED symbols are set up in typeEnv. 


makeSizeExpressions <- function(sizeVec, name) {
    sizeExprs <- vector('list', length(sizeVec))
    if(is.character(name)) name <- as.name(name)
    for(i in seq_along(sizeVec)) {
        sizeExprs[[i]] <- if(is.na(sizeVec[i])) substitute(dim(A)[I], list(A = name, I = i)) else sizeVec[i]
    }
    sizeExprs
}

addToTypeEnv <- function(sym, typeEnv, name) {
    ## If it came from symbolTable:
    if(inherits(sym, 'symbolBase')) {
        if(inherits(sym, 'symbolBasic')) {
            sizeVec <- sym$size
            type <- sym$type
            nDim <- sym$nDim
        } 
        else if(inherits(sym, 'symbolNimbleList')){ ## if return type is a nl, need to store info in typeEnv for use in sizeProcessing
          sizeVec <- 0
          type <- sym$type
          nDim <- 0
        }  
        else {
            return(typeEnv) ## symbol exists but it is something without numeric type info
        }
    } else {
        ## If it came from sizeList
        sizeVec <- sym[[2]]
        type <- sym[[1]]
        nDim <- length(sizeVec)
    }
    
    if(nDim == 0)  {
      if(inherits(sym, 'symbolNimbleList')){
        sizeExprs <- sym
      }
      else{
        ## If nDim == 0, set sizes to 1
        sizeExprs <- list()
      }
    } else {
        ## Otherwise iterate over sizes and build a list
        ## with expressions like dim(A)[2] for unknown second size
        ## or numbers for known sizes.
        if(length(sizeVec)==0) sizeVec <- rep(NA, nDim)
        sizeExprs <- makeSizeExpressions(sizeVec, name)
    }
    
    ## Put the new object in the typeEnv
    assign(name, exprTypeInfoClass$new(nDim = nDim, sizeExprs = sizeExprs, type = type), envir = typeEnv)
    return(typeEnv)
}

exprClasses_initSizes <- function(code, symTab = NULL, sizeList = NULL, typeEnv = new.env(), returnSymbol = NULL) {
    ## If it is a name
    if(!is.null(returnSymbol)) {
        addToTypeEnv(returnSymbol, typeEnv, "return")
    }

    if(code$isName) {
        ## If it is not already in typeEnv
        if(code$name != "") { ## in A[i,], the second index has name == "".  exists() will not accept ""
            if(!exists(code$name, typeEnv, inherits = FALSE)) {
                ## Try finding in symTab
                sym <- if(!is.null(symTab)) symTab$getSymbolObject(code$name, inherits = TRUE) else NULL
                ## If not there, try finding in sizeList
                if(is.null(sym)) sym <- sizeList[[code$name]]
                
                ## If not there, nothing to do (it may be in need of type inference) 
                if(is.null(sym)) return(typeEnv) ## No size info available

                return(addToTypeEnv(sym, typeEnv, code$name))
            }
        }
    }
    ## If it is a call, go recursively into non-constant arguments
    if(code$isCall) {
        for(i in seq_along(code$args)) {
            if(inherits(code$args[[i]], 'exprClass'))
                exprClasses_initSizes(code$args[[i]], symTab, sizeList, typeEnv)
        }
    }
    return(typeEnv)
}
