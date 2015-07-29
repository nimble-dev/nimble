## CASE 2: modelVariableAccessorVector
## copy(model1, xxx, nodes) becomes:
## model1_nodes_accessors <- modelVariableAccessorVector(model1, nodes)
## copy(model1_nodes_accessors, xxx)
## model1_nodes_accessors$getAccessors() returns a list of the modelVariableAccessor objects

## I think this is not needed anymore
## modelVariableAccessor <- setRefClass(
##     Class = 'modelVariableAccessor',
##     fields = list(model = 'ANY',
##                   var   = 'ANY', 		# 'character',
##                   first = 'ANY', 		#'numeric',
##                   last  = 'ANY', 		#'numeric',
##                   length = 'ANY' 		#'numeric'
##     ),
##     methods = list(toStr = function() paste0(var, '[', first, ':', last, ']'),
##                    show  = function() cat(paste0(toStr(), '\n'))
##     )
## )

modelVariableAccessorVector <- function(mv, nodeNames, logProb = FALSE) {
    ans <- list(mv, substitute(nodeNames), logProb, parent.frame())
    class(ans) <- c('modelVariableAccessorVector', 'valuesAccessorVector')
    ans
}

## getValues.modelVariableAccessorVector <- function(accessorVector, i) {
##     nodeName <- eval(accessorVector[[2]], envir = accessorVector[[4]])[i]
##     nodeCode <- parse(text = nodeName, keep.source = FALSE)[[1]]
##     sourceObject <- accessorVector[[1]]
##     if(is.name(nodeCode)) return(sourceObject[[nodeName]])
##     nodeCode[[2]] <- substitute(sourceObject[[VN]], list(VN = as.character(nodeCode[[2]]) ) )
##     eval(nodeCode)
## }
## setValues.modelVariableAccessorVector <- function(accessorVector, i, vals, row = NA) {
##     if(is.na(row)) stop('Can no longer use row = NA in setValues.modelVariableAccessorVector')
##     nodeName <- eval(accessorVector[[2]], envir = accessorVector[[4]])[i]
##     nodeCode <- parse(text = nodeName, keep.source = FALSE)[[1]]
##     sourceObject <- accessorVector[[1]]
##     if(is.name(nodeCode)) nodeCode <- quote(sourceObject[[nodeName]])
##     nodeCode[[2]] <- substitute(sourceObject[[VN]], list(VN = as.character(nodeCode[[2]]) ) )
##     eval(substitute(A <- vals), list(A = nodeCode))
## }


## modelVariableAccessorVector<- setRefClass( ## new implementation
##     Class = 'modelVariableAccessorVector',
##     contains = 'valuesAccessorVector',
##     methods = list(
##         initialize = function(...) {
##             callSuper(...)
##             makeAccessAndSetCode()
##         },
##         makeAccessAndSetCode = function() {
##             accessCode <<- lapply(code, function(temp) {
##                 if(is.name(temp)) return(substitute(sourceObject$B, list(B = temp)))
##                 temp[[2]] <- substitute(sourceObject$B, list(B = temp[[2]]))
##                 temp
##             })
##             setCode <<- lapply(accessCode, function(x) substitute(A <- vals, list(A = x)))
##         },
##         getValues = function(i) {
##             eval(accessCode[[i]])
##         },
##         setValues = function(i, vals) {
##             eval(setCode[[i]])
##         }
##     ))


## # This function allows you to just call "length(access)", rather than access$length
## # Motivation for this is to make it easier to generically use accessors
## length.modelVariableAccessorVector <- function(access)
##     return(access$length)
