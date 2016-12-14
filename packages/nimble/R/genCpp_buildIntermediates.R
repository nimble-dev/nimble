IntermLabelMaker <- labelFunctionCreator('Interm')

## removing nimArrayGeneral here because now we allow matrix(non-scalar,...) amid an expression
buildIntermCalls <- c(makeCallList(c('eigen', 'chol','run.time'), 'buildSimpleIntermCall')
                      )


## exprClasses_buildInterms
##
## This works recursively through the parse tree (an exprClass object)
## It checks every call for whether the name of the call is in buildIntermCalls
## If so it calls the function named by that element of buildIntermCalls and replaces the original
## line with a bracket expression (`{`) containing the intermediate calls followed by the modified original
## This turned out not to be needed much in the compilation system, but it is handy in a few cases.
## Many intermediates are created during the size processing and eigenization.
exprClasses_buildInterms <- function(code) {
    tempExprs <- list()
    if(code$isCall) {
        if(code$name == '{') {
            for(i in seq_along(code$args)) {
                if(!(code$args[[i]]$isCall)) {
                    writeLines('Warning in buildInterms: there is a line of code that is not a call.')
                    return(NULL)
                } else {
                    intermCalls <- unlist(exprClasses_buildInterms(code$args[[i]]))
                    if(length(intermCalls) > 0) {
                        newExpr <- newBracketExpr(args = c(intermCalls, code$args[i]))
                        setArg(code, i, newExpr)
                    }
                }
            }
            return(NULL)
        }

        ## first get temps from any arguments
        for(i in seq_along(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
                if(code$args[[i]]$isCall) {
                    tempExprs[[length(tempExprs) + 1]] <- exprClasses_buildInterms(code$args[[i]])
                }
            }
        }
     
        ## second generate any temps for this call itself
        intermCall <- buildIntermCalls[[code$name]]
        if(!is.null(intermCall)) {
            newpiece <- eval(call(buildIntermCalls[[code$name]], code))
            tempExprs[[length(tempExprs) + 1 ]] <- newpiece
        }
    }
    if(length(tempExprs) == 0) NULL else tempExprs
}

buildIntermNone <- function(code) {
    return(NULL)
}
