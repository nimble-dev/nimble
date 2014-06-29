
## This recurses through the exprClass parse tree and pulls out any maps
## That means it generates the new varianble and the setMap call.
## This uses the assertions field, which will be cleared from each line first
## This must be done after eigenization, since that will grab some of the maps and make them Eigen Map objects.
## What is left will be maps to be handled with NimArrs
exprClasses_liftMaps <- function(code, symTab, typeEnv) {
    if(code$isCall) {
        if(code$name == '{') {
            for(i in seq_along(code$args)) {
                if(inherits(code$args[[i]], 'exprClass')) {
                    newAsserts <- exprClasses_liftMaps(code$args[[i]], symTab, typeEnv)
                    newExpr <- newBracketExpr(args = c(lapply(newAsserts, function(x) if(inherits(x, 'exprClass')) x else RparseTree2ExprClasses(x)), code$args[i]))
                    setArg(code, i, newExpr)
                }
            }
            return(NULL)
        }
        if(code$name == 'map') {
            forceAssign <- code$caller$name == 'return' ## this case should be moot: sizeReturn was modified so any maps should be lifted earlier
            newAssert <- sizeInsertIntermediate(code$caller, code$callerArgID, symTab, typeEnv, forceAssign = forceAssign)
            return(newAssert)
        }
        asserts <- list()
        for(i in seq_along(code$args)) {
            if(inherits(code$args[[i]], 'exprClass')) {
                if(code$args[[i]]$isCall) {
                    newAsserts <- exprClasses_liftMaps(code$args[[i]], symTab, typeEnv)
                    asserts <- c(asserts, newAsserts)
                }
            }
        }
        return(if(length(asserts)==0) NULL else asserts)
    }
}
