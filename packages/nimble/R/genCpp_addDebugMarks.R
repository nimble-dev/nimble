###########################
## Add debugging markers ##
###########################

exprClasses_addDebugMarks <- function(code, label = '', counterEnv = NULL) {
    if(is.null(counterEnv)) {counterEnv <- new.env(); counterEnv$count <- 1}
    if(code$isCall) {
        if(code$name == '{') {
            for(i in seq_along(code$args)) {
                ## replace with recursed result
                setArg(code, i, exprClasses_addDebugMarks(code$args[[i]], label, counterEnv))
            }
            return(code)
        }
        if(code$name == 'for') {
            setArg(code, 3, exprClasses_addDebugMarks(code$args[[3]], label, counterEnv))
            return(code)
        }
        if(code$name %in% ifOrWhile) {
            setArg(code, 2, exprClasses_addDebugMarks(code$args[[2]], label, counterEnv))
            if(length(code$args) == 3) setArg(code, 3, exprClasses_addDebugMarks(code$args[[3]], label, counterEnv))
            return(code)
        }
        ## ditto for if
        ## return bracketed with PRINTF
        newLine <- RparseTree2ExprClasses(substitute(Rprintf(MSG), list(MSG = paste(label, counterEnv$count, "\\n"))))
        ans <- RparseTree2ExprClasses(quote({A; B}))
        setArg(ans, 1, newLine)
        setArg(ans, 2, code)
        counterEnv$count <- counterEnv$count + 1
        return(ans)
    }
    return(code)
}
