
##################################
## Section for insertAssertions ##
##################################

## exprClasses_insertAssertions
## Takes assertions, which were built into exprClass objects (in the assertions field)
## during exprClasses_setSizes, and turns them into code inserted into the parse tree (exprClass)
## This works by heavy use of '{' expressions to expand code
## E.g. if there is a line A <- B, for which the exprClass object has an assertion to resize A to the dimensions of B
## The new code will be { A.resize(dim(B)[1], dim(B)[2]); A <- B }
## Since '{' is a call, with each line an argument, this is all just one call and so can replace the oringal A<-B call
## During C++ generation, such extra brackets are omitted 
exprClasses_insertAssertions <- function(code) {
    if(code$name != '{') return(invisible(NULL)) ## This function only iterates over lines enclosed in '{'
    for(i in seq_along(code$args)) {
        if(code$args[[i]]$isCall) {
            if(code$args[[i]]$name == 'for') {
                exprClasses_insertAssertions(code$args[[i]]$args[[3]])
            }
            if(code$args[[i]]$name %in% ifOrWhile) {
                exprClasses_insertAssertions(code$args[[i]]$args[[2]])
                if(length(code$args[[i]]$args) == 3) exprClasses_insertAssertions(code$args[[i]]$args[[3]])
            }
            if(code$args[[i]]$name == '{') {
                exprClasses_insertAssertions(code$args[[i]])
            }
            if(code$args[[i]]$name == 'nimSwitch') {
                if(length(code$args[[i]]$args) > 2)
                    for(j in 3:length(code$args[[i]]$args)) {
                        exprClasses_insertAssertions(code$args[[i]]$args[[j]])
                    }
            }
            if(length(code$args[[i]]$assertions) > 0) {
                toInsert <- lapply(code$args[[i]]$assertions, function(x) if(inherits(x, 'exprClass')) x else RparseTree2ExprClasses(x))
                before <- unlist(lapply(toInsert, function(x) {if(x$name == 'after') FALSE else TRUE}))
                if(nimbleOptions('useRefactoredSizeProcessing')) code$args[[i]]$assertions <- list()  ## Clear assertions field so it can be used by later compiler stages
                newExpr <- newBracketExpr(args = c(lapply(toInsert[before], ## assertions will be inserted recursively IF they are in {}
                                                          function(z) {exprClasses_insertAssertions(z); z}),
                                                   code$args[i],
                                                   lapply(toInsert[!before],
                                                          function(x) {lapply(x$args[1],
                                                                              function(z) {exprClasses_insertAssertions(z)}
                                                                              ); x$args[[1]]}
                                                          )))
                setArg(code, i, newExpr)
            }
        }
    }
    invisible(NULL)
}
