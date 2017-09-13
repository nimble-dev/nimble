
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
            handleWhile <- FALSE
            if(code$args[[i]]$name %in% ifOrWhile) {
                exprClasses_insertAssertions(code$args[[i]]$args[[2]])
                if(length(code$args[[i]]$args) == 3) exprClasses_insertAssertions(code$args[[i]]$args[[3]])
                handleWhile <- code$args[[i]]$name == 'while'
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
            beforeContents <- afterContents <- list()
            if(length(code$args[[i]]$assertions) > 0) {
                toInsert <- lapply(code$args[[i]]$assertions, function(x) if(inherits(x, 'exprClass')) x else RparseTree2ExprClasses(x))
                before <- unlist(lapply(toInsert, function(x) {if(x$name == 'after') FALSE else TRUE}))
                if(nimbleOptions('experimentalNewSizeProcessing')) code$args[[i]]$assertions <- list()  ## Clear assertions field so it can be used by later compiler stages
                beforeContents <- lapply(toInsert[before], ## assertions will be inserted recursively IF they are in {}
                                         function(z) {
                                             exprClasses_insertAssertions(z)
                                             z
                                         })
                afterContents <- lapply(toInsert[!before],
                                        function(x) {
                                            lapply(x$args[1],
                                                   function(z) {
                                                       exprClasses_insertAssertions(z)
                                                   }
                                                   )
                                            x$args[[1]]
                                        }
                                        )
                newExpr <- newBracketExpr(args = c(beforeContents,
                                                   code$args[i],
                                                   afterContents))
                setArg(code, i, newExpr)
                newExpr <- NULL
            }
            ## For 'if' or 'while', recurse into the if/while clause and (possibly) else clause.
            ## 'while' causes a special challenge.  Any assertions must go before the condition
            ## and also at the end of the clause for execution.  Otherwise the assertions will be
            ## done only once on the first time through.
            if(handleWhile) {
                if(length(afterContents) > 0)
                    stop('Invalid conditional in while()')
                if(length(beforeContents) > 0) {
                    ## The "while" is now embedded in newExpr, so it is in a different place in the AST.
                    newWhileBody <- newBracketExpr(
                        args = c(code$
                                 args[[i]]$ ## entry in original {
                                 args[[length(beforeContents)+1]]$ ## entry in newExpr {
                                 args[2],  ## body of while clause, as list
                                 lapply(beforeContents, copyExprClass)) ## new piece
                    )
                    setArg(
                        code$args[[i]]$args[[length(beforeContents)+1]],
                        2,
                        newWhileBody)
                }
            }
        }
    }
    invisible(NULL)
}
