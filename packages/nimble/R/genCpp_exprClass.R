## exprClass is a class for parse trees (expressions) that enhances
## the information in the regular R parse tree.  It makes it so we can
## navigate the parse tree more easily, keeping track of what we need, and
## not being awkward with the R parse tree all the time.  In particular, it has fields
## flagging whether the expression is a name, a call, and/or an assignment;
## it has information on the types, number of dimensions, and size expressions, when appropriate;
## it knows what exprClass object it is an argument to, and what number argument it is;
## and it has a field eigMatrix used for tracking matrix vs. array status of Eigen expressions
##
## In places it is still convenient to use R parse trees, or to generate
## code in an R parse tree and then call RparseTree2ExprClasses to convert it
exprClass <- R6::R6Class(
    'exprClass',
    portable = FALSE,
    public = list(
        expr = NULL, ## optional for the original R expr. Not used much except at setup
        isName = NULL,		#'logical', ## is it a name
        isCall = NULL,		# 'logical', ## is it a call
        isAssign =  NULL,		#'logical', ## it is an assignment (all assignments are also calls)
        name =  NULL,		#'character', ## what is the name of the call or the object (e.g. 'a' or '+')
        nDim =  NULL,		#'numeric', ## how many dimensions
        sizeExprs =  list(),	#'list', ## a list of size expressions (using R parse trees for each non-numeric expression
        type =  NULL,		#'character', ## type label
        args =  list(),		#'list', ## list of exprClass objects for the arguments
        eigMatrix =  logical(),	#'logical', ## vs. Array.  Only used for Eigenized expressions
        toEigenize =  'unknown',	#'character', ##'yes', 'no', or 'maybe'
        caller = NULL,           # exprClass object for the call to which this is an argument (if any)
        callerArgID =  NULL,	#'numeric', ## index in the calling object's args list for this object.
        assertions =  list(),	#'list'
        cppADCode = FALSE,         #'logical' ## is expr in code generated for cppad?
        aux = NULL,               # anything needed for specific operators
        initialize = function(...) {
            dotsList <- list(...)
            for(v in names(dotsList))
                self[[v]] <- dotsList[[v]]
        },
        ## This displays the parse tree using indentation on multiple rows of output
        ## It also checks that the caller and callerArgID fields are all correct
        ## For deparsing, call nimDeparse
        print = function(indent = '', showType = FALSE, showAssertions = FALSE, showToEigenize = FALSE) {
            ## optionally include size information
            sizeDisp <- if(showType) { 
                if(is.null(type)) {
                    paste0(' \t| type = (NULL)(', length(sizeExprs), ', c(', paste0(unlist(lapply(sizeExprs, function(x) if(inherits(x, 'exprClass')) nimDeparse(x) else deparse(x))), collapse = ', '), '))')
                } else {
                    paste0(' \t| type = ', type,'(', length(sizeExprs), ', c(', paste0(unlist(lapply(sizeExprs, function(x) if(inherits(x, 'exprClass')) nimDeparse(x) else deparse(x))), collapse = ', '), '))')
                }
            } else
                character()
    
            assertDisp <- if(showAssertions & length(assertions) > 0) {
                paste(' \t| assert:',paste(unlist(lapply(assertions, deparse)), collapse = ', '))
            } else character()
    
            toEigenizeDisp <- if(showToEigenize) {
                paste(' \t| ', toEigenize)
            } else character()
            ## write self label
            writeLines(paste0(indent, name, sizeDisp, assertDisp, toEigenizeDisp))
            ## Iterate through arguments and show them with more indenting
            for(i in seq_along(args)) {
                if(inherits(args[[i]], 'exprClass')) {
                    args[[i]]$print(paste0(indent, '  '), showType, showAssertions, showToEigenize)
                    ## check caller and callerArgID validity
                    if(!identical(args[[i]]$caller, self) |
                       args[[i]]$callerArgID != i)
                        writeLines(paste0(indent, '****',
                                          'Warning: caller and/or callerID are not set correctly.'))
                } else {
                    writeLines(paste0(indent, '  ', if(is.null(args[[i]])) 'NULL' else args[[i]]))
                }
            }
            ## close brackets
            if(name=='{') writeLines(paste0(indent,'}'))
        }
    )
)

## This class is for keeping track of symbolic dimension and size information
## for a variable name, including any known changes as code processing progresses.
## This differs from the information in the symbol table in that the latter is static: it is not modified as a result of code processing.
exprTypeInfoClass <- setRefClass('exprTypeInfoClass',
    fields = list(
        nDim =  'ANY',		#'numeric',
        sizeExprs =  'ANY',	#'list',
        type =  'ANY'		#'character'
    ),
    methods = list(
    	initialize = function(...) {
            sizeExprs <<- list()
            callSuper(...)
        },
        show = function() {
            writeLines(paste0('exprTypeInfoClass: nDim = ', nDim, '. type = ', type,'. sizeExprs = ', paste0( lapply(sizeExprs, deparse), collapse = ',')))
        }
    )
)

copyExprClass <- function(original) {
    result <- original$clone(deep = FALSE)
    ## shallow=FALSE does not deep-copy on list elements, so it is
    ## useless for args list.  Another reason for shallow = TRUE
    ## is we do not want to deep copy 'caller' here.  Instead it is
    ## re-assigned below.
    for(i in seq_along(result$args)) {
        if(inherits(result$args[[i]], 'exprClass')) {
            result$args[[i]] <- copyExprClass(result$args[[i]])
            result$args[[i]]$caller <- result
        }
    }
    result
}

## Add indendation to every character() element in a list, doing so recursively for any nested list()
addIndentToList <- function(x, indent) {
    if(is.list(x)) return(lapply(x, addIndentToList, indent))
    paste0(indent, x)
}

### Deparse from exprClass back to R code: not guaranteed to be identical, but valid.
nimDeparse <- function(code, indent = '') {
    ## numeric case
    if(is.numeric(code) | is.logical(code)) return(code)
    if(is.character(code)) return(paste0('\"', code, '\"'))
    if(is.null(code)) return('NULL')
    ## name
    if(!inherits(code, 'exprClass')) return(class(code))
    if(code$isName) return(code$name)
    ## { : iterate through contents, deparsing and adding indentation
    if(code$name == '{') {
        ans <- vector('list', length(code$args) + 2)
        ans[[1]] <- paste0(indent, '{')
        for(i in seq_along(code$args)) {
            newInd <- paste0(indent, '  ')
            ans[[i+1]] <- addIndentToList(nimDeparse(code$args[[i]], newInd), newInd)
        }
        ans[[length(code$args) + 2]] <- paste0(indent, '}')
        return(ans)
    }
    ## for: 
    if(code$name == 'for') {
        part1 <- paste0('for(', nimDeparse(code$args[[1]]), ' in ', nimDeparse(code$args[[2]]), ')')
        part2 <- nimDeparse(code$args[[3]]) ## body of loop
        if(is.list(part2)) { ## should be TRUE if the created by RparseTree2ExprClasses, which ensures that for bodys are embedded in {
            ## put the { on the same line as the for()
            part2[[1]] <- paste(part1, part2[[1]])
            part2 <- addIndentToList(part2, indent)
            return(part2)
        } else {
            ## If the body is not in a {, put it on the same line as for()
            return(paste(part1, part2))
        }
    }
    ## if:
    if(code$name %in% ifOrWhile) {
        part1 <- paste0('if(', nimDeparse(code$args[[1]]),')')
        ## body of 'then' or 'while' clause:
        part2 <- nimDeparse(code$args[[2]])
        if(is.list(part2)) { ## If it is in a {, put { on same line as if()
            part2[[1]] <- paste(part1, part2[[1]])
            part2 <- addIndentToList(part2, indent)
        } else {
            ## If it is not in a {, put on same line as if()
            part2 <- paste(part1, part2)
        }
        ## If there is no 'else' clause, return
        if(length(code$args) == 2) return(part2) ## no else clause
        ## body of 'else' clause:
        part3 <- nimDeparse(code$args[[3]])
        if(is.list(part3)) { ## If it is in {, merge lines to get '} else {'
            part3[[1]] <- paste(part2[[length(part2)]], 'else', part3[[1]])
            part3[-1] <- addIndentToList(part3[-1], indent)
            part3 <- c(part2[-length(part2)], part3)
            return(part3)
        } else {
            ## If it is not in {, merge lines to get '} else else-body'
            part2[[length(part2)]] <- paste(part2[[length(part2)]], 'else', part3)
            return(part2)
        }
    }
    ## for operators like '+' that go after the first argument:
    if(code$name %in% names(midOperators)) {
        if(length(code$args) == 2) {
            if(is.null(code$caller)) useParens <- FALSE
            else {
                thisRank <- operatorRank[[code$name]]
                callingRank <- if(!is.null(code$caller)) operatorRank[[code$caller$name]] else NULL
                useParens <- FALSE ## default to FALSE
                if(!is.null(callingRank)) {
                    if(!is.null(thisRank)) {
                        if(callingRank < thisRank) useParens <- TRUE
                    }
                }
            }
            if(useParens)
                return( paste0('(', nimDeparse(code$args[[1]]), midOperators[[code$name]], paste0(unlist(lapply(code$args[-1], nimDeparse) ), collapse = ', '), ')' ) )
            else
                return( paste0(nimDeparse(code$args[[1]]), midOperators[[code$name]], paste0(unlist(lapply(code$args[-1], nimDeparse) ), collapse = ', ') ) )
        ## unary case will go to end, so -A will become -(A)
        }
    }
    ## for operators like '[' that go after the first argument and then have a closing part like ']'
    if(code$name %in% names(brackOperators)) {
        return( paste0(nimDeparse(code$args[[1]]), brackOperators[[code$name]][1], paste0(unlist(lapply(code$args[-1], nimDeparse) ), collapse = ', '), brackOperators[[code$name]][2] ) )
    }
    if(code$name == '(') {
        return(paste0('(', paste0(unlist(lapply(code$args, nimDeparse) ), collapse = ', '), ')'))
    }
    if(code$name == 'chainedCall') {
        return(paste0(nimDeparse(code$args[[1]]),
                      '(', paste0(unlist(lapply(code$args[-1], nimDeparse) ), collapse = ', '), ')' ) )
    }
    ## for a general function call. Modified to preseve any argument names
    deparsedArguments <- lapply(code$args, nimDeparse)
    argumentText <- if(!is.null(names(deparsedArguments))) paste(names(deparsedArguments), deparsedArguments, sep = '=') else unlist(deparsedArguments)
    return(paste0(code$name, '(', paste0(argumentText, collapse = ','), ')' ))
}

## error trapping utilities to be used from the various processing steps
exprClassProcessingErrorMsg <- function(code, msg) {
    contextCode <- if(!is.null(code$caller)) paste(unlist(nimDeparse(code$caller)), collapse = '\n') else character()
    ans <- paste0(msg, '\n This occurred for: ', nimDeparse(code),'\n', collapse = '')
    if(!is.null(contextCode)) ans <- paste(ans, 'This was part of the call: ', contextCode, collapse = '')
    ans
}

## some utilities

## return a list of size expressions with any "1"s dropped
## E.g. (dim(A)[1], 1, 2) returns (dim(A)[1], 2)
## This would be useful for determining e.g. if (dim(B)[1], 2) and (dim(A)[1], 1, 2) can be added as like types
## For now, the only case we use it for is determining is something can be treated like a scalar, e.g. with sizes of 1 in every dimension
dropSingleSizes <- function(sizeList) {
    drop <- logical(length(sizeList))
    for(i in seq_along(sizeList)) {
        drop[i] <- if(is.numeric(sizeList[[i]])) sizeList[[i]] == 1 else FALSE 
    }
    sizeList[drop] <- NULL
    list(sizeExprs = sizeList, drop = drop)
}

insertExprClassLayer <- function(code, argID, funName, isName = FALSE, isCall = TRUE, isAssign = FALSE, ... ) {
    newExpr <- exprClass$new(name = funName, isName = isName, isCall = isCall, isAssign = isAssign,
                             args = list(code$args[[argID]]), ...)
    setCaller(code$args[[argID]], newExpr, 1)
    setArg(code, argID, newExpr)
    newExpr
}

insertIndexingBracket <- function(code, argID, index) {
    insertExprClassLayer(code, argID, '[')
    setArg(code$args[[argID]], 2, index)
}

## make a new bracket { expression in exprClass objects
newBracketExpr <- function(args = list()) {
    ans <- exprClass$new(isCall = TRUE, isName = FALSE, isAssign = FALSE,
                  name = '{', args = args)
    for(i in seq_along(args)) {
        setCaller(ans$args[[i]], ans, i)
    }
    ans
}

removeExprClassLayer <- function(code, argID = 1) {
    setArg(code$caller, code$callerArgID, if(length(code$args) >= argID) code$args[[argID]] else NULL)
}

setCaller <- function(value, expr, ID) {
    value$caller <- expr
    value$callerArgID <- ID
    invisible(value)
}

setArg <- function(expr, ID, value) {
    expr$args[[ID]] <- value
    if(inherits(value, 'exprClass')) setCaller(value, expr, ID)
    invisible(value)
}

newAssignmentExpression <- function() {
    exprClass$new(isName = FALSE, isCall = TRUE, isAssign = TRUE, name = '<-')
}

## This modifies the code$caller in place
## and generates the temp expr
buildSimpleIntermCall <- function(code) {
    if(code$caller$isAssign) return(NULL)
    newName <- IntermLabelMaker()

    ## change my argument from the caller to the new temp
    setArg(code$caller, code$callerArgID, RparseTree2ExprClasses(as.name(newName)))
    
    newExpr <- newAssignmentExpression()
    setArg(newExpr, 1, RparseTree2ExprClasses(as.name(newName))) 
    setArg(newExpr, 2, code) ## The setArg function should set code$caller (to newExpr) and code$callerArgID (to 3)

    return(newExpr)
}

isCodeScalar <- function(code) {
    for(i in seq_along(code$sizeExprs))
        if(!identical(code$sizeExprs[[i]], 1)) return(FALSE)
    TRUE
}

anyNonScalar <- function(code) {
    if(!inherits(code, 'exprClass')) return(FALSE)
    if(code$name == 'map') return(TRUE)
    if(is.character(code$type))
        if(code$type[1] == 'nimbleList') return(FALSE)
    if(code$isName) {
        return(!isCodeScalar(code))
    }
    if(code$isCall) {
        if(code$name == 'nfVar') ## don't recurse just for nested member access
            return(!isCodeScalar(code))
        skipFirst <- FALSE
        if(code$name == '[') skipFirst <- TRUE
        if(code$name == 'size') skipFirst <- TRUE
        for(i in seq_along(code$args)) {
            if(!(skipFirst & i == 1)) {
                if(anyNonScalar(code$args[[i]])) return(TRUE)
            }
        }
    }
    FALSE
}
