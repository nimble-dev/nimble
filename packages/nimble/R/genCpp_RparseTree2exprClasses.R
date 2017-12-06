isRbracket <- function(code) {
    if(!is.call(code)) return(FALSE)
    return(code[[1]] == '{')
}

## put a single exprObject within a '{'
embedInRbracket <- function(code) {
    template <- quote({A})
    template[[2]] <- code
    template
}

embedListInRbracket <- function(code) {
    if(!is.list(code)) stop('Error: embedListInRbracket called with code that is not a list')
    as.call(c(list(as.name('{')), code))
}


## build exprClasses from an R parse tree.
## caller and callerArgID are for recursion, not to be used on first entry
RparseTree2ExprClasses <- function(code, caller = NULL, callerArgID = numeric()) { ## input code is R parse tree
    ## name:
    if(is.name(code)) return(exprClass$new(expr = code, isName = TRUE, isCall = FALSE, isAssign = FALSE, name = as.character(code), caller = caller, callerArgID = callerArgID))
    ## call
    if(is.call(code)) {
        if(is.call(code[[1]])) { ## chained calls like a(b)(c) or a[[b]](c).  We wrap these as chainedCall(a(b), c) or chainedCall(a[[b]], c)
            code <- as.call(c(list(as.name('chainedCall')), as.list(code))) 
        }
        name <- as.character(code[[1]])
        isAssign <- name %in% c('<-','=','<<-')
        args <- vector('list', length = length(code)-1)
        ## build the object
        ans <- exprClass$new(expr = code, isName = FALSE, isCall = TRUE, isAssign = isAssign, name = name, args = args, caller = caller, callerArgID = callerArgID)

        ## Ensure that bodies of for and if are in { expressions.  Makes for less special-case checking in later processing
        if(name == 'for') {
            if(!isRbracket(code[[4]])) code[[4]] <- embedInRbracket(code[[4]])
        }
        if(name %in% ifOrWhile) {
            ## 'then' or 'while' clause
            if(!isRbracket(code[[3]])) code[[3]] <- embedInRbracket(code[[3]])
            ## 'else' clause
            if(length(code)==4) {
                if(!isRbracket(code[[4]])) code[[4]] <- embedInRbracket(code[[4]])
            }
        }
        if(name == 'nimSwitch') {
            if(length(code) > 3)
                for(iSwitch in 4:length(code))
                    if(!isRbracket(code[[iSwitch]])) code[[iSwitch]] <- embedInRbracket(code[[iSwitch]])
        }
        if(name == 'run.time') {
            if(!isRbracket(code[[2]])) code[[2]] <- embedInRbracket(code[[2]])
        }
        if(name == 'map') { ## special treatment. just stick the remaining arguments in as a list
            ans$args <- as.list(code[-1])
            return(ans)
        }
        
        ## populate args with recursive calls
        if(length(code) > 1) {
            if(!is.null(names(code))) names(ans$args) <- names(code)[-1] ## entries like "" on the RHS leave no name on LHS --> good.
            for(i in 2:length(code)) ## Note for NULL this removes the list entry.  Not very general, but handles return(invisible(NULL))
                if(is.logical(code[[i]])) {
                    if(name == '[' & i == length(code))
                        ans$args[[i-1]] <- code[[i]]
                    else ## cast logical to numeric unless it is last arg of a [, in which case it could be for drop.  This is not a very logical place for this step, but it works.
                        ans$args[[i-1]] <- code[[i]] ## try keeping logicals instead of casting to numeric as.numeric(code[[i]])
                } else {
                    ans$args[[i-1]] <- if(is.numeric(code[[i]]) | is.character(code[[i]]) | is.null(code[[i]])) code[[i]] else RparseTree2ExprClasses(code[[i]], caller = ans, callerArgID = i-1)
                }
        }
        ans
    }
}
