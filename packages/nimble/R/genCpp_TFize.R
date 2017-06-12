TFrunnerLabelGenerator <- labelFunctionCreator("TFrunner")

exprClasses_TFize <- function(code, symTab, typeEnv, workEnv = new.env()) {
    ## This imitates exprClasses_eigenize
    ## It is possible they could later be combined and simply use different handler tables
    setupExprs <- list()
    if(code$isCall) {
        if(code$name == '{') {
            for(i in seq_along(code$args)) {
                recurse <- FALSE                         ## Decide whether to recurse
                if(code$args[[i]]$name == 'eigenize') {  ## No "TFize" label yet
                    removeExprClassLayer(code$args[[i]]) ## strip the eigenize()
                    recurse <- TRUE
                }
                if(code$args[[i]]$name %in% c('for', ifOrWhile, '{', 'nimSwitch')) recurse <- TRUE
                if(recurse) {
                    ## return object collects setup calls...
                    setupCalls <- unlist(exprClasses_TFize(code$args[[i]], symTab, typeEnv, workEnv = new.env())) ## a new line
                    ## ... which are immediately inserted 
                    if(length(setupCalls) > 0) {
                        newExpr <- newBracketExpr(args = c(setupCalls, code$args[[i]]))
                        setArg(code, i, newExpr)
                    }
                }
            }
            return(invisible(NULL))
        }
        if(code$name == 'for') {
            exprClasses_TFize(code$args[[3]], symTab, typeEnv); return(invisible(NULL))
        }
        if(code$name %in% ifOrWhile) {
            exprClasses_TFize(code$args[[2]], symTab, typeEnv)
            if(length(code$args)==3) exprClasses_TFize(code$args[[3]], symTab, typeEnv)
            return(invisible(NULL))
        }
        if(code$name == 'nimSwitch') {
            message('need to set up TFize for nimSwitch')
        }
        if(code$name == 'map') {
            message('need to set up TFize (or determine unnecessary) for map')
        }
        if(code$name == '[') {
            message('need to set up TFize (or determine unnecessary) for `[`')
        }
        TFize_oneStatement(code, symTab, typeEnv, workEnv) ## A single handler
        ## various further steps in exprClasses_eigenize may not be relevant - TBD.
    }
    return(if(length(setupExprs) == 0) NULL else setupExprs)
}

## Manage tensorflow-ization of code (e.g. one statement)
TFize_oneStatement <- function(code, symTab, typeEnv, workEnv) {
    TFcontent <- exprClasses2serializedTF(code, symTab)
    TFrunnerName <- TFrunnerLabelGenerator()
    TFrunnerSym <- symbolTensorFlowRunner()
}



exprClasses2serializedTF <- function(code, symTab) {
    ## 

    ##
    ans <- nimDeparse(code)
    ans
}

TFoperationNames = list(
    '*' = 'multiply',
    '+' = 'add'
)

## build up an object 
exprClasses2TFgraph <- function(code, symTab, tfProxy) {
    if(!inherits(code, 'exprClass')) return; ## could be a numeric constant
    if(code$isName) {
        tfProxy$addVar(code$name, code$type, code$nDim)
        return;
    }
    if(code$isCall) {
        TFop <- TFoperationNames[[code$name]]
        if(is.null(TFop)) stop(paste0('Function ', code$name, ' not supported for use with TF.'))
        names <- unique(unlist(lapply(code$args,
                                      function(x) {
                                          exprClasses2TFgraph(x, symTab, tfProxy)
                                          x$aux$names
                                      })))
    }
}

## This is a proxy for (or may evolved to contain) the tf graph
tfproxy <- setRefClass(
    fields = list(
        inputNames = 'ANY', ## character vector of arguments *representing canonical order*
        outputNames= 'ANY', ## initially, a single name
        serializedGraph = 'ANY'),
    methods = list(
        addVar = function(name, type, nDim, output = FALSE) { ## not sure if nDim is needed
            if(!(name %in% inputNames)) {
                inputNames <- append(inputNames, name)
                
            }
        }
    ))
