TFrunnerLabelGenerator <- labelFunctionCreator("TFrunner")

## This imitates exprClasses_eigenize.
## It is possible they could later be combined and simply use different handler tables.
exprClasses_TFize <-
    function(code, symTab, typeEnv, workEnv = new.env()) {
        setupExprs <- list()
        if (code$isCall) {
            if (code$name == '{') {
                for (i in seq_along(code$args)) {
                    ## Decide whether to recurse
                    recurse <- FALSE
                    if (code$args[[i]]$name == 'eigenize') {
                        ## No "TFize" label yet.
                        removeExprClassLayer(code$args[[i]]) ## Strip the eigenize().
                        recurse <- TRUE
                    }
                    if (code$args[[i]]$name %in% c('for', ifOrWhile, '{', 'nimSwitch'))
                        recurse <- TRUE
                    if (recurse) {
                        ## Return object collects setup calls...
                        setupCalls <-
                            unlist(
                                exprClasses_TFize(code$args[[i]], symTab, typeEnv, workEnv = new.env())
                            ) ## a new line
                        ## ... which are immediately inserted.
                        if (length(setupCalls) > 0) {
                            newExpr <- newBracketExpr(args = c(setupCalls, code$args[[i]]))
                            setArg(code, i, newExpr)
                        }
                    }
                }
                return(invisible(NULL))
            }
            if (code$name == 'for') {
                exprClasses_TFize(code$args[[3]], symTab, typeEnv)
                return(invisible(NULL))
            }
            if (code$name %in% ifOrWhile) {
                exprClasses_TFize(code$args[[2]], symTab, typeEnv)
                if (length(code$args) == 3)
                    exprClasses_TFize(code$args[[3]], symTab, typeEnv)
                return(invisible(NULL))
            }
            if (code$name == 'nimSwitch') {
                message('need to set up TFize for nimSwitch')
            }
            if (code$name == 'map') {
                message('need to set up TFize (or determine unnecessary) for map')
            }
            if (code$name == '[') {
                message('need to set up TFize (or determine unnecessary) for `[`')
            }
            ## All other cases are handled by a single handler.
            setupExprs <-
                c(setupExprs,
                  TFize_oneStatement(code, symTab, typeEnv, workEnv))
            ## Various further steps in exprClasses_eigenize may not be relevant - TBD.
        }
        return(if (length(setupExprs) == 0)
            NULL
            else
                setupExprs)
    }

## Manage tensorflow-ization of code (e.g. one statement).
TFize_oneStatement <- function(code, symTab, typeEnv, workEnv) {
    TfBuilder <- exprClasses2serializedTF(code, symTab)
    ## This code-generation problem is different than others we've supported:
    ## We have a constructor that must happen as a constructor because it is for a
    ## static pointer.  But the cppVarFull::generate takes constructor as simple text
    ## and does not recursively generate constructor for exprClasses.  Hence
    ## for a first step now we will just paste together the constructor code
    TFconstructor <- TfBuilder$generateConstructor()
    TFrunnerName <- TFrunnerLabelGenerator()
    TFrunnerSym <-
        symbolTensorflowRunner(name = TFrunnerName,
                               constructor = TFconstructor,
                               type = "symbolTensorflowRunner")
    symTab$addSymbol(TFrunnerSym)
    TFsetupExprs <- makeTFsetupExprs(TfBuilder, TFrunnerName)
    ## Convert original code line to a comment.
    original <- nimDeparse(code)
    code$name <- "cppComment"
    code$args <-
        list(paste0("End: Tensorflow implementation for: ", original))
    TFsetupExprs
}

makeTFsetupExprs <- function(TfBuilder, TFrunnerName) {
    ## This uses the prefix NimTf_ so that we can assume that internal
    ## symbols will not conflict with, say, a method written by a user during
    ## generateCpp step (which operates by names only, not smart lookup
    ## in classes etc).
    setupExprs <- list()
    for (v in TfBuilder$inputNames) {
        setupExprs[[length(setupExprs) + 1]] <-
            RparseTree2ExprClasses(substitute(
                cppMemberFunction(NimTf_setInput(TFRN, V)),
                list(TFRN = as.name(TFrunnerName), V = as.name(v))
            ))
    }
    setupExprs[[length(setupExprs) + 1]] <-
        RparseTree2ExprClasses(substitute(cppMemberFunction(NimTf_run(TFRN)),
                                          list(TFRN = as.name(TFrunnerName))))
    for (v in TfBuilder$outputNames) {
        setupExprs[[length(setupExprs) + 1]] <-
            RparseTree2ExprClasses(substitute(
                cppMemberFunction(NimTf_getOutput(TFRN, V)),
                list(TFRN = as.name(TFrunnerName), V = as.name(v))
            ))
    }
    setupExprs
}

## Recursively collects all names in an exprClass and create tf$placeholders for each.
TfCollectPlaceholders <- function(code, symTab, placeholders = NULL) {
    if (is.null(placeholders)) {
        placeholders = new.env()
    }
    if (code$isName) {
        if (is.null(placeholders[[code$name]])) {
            nimType2TfDtype = list('double' = tf$float64)
            sym <- symTab$getSymbolObject(code$name)
            dtype <- nimType2TfDtype[[sym$type]]
            size <- as.list(rev(sym$size))  ## Note the transpose
            for (i in 1:length(size)) {
                if (is.na(size[i])) {
                    size[i] <- NULL
                }
            }
            placeholders[[code$name]] = tf$placeholder(name = code$name,
                                                       dtype = dtype,
                                                       shape = size)
        }
        return(placeholders)
    }
    if (code$isAssign) {
        return(TfCollectPlaceholders(code$args[[2]], symTab, placeholders))
    }
    for (arg in code$args) {
        if (class(arg) == 'exprClass') {
            TfCollectPlaceholders(arg, symTab, placeholders)
        }
    }
    return(placeholders)
}

TfTensorize <- function(code, placeholders) {
    if (!is.environment(placeholders)) stop()

    if (code$isName) {
        return(placeholders[[code$name]])
    }
    if (code$name %in% c('+', '-', '*', '/')) {
        args <- lapply(code$args, TfTensorize, placeholders)
        return(do.call(code$name, args))
    }
    if (code$name == '%*%') {
        lhs <- TfTensorize(code$args[[1]], placeholders)
        rhs <- TfTensorize(code$args[[2]], placeholders)
        return(tf$matmul(lhs, rhs))
    }
    stop(paste('Not implemented:', code))
}

exprClasses2serializedTF <- function(code, symTab) {
    if (!require(tensorflow)) {
        stop('Failed to load tensorflow package')
    }
    if (!require(reticulate)) {
        stop('Failed to load reticulate package')
    }
    
    if (code$name != '<-') {
        stop(paste('Not implemented:', code$name))
    }
    target <- code$args[[1]]
    expr <- code$args[[2]]

    ## Construct a tensorflow graph.
    tf$reset_default_graph()
    placeholders <- TfCollectPlaceholders(expr, symTab)
    tensor <- TfTensorize(expr, placeholders)
    tensor <- tf$identity(tensor, name = target$name)

    ## Serialize the graph as a string.
    base64 <- reticulate::import('base64')
    graph <- base64$b64encode(tensor$graph$as_graph_def()$SerializeToString())

    ## Create NimTf_Builder.
    tfBuilder <- TfBuilder()
    tfBuilder$setSerializedGraph(graph)
    for (name in names(placeholders)) {
        tfBuilder$addInputVar(name)
    }
    tfBuilder$addOutputVar(target$name)
    return(tfBuilder)
}

## This is a proxy for (or may evolved to contain) the tf graph
TfBuilder <- setRefClass(
    Class = "TfBuilder",
    fields = list(
        inputNames = 'ANY',
        ## character vector of arguments *representing canonical order*
        outputNames = 'ANY',
        ## initially, a single name
        serializedGraph = 'ANY'
    ),
    methods = list(
        initialize = function() {
            inputNames <<- outputNames <<- character()
            serializedGraph <<- character()
        },
        addInputVar = function(name, type = 'double', output = FALSE) {
            ## not sure if nDim is needed
            if (!(name %in% inputNames)) {
                inputNames <<- append(inputNames, name)
            }
        },
        addOutputVar = function(name, type = 'double', output = FALSE) {
            ## not sure if nDim is needed
            if (!(name %in% outputNames)) {
                outputNames <<- append(outputNames, name)
            }
        },
        setSerializedGraph = function(graphText) {
            serializedGraph <<- graphText
        },
        generateConstructor = function() {
            ## Paste code, instead of creating a parse tree.
            paste(
                paste0(' = *NimTf_Builder("', serializedGraph, '")'),
                paste0(
                    paste0('.NimTf_withInput("', inputNames, '")'),
                    collapse = "\n"
                ),
                paste0(
                    paste0('.NimTf_withOutput("', outputNames, '")'),
                    collapse = "\n"
                ),
                '.NimTf_build()',
                sep = "\n"
            )
        }
    )
)
