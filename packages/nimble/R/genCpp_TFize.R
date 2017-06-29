TFrunnerLabelGenerator <- labelFunctionCreator("TFrunner")

exprClasses_TFize <-
    function(code, symTab, typeEnv, workEnv = new.env()) {
        ## This imitates exprClasses_eigenize
        ## It is possible they could later be combined and simply use different handler tables
        setupExprs <- list()
        if (code$isCall) {
            if (code$name == '{') {
                for (i in seq_along(code$args)) {
                    recurse <-
                        FALSE                         ## Decide whether to recurse
                    if (code$args[[i]]$name == 'eigenize') {
                        ## No "TFize" label yet
                        removeExprClassLayer(code$args[[i]]) ## strip the eigenize()
                        recurse <- TRUE
                    }
                    if (code$args[[i]]$name %in% c('for', ifOrWhile, '{', 'nimSwitch'))
                        recurse <- TRUE
                    if (recurse) {
                        ## return object collects setup calls...
                        setupCalls <-
                            unlist(
                                exprClasses_TFize(code$args[[i]], symTab, typeEnv, workEnv = new.env())
                            ) ## a new line
                        ## ... which are immediately inserted
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
            setupExprs <-
                c(setupExprs,
                  TFize_oneStatement(code, symTab, typeEnv, workEnv)) ## A single handler
            ## various further steps in exprClasses_eigenize may not be relevant - TBD.
        }
        return(if (length(setupExprs) == 0)
            NULL
            else
                setupExprs)
    }

## Manage tensorflow-ization of code (e.g. one statement)
TFize_oneStatement <- function(code, symTab, typeEnv, workEnv) {
    TFcontent <- exprClasses2serializedTF(code, symTab)
    ## This code-generation problem is different than others we've supported:
    ## We have a constructor that must happen as a constructor because it is for a
    ## static pointer.  But the cppVarFull::generate takes constructor as simple text
    ## and does not recursively generate constructor for exprClasses.  Hence
    ## for a first step now we will just paste together the constructor code
    TFconstructor <- makeTFconstruction(TFcontent)
    TFrunnerName <- TFrunnerLabelGenerator()
    TFrunnerSym <-
        symbolTensorflowRunner(name = TFrunnerName,
                               constructor = TFconstructor,
                               type = "symbolTensorflowRunner")
    symTab$addSymbol(TFrunnerSym)
    TFsetupExprs <- makeTFsetupExprs(TFcontent, TFrunnerName)
    ## convert original code line to a comment
    original <- nimDeparse(code)
    code$name <- "cppComment"
    code$args <-
        list(paste0("End: Tensorflow implementation for: ", original))
    TFsetupExprs
}

makeTFsetupExprs <- function(TFcontent, TFrunnerName) {
    ## This uses the prefix NimTf_ so that we can assume that internal
    ## symbols conflict with say a method written by a user during
    ## generateCpp step (which operates by names only, not smart lookup
    ## in classes etc).
    setupExprs <- list()
    for (v in TFcontent$inputNames) {
        setupExprs[[length(setupExprs) + 1]] <-
            RparseTree2ExprClasses(substitute(
                cppMemberFunction(NimTf_setInput(TFRN, V)),
                list(TFRN = as.name(TFrunnerName), V = as.name(v))
            ))
    }
    setupExprs[[length(setupExprs) + 1]] <-
        RparseTree2ExprClasses(substitute(cppMemberFunction(NimTf_run(TFRN)),
                                          list(TFRN = as.name(TFrunnerName))))
    for (v in TFcontent$outputNames) {
        setupExprs[[length(setupExprs) + 1]] <-
            RparseTree2ExprClasses(substitute(
                cppMemberFunction(NimTf_getOutput(TFRN, V)),
                list(TFRN = as.name(TFrunnerName), V = as.name(v))
            ))
    }
    setupExprs
}

makeTFconstruction <- function(TFcontent) {
    ## pasting code instead of creating a parse tree
    paste(
        paste0(' = *NimTf_Builder("', TFcontent$serializedGraph, '")'),
        paste0(
            paste0('.NimTf_withInput("', TFcontent$inputNames, '")'),
            collapse = "\n"
        ),
        paste0(
            paste0('.NimTf_withOutput("', TFcontent$outputNames, '")'),
            collapse = "\n"
        ),
        '.NimTf_build()',
        sep = "\n"
    )
}



exprClasses2serializedTF <- function(code, symTab) {
    ## This prototype assumes arguments are 'ARG1_a_','ARG1_x_', and 'ARG1_y_' and output is 'z'.
    ## TODO Determine intputs, outputs, and tensorflow graph from `code`.
    
    ## Construct a tensorflow graph.
    if (!require(tensorflow))
        stop('Failed to load tensorflow package')
    tf$reset_default_graph()
    a <-
        tf$placeholder(name = 'ARG1_a_',
                       dtype = tf$float64,
                       shape = list())
    x <-
        tf$placeholder(name = 'ARG2_x_',
                       dtype = tf$float64,
                       shape = list(NULL))
    y <-
        tf$placeholder(name = 'ARG3_y_',
                       dtype = tf$float64,
                       shape = list(NULL))
    z <- tf$add(tf$multiply(a, x), y, name = 'z')
    
    ## Serialize the graph as a string.
    if (!require(reticulate))
        stop('Failed to load reticulate package')
    base64 <- import('base64')
    graph <-
        base64$b64encode(z$graph$as_graph_def()$SerializeToString())
    
    ## Create invocation of NimTf_Builder().
    tfProxy <- tfProxyClass()
    for (var in c('ARG1_a_', 'ARG2_x_', 'ARG3_y_'))
        tfProxy$addInputVar(var)
    tfProxy$addOutputVar('z')
    tfProxy$setSerializedGraph(graph)
    
    return(tfProxy)
}

## This is a proxy for (or may evolved to contain) the tf graph
tfProxyClass <- setRefClass(
    Class = "tfProxyClass",
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
        }
        
    )
)
