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
    TFconstructor <- TfBuilder$generateConstructor(code$cppADCode)
    TFrunnerName <- TFrunnerLabelGenerator()
    if (code$cppADCode) {
        TFrunnerSym <-
            symbolTensorflowOp(name = TFrunnerName,
                               constructor = TFconstructor,
                               type = "symbolTensorflowOp")
    } else {
        TFrunnerSym <-
            symbolTensorflowRunner(name = TFrunnerName,
                                   constructor = TFconstructor,
                                   type = "symbolTensorflowRunner")
    }
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
## Returns an environment mapping placeholder names to tensorflow tensors.
TfCollectPlaceholders <- function(code, symTab, placeholders = NULL) {
    if (is.null(placeholders)) {
        placeholders = new.env()
    }
    if (code$isName) {
        if (is.null(placeholders[[code$name]])) {
            nimType2TfDtype = list('double' = tensorflow::tf$float64)
            sym <- symTab$getSymbolObject(code$name)
            if (is.null(sym)) {
                ## TODO this case may be unnecessar.
                ## We make a best guess to push past a bug.
                dtype <- tensorflow::tf$float64
                size <- list()
            } else {
                dtype <- nimType2TfDtype[[sym$type]]
                if (class(sym$size) == 'uninitializedField') {
                    size <- as.list(rep(NA, sym$nDim))
                    for (i in 1:sym$nDim) {
                        size[i] <- list(NULL)
                    }
                } else {
                    size <- as.list(rev(sym$size))  ## Note the transpose
                    for (i in 1:length(size)) {
                        if (is.na(size[i])) {
                            size[i] <- list(NULL)
                        }
                    }
                }
            }
            placeholders[[code$name]] = tensorflow::tf$placeholder(name = code$name,
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

## This constructs the tfTranslate list lazily,
## since the tensorflow package may not be loaded.
.tfLazyData <- new.env()
tfTranslate <- function(name) {
    developing <- FALSE  ## Switch to TRUE while developing.
    if (!developing && !is.null(.tfLazyData$tfTranslate)) {
        return(.tfLazyData$tfTranslate[[name]])
    }
    ## This translates R's functions to tensorflow functions, which are
    ## implemented in Python and wrapped by the reticulate package.
    ## For documentation, see the tensorflow python api docs at:
    ## https://www.tensorflow.org/api_docs/python/tf
    .tfLazyData$tfTranslate <- list(
        '+' = '+',
        '-' = '-',
        '*' = '*',
        '/' = '/',
        'abs' = tensorflow::tf$abs,
        'ceil' = tensorflow::tf$ceil,
        'floor' = tensorflow::tf$floor,
        'round' = tensorflow::tf$round,
        'nimRound' = tensorflow::tf$round,
        'trunc' = function(x) {
            if (x$dtype == tensorflow::tf$int64) {
                tensorflow::tf$truncatediv(x, tensorflow::tf$constant(1, tensorflow::tf$int64))
            } else {
                tensorflow::tf$cast(tensorflow::tf$cast(x, tensorflow::tf$int64), tensorflow::tf$float64)
            }
        },
        'ftrunc' = function(x) tensorflow::tf$cast(tensorflow::tf$cast(x, tensorflow::tf$int64), tensorflow::tf$float64),
        'cos' = tensorflow::tf$cos,
        'sin' = tensorflow::tf$sin,
        'tan' = tensorflow::tf$tan,
        'acos' = tensorflow::tf$acos,
        'asin' = tensorflow::tf$asin,
        'atan' = tensorflow::tf$atan,
        'cosh' = tensorflow::tf$cosh,
        'sinh' = tensorflow::tf$sinh,
        'tanh' = tensorflow::tf$tanh,
        'acosh' = tensorflow::tf$acosh,
        'asinh' = tensorflow::tf$asinh,
        'atanh' = tensorflow::tf$atanh,
        'exp' = tensorflow::tf$exp,
        'log' = tensorflow::tf$log,
        'log1p' = tensorflow::tf$log1p,
        'pow' = tensorflow::tf$pow,
        'sqrt' = tensorflow::tf$sqrt,
        'square' = tensorflow::tf$square,
        'cube' = function(x) tensorflow::tf$pow(x, 3L),
        'gamma' = function(x) tensorflow::tf$exp(tensorflow::tf$lgamma(x)),
        'gammafn' = function(x) tensorflow::tf$exp(tensorflow::tf$lgamma(x)),
        'lgamma' = tensorflow::tf$lgamma,
        'lgammafn' = tensorflow::tf$lgamma,
        'loggam' = tensorflow::tf$lgamma,
        'factorial' = function(x) tensorflow::tf$exp(tensorflow::tf$lgamma(tensorflow::tf$constant(1, tensorflow::tf$float64) + x)),
        'lfactorial' = function(x) tensorflow::tf$lgamma(tensorflow::tf$constant(1, tensorflow::tf$float64) + x),
        'cloglog' = function(x) tensorflow::tf$log(-tensorflow::tf$log1p(-x)),
        'icloglog' = function(x) tensorflow::tf$constant(1, tensorflow::tf$float64) - tensorflow::tf$exp(-tensorflow::tf$exp(x)),
        'logit' = function(x) tensorflow::tf$log(x / (tensorflow::tf$constant(1, tensorflow::tf$float64) - x)),
        'ilogit' = tensorflow::tf$sigmoid,
        'expit' = tensorflow::tf$sigmoid,
        'probit' = function(x) {
            one <- tensorflow::tf$constant(1.0, tensorflow::tf$float64)
            two <- tensorflow::tf$constant(2.0, tensorflow::tf$float64)
            ## tensorflow::tf$erfinv is implemented in C++ but not yet exposed in python as of tensorflow 1.3.0.
            tensorflow::tf$sqrt(two) * tensorflow::tf$erfinv(two * x - one)
        },
        'iprobit' = function(x) {
            one <- tensorflow::tf$constant(1.0, tensorflow::tf$float64)
            half <- tensorflow::tf$constant(0.5, tensorflow::tf$float64)
            (one + tensorflow::tf$erf(x * tensorflow::tf$sqrt(half))) * half
        },
        'inprod' = function(lhs, rhs) tensorflow::tf$reduce_sum(lhs * rhs),
        'min' = tensorflow::tf$reduce_min,
        'max' = tensorflow::tf$reduce_max,
        'pmin' = tensorflow::tf$minimum,
        'pmax' = tensorflow::tf$maximum,
        't' = function(mat) tensorflow::tf$transpose(mat, c(1L, 0L)),
        'asCol' = function(x) x,  ## TODO Does this suffice?
        'diagonal' = tensorflow::tf$matrix_diag,  ## TODO Decide between diag() and diag_part().
        'det' = tensorflow::tf$matrix_determinant,
        'logdet' = function(x) tensorflow::tf$log(tensorflow::tf$matrix_determinant(x)),
        'inverse' = tensorflow::tf$matrix_inverse,
        'sum' = tensorflow::tf$reduce_sum,
        'prod' = tensorflow::tf$reduce_prod,
        'mean' = tensorflow::tf$reduce_mean,
        'var' = function(x) {
            n_minus_one <- tensorflow::tf$cast(tensorflow::tf$size(x), tensorflow::tf$float64) - tensorflow::tf$constant(1, tensorflow::tf$float64)
            tensorflow::tf$reduce_sum(tensorflow::tf$square(x - tensorflow::tf$reduce_mean(x))) / n_minus_one
        },
        'sd' = function(x) {
            n_minus_one <- tensorflow::tf$cast(tensorflow::tf$size(x), tensorflow::tf$float64) - tensorflow::tf$constant(1, tensorflow::tf$float64)
            tensorflow::tf$norm(x - tensorflow::tf$reduce_mean(x), ord = 2) / tensorflow::tf$sqrt(n_minus_one)
        },
        '%*%' = function(x, y) {
            while (length(x$shape$as_list()) < 2) {
                x <- tensorflow::tf$expand_dims(x, 1L)
            }
            while (length(y$shape$as_list()) < 2) {
                y <- tensorflow::tf$expand_dims(y, 0L)
            }
            tensorflow::tf$matmul(y, x)  ## Note the transpose.
        },
        'chol' = tensorflow::tf$cholesky,
        'backsolve' = function(r, x) {
            ## TODO Decide whether to transpose.
            tensorflow::tf$matrix_triangular_solve(r, x, lower = TRUE)
        },
        'forwardsolve' = function(l, x) {
            ## TODO Decide whether to transpose.
            tensorflow::tf$matrix_triangular_solve(l, x, lower = FALSE)
        },
        ## These may be unnecessary. They were added to push past a possible bug
        ## in AD processing of tensorflow graphs.
        'nimCd' = tensorflow::tf$stack,
        'concatenateTemp' = function(x) x,
        'eigenBlock' = function(x) x
    )

    return(.tfLazyData$tfTranslate[[name]])
}

## Translates code from exprClasses format to tensorflow graph format.
## Returns a tensorflow Tensor that is the final result of the code.
TfTensorizeExpr <- function(code, placeholders) {
    if (!is.environment(placeholders)) stop()
    if (class(code) == 'exprClass') {
        while (code$name == '(') {
            code <- code$args[[1]]
        }
        if (code$isName) {
            return(placeholders[[code$name]])
        }
        translated <- tfTranslate(code$name)
        if (!is.null(translated)) {
            args <- lapply(code$args, TfTensorizeExpr, placeholders)
            return(do.call(translated, args))
        }
        stop(paste('Not implemented:', code$name))
    }
    if (class(code) == 'numeric') {
        return(tensorflow::tf$constant(code, dtype = tensorflow::tf$float64))
    }
    if (class(code) == 'integer') {
        return(tensorflow::tf$constant(code, dtype = tensorflow::tf$int64))
    }
    if (class(code) == 'logical') {
        return(tensorflow::tf$constant(code, dtype = tensorflow::tf$bool))
    }
    stop(paste('Not implemented:', class(code)))
}

## Creates a TfBuilder object representing the given code.
## Also computes gradients wrt each input, if possible.
## Setting threads=0 lets tensorflow choose the number of threads.
exprClasses2serializedTF <- function(code, symTab, threads = 0L) {
    if (!requireNamespace('tensorflow', 'quietly' = TRUE)) {
        stop('Failed to load tensorflow package')
    }
    if (!requireNamespace('reticulate', 'quietly' = TRUE)) {
        stop('Failed to load reticulate package')
    }
    
    if (code$name != '<-') {
        stop(paste('Not implemented:', code$name))
    }
    target <- code$args[[1]]
    expr <- code$args[[2]]

    ## Construct a tensorflow graph.
    tensorflow::tf$reset_default_graph()
    placeholders <- TfCollectPlaceholders(expr, symTab)
    tensor <- TfTensorizeExpr(expr, placeholders)
    tensor <- tensorflow::tf$identity(tensor, name = target$name)

    ## Optionally add a gradient for each input arg.
    ## Note that the gradients only make sense when the target tensor has size 1
    ## (i.e. when it is either a scalar or a tensor with a single entry).
    ## When the target has size greater than 1, then the gradients are wrt the
    ## sum of all target entries. For details on tf.gradients, see
    ## https://www.tensorflow.org/api_docs/python/tf/gradients
    if (code$cppADCode) {
        placeholders <- as.list(placeholders)
        placeholderValues <- as.list(placeholders)
        placeholderNames <- names(placeholderValues)
        names(placeholderValues) <- NULL
        gradients <- tensorflow::tf$gradients(tensor, placeholderValues)
        if (length(gradients) != length(placeholderNames)) {
            stop('programmer error')
        }
        for (i in seq_along(gradients)) {
            ## This suffixed name is expected by include/nimble/tensorflow.hpp:
            name <- paste0(placeholderNames[i], '_gradient')
            tensorflow::tf$identity(gradients[i], name = name)
        }
    }

    ## Configure the tensorflow session.
    configProto <- tensorflow::tf$ConfigProto()
    configProto$inter_op_parallelism_threads <- threads;
    configProto$intra_op_parallelism_threads <- threads;
    ## This tryCatch works around the following bug in reticulate:
    ## https://github.com/rstudio/reticulate/issues/92
    tryCatch({
        configProto$gpu_options$allow_growth <- TRUE
    }, error = function(e) invisible(NULL))

    ## Serialize protos to strings.
    base64 <- reticulate::import('base64')
    graph <- tensorflow::tf$get_default_graph()
    graph <- base64$b64encode(graph$as_graph_def()$SerializeToString())
    config <- base64$b64encode(configProto$SerializeToString())

    ## Create NimTf_Builder.
    tfBuilder <- TfBuilder()
    tfBuilder$setSerializedGraph(graph)
    tfBuilder$setSerializedConfig(config)
    for (name in names(placeholders)) {
        tfBuilder$addInputVar(name)
    }
    tfBuilder$addOutputVar(target$name)
    return(tfBuilder)
}

## This is a wrapper class around the C++ class TfBuilder.
## The purpose is to build a tensorflow graph and session.
TfBuilder <- setRefClass(
    Class = "TfBuilder",
    fields = list(
        serializedGraph = 'ANY',
        serializedConfig = 'ANY',
        ## Character vector of arguments in canonical order.
        inputNames = 'ANY',
        ## Character vector of outputs in canonical order.
        outputNames = 'ANY'
    ),
    methods = list(
        initialize = function() {
            inputNames <<- outputNames <<- character()
            serializedGraph <<- character()
        },
        setSerializedGraph = function(graphText) {
            serializedGraph <<- graphText
        },
        setSerializedConfig = function(configText) {
            serializedConfig <<- configText
        },
        addInputVar = function(name) {
            inputNames <<- append(inputNames, name)
        },
        addOutputVar = function(name) {
            outputNames <<- append(outputNames, name)
        },
        generateConstructor = function(cppADCode) {
            ## Paste code, instead of creating a parse tree.
            ## Note that we have decided to use the static pointer trick to avoid
            ## crashes due to out-of-order destructors or unlinking. This
            ## solution leaks memory (including GPU memory), but it is tricky and
            ## error-prone to implement a non-leaking solution.
            paste(
                paste0(' = *NimTf_Builder("', serializedGraph, '")'),
                paste0('.NimTf_withConfig("', serializedConfig, '")'),
                paste0(
                    paste0('.NimTf_withInput("', inputNames, '")'),
                    collapse = "\n"
                ),
                paste0(
                    paste0('.NimTf_withOutput("', outputNames, '")'),
                    collapse = "\n"
                ),
                if (cppADCode) '.NimTf_op()' else '.NimTf_runner()',
                sep = "\n"
            )
        }
    )
)
