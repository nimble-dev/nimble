# for use in DSL code check:
otherDSLcalls <- c("{", "[[", "$", "resize", "declare", "returnType", "seq_along", "double", "rankSample")

nimKeyWords <- list(copy = 'nimCopy',
                    print = 'nimPrint',
                    cat = 'nimCat',
                    step = 'nimStep',
                    equals = 'nimEquals',
                    dim = 'nimDim',
                    stop = 'nimStop',
                    numeric = 'nimNumeric',
                    integer = 'nimInteger',
                    matrix = 'nimMatrix',
                    array = 'nimArray',
                    round = 'nimRound',
                    c = 'nimC',
                    rep = 'nimRep',
                    seq = 'nimSeq')

nfMethodRC <- 
    setRefClass(Class   = 'nfMethodRC',
                fields  = list(
                    argInfo    = 'ANY',
                    arguments  = 'ANY',
                    template   = 'ANY',
                    code       = 'ANY',
                    returnType = 'ANY',
                    uniqueName = 'character',
                    neededRCfuns = 'ANY'		#list
                ),
                methods = list(
                    initialize = function(method, name, check = FALSE) {
                    	
                        if(!missing(name)) uniqueName <<- name ## only needed for a pure RC function. Not needed for a nimbleFunction method
                        neededRCfuns <<- list()	
                        argInfo <<- formals(method)
                        code <<- nf_changeNimKeywords(body(method))  ## changes all nimble keywords, e.g. 'print' to 'nimPrint'; see 'nimKeyWords' list at bottom
                        if(code[[1]] != '{')  code <<- substitute({CODE}, list(CODE=code))
                        if(check && "package:nimble" %in% search()) # don't check nimble package nimbleFunctions
                            nf_checkDSLcode(code)
                        generateArgs()
                        generateTemplate() ## used for argument matching
                        removeAndSetReturnType()
                    },
                    generateArgs = function() {
                        argsList <- nf_createAList(names(argInfo))
                        for(i in seq_along(argsList)) { if('default' %in% names(argInfo[[i]]))      argsList[[i]] <- argInfo[[i]]$default }
                        arguments <<- as.pairlist(argsList)
                    },
                    generateTemplate = function() {
                        functionAsList <- list(as.name('function'))
                        functionAsList[2] <- list(NULL)
                        if(!is.null(args)) functionAsList[[2]] <- arguments
                        functionAsList[[3]] <- quote({})
                        template <<- eval(as.call(functionAsList))
                    },
                    removeAndSetReturnType = function() {
                        returnLineNum <- 0
                        for(i in seq_along(code)) {
                            if(length(code[[i]]) > 1) {
                                if(is.name(code[[i]][[1]])) {
                                    if(code[[i]][[1]] == 'returnType') {
                                        returnLineNum <- i
                                        break;
                                    }
                                }
                            }
                        }
                        if(sum(all.names(code) == 'returnType') > 1) stop('multiple returnType() declarations in nimbleFunction method; only one allowed')
                        if(returnLineNum == 0) {   ## no returnType() declaration found; default behavior
                            returnTypeDeclaration <- quote(void())
                        } else {                   ## returnType() declaration was found
                            returnTypeDeclaration <- code[[returnLineNum]][[2]]
                            code[returnLineNum] <<- NULL
                        }
                        ## a very patchy solution: switch nimInteger back to integer
                        if(as.character(returnTypeDeclaration[[1]]) == 'nimInteger') returnTypeDeclaration[[1]] <- as.name('integer')
                        returnType <<- returnTypeDeclaration
                    },
                    generateFunctionObject = function(keep.nfMethodRC = FALSE) {
                        functionAsList <- list(as.name('function'))
                        functionAsList[2] <- list(NULL)
                        if(!is.null(args)) functionAsList[[2]] <- arguments
                        functionAsList[[3]] <- code
                        ans <- eval(parse(text=deparse(as.call(functionAsList)), keep.source = FALSE)[[1]])
                        environment(ans) <- new.env()
                        if(keep.nfMethodRC) {environment(ans)$nfMethodRCobject <- .self}
                        ans
                    },
                    getArgInfo       = function() { return(argInfo)    },
                    getReturnType    = function() { return(returnType) })
    )



nf_checkDSLcode <- function(code) {  
    dslCalls <- c(names(sizeCalls), otherDSLcalls, names(specificCallReplacements), nimKeyWords)
    calls <- setdiff(all.names(code), all.vars(code))
    # find cases of x$y() and x[]$y and x[[]]$y (this also unnecessarily finds x$y)
    names <- all.names(code)
    dollars <- which(names == "$")
    if(length(dollars)) {
        lhs <- dollars+1
        while(sum(names[lhs] %in% c('[', '[[')))
            lhs[names[lhs] %in% c('[', '[[')] <- lhs[names[lhs] %in% c('[', '[[')] + 1
        lhs <- lhs[lhs <= length(names)]
        otherVars <- unique(names[lhs])
    } else otherVars <- NULL
    
    nonDSLcalls <- calls[!(calls %in% c(dslCalls, otherVars))]
    if(length(nonDSLcalls)) {
        objInR <- sapply(nonDSLcalls, exists)
        nonDSLnonR <- nonDSLcalls[!objInR]
        nonDSLinR <- nonDSLcalls[objInR]
        if(length(nonDSLinR)) {
            # problem with passing inputIsName when run through roxygen...
            nonDSLinR <- nonDSLinR[!(sapply(nonDSLinR, function(x) is.nf(x, inputIsName = TRUE)) |
                                     sapply(nonDSLinR, function(x) is.rcf(x, inputIsName = TRUE)))]
            warning(paste0("Detected possible use of R functions in nimbleFunction run code. For this nimbleFunction to compile, these functions must defined as nimbleFunctions or nimbleFunction methods: ", paste(nonDSLinR, collapse = ', '), "."))
            if("c" %in% nonDSLinR) warning("Note that until version 0.6-3 of NIMBLE, c() cannot be used as a stand-alone function, but its use to create vector arguments to a function may be valid.")
        }
        if(length(nonDSLnonR))
            warning(paste0("For this nimbleFunction to compile, these functions must be defined as nimbleFunctions or nimbleFunction methods: ", paste(nonDSLnonR, collapse = ', '), "."))
    }
    return(0)
}

nf_changeNimKeywords <- function(code){
    if(length(code) > 0){
        for(i in seq_along(code) ) {
            if(is.call(code) ) {
                if(!is.null(code[[i]]) ) {
                    code[[i]] <- nf_changeNimKeywordsOne(code[[i]])
                }
            }
        }
    }
    return(code)
}


nf_changeNimKeywordsOne <- function(code, first = FALSE){
    if(length(code) == 1){
        if(as.character(code) %in% names(nimKeyWords)) {
            if(is.call(code)) {
                code[[1]] <- as.name( nimKeyWords[[as.character(code)]] )
            } else {
                if(!is.character(code) & first)
                    code <- as.name( nimKeyWords[[as.character(code)]] )
            }
        }
    }
    else if(length(code) > 1){ 
        for(i in seq_along(code) ) {
            if(!is.null(code[[i]]) )
                code[[i]] <- nf_changeNimKeywordsOne(code[[i]], first = i == 1)
        }
    }
    return(code)
}
