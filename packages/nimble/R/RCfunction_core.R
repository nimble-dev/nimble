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
                    initialize = function(method, name) {
                    	
                        if(!missing(name)) uniqueName <<- name ## only needed for a pure RC function. Not needed for a nimbleFunction method
                        neededRCfuns <<- list()	
                        argInfo <<- formals(method)
                        code <<- nf_changeNimKeywords(body(method))  ## changes all nimble keywords, e.g. 'print' to 'nimPrint'; see 'nimKeyWords' list at bottom
                        if(code[[1]] != '{')  code <<- substitute({CODE}, list(CODE=code))
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
