##############################################################
## Section for outputting C++ code from an exprClass object ##
##############################################################

cppOutputCalls <- c(makeCallList(recyclingRuleOperatorsAD, 'cppOutputRecyclingRuleADFunction'),
                    makeCallList(nimDerivsPrependTypeOperators, 'cppOutputNimDerivsPrependType'),
                    makeCallList(binaryMidOperators, 'cppOutputMidOperator'),
                    makeCallList(binaryMidLogicalOperators, 'cppOutputMidOperator'),
                    makeCallList(binaryOrUnaryOperators, 'cppOutputBinaryOrUnary'),
                    makeCallList(assignmentOperators, 'cppOutputMidOperator'),
                    makeCallList(nonNativeEigenProxyCalls, 'cppOutputNonNativeEigen'),
                    makeCallList(eigProxyCalls, 'cppOutputEigMemberFunction'),
                    makeCallList(eigCalls, 'cppOutputMemberFunction'),
                    makeCallList(c('setSize', 'initialize', 'getPtr', 'dim', 'getOffset', 'strides', 'isMap', 'mapCopy', 'setMap'), 'cppOutputMemberFunction'),
                    makeCallList(eigOtherMemberFunctionCalls, 'cppOutputEigMemberFunctionNoTranslate'),
                    makeCallList(eigProxyCallsExternalUnary, 'cppOutputEigExternalUnaryFunction'),
                    makeCallList(c('startNimbleTimer','endNimbleTimer','push_back'), 'cppOutputMemberFunction'),
                    makeCallList(c('nimSeqBy','nimSeqLen', 'nimSeqByLen'), 'cppOutputCallAsIs'),
                    list(nimDerivs_dummy = 'cppOutputNimDerivs'),
                    makeCallList(nimbleListReturningOperators, 'cppNimbleListReturningOperator'),
                    makeCallList(c("TFsetInput_", "TFgetOutput_", "TFrun_"), 'cppOutputMemberFunctionDeref'),
                    list(
                        'next' = 'cppOutputNext',
                        cppComment = 'cppOutputComment',
                        eigenCast = 'cppOutputEigenCast',
                        memberData = 'cppOutputMemberData',
                        fill = 'cppOutputEigMemberFunctionNoTranslate',
                        MAKE_FIXED_VECTOR = 'cppOutputMakeFixedVector',
                        concatenateTemp = 'cppOutputEigBlank',
                        ':' = 'cppOutputColon',
                        '::' = 'cppOutputMidOperator',
                        size = 'cppOutputSize',
                         'for' = 'cppOutputFor',
                         'if' = 'cppOutputIfWhile',
                         'while' = 'cppOutputIfWhile',
                         '[' = 'cppOutputBracket',
                         mvAccessRow = 'cppOutputBracket',
                         nimSwitch = 'cppOutputNimSwitch',
                         getParam = 'cppOutputGetParam',
                         getBound = 'cppOutputGetBound',
                         nimDerivs_calculate = 'cppOutputNimDerivs_calculate',
                         '(' = 'cppOutputParen',
                         resize = 'cppOutputMemberFunctionDeref',
                         nfMethod = 'cppOutputNFmethod',
                         nfVar = 'cppOutputNFvar',
                         makeNewNimbleListObject = "cppNewNimbleList",
                         getsize = 'cppOutputMemberFunctionDeref',
                         getNodeFunctionIndexedInfo = 'cppOutputGetNodeFunctionIndexedInfo',
                         resizeNoPtr = 'cppOutputMemberFunction',
                         cppMemberFunction = 'cppOutputMemberFunctionGeneric',
                         AssignEigenMap = 'cppOutputEigenMapAssign',
                         chainedCall = 'cppOutputChainedCall',
                         template = 'cppOutputTemplate',
                         nimPrint = 'cppOutputCout',
                         nimCat = 'cppOutputCoutNoNewline',
                         return = 'cppOutputReturn',
                         cppPtrType = 'cppOutputPtrType', ## mytype* (needed in templates like myfun<a*>(b)
                         cppDereference = 'cppOutputDereference', ## *(arg)
                         cppPointerDereference = 'cppOutputPointerDereference',  ##(*arg)
                         cppMemberDereference = 'cppOutputMidOperator', ## arg1->arg2
                         cppLiteral = 'cppOutputLiteral',
                         '[[' = 'cppOutputDoubleBracket',
                         as.integer = 'cppOutputCast',
                         as.numeric = 'cppOutputCast',
                         numListAccess = 'cppOutputNumList',
                        blank = 'cppOutputBlank',
                        nimVerbatim = 'cppOutputSkip',
                         callC = 'cppOutputEigBlank', ## not really eigen, but this just jumps over a layer in the parse tree
                         eigBlank = 'cppOutputEigBlank',
                         voidPtr = 'cppOutputVoidPtr',
                        cppLiteral = 'cppOutputLiteral'                        )
                    )
cppOutputCalls[['pow']] <-  'cppOutputPow'
cppMidOperators <- midOperators
cppMidOperators[['::']] <- '::'
cppMidOperators[['%*%']] <- ' * '
cppMidOperators[['cppMemberDereference']] <- '->'
cppMidOperators[['nfVar']] <- '->'
cppMidOperators[['&']] <- ' && '
cppMidOperators[['|']] <- ' || '
for(v in c('$', ':')) cppMidOperators[[v]] <- NULL
for(v in assignmentOperators) cppMidOperators[[v]] <- ' = '

nimCppKeywordsThatFillSemicolon <- c(
    '{',
    'for',
    ifOrWhile,
    'nimSwitch',
    'cppLiteral',
    'cppComment')

## In the following list, the names are names in the parse tree (i.e. the name field in an exprClass object)
## and the elements are the name of the function to use to generate output for that name
## e.g. if the parse tree has "dim(A)" (meaning there is an exprClass object with name = "dim", isCall = TRUE, and a single
## argument that is an exprClass object with name = "A" and isName = TRUE)
## then cppOutputMemberFunction will be called, which will generate A.dim()

## Main function for generating C++ output
nimGenerateCpp <- function(code, symTab = NULL, indent = '', showBracket = TRUE, asArg = FALSE) {
    if(is.numeric(code)) return(if(is.nan(code)) "(nimble_NaN())" else code)
    if(is.character(code)) return(paste0('\"', gsub("\\n","\\\\n", code), '\"'))
    if(is.null(code)) return('R_NilValue')
    if(is.logical(code) ) return(if(code) 'true' else 'false')
    if(is.list(code) ) stop("Error generating C++ code, there is a list where there shouldn't be one.  It is probably inside map information.", call. = FALSE)

    if(length(code$isName) == 0) stop("Error generating C++ code, length(code$isName) == 0.", call. = FALSE)
    if(code$isName) return(exprName2Cpp(code, symTab, asArg))
    if(code$name == '{') {
        iOffset <- as.integer(showBracket)
        ans <- vector('list', length(code$args) + 2*iOffset)
        if(showBracket) ans[[1]] <- paste0(indent, '{')
        newInd <- if(showBracket) paste0(indent, ' ') else indent
        for(i in seq_along(code$args)) {
            oneEntry <- nimGenerateCpp(code$args[[i]], symTab, newInd, FALSE)
            if(code$args[[i]]$isCall) if(!(code$args[[i]]$name %in% nimCppKeywordsThatFillSemicolon)) oneEntry <- pasteSemicolon(oneEntry)
            ans[[i + iOffset]] <- if(showBracket) addIndentToList(oneEntry, newInd) else oneEntry
        }
        if(showBracket) ans[[length(code$args) + 2]] <- paste0(indent, '}')
        return(ans)
    }
    operatorCall <- cppOutputCalls[[code$name]]
    if(!is.null(operatorCall)) return(eval(call(operatorCall, code, symTab)))
    return(cppOutputCallAsIs(code, symTab))
}

exprName2Cpp <- function(code, symTab, asArg = FALSE) {
    if(!is.null(symTab)) {
        sym <- symTab$getSymbolObject(code$name, inherits = TRUE)
        if(!is.null(sym)) return(sym$generateUse(asArg = asArg))
        return(code$name)
    } else {
        return(code$name)
    }
}

cppOutputNext <- function(code, symTab) 'continue'

cppOutputEigenCast <- function(code, symTab) {
    paste0( '(',nimGenerateCpp(code$args[[1]], symTab), ').cast<', code$args[[2]], '>()')
}

cppOutputMakeFixedVector <- function(code, symTab) {
    type <- code$args[[5]]
    fixedName <- code$args[[2]]
    vecName <- code$args[[1]]
    len <- code$args[[3]]
    fixedValues <- paste0(unlist(lapply(code$args[[4]]$args, nimGenerateCpp, symTab) ), collapse = ',')
    paste0(type, ' ', fixedName,'[] = {', fixedValues, '}; std::vector<', type, '> ', vecName,'(', fixedName, ',', fixedName, ' + ', len, ')')
}

cppOutputVoidPtr <- function(code, symTab) {
    paste('static_cast<void *>(&',code$args[[1]]$name,')')
}

cppOutputBlank <- function(code, symTab) NULL

cppOutputSkip <- function(code, symTab) nimGenerateCpp(code$args[[1]], symTab)

cppOutputEigBlank <- function(code, symTab) {
    paste0('(', nimGenerateCpp(code$args[[1]], symTab), ')')
}

cppOutputRecyclingRuleADFunction <- function(code, symTab) {
  if(identical(nimbleUserNamespace$cppADCode, 2L)) {
    code$name <- paste0('nimDerivs_', code$name)
    code$name <- gsub('::', '::nimDerivs_', code$name)
  }
  cppOutputCallAsIs(code, symTab)
}

cppOutputNimDerivsPrependType <- function(code, symTab){
  ## if(isTRUE(nimbleUserNamespace$cppADCode)){
  ##   paste0('nimDerivs_', code$name, '(',
  ##          paste0(unlist(lapply(code$args, function(x){
  ##            if(is.numeric(x) || is.logical(x)){
  ##              return(paste0('TYPE_(', nimGenerateCpp(x, symTab, asArg = TRUE), ')'))
  ##            }
  ##            return(nimGenerateCpp(x, symTab, asArg = TRUE))
  ##          })), collapse = ', '), ')')
  ## }
  ## else
      if(identical(nimbleUserNamespace$cppADCode, 2L)) {
    argList <- list()
    logFixedString <- ''
    for(i in seq_along(code$args)){
      iArg <- code$args[[i]]
      iName <- names(code$args)[i]
      if(iName == 'log' && (is.numeric(iArg) | is.logical(iArg))){
        logFixedString <- '_logFixed'
        argList[[length(argList) + 1]] <- nimGenerateCpp(iArg, symTab,
                                                         asArg = TRUE)
      }
      else if(is.numeric(iArg) || is.logical(iArg)){
        argList[[length(argList) + 1]] <- paste0('CppAD::AD<double>(', 
                                                 nimGenerateCpp(iArg, symTab,
                                                                asArg = TRUE),
                                                 ')')
      }
      else{
        argList[[length(argList) + 1]] <- nimGenerateCpp(iArg, symTab,
                                                         asArg = TRUE)
      }
    }
    paste0('nimDerivs_', code$name, logFixedString, '(', 
           paste0(unlist(argList), collapse = ', '), ')')
  } else {
     paste0(code$name, '(',
           paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), ')')
  }
}

cppOutputNumList <- function(code, symTab) {
    paste0( nimGenerateCpp(code$args[[1]], symTab))
}

cppOutputGetNodeFunctionIndexedInfo <- function(code, symTab) {
    paste0(code$args[[1]]$name,'.info[', code$args[[2]] - 1, ']')
}

cppOutputReturn <- function(code, symTab) {
    if(length(code$args) == 0) {
        return('return')
    }
    cppOutputCallAsIs(code, symTab)
}

cppOutputCout <- function(code, symTab) {
    paste0('_nimble_global_output <<', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = '<<'), '<<\"\\n\"; nimble_print_to_R(_nimble_global_output)')
}

cppOutputCoutNoNewline <- function(code, symTab) {
    paste0('_nimble_global_output <<', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = '<<'), '; nimble_print_to_R(_nimble_global_output)')
}

cppOutputChainedCall <- function(code, symTab) {
    firstCall <- nimGenerateCpp(code$args[[1]], symTab)
    ## now similar to cppOutputCallAsIs
    paste0(firstCall, '(', paste0(unlist(lapply(code$args[-1], nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), ')' )
}

cppOutputNimDerivs <- function(code, symTab) {
    origName <- code$args[[1]]$name
    derivName <- paste0(origName, '_deriv_')
    innerArgs <- code$args[[1]]$args
    iOmit <- which(names(code$args) == 'dropArgs')
    outerArgs <- code$args[-c(1, iOmit)]
    ADinfoArg <- code$aux$ADinfoName
    updateNodesName <- code$aux[['updateNodesName']]
    if(!is.null(updateNodesName)) {
      ADinfoArg <- paste0( ADinfoArg, ".setUpdaterNV(", updateNodesName, ")" )
    }
    paste0(derivName, '(', paste0(unlist(lapply( c(innerArgs, outerArgs), nimGenerateCpp, symTab, asArg = TRUE ) ), collapse = ', '), ',', ADinfoArg, ')')
}

cppOutputNimDerivs_calculate <- function(code, symTab) {
  ## Is there  a dropArgs option for this?
  ADinfoArg <- code$aux$ADinfoName
  paste0('nimDerivs_calculate', '(', ADinfoArg, ',', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE ) ), collapse = ', '), ')')
}

cppNewNimbleList <- function(code, symTab) {
    ## This won't work for something like A$B <- nl$new()
    ## because the first arg of the caller is A$B, not a simple name
    ## But the generator info is embedded in sizeExprs
    listType <- code$sizeExprs$nlProc$cppDef$name
    paste0("new ", listType)
}

cppNimbleListReturningOperator <- function(code, symTab) {
  code$name <- nimbleListReturningFunctionList[[code$name]]$cppName
  cppOutputCallAsIs(code, symTab)
}

cppOutputFor <- function(code, symTab) {
    if(code$args[[2]]$name != ':') stop('Error: for now for loop ranges must be defined with :')
    begin <- nimGenerateCpp(code$args[[2]]$args[[1]], symTab)
    end <- nimGenerateCpp(code$args[[2]]$args[[2]], symTab)
    iterVar <- nimGenerateCpp(code$args[[1]], symTab)
    part1 <- paste0('for(', iterVar ,'=', begin,'; ', iterVar, '<= static_cast<int>(', end, '); ++', iterVar,')')
    part2 <- nimGenerateCpp(code$args[[3]], symTab)
    if(is.list(part2)) {
        part2[[1]] <- paste(part1, part2[[1]])
        return(part2)
    } else {
        return(paste(part1, part2))
    }
}

cppOutputIfWhile <- function(code, symTab) {
    ## if(identical(nimbleUserNamespace$cppADCode, 2L))
    ##     stop("Cannot use 'if' or 'while' statements in derivative-enabled nimbleFunctions.")
    part1 <- paste0(code$name,'(', nimGenerateCpp(code$args[[1]], symTab), ')')
    part2 <- nimGenerateCpp(code$args[[2]], symTab)
    if(is.list(part2)) {
        part2[[1]] <- paste(part1, part2[[1]])
    } else {
        part2 <- list(paste(part1, part2))
    }
    if(length(code$args)==2) return(part2)

    part3 <- nimGenerateCpp(code$args[[3]], symTab)
    if(is.list(part3)) {
        part2[[length(part2)]] <- paste(part2[[length(part2)]], 'else', part3[[1]])
        part3 <- c(part2, part3[-1])
        return(part3)
    } else {
        part2[[length(part2)]] <- paste(part2[[length(part2)]], 'else', part3)
        return(part2)
    }
    stop('Error in cppOutputIf')
}

cppOutputNimSwitch <- function(code, symTab) {
    numChoices <- length(code$args)-2
    if(numChoices <= 0) return('')
    choicesCode <- vector('list', numChoices)
    choiceValues <- eval(nimbleGeneralParseDeparse(code$args[[2]]))
    if(length(choiceValues) != numChoices) stop(paste0('number of switch choices does not match number of indices for ',nimDeparse(code)))
    for(i in 1:numChoices) {
        if(code$args[[i+2]]$name != '{')
            bracketedCode <- insertExprClassLayer(code, i+2, '{')
        choicesCode[[i]] <- list(paste0('case ',choiceValues[i],':'), nimGenerateCpp(code$args[[i+2]], symTab, showBracket = FALSE), 'break;')
    }
    ans <- list(paste('switch(',code$args[[1]]$name,') {'), choicesCode, '}')
    ans
}

cppOutputGetParam <- function(code, symTab) {
    if(length(code$args) < 4) {  ## code$args[[3]] is used for the paramInfo that is only used in size processing
        ans <- paste0('getParam_',code$nDim,'D_',code$type,'(',code$args[[2]]$name,',',code$args[[1]]$name,'.getInstructions()[0])')
    } else {
        iNodeFunction <- paste(cppMinusOne(nimDeparse(code$args[[4]])))
        ans <- paste0('getParam_',code$nDim,'D_',code$type,'(',code$args[[2]]$name,',',code$args[[1]]$name,
                      '.getInstructions()[',iNodeFunction,'],', iNodeFunction ,')')
    }
    return(ans)
}

cppOutputGetBound <- function(code, symTab) {
    if(length(code$args) < 4) {  ## code$args[[3]] is used for the paramInfo that is only used in size processing
        ans <- paste0('getBound_',code$nDim,'D_',code$type,'(',code$args[[2]]$name,',',code$args[[1]]$name,'.getInstructions()[0])')
    } else {
        iNodeFunction <- paste(cppMinusOne(nimDeparse(code$args[[4]])))
        ans <- paste0('getBound_',code$nDim,'D_',code$type,'(',code$args[[2]]$name,',',code$args[[1]]$name,
                      '.getInstructions()[',iNodeFunction,'],', iNodeFunction ,')')
    }
    return(ans)
}

cppOutputEigenMapAssign <- function(code, symTab) {
    useStrides <- length(code$args) > 5
    strideTemplateDec <- if(useStrides) {
        if(!(is.numeric(code$args[[6]]) & is.numeric(code$args[[7]]) ) ) {
            bothStridesDyn <- TRUE
            paste0('EigStrDyn')
        } else {
            bothStridesDyn <- FALSE
            paste0(', Stride<', if(is.numeric(code$args[[6]])) code$args[[6]] else 'Dynamic', ', ', if(is.numeric(code$args[[7]])) code$args[[7]] else 'Dynamic','>')
        }
    } else character()
    strideConstructor <- if(useStrides) {
        paste0(strideTemplateDec, '(', nimGenerateCpp(code$args[[6]], symTab),', ', nimGenerateCpp(code$args[[7]], symTab), ')')
    } else character()
    MapType <- if(!useStrides) {
        paste0('Map< ', code$args[[3]]$name,' >')
    } else {
        if(bothStridesDyn) {
            symTab$getSymbolObject(nimDeparse(code$args[[1]]))$baseType
        } ##'EigenMapStr'
        else paste0('Map< ', code$args[[3]]$name, ', Unaligned, ', strideTemplateDec,' >')
    }
    paste0('new (&', nimGenerateCpp(code$args[[1]], symTab),') ', MapType, '(', paste(c(nimGenerateCpp(code$args[[2]], symTab), nimGenerateCpp(code$args[[4]], symTab), nimGenerateCpp(code$args[[5]], symTab), strideConstructor), collapse = ','), ')')
}

cppOutputSize <- function(code, symTab) {
    ## Windows compiler will give warnings if something.size(), which returns unsigned int, is compared to an int.  Since R has no unsigned int, we cast .size() to int.
    paste0('static_cast<int>(', nimGenerateCpp(code$args[[1]], symTab), '.size())')
}

cppOutputMemberData <- function(code, symTab) {
    paste0( nimGenerateCpp(code$args[[1]], symTab), '.', code$args[[2]]$name)
}

cppOutputMemberFunction <- function(code, symTab) {
    paste0( nimGenerateCpp(code$args[[1]], symTab), '.', paste0(code$name, '(', paste0(unlist(lapply(code$args[-1], nimGenerateCpp, symTab) ), collapse = ', '), ')' ))
}

cppOutputMemberFunctionGeneric <- function(code, symTab) { ##cppMemberFunction(myfun(myobj, arg1, arg2)) will generate myobj.myfun(arg1, arg2)
    cppOutputMemberFunction(code$args[[1]], symTab)
}

cppOutputEigExternalUnaryFunction <- function(code, symTab) {
    info <-  eigProxyTranslateExternalUnary[[code$name]]
    if(length(info) < 3) stop(paste0("Invalid information entry for outputting eigen version of ", code$name), call. = FALSE)
    info1 <- info[1]
    info2 <- info[2]
    info3 <- info[3]
    if(identical(nimbleUserNamespace$cppADCode, 2L)) {
      info1 <- paste0('nimDerivs_', info1)
      info2 <- paste0('CppAD::AD<', info2, '>')
      info3 <- paste0('CppAD::AD<', info3, '>')
    }
    paste0(
      '(', nimGenerateCpp(code$args[[1]], symTab),
      ').unaryExpr(std::ptr_fun<',info2,', ',info3,'>(', info1, '))'
    )
}

## like cppOutputCallAsIs but using eigProxyTranslate on the name
cppOutputNonNativeEigen <- function(code, symTab) {
    ans <- paste0(eigProxyTranslate[code$name], '(', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), ')' )
    if(identical(nimbleUserNamespace$cppADCode, 2L))
        ans <- paste0('nimDerivs_', ans)
    ans
}

cppOutputEigMemberFunction <- function(code, symTab) {
    paste0( '(', nimGenerateCpp(code$args[[1]], symTab), ').', paste0(eigProxyTranslate[code$name], '(', paste0(unlist(lapply(code$args[-1], nimGenerateCpp, symTab) ), collapse = ', '), ')' ))
}

cppOutputEigMemberFunctionNoTranslate <- function(code, symTab) {
    paste0( '(', nimGenerateCpp(code$args[[1]], symTab), ').', paste0(code$name, '(', paste0(unlist(lapply(code$args[-1], nimGenerateCpp, symTab) ), collapse = ', '), ')' ))
}

cppOutputMemberFunctionDeref <- function(code, symTab) {
    paste0( nimGenerateCpp(code$args[[1]], symTab), '->', paste0(code$name, '(', paste0(unlist(lapply(code$args[-1], nimGenerateCpp, symTab) ), collapse = ', '), ')' ))
}

cppOutputNFvar <- function(code, symTab) {
    if(length(code$args) != 2) stop('Error: expecting 2 arguments for operator ',code$name)
    paste0( nimGenerateCpp(code$args[[1]], symTab), '.', code$args[[2]] ) ## No nimGenerateCpp on code$args[[2]] because it should be a string
}

cppOutputNFmethod <- function(code, symTab) {
    if(length(code$args) < 2) stop('Error: expecting at least 2 arguments for operator ',code$name)
    if(code$caller$name == 'chainedCall') {
        paste0( nimGenerateCpp(code$args[[1]], symTab), '.', code$args[[2]]) ##, ## No nimGenerateCpp on code$args[[2]] because it should be a string
    } else {
        ## This bound method obj$run is not part of a call, so we transform it to NimBoundMethod<T>(&T::run, obj).
        
        objectName <- code$args[[1]]$name
        if(objectName == "this") {
            typeName <- code$args[[1]]$type
        } else {
            typeName <- symTab$getSymbolObject(objectName, inherits = TRUE)$baseType
        }
        methodName <- code$args[[2]]
        #paste0('NimBoundMethod<', typeName, '>(&', typeName, '::', methodName, ', ', objectName, ')')
        paste0('NimBind(&', typeName, '::', methodName, ', ', objectName, ')')
    }
    ## This used to take method args in this argList.  But now they are in a chainedCall
}

cppOutputColon <- function(code, symTab) {
    if(length(code$args) != 2) stop('Error: expecting 2 arguments for operator ',code$name)
    paste0( 'nimSeqByD(', paste(nimGenerateCpp(code$args[[1]], symTab), nimGenerateCpp(code$args[[2]], symTab), 1, 0, sep = ','),')');
}

cppOutputMidOperator <- function(code, symTab) {
    if(length(code$args) != 2) stop('Error: expecting 2 arguments for operator ',code$name)
    if(is.null(code$caller)) useParens <- FALSE
    else {
        thisRank <- operatorRank[[code$name]]
        callingRank <- if(!is.null(code$caller)) operatorRank[[code$caller$name]] else NULL
        useParens <- FALSE ## default to FALSE - possibly dangerous if we've missed a case
        if(!is.null(callingRank)) {
            if(!is.null(thisRank)) {
                if(callingRank <= thisRank) useParens <- TRUE
            }
        }
    }

    useDoubleCast <- FALSE
##   if(!isTRUE(nimbleUserNamespace$cppADCode)){
      if(code$name == '/') ## cast the denominator to double if it is any numeric or if it is an scalar integer expression
          if(is.numeric(code$args[[2]]) ) useDoubleCast <- TRUE
          else ## We have cases where a integer ends up with type type 'double' during compilation but should be cast to double for C++, so we shouldn't filter on 'integer' types here
              if(identical(code$args[[2]]$nDim, 0)) useDoubleCast <- TRUE
##    }

    secondPart <- nimGenerateCpp(code$args[[2]], symTab)
    if(useDoubleCast) {
        static_cast <- 'static_cast<double>('
        if(identical(nimbleUserNamespace$cppADCode, 2L))
            static_cast <- 'CppAD::AD<double>('
        secondPart <- paste0(static_cast, secondPart, ')')
    }
    if(useParens)
        paste0( '(',nimGenerateCpp(code$args[[1]], symTab), cppMidOperators[[code$name]],secondPart,')' )
    else
        paste0( nimGenerateCpp(code$args[[1]], symTab), cppMidOperators[[code$name]],secondPart)
}

cppOutputBinaryOrUnary <- function(code, symTab) {
    if(length(code$args) == 2) return(cppOutputMidOperator(code, symTab))
    cppOutputCallAsIs(code, symTab)
}

cppMinusOne <- function(x) {
    if(is.numeric(x)) return(x-1)
    paste0('(',x,') - 1')
}

cppOutputBracket <- function(code, symTab) {
    brackets <- if(length(code$args) <= 2) c('[',']') else c('(',')')
    paste0( nimGenerateCpp(code$args[[1]], symTab), brackets[1], paste0(unlist(lapply(code$args[-1], function(x) cppMinusOne(nimGenerateCpp(x, symTab) ) ) ), collapse = ', '), brackets[2] )
}

cppOutputDoubleBracket <- function(code, symTab) {
    ## right now the only case is from a functionList, so we'll go straight there.
    cppOutputNimFunListAccess(code, symTab)
}

cppOutputNimFunListAccess <- function(code, symTab) {
    paste0( '(*', nimGenerateCpp(code$args[[1]], symTab), '[', cppMinusOne(nimGenerateCpp(code$args[[2]], symTab) ), '])' )
}

cppOutputParen <- function(code, symTab) {
    paste0('(', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab) ), collapse = ', '), ')' )
}

cppOutputCall <- function(code, symTab) {
    paste0(cppOutputCalls[[code$name]], '(', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab ) ), collapse = ', '), ')' )
}

cppOutputPow <- function(code, symTab) {
    useStaticCast <-
        ## if(isTRUE(nimbleUserNamespace$cppADCode))
        ##     FALSE
        ## else
            if(is.numeric(code$args[[2]]) )
            TRUE
        else identical(code$args[[2]]$nDim, 0)
    if(useStaticCast) {
        static_cast <- '( static_cast<double>('
        nimDerivs_text <- ''
          if(identical(nimbleUserNamespace$cppADCode, 2L)){
            static_cast <- '( CppAD::AD<double>('
            nimDerivs_text <- 'nimDerivs_'
        }
        paste0(nimDerivs_text, exprName2Cpp(code, symTab), static_cast,
               nimGenerateCpp(code$args[[1]], symTab, asArg = TRUE),'),', 
               nimGenerateCpp(code$args[[2]], symTab, asArg = TRUE),')')
      } else{
      if(identical(nimbleUserNamespace$cppADCode, 2L)){
        paste0('nimDerivs_', exprName2Cpp(code, symTab), 
               '(', nimGenerateCpp(code$args[[1]], symTab, asArg = TRUE),',',
               nimGenerateCpp(code$args[[2]], symTab, asArg = TRUE),')')
      }
      else{
        paste0(exprName2Cpp(code, symTab), '(', 
               nimGenerateCpp(code$args[[1]], symTab, asArg = TRUE),',',
               nimGenerateCpp(code$args[[2]], symTab, asArg = TRUE),')')
      }
    }
}

cppOutputCallAsIs <- function(code, symTab) {
    paste0(exprName2Cpp(code, symTab), '(', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), ')' )
}

cppOutputCast <- function(code, symTab) {
    paste0('static_cast<', cppCasts[[code$name]], '>(', nimGenerateCpp(code$args[[1]], symTab), ')')
}

cppOutputCallWithName <- function(code, symTab, name) {
    paste0(name, '(', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), ')' )
}

cppOutputPtrType <- function(code, symTab) {
    paste0( nimGenerateCpp(code$args[[1]], symTab), '*' )
}

cppOutputDereference <- function(code, symTab) {
    cppOutputCallWithName(code, symTab, '*')
}

cppOutputPointerDereference <- function(code, symTab) {
  paste0('(*', paste0(unlist(lapply(code$args, nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), ')' )
}

cppOutputLiteral <- function(code, symTab) {
  return(eval(code$args[[1]]))
}

cppOutputTemplate <- function(code, symTab) {
    paste0(code$args[[1]]$name, '<', paste0(unlist(lapply(code$args[-1], nimGenerateCpp, symTab, asArg = TRUE) ), collapse = ', '), '>' )
}

cppOutputLiteral <- function(code, symTab) code$args[[1]]

cppOutputComment <- function(code, symTab) paste0("/* ", code$args[[1]], " */")

