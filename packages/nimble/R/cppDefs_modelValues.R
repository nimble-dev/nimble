cppModelValuesClass <- setRefClass('cppModelValuesClass',
                                   contains = 'cppNamedObjectsClass',
                                   fields = list(
                                       vars = 'ANY'
                                       ),
                                   methods = list(
                                       initialize = function(...) {
                                           Hincludes <<- c(Hincludes, nimbleIncludeFile("NimArr.h"), nimbleIncludeFile("Values.h"))
                                           callSuper(...)
                                           inheritance <<- inheritance[inheritance != 'NamedObjects']
                                           inheritance <<- c(inheritance, 'Values') ## and Values inherits from NamedObjects
                                       },
                                       makeCppNames = function() {
                                           if(is.list(vars)) {
                                               Rnames2CppNames <<- as.list(makeVecName(Rname2CppName(names(vars))))
                                               names(Rnames2CppNames) <<- names(vars)
                                               return(invisible(NULL))
                                           }
                                           if(inherits(vars, 'symbolTable')) {
                                               Rnames2CppNames <<- as.list(makeVecName(Rname2CppName(vars$getSymbolNames())))
                                               names(Rnames2CppNames) <<- vars$getSymbolNames()
                                               return(invisible(NULL))
                                           }
                                           stop('Invalid format for vars in cppModelValuesClass.  Must be a list or a symbolTable')
                                       },
                                       buildVars = function(forAD = FALSE) {
                                           cppDef_MV_buildVars_impl(.self, forAD)
                                       },
                                       buildConstructorFunctionDef = function() {
                                           lastLine <- cppLiteral(c("resize(1);",
                                                                    paste0("buildName = ", paste0('"new_',name,'";'))))
                                           functionDefs[['constructor']] <<- cppFunctionDef(name = name,
                                                                                            returnType = emptyTypeInfo(),
                                                                                            code = cppCodeBlock(code = putCodeLinesInBrackets(list(namedObjectsConstructorCodeBlock(), lastLine)), skipBrackets = TRUE))
                                       },
                                       buildResizeFunctionDef = function() {
                                           template <- quote(VAR %.% resize(nrow))
                                           codeLines <- lapply(Rnames2CppNames, function(x) codeSubstitute(template, list(VAR = as.name(x))))
                                           resizeCodeBlock <- cppCodeBlock(code = putCodeLinesInBrackets(codeLines), skipBrackets = TRUE)

                                           template2 <- quote(VAR[i-1]%.%SETSIZECALL)
                                           codeLines2 <- list()
                                           if(is.list(vars)) {
                                               for(v in names(Rnames2CppNames)) {
                                                   ssc <- as.call(c(list(as.name('setSize')), as.list(vars[[v]])))
                                                   codeLines2[[v]] <- codeSubstitute(template2, list(SETSIZECALL = ssc, VAR = as.name(Rnames2CppNames[[v]])))
                                                   
                                               }
                                           }
                                           if(inherits(vars, 'symbolTable')) {
                                               for(v in names(Rnames2CppNames)) {
                                                   thisSize <- vars$getSymbolObject(v)$size
                                                   if(length(thisSize)==0) thisSize <- 1
                                                   ssc <- as.call(c(list(as.name('setSize')), as.list(thisSize)))
                                                   codeLines2[[v]] <- codeSubstitute(template2, list(SETSIZECALL = ssc, VAR = as.name(Rnames2CppNames[[v]])))
                                               }
                                           }
                                           setSizeForLoop <- substitute(cppFor(i, 1, nrow, CONTENTS), list(CONTENTS = cppCodeBlock(code = putCodeLinesInBrackets(codeLines2), skipBrackets = TRUE)))
                                           
                                           
                                           setSizeCodeBlock <- cppCodeBlock(code = setSizeForLoop, skipBrackets = TRUE)
                                           
                                           functionDefs[['resize']] <<- cppFunctionDef(name = 'resize',
                                                                                       args = list(cppInt('nrow')),
                                                                                       returnType = cppVoid(),
                                                                                       virtual = TRUE,
                                                                                       code = cppCodeBlock(code = putCodeLinesInBrackets(list(resizeCodeBlock, setSizeCodeBlock, cppLiteral("numRows = nrow;"))), objectDefs = list(cppInt('i')), skipBrackets = TRUE))
                                       },
                                       buildAll = function(forAD = FALSE) {
                                           makeCppNames()
                                           buildVars(forAD = forAD)
                                           buildConstructorFunctionDef()
                                           buildResizeFunctionDef()
                                           buildSEXPgenerator(finalizer = 'namedObjects_Finalizer')
                                       }
                                       )
                                   )

cppDef_MV_buildVars_impl <- function(.self, forAD) {
    vars <- .self$vars
    Rnames2CppNames <- .self$Rnames2CppNames
    if(is.list(vars)) {
        for(v in names(Rnames2CppNames)) {
            cName <- Rnames2CppNames[[v]]
            nDim <- max(length(vars[[v]]), 1)
            if(!forAD)
                .self$addObject(cName,
                                cppVecNimArr(name = cName,
                                             nDim = nDim,
                                             type = 'double'))
            else
                .self$addObject(cName,
                                cppVecNimArr(name = cName,
                                             nDim = nDim,
                                             type = 'CppAD::AD<double> '))
        }
        return(invisible(NULL))
    }
    if(inherits(vars, 'symbolTable')) {
        ## Actually in this case we should be able to use symbolTable2cppVars
        for(v in names(Rnames2CppNames)) {
            cName <- Rnames2CppNames[[v]]
            thisSym <- vars$getSymbolObject(v)
            nDim <- max(thisSym$nDim, 1)
            type <- thisSym$type
            if(type == 'integer') type <- 'int'
            if(!forAD) {
                .self$addObject(cName,
                                cppVecNimArr(name = cName,
                                             nDim = nDim,
                                             type = type))
            } else {
                .self$addObject(cName,
                                cppVecNimArr(name = cName,
                                             nDim = nDim,
                                             type = paste0("CppAD::AD<",type,"> ")))
            }
        }
        return(invisible(NULL))
    }
}
