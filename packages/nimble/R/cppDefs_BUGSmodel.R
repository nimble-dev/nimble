## This class represents a C++ BUGS model class
cppBUGSmodelClass <- setRefClass('cppBUGSmodelClass',
                                 contains = 'cppNamedObjectsClass',
                                 fields = list(
                                     modelDef = 'ANY', ## a modelDefClass, but saying 'ANY' to avoid hassles,
                                     model = 'ANY', ## an RmodelBaseClass, needed for the nodeFunctions
                                     nodeFuns = 'ANY',		#list		 ## list of cppNimbleFunctionClass objects
                                     CmodelValuesClassName = 'ANY'		#character
                                     ),
                                 methods = list(
                                     initialize = function(...) {
                                     	nodeFuns <<- list()
                                     	CmodelValuesClassName <<- character()
                                         Hincludes <<- c(Hincludes, nimbleIncludeFile("NimArr.h"),
                                                                    nimbleIncludeFile("ModelClassUtils.h"))
                                         callSuper(...)
                                         CmodelValuesClassName <<- Rname2CppName(modelDef$modelValuesClassName) 
                                         inheritance <<- inheritance[inheritance != 'NamedObjects']
                                         inheritance <<- c(inheritance, 'ModelBase')
                                     },
                                     ## This returns a cppModelValuesClass object
                                     ## so we have the option of putting it in the same file or different file
                                     genModelValuesCppClass = function() {
                                         mvc <- nimbleProject$getModelValuesCppDef(modelDef$modelValuesClass, NULLok = TRUE)
                                         if(is.null(mvc)) {
                                             mvc <- nimbleProject$needModelValuesCppClass(modelDef$modelValuesClass, fromModel = TRUE)
                                         }
                                         mvc
                                     },
                                     ## makeCppNames populates the Rnames2CppNames, which is inherited from the cppNamedObjectsClass
                                     ## This is the list used to generated namedObjects assignments in the constructor
                                     makeCppNames = function() {
                                         ## also need log prob names
                                         Rnames2CppNames <<- as.list(c(Rname2CppName(names(modelDef$varInfo)), Rname2CppName(names(modelDef$logProbVarInfo))))
                                         names(Rnames2CppNames) <<- c(names(modelDef$varInfo), names(modelDef$logProbVarInfo))
                                     },
                                     ## buildVars creates the cppVar objects for each variable needed for the class. Currently those are NimArr<> objects for the variables inthe model
                                     buildVars = function() {
                                         for(v in names(modelDef$varInfo)) {
                                             cname <- Rnames2CppNames[[v]]
                                             nDim <- max(modelDef$varInfo[[v]]$nDim, 1) ## set scalars (nDim = 0) to vectors (dim = 1)
                                             addObject(cname, cppNimArrPtr(name = cname, nDim = nDim, type = 'double'))
                                         }
                                         for(v in names(modelDef$logProbVarInfo)) {
                                             cname <- Rnames2CppNames[[v]]
                                             nDim <- max(modelDef$logProbVarInfo[[v]]$nDim, 1) 
                                             addObject(cname, cppNimArrPtr(name = cname, nDim = nDim, type = 'double'))
                                         }
                                         addObject('_defaultModelValues', cppVar(name = '_defaultModelValues', baseType = CmodelValuesClassName))
                                     },
                                     ## buildConstructorFunctionDef adds a cppFunctionDef object for the constructor to the functionDefs field (inherited from cppNamespaceClass via cppNamedObjectsClass).  Most of this comprises the namedObjects assignments generated from  namedObjectsConstructorCodeBlock() inherited from cppNamedObjectsClass
                                     buildConstructorFunctionDef = function() {

                                         lines12 <- cppLiteral(c("_defaultModelValues.resize(1);",
                                                                 "pointAtAll(&_defaultModelValues, 0);",
                                                                 "_modelValues = static_cast<Values *>(&_defaultModelValues);")) 
                                         
                                         code <- putCodeLinesInBrackets(list(lines12, namedObjectsConstructorCodeBlock())) 
                                         conFunDef <- cppFunctionDef(name = name,
                                                                     returnType = emptyTypeInfo(),
                                                                     code = cppCodeBlock(code = code, skipBrackets = TRUE))
                                         functionDefs[['constructor']] <<- conFunDef
                                     },
                                     ## buildPointAtAll adds a cppFunctionDef for pointAtAll to the functionDefs field.  This consists of lines like one__over_sqrt_tau__between = &(values->one__over_sqrt_tau__between_Vec[i]);  These point each model variable at the corresponding element of a modelValues object.
                                     buildPointAtAll = function() {
                                         template = quote( {NAME = cppReference( values %->% VECNAME[i] ) } )[[2]]

                                         numVars <- length(Rnames2CppNames)
                                         codeLines <- vector('list', length = numVars)
                                         if(numVars > 0) {
                                             
                                             for(iNN in 1:numVars) {
                                                 codeLines[[iNN]] <- codeSubstitute(template, list(NAME = as.name(Rnames2CppNames[[iNN]]),
                                                                                                   VECNAME = as.name(makeVecName(Rnames2CppNames[[iNN]])) ) )
                                             }
                                         }
                                         newFun <- cppFunctionDef(name = 'pointAtAll',
                                                                  returnType = cppVoid(),
                                                               
                                                                  args = list(cppVar(ptr = 1, baseType = CmodelValuesClassName, name = 'values'),
                                                                      cppInt('i')),
                                                                  code = cppCodeBlock( code = putCodeLinesInBrackets(codeLines)))
                                         functionDefs[['pointAtAll']] <<- newFun
                                     },

                                     buildNodes = function(where = globalenv(), debugCpp = FALSE) {
                                         nimbleProject$addNimbleFunctionMulti(model$nodeFunctions, fromModel = TRUE, model$nodeFunctionGeneratorNames)
                                         
                                         nodeFuns <<- nimbleProject$compileNimbleFunctionMulti(model$nodeFunctions, isNode = TRUE,
                                                                                               returnCppClass = TRUE,
                                                                                               fromModel = TRUE,
                                                                                               generatorFunNames = model$nodeFunctionGeneratorNames,
                                                                                               alreadyAdded = TRUE) ## fromModel is redundant here
                                         
                                         browser()
                                         
                                         nodeFuns$y_L3_UID_7$addADclassContent()
                                         addADclassContent()
                                         
                                     },
                                     addADclassContentOneNode = function(nodeFun) {
                                       browser()
                                       funName <- 'calculate_AD_'
                                       independentVarNames <- names(formals(nodeFun$calculate_AD_))
                                       independentVarNames <- independentVarNames[-1] ## For the time being, INDEXNODEINFO_ is still the first argument, but don't want to add it to tape!
                                       ## don't check args for calculate_AD_ currently.  Will we ever?
                                       # for(iArg in seq_along(functionDefs[[funName]]$args$symbols)){
                                       #   arg <- functionDefs[[funName]]$args$symbols[[iArg]]
                                       #   argSym <- nfProc$RCfunProcs[[funName]]$compileInfo$origLocalSymTab$getSymbolObject(arg$name)
                                       #   argName <- names(nfProc$RCfunProcs[[funName]]$nameSubList)[iArg]
                                       #   checkADargument(funName, argSym, argName = argName)
                                       #}
                                       addTypeTemplateFunction(funName)
                                       addADtapingFunction(funName, independentVarNames = independentVarNames, dependentVarNames = 'ANS_' )
                                       addADargumentTransferFunction(funName, independentVarNames = independentVarNames)
                                     },
                                     addADtapingFunction = function( funName, independentVarNames, dependentVarNames, nodeFun) {
                                       ADfunName <- paste0(funName, '_AD_')
                                       regularFun <- nodeFun[[funName]]
                                       newFunName <- paste0(funName, '_callForADtaping_')
                                       tst <- nodeFun$model$modelDef$symTab$copy()
                                       for(i in seq_along(tst$symbols)){
                                         if(!(names(tst$symbols)[i] %in% independentVarNames)){
                                           tst$removeSymbol(names(tst$symbols)[i])
                                         }
                                       }
                                       tst$removeSymbol('y')
                                       functionDefs[[newFunName]] <<- makeADtapingFunction(newFunName, regularFun, ADfunName, independentVarNames, dependentVarNames,
                                                                                           calculateFunc = T, calculateFuncSymTab = )
                                       invisible(NULL)
                                     },
                                     addTypeTemplateFunction = function( funName ) {
                                       browser()
                                       newFunName <- paste0(funName, '_AD_')
                                       regularFun <- RCfunDefs[[funName]]
                                       functionDefs[[newFunName]] <<- makeTypeTemplateFunction(newFunName, regularFun)
                                       invisible(NULL)
                                     },
                                     addADclassContent = function() {
                                       CPPincludes <<- c("<cppad/cppad.hpp>", CPPincludes)
                                       Hincludes <<- c("<cppad/cppad.hpp>", nimbleIncludeFile("nimbleCppAD.h"), Hincludes)
                                       addInheritance("nimbleFunctionCppADbase")
                                       ## cppClass$objectDefs$addSymbol(cppVarFull(name = 'allADtapePtrs_', static = TRUE, baseType = 'vector', templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun', templateArgs = list('double'), ptr = 1))))
                                       ##objectDefs[['vectorADtapePtrs']] <<- cppVarFull(baseType = 'vector', templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun', templateArgs = list('double'), ptr = 1)), static = TRUE, name = 'allADtapePtrs_')
                                       objectDefs[['allADtapePtrs_']] <<- cppVarFull(baseType = 'vector', templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun', templateArgs = list('double'), ptr = 1)), static = TRUE, name = 'allADtapePtrs_')
                                       # ##cppClass$objectDefs$addSymbol(cppVarFull(name = 'ADtapeSetup', baseType = 'nimbleCppADinfoClass'))
                                       objectDefs[['ADtapeSetup']] <<- cppVarFull(name = 'ADtapeSetup', baseType = 'nimbleCppADinfoClass')

                                       for(nodeFun in nodeFuns){
                                         addADclassContentOneNode(nodeFun)
                                       }
                                       # ## static declaration in the class definition
                                       # 
                                       # ## globals to hold the global static definition
                                       globals <- cppGlobalObjects(name = 'staticGlobals', staticMembers = TRUE)
                                       globals$objectDefs[['staticGlobalTape']] <- cppVarFull(baseType = 'vector', templateArgs = list(cppVarFull(baseType = 'CppAD::ADFun', templateArgs = list('double'), ptr = 1)), name = paste0(name,'::allADtapePtrs_'))
                                       # ##globalObjectsDefs[['allADtapePtrs_']] <<- globals
                                       neededTypeDefs[['allADtapePtrs_']] <<- globals
                                       # 
                                       # addStaticInitClass()
                                       # 
                                       # invisible(NULL)
                                     },
                                     
                                     buildAll = function(buildNodeDefs = TRUE, where = globalenv(), ...) {
                                         makeCppNames() 
                                         buildVars()
                                         buildConstructorFunctionDef()
                                         buildSEXPgenerator(finalizer = 'namedObjects_Finalizer')
                                         ##buildSEXPgenerator()
                                         ##buildSEXPfinalizer()
                                         buildPointAtAll()
                                         if(buildNodeDefs) buildNodes(where = where, debugCpp = debugCpp)
                                     }
                                     )
                                 )

compileBUGSmodel <- function(model, name, fileName, dirName, compileNodes = TRUE, writeFiles = !(model$cWritten), 
                             compileCpp = !(model$compiled), loadSO = !(model$loaded), buildInterface = TRUE, 
                             createModel = TRUE, returnInternals = FALSE, where = globalenv(), debugCpp = FALSE) {
    ## Start assuming model is a model (not just a modelDef)
    if(missing(name)) name <- model$getModelDef()$name
    Cname <- Rname2CppName(name)
    if(missing(fileName)) fileName <- Cname
    if(missing(dirName))    dirName <- makeDefaultDirName()

    if(inherits(model, 'RmodelBaseClass')) {
        modelDef <- model$getModelDef()
        modelDefCpp <- cppBUGSmodelClass$new(modelDef = modelDef, model = model, name = Cname) ## model is needed for nodeFunctions
        modelDefCpp$buildAll(buildNodeDefs = compileNodes, where = where, debugCpp = debugCpp)
    }
    if(inherits(model, 'cppBUGSmodelClass')) {
        modelDefCpp <- model
    }
    if(!inherits(model, 'cppProjectClass')) {
        cppProj <- cppProjectClass$new(dirName = dirName)
        mvc <- modelDefCpp$genModelValuesCppClass() 
        cppProj$addClass(mvc, filename = Cname)
        cppProj$addClass(modelDefCpp, Cname)
        if(compileNodes) {
            nfFileName <- paste0(Cname,'_nfCode')
            for(i in names(modelDefCpp$nodeFuns)) {
                cppProj$addClass(modelDefCpp$nodeFuns[[i]], filename = nfFileName)
            }
        }
    } else {
        cppProj <- model
        Cname <- names(cppProj$cppDefs)[2]
        if(compileNodes) nfFileName <- names(cppProj$cppDefs)[3] 
    }

    if(writeFiles) {
        cppProj$writeFiles(Cname)
        if(compileNodes) cppProj$writeFiles(nfFileName)
        model$cWritten = TRUE
    }
    if(compileCpp) {
        compileList <- Cname
        if(compileNodes) compileList <- c(compileList, nfFileName)
        cppProj$compileFile(compileList)
        model$compiled = TRUE
    }
    if(loadSO) {
        cppProj$loadSO(Cname)
        model$loaded = TRUE
    }
    if(returnInternals) {
        return(cppProj)
    } else {
        if(buildInterface) {
            interfaceName <- paste0('C', Cname)
            compiledModel <- cppProj$cppDefs[[2]] 
            newCall <- paste0('new_',Cname) 
            ans <- buildModelInterface(interfaceName, compiledModel, newCall, where = where, dll = cppProj$dll)
            if(!createModel) return(ans) else return(ans(model, where, dll = cppProj$dll))
        }
        return(NULL)
    }
}
