
## proxy model for model_AD
## This class needs just enough pieces to be used like a model
## for purposes of nodeFunction compilation.
## The model will contain an ADproxyModel
## and then the nodeFunction setup code will extract it.
## The model interface will population the proxy model's CobjectInterface
ADproxyModelClass <- setRefClass(
    'ADproxyModelClass',
    fields = list(
        CobjectInterface = 'ANY' ## needs .basePtr
      , model = 'ANY'
    ),
    methods = list(
        getVarInfo = function(...) {model$getVarInfo(...)},
        initialize = function(Rmodel) {
            model <<- Rmodel
        }
        ## symbolNimArrDoublePtr
    )
)

## This is not a reference class (or other class)
## definition because if it was, then when it is used
## in a nodeFunction, R would check to see that variable
## names exist as fields in the class, and that is
## more trouble than it is worth.  A simple environment
## will pass muster here.
## ADproxyModelClass <- function(Rmodel) {
##     ans <- new.env()
##     model <- Rmodel ## for getVarInfo
##     ans$getVarInfo <- function(...) model$getVarInfo(...)
##     ans$CobjectInterface <- NULL
##     ans
## }

## Convert one symbol object for a C++ var into a symbol object for C++ templated CppAD code
# See symbolTable2templateTypeSymbolTable
cppVarSym2templateTypeCppVarSym <- function(oldSym,
                                            addRef = FALSE,
                                            clearRef = FALSE,
                                            replacementBaseType = 'TYPE_',
                                            replacementTemplateArgs = list()) {
    if(oldSym$baseType == 'double') {
        newSym <- cppVarFull(name = oldSym$name, baseType = replacementBaseType, ref = addRef, templateArgs = replacementTemplateArgs)
        return(newSym)
    }

    newSym <- oldSym$copy()
    if(newSym$baseType == 'NimArr') {
        if(newSym$templateArgs[[2]] == 'double') {
            if(length(replacementTemplateArgs)==0)
                newSym$templateArgs[[2]] <- replacementBaseType
            else
                newSym$templateArgs[[2]] <- cppVarFull(name='', baseType = replacementBaseType, templateArgs = replacementTemplateArgs)
            if(clearRef)
                newSym$ref <- FALSE
        }
    }
    ## Next we replace nimSmartPtr<NIMBLE_ADCLASS> with nimSmartPtr<NIMBLE_ADCLASS_META>
    ## It is somewhat ad hoc to do it here.
    ## At the time of writing this, it always makes sense to do it here.
    ## If there are more general use cases, it might become necessary to split this off to a
    ## separate step
    if(newSym$baseType == "nimSmartPtr")
        if(newSym$templateArgs[[1]] == "NIMBLE_ADCLASS")
            newSym$templateArgs[[1]] <- "NIMBLE_ADCLASS_META"
    
    newSym
}

## Convert a symbol table for C++ vars into a symbol table for C++ for templated CppAD code
## For CppAD, we wrap C++ code in template<class TYPE_> 
## and replace any double with TYPE_
## This includes NimArr<nDim, double> with NimArr<nDim, TYPE_>
## and similar treatmnt for Eigen templated types.
symbolTable2templateTypeSymbolTable <- function(symTab,
                                                addRef = FALSE,
                                                clearRef = FALSE,
                                                replacementBaseType = 'TYPE_',
                                                replacementTemplateArgs = list(),
                                                ignore = character()) {
  newSymTab <- symbolTable()
  symNames <- symTab$getSymbolNames()
  for(sn in symNames) {
    oldSym <- symTab$getSymbolObject(sn)
    inIgnore <- any(unlist(lapply(ignore, function(x) grepl(x, sn))))
    if(inIgnore)
        newSym <- oldSym$copy()
    else
        newSym <- cppVarSym2templateTypeCppVarSym(oldSym,
                                                  addRef = addRef,
                                                  clearRef = clearRef,
                                                  replacementBaseType = replacementBaseType,
                                                  replacementTemplateArgs = replacementTemplateArgs)
    newSymTab$addSymbol(newSym)
  }
  newSymTab
}

## This makes a Cpp function definition object wrapped in template<class TYPE_> and with
## doubles converted to TYPE_s (including in templated use if NimArr and Eigen).
## This is called from an existing version of the cppFunctionDef and returns a separate one
makeTypeTemplateFunction <- function(newName,
                                     .self,
                                     useRecordingInfo = FALSE,
                                     derivControl = list()) {
    newCppFunDef <- RCfunctionDef$new()
    ## use typedefs to change nimble's general typedefs for Eigen locally
    typeDefs <- symbolTable()
    typeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<TYPE_>::typeEigenMapStrd", name = "EigenMapStrd") ) ## these coerces the cppVar system to generate a line of typedef code for us
    typeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<TYPE_>::typeMatrixXd", name = "MatrixXd") )
    newCppFunDef$name <- newName
    newCppFunDef$template <- cppVarFull(name = character(), baseType = 'template', templateArgs = list('class TYPE_'))
    ignore <- derivControl[['nonTemplateArgs']]
    if(is.null(ignore)) ignore <- character()
    newCppFunDef$args <- symbolTable2templateTypeSymbolTable(.self$args, addRef = FALSE, ignore = ignore) ## addRef = TRUE breaks if a literal number is passed.
    if(useRecordingInfo) {
      recordingInfoArg <- cppVarFull(baseType = "nimbleCppADrecordingInfoClass", name = "recordingInfo_")
      newCppFunDef$args$addSymbol(recordingInfoArg)
    }
    localArgs <- symbolTable2templateTypeSymbolTable(.self$code$objectDefs)
    localArgs$setParentST( .self$code$objectDefs$getParentST() ) ## this is the argument symTab, but it should be ok b/c it's only used for names
    newCppFunDef$returnType <- cppVarSym2templateTypeCppVarSym(.self$returnType)
    newCode <- copyExprClass(.self$code$code)
    workEnv <- new.env()
    workEnv$RsymTab <- .self$RCfunProc$compileInfo$newLocalSymTab
    ## Access to symbols by workEnv$RCfunProc$compileInfo$newLocalSymTab$getSymbolObject("run", TRUE)
    ## or workEnv$RCfunProc$compileInfo$newLocalSymTab$symbolExists("run", TRUE)
    ## Access to info about a method by workEnv$RCfunProc$compileInfo$newLocalSymTab$getSymbolObject("run", TRUE)$nfMethodRCobj$enableDerivs
    ## 
    exprClasses_modifyForAD(newCode, localArgs, workEnv = workEnv)
    workEnv$RCfunDef <- NULL
    newCppFunDef$code <- cppCodeBlock(code = newCode, objectDefs = localArgs, typeDefs = typeDefs, 
                                      generatorSymTab = .self$code$objectDefs, cppADCode = 2L)
    list(fun = newCppFunDef,
         nodeFxnVector_name = workEnv[['nodeFxnVector_name']])
}

make_deriv_function <- function(origFun,
                                newFunName,
                                argTransferFunName,
                                meta = FALSE) {
  ADNimbleListName <- nl.getDefinitionContent(ADNimbleList, 'name')
  newFun <- RCfunctionDef$new()
  newFun$name <- newFunName
  typeDefs <- symbolTable()
  if(meta) {
    typeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<TYPE_>::typeEigenMapStrd", name = "EigenMapStrd") ) ## these coerces the cppVar system to generate a line of typedef code for us
    typeDefs$addSymbol(cppVarFull(baseType = "typedef typename EigenTemplateTypes<TYPE_>::typeMatrixXd", name = "MatrixXd") )
    newFun$template <- cppVarFull(name = character(), baseType = 'template', templateArgs = list('class TYPE_'))
  }
  newFun$returnType <- cppVarFull(baseType = 'nimSmartPtr',
                                  templateArgs = ADNimbleListName,
                                  name = 'RETURN_OBJ')
  if(meta)
    newFun$returnType <- cppVarSym2templateTypeCppVarSym(newFun$returnType)
  ## 0. Add argument copies.
  newFun$args <- origFun$args$copy()
  if(meta)
    newFun$args <- symbolTable2templateTypeSymbolTable(newFun$args)
  if(meta)
    newFun$args$addSymbol(cppVarFull(baseType = "nimbleCppADrecordingInfoClass", name = "recordingInfo_"))
  newFun$args$addSymbol(cppNimArr(name = 'ARGZ_nimDerivsOrders_',
                                  nDim = 1, type = 'double', ref = TRUE, const = TRUE))
  newFun$args$addSymbol(cppNimArr(name = 'ARGZ_wrtVector_',
                                  nDim = 1, type = 'double', ref = TRUE, const = TRUE))
  newFun$args$addSymbol(cppVar(name = 'ARGZ_ADinfo_',
                               ref = TRUE,
                               baseType = "nimbleCppADinfoClass"))

  ## 0b. add orders, wrtVector, and ADinfo arguments
  ## 1. add ansList to local symTab
  localVars <- symbolTable()
  returnSym <- cppVarFull(baseType = 'nimSmartPtr',
                                 templateArgs = ADNimbleListName,
                          name = 'returnList_')
  if(meta)
    returnSym <- cppVarSym2templateTypeCppVarSym(returnSym)
  localVars$addSymbol(returnSym)
  ## 2. Create getDerivs_wrapper line
  innerRcall <- do.call('call',
                        c(list(argTransferFunName),
                          lapply(origFun$args$getSymbolNames(), as.name),
                          list(as.name("ARGZ_ADinfo_"))),
                        quote = TRUE
                        )

  getDerivs_wrapper <- if(!meta) 'getDerivs_wrapper' else 'getDerivs_wrapper_meta'
  getDerivsRcall <- substitute(returnList_ <- GETDERIVS_WRAPPER( INNERCALL,
                                                                ARGZ_nimDerivsOrders_,
                                                                ARGZ_wrtVector_ ),
                               list(INNERCALL = innerRcall,
                                    GETDERIVS_WRAPPER = as.name(getDerivs_wrapper)))

  ## 3. create return list
  returnCall <- cppLiteral("return(returnList_);")

  allRCode <- do.call('call',c(list('{'),
                               list(getDerivsRcall,
                                    returnCall)),
                      quote = TRUE)

  allCode <- RparseTree2ExprClasses(allRCode)
  
  newFun$code <- cppCodeBlock(code = allCode,
                              typeDefs = typeDefs,
                              objectDefs = localVars)
  newFun
}

## This makes the function to be called once for CppAD taping
## It sets up AD variables, copies from regular variables into them
## calls the templated version of the member function
## copies the results back out.
## Not that values in the regular variables are not really important during taping.
## Currently those values are intialized to 0.5, which should satisfy needs for (-inf, inf), [0, inf) and [0, 1].
## Ending the tape is not done here.  That is done from the calling function
## (which is in permanent C++, not generated from R)
## We do not assume that in the target function the arguments are independent variables and the
## returned value is the dependent variable.  Those are set by the independentVarNames and dependentVarNames
## makeADtapingFunction <- function(newFunName = 'callForADtaping', targetFunDef, ADfunName, independentVarNames, dependentVarNames, isNode, className = "className") {
##     ## Make new function definition to call for taping (CFT)
##     CFT <- RCfunctionDef$new(static = TRUE)
##     CFT$returnType <- cppVarFull(baseType = "CppAD::ADFun", templateArgs = list('double'), ptr = 1, name = 'RETURN_TAPE_') ##cppVoid()
##     CFT$name <- newFunName
    
##     ## args will always be same; these do not depend on case.  actually now these will be empty.
##     CFT$args <- symbolTable()
##     ## create vector< CppAD::AD<double> > ADindependentVars
##     ADindependentVarsSym <- cppVarFull(name = 'ADindependentVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) ## was ref = TRUE if taking as argument
##     ## create vector< CppAD::AD<double> ADresponseVars
##     ADresponseVarsSym <- cppVarFull(name = 'ADresponseVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) ## ditto
##     ## Add them to arguments symbol table ## switch design and make these local
## ##    CFT$args$addSymbol( ADindependentVarsSym )
## ##    CFT$args$addSymbol( ADresponseVarsSym )
##     ## Make local AD variables for all function inputs and outputs
##     ## e.g. if the original targetFun takes NimArr<1, double>, it's templated CppAD version will take NimArr<1, TYPE_>
##     ## Next line creates local variables for passing to that templated CppAD version
##     localVars <- symbolTable2templateTypeSymbolTable(targetFunDef$args, clearRef = TRUE, replacementBaseType = 'CppAD::AD', replacementTemplateArgs = list('double'))
##     if(isNode){
##       localVars$removeSymbol('ARG1_INDEXEDNODEINFO__')
##       indexNodeInfoSymbol <- symbolInternalType(name = 'ARG1_INDEXEDNODEINFO__', argList = list('indexedNodeInfoClass'))
##     }
    
##     ## and similar for the return variable
##     initADptrCode <- cppLiteral("RETURN_TAPE_ = new CppAD::ADFun<double>;")
##     ansSym <- cppVarSym2templateTypeCppVarSym(targetFunDef$returnType, clearRef = TRUE, replacementBaseType = 'CppAD::AD', replacementTemplateArgs = list('double'))
##     ansSym$name <- 'ANS_'
##     localVars$addSymbol(ansSym)
##     symNames <- localVars$getSymbolNames()
##     ## set up a set of index variables for copying code, up to six to be arbitrary (allowing up to 6-dimensional nimble objects to be handled)
##     indexVarNames <- paste0(letters[9:14],'_')
##     ## set any sizes, which must be known
##     nimbleSymTab <- targetFunDef$RCfunProc$compileInfo$newLocalSymTab

##     ## This creates lines like setSize(z, 2 3)
##     ## which the C++ output generator turns into something like z.resize(2, 3)
##     setSizeLines <- vector('list', length(symNames) + 2) ## extra 2 are for the ADindependentVars and ADresponseVars
##     iNextLine <- 1
    
##     for(iSym in seq_along(symNames)) {
##         thisSymName <- symNames[iSym]
##         if(thisSymName == 'ANS_') {
##             thisSym <- targetFunDef$RCfunProc$compileInfo$returnSymbol
##         } else {
##             thisSym <- nimbleSymTab$getSymbolObject(thisSymName)
##         }
##         if(thisSym$nDim > 0) {
##             setSizeCall <- do.call('call',c(list('setSize', quote(as.name(thisSymName))), as.list(thisSym$size))) 
##             setSizeLines[[iNextLine]] <- setSizeCall ##RparseTree2ExprClasses(setSizeCall)
##             iNextLine <- iNextLine + 1
##         } else {
##             setSizeLines[[iNextLine]] <- NULL
##         }
##     }

##     localVars$addSymbol( ADindependentVarsSym )
##     localVars$addSymbol( ADresponseVarsSym )
##     localVars$addSymbol( CFT$returnType )

##     ## call CppAD::Independent(ADindependentVars)
##     ## This starts CppADs taping system
##     CppADindependentCode <- quote(`CppAD::Independent`(ADindependentVars)) ##nimble:::RparseTree2ExprClasses(quote(`CppAD::Independent`(ADindependentVars)))

##     ## make copying blocks into independent vars
##     ## This looks like e.g.
##     ## for(i_ in 1:3) {ADindependentVars[netIncrement_] = x[i]; netIncrement_ <- netIncrement + 1;}
##     numIndependentVars <- length(independentVarNames)
##     copyIntoIndepVarCode <- vector('list', numIndependentVars+1)
##     ## create the netIncrement_ variable and code to initialize it to 1
##     localVars$addSymbol( cppVar(name = 'netIncrement_', baseType = 'int') )
##     copyIntoIndepVarCode[[1]] <- quote(netIncrement_ <- 1) 
##     ## getting the sizes is going to be trickier when an independent var is really an expression, in particular with indexing, like model$x[3]
##     ## for now let's assume only cleanly defined vars.
##     ## one approach would be intermediate variables
##     totalIndependentLength <- 0
##     maxSize <- 1
##     for(ivn in seq_along(independentVarNames)) {
##         thisName <- independentVarNames[ivn]
##         thisSym <- nimbleSymTab$getSymbolObject(thisName)
##         if(thisSym$nDim > 0) {
##             thisSizes <- thisSym$size
##             sizeList <- lapply(thisSizes, function(x) c(1, x))
##             names(sizeList) <- indexVarNames[1:length(sizeList)]
##             if(length(sizeList) > maxSize) maxSize <- length(sizeList)
##             newRcode <- makeCopyingCodeBlock(as.name(thisName), quote(ADindependentVars), sizeList, indicesRHS = FALSE, incrementIndex = quote(netIncrement_), isNode)
##             copyIntoIndepVarCode[[ivn+1]] <- newRcode 
##             totalIndependentLength <- totalIndependentLength + prod(thisSizes)
##         } else {
##             copyIntoIndepVarCode[[ivn+1]] <- substitute({LHS <- ADindependentVars[netIncrement_]; netIncrement_ <- netIncrement_ + 1}, list(LHS = as.name(thisName))) 
##             totalIndependentLength <- totalIndependentLength + 1
##         }
##     }

##     ## put dummy values in ADindependentVars
##     dummyValueRcode <- substitute(for(III in 1:TOTLENGTH) ADindependentVars[III] = 1, list(III = as.name(indexVarNames[1]), TOTLENGTH = totalIndependentLength))
    
##     if(isNode){
##       dummyIndexNodeInfoCode <- list(cppLiteral('indexedNodeInfo ARG1_INDEXEDNODEINFO__ = generateDummyIndexedNodeInfo();'))
##     }
##     else   dummyIndexNodeInfoCode <- list()
##     ## call the taping function
##     TCFcall <- do.call('call', c(list(ADfunName), lapply(targetFunDef$args$getSymbolNames(), as.name)), quote = TRUE)
##     tapingCallRCode <- substitute(ANS_ <- TCF, list(TCF = TCFcall))
    
##     ## make copying blocks from dependent vars
##     numDependentVars <- length(dependentVarNames)
##     copyFromDepVarCode <- vector('list', numDependentVars+1)
##     copyFromDepVarCode[[1]] <- quote(netIncrement_ <- 1) 
##     totalDepLength <- 0;
##     for(ivn in seq_along(dependentVarNames)) {
##         thisName <- dependentVarNames[ivn]
##         if(thisName == 'ANS_') {
##             thisSym <- targetFunDef$RCfunProc$compileInfo$returnSymbol
##         } else {
##             thisSym <- nimbleSymTab$getSymbolObject(thisName)
##         }
##         if(thisSym$nDim > 0) {
##             thisSizes <- thisSym$size
##             sizeList <- lapply(thisSizes, function(x) c(1, x))
##             names(sizeList) <- indexVarNames[1:length(sizeList)]
##             if(length(sizeList) > maxSize) maxSize <- length(sizeList)
##             newRcode <- makeCopyingCodeBlock(quote(ADresponseVars), as.name(thisName), sizeList, indicesRHS = TRUE, incrementIndex = quote(netIncrement_))
##             copyFromDepVarCode[[ivn+1]] <- newRcode 
##             totalDepLength <- totalDepLength + prod(thisSizes)
##         } else {
##             copyFromDepVarCode[[ivn+1]] <- substitute({ADresponseVars[netIncrement_] <- RHS; netIncrement_ <- netIncrement_ + 1}, list(RHS = as.name(thisName))) 
##             totalDepLength <- totalDepLength + 1
##         }
##     }

##     for(ivn in 1:maxSize)
##       localVars$addSymbol( cppVar(name = indexVarNames[ivn], baseType = 'int') )
    
    
##     ## Now that we know how big ADindependenVars and ADresponseVars should be, 
##     ## we can make two more entries to setSizeCalls for them
##     ## Note that code for these will appear above code that uses them.
##     setSizeLines[[iNextLine]] <- substitute(cppMemberFunction(resize(ADindependentVars, TIL)), list(TIL = totalIndependentLength))
##     iNextLine <- iNextLine + 1
##     setSizeLines[[iNextLine]] <- substitute(cppMemberFunction(resize(ADresponseVars, TDL)), list(TDL = totalDepLength))

##     ## line to finish taping
##     finishTapingCall <- cppLiteral('RETURN_TAPE_->Dependent(ADindependentVars, ADresponseVars);')

##     ADoptimizeCalls <- list(
##         # cppLiteral(paste0("std::cout<<\"about to optimize for ", className,"\"<<std::endl;")),
##         # cppLiteral("std::cout<<\"size before optimize = \"<< RETURN_TAPE_->size_var() <<\"\\n\";"),
##                             cppLiteral("RETURN_TAPE_->optimize();"))
##         #                     cppLiteral("std::cout<<\"size after optimize = \"<< RETURN_TAPE_->size_var() <<\"\\n\";"))

##     returnCall <- cppLiteral("return(RETURN_TAPE_);")
    
##     ## Finally put together all the code, parse it into the nimble exprClass system,
##     ## and add it to the result (CFT)
##     allRcode <- do.call('call', c(list('{'), setSizeLines, dummyIndexNodeInfoCode, list(initADptrCode, dummyValueRcode, CppADindependentCode), copyIntoIndepVarCode, list(tapingCallRCode), copyFromDepVarCode, list(finishTapingCall), ADoptimizeCalls, list(returnCall)), quote=TRUE)
##     allCode <- RparseTree2ExprClasses(allRcode)
##     CFT$code <- cppCodeBlock(code = allCode, objectDefs = localVars)
##     CFT
## }

makeADtapingFunction2 <- function(newFunName = 'callForADtaping',
                                  targetFunDef,
                                  ADfunName,
                                  independentVarNames,
                                  dependentVarNames,
                                  isNode,
                                  className = "className",
                                  useModelInfo = list()) {
  nodeFxnVector_name <- useModelInfo[['nodeFxnVector_name']]
  usesModelCalculate <- length(nodeFxnVector_name) > 0
  ## Make new function definition to call for taping (CFT)
  if(isNode) warning("makeADtapingFunction2 has not been updated for isNode==TRUE")
  CFT <- RCfunctionDef$new(static = FALSE)
  CFT$returnType <- cppVarFull(baseType = "CppAD::ADFun", templateArgs = list('double'), ptr = 1, name = 'RETURN_TAPE_') ##cppVoid()
  CFT$name <- newFunName
  CFT$args <- targetFunDef$args$copy()

    ## create vector< CppAD::AD<double> > ADindependentVars
  ADindependentVarsSym <- cppVarFull(name = 'ADindependentVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) 
  ## create vector< CppAD::AD<double> > ADdynamicVars (needed only when a model will be used)
  ADdynamicVarsSym <- cppVarFull(name = 'ADdynamicVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) 
  ## create vector< CppAD::AD<double> ADresponseVars
  ADresponseVarsSym <- cppVarFull(name = 'ADresponseVars', baseType = 'vector', templateArgs = list( cppVarFull(baseType = 'CppAD::AD', templateArgs = 'double', name = character()) ), ref = FALSE) 
  ## Add them to arguments symbol table ## switch design and make these local
##    CFT$args$addSymbol( ADindependentVarsSym )
##    CFT$args$addSymbol( ADresponseVarsSym )
    ## Make local AD variables for all function inputs and outputs
    ## e.g. if the original targetFun takes NimArr<1, double>, it's templated CppAD version will take NimArr<1, TYPE_>
    ## Next line creates local variables for passing to that templated CppAD version
  localVars <- symbolTable2templateTypeSymbolTable(targetFunDef$args,
                                                   clearRef = TRUE,
                                                   replacementBaseType = 'CppAD::AD',
                                                   replacementTemplateArgs = list('double'))
  makeADname <- function(x) paste0(x, "AD_")
  for(varName in CFT$args$getSymbolNames()) {
    ## Move this to a symbolTable$changeSymbolName method
    iName <- which(names(localVars$symbols)==varName)
    newName <- makeADname(varName)
    if(length(iName) > 0) {
      localVars$symbols[[iName[1] ]]$name <- newName
      names(localVars$symbols)[iName[1] ] <- newName
    }
      
  }

    if(isNode){
      localVars$removeSymbol('ARG1_INDEXEDNODEINFO__')
      indexNodeInfoSymbol <- symbolInternalType(name = 'ARG1_INDEXEDNODEINFO__', argList = list('indexedNodeInfoClass'))
    }
    
    ## and similar for the return variable
    initADptrCode <- cppLiteral("RETURN_TAPE_ = new CppAD::ADFun<double>;")
  ansSym <- cppVarSym2templateTypeCppVarSym(targetFunDef$returnType,
                                            clearRef = TRUE,
                                            replacementBaseType = 'CppAD::AD',
                                            replacementTemplateArgs = list('double'))
    ansSym$name <- 'ANS_'
  localVars$addSymbol(ansSym)


    ## Add a model initialization step if a model is used.
    ## Arguably this should go in the TypeTemplateFunction, using CppAD::Value() to copy values without recording in the tape.
    modelInitCode <- quote(blank())
    if(usesModelCalculate) {
        modelInitCode <- cppLiteral("initialize_AD_model_before_recording(*ADinfo.updaterNV());")
##            substitute(initialize_AD_model_before_recording(NV),
##                                    list(NV = as.name(nodeFxnVector_name[1])))
    }

    ## Make a separate resize line for ADresponseVars.  ANS_ does not need a resize line.
    symNames <- CFT$args$getSymbolNames()
    ## set up a set of index variables for copying code, up to six to be arbitrary (allowing up to 6-dimensional nimble objects to be handled)
    indexVarNames <- paste0(letters[9:14],'_')
    ## set any sizes, which must be known
    nimbleSymTab <- targetFunDef$RCfunProc$compileInfo$newLocalSymTab

    ## This creates lines like setSize(z, 2 3)
    ## which the C++ output generator turns into something like z.resize(2, 3)
    setSizeLines <- vector('list', length(symNames) + 1) ## extra 1 is for the ADindependentVars. ADresponseVars is one at a later step
    iNextLine <- 1
    
    for(iSym in seq_along(symNames)) {
      thisSymName <- symNames[iSym]
      thisADname <- makeADname(thisSymName)
        ## if(thisSymName == 'ANS_') {
        ##     thisSym <- targetFunDef$RCfunProc$compileInfo$returnSymbol
        ## } else {
        ##     thisSym <- nimbleSymTab$getSymbolObject(thisSymName)
        ## }
        thisSym <- nimbleSymTab$getSymbolObject(thisSymName)
        if(thisSym$nDim > 0) {
            dimExprs <- lapply(1:thisSym$nDim,
                               function(i)
                                   substitute(dim(NAME)[i],
                                              list(NAME = as.name(thisSymName),
                                                   i = i)))
            setSizeCall <- do.call('call',c(list('setSize', as.name(thisADname)), dimExprs), quote = TRUE) 
            setSizeLines[[iNextLine]] <- setSizeCall 
            iNextLine <- iNextLine + 1
        } else {
            setSizeLines[[iNextLine]] <- NULL
        }
    }

    localVars$addSymbol( ADindependentVarsSym )
    localVars$addSymbol( ADdynamicVarsSym )
    localVars$addSymbol( ADresponseVarsSym )
    localVars$addSymbol( CFT$returnType )

  recordingInfoSym <- cppVarFull(name = "recordingInfo_", baseType = "nimbleCppADrecordingInfoClass",
                                 constructor = "(CppAD::AD<double>::get_tape_id_nimble(), CppAD::AD<double>::get_tape_handle_nimble())")
  localVars$addSymbol(recordingInfoSym)
  setRecordingFalseLine <- cppLiteral("recordingInfo_.recording()=false;")
  setRecordingTrueLine <- cppLiteral("recordingInfo_.recording()=true;")
  
    ## call CppAD::Independent(ADindependentVars)
    ## This starts CppADs taping system
    CppADindependentCode <- if(usesModelCalculate)
                                quote(`CppAD::Independent`(ADindependentVars, 0, true, ADdynamicVars)) ## consider switching to false for speed?
                            else
                                quote(`CppAD::Independent`(ADindependentVars))
    ## make copying blocks into independent vars
    ## This looks like e.g.
    ## for(i_ in 1:3) {ADindependentVars[netIncrement_] = x[i]; netIncrement_ <- netIncrement + 1;}
    numIndependentVars <- length(independentVarNames)
    copyIntoIndepVarCode <- vector('list', numIndependentVars+1)
    ## create the netIncrement_ variable and code to initialize it to 1
    localVars$addSymbol( cppVar(name = 'netIncrement_', baseType = 'int') )
    copyIntoIndepVarCode[[1]] <- quote(netIncrement_ <- 1) 
    ## getting the sizes is going to be trickier when an independent var is really an expression, in particular with indexing, like model$x[3]
    ## for now let's assume only cleanly defined vars.
    ## one approach would be intermediate variables
    ## totalIndependentLength <- 0
    initADindependentVarsCode <- vector('list', numIndependentVars+1)
    initADindependentVarsCode[[1]] <- quote(netIncrement_ <- 1)
    
    maxSize <- 1
    for(ivn in seq_along(independentVarNames)) {
      thisName <- independentVarNames[ivn]
      thisNameAD <- paste0(thisName, "AD_")
        thisSym <- nimbleSymTab$getSymbolObject(thisName)
        if(thisSym$nDim > 0) {
            thisSizes <- thisSym$size
            sizeList <- lapply(1:thisSym$nDim,
                                 function(x) list(1,
                                                  substitute(dim(RHS)[INDEX],
                                                             list(RHS = as.name(thisName),
                                                                  INDEX = x))))
            names(sizeList) <- indexVarNames[1:length(sizeList)]
            if(length(sizeList) > maxSize) maxSize <- length(sizeList)
            
            newRcode <- makeCopyingCodeBlock(as.name(thisNameAD), quote(ADindependentVars), sizeList, indicesRHS = FALSE, incrementIndex = quote(netIncrement_), isNode)
            copyIntoIndepVarCode[[ivn+1]] <- newRcode
            newRcode <- makeCopyingCodeBlock(quote(ADindependentVars),
                                             as.name(thisName),
                                             sizeList,
                                             indicesRHS = TRUE,
                                             incrementIndex = quote(netIncrement_),
                                             isNode)
            initADindependentVarsCode[[ivn+1]] <- newRcode
      ##      totalIndependentLength <- totalIndependentLength + prod(thisSizes)
        } else {
            copyIntoIndepVarCode[[ivn+1]] <- substitute({LHS <- ADindependentVars[netIncrement_]; netIncrement_ <- netIncrement_ + 1},
                                                        list(LHS = as.name(thisNameAD)))
            initADindependentVarsCode[[ivn+1]] <- substitute({ADindependentVars[netIncrement_] <- RHS; netIncrement_ <- netIncrement_ + 1},
                                                             list(RHS = as.name(thisName)))
      ##      totalIndependentLength <- totalIndependentLength + 1
        }
    }

    ## put dummy values in ADindependentVars
    ## dummyValueRcode <- substitute(for(III in 1:TOTLENGTH) ADindependentVars[III] = 1, list(III = as.name(indexVarNames[1]), TOTLENGTH = totalIndependentLength))
    
    if(isNode){
      dummyIndexNodeInfoCode <- list(cppLiteral('indexedNodeInfo ARG1_INDEXEDNODEINFO__ = generateDummyIndexedNodeInfo();'))
    }
    else   dummyIndexNodeInfoCode <- list()
    ## call the taping function
  TCFcall <- do.call('call', c(list(ADfunName),
                               lapply(targetFunDef$args$getSymbolNames(),
                                      function(x) as.name(makeADname(x))),
                               list(as.name(recordingInfoSym$name))),
                     quote = TRUE)
    tapingCallRCode <- substitute(ANS_ <- TCF, list(TCF = TCFcall))
    
    ## make copying blocks from dependent vars
    numDependentVars <- length(dependentVarNames)
    copyFromDepVarCode <- vector('list', numDependentVars+1)
    copyFromDepVarCode[[1]] <- quote(netIncrement_ <- 1) 
    totalDepLength <- 0;
    for(ivn in seq_along(dependentVarNames)) { ## typically this will just be "ANS_"  It could also be an argument that will be used to get results.
        thisName <- dependentVarNames[ivn]
        if(thisName == 'ANS_') {
            thisSym <- targetFunDef$RCfunProc$compileInfo$returnSymbol
        } else {
            thisSym <- nimbleSymTab$getSymbolObject(thisName)
        }
        if(thisSym$nDim > 0) {
            thisSizes <- thisSym$size
            sizeList <- lapply(1:thisSym$nDim,
                               function(x) list(1,
                                                substitute(dim(RHS)[INDEX],
                                                           list(RHS = as.name(thisName),
                                                                INDEX = x))))
            names(sizeList) <- indexVarNames[1:length(sizeList)]
            if(length(sizeList) > maxSize) maxSize <- length(sizeList)
            newRcode <- makeCopyingCodeBlock(quote(ADresponseVars), as.name(thisName), sizeList, indicesRHS = TRUE, incrementIndex = quote(netIncrement_))
            copyFromDepVarCode[[ivn+1]] <- newRcode 
            totalDepLength <- totalDepLength + prod(thisSizes)
        } else {
            copyFromDepVarCode[[ivn+1]] <- substitute({ADresponseVars[netIncrement_] <- RHS; netIncrement_ <- netIncrement_ + 1}, list(RHS = as.name(thisName))) 
            totalDepLength <- totalDepLength + 1
        }
    }

    for(ivn in 1:maxSize)
      localVars$addSymbol( cppVar(name = indexVarNames[ivn], baseType = 'int') )
    
    
    ## Now that we know how big ADindependenVars and ADresponseVars should be, 
    ## we can make two more entries to setSizeCalls for them
    ## Note that code for these will appear above code that uses them.
    localVars$addSymbol(cppInt(name = "totalIndependentLength_"))
    localVars$addSymbol(cppInt(name = "totalDepLength_"))
    calcTotalLengthCode <- makeCalcTotalLengthBlock(independentVarNames,
                                                    nimbleSymTab,
                                                    "totalIndependentLength_",
                                                    0) #if(usesModelCalculate) 1 else 0) ## for extraInputDummy
  ansNDim <- targetFunDef$RCfunProc$compileInfo$returnSymbol$nDim
  if(ansNDim > 0)
    calcTotalResponseLengthCode <- quote(cppLiteral("totalDepLength_ = ANS_.size();"))
  else
    calcTotalResponseLengthCode <- quote(cppLiteral("totalDepLength_ = 1;"))
      ## calcTotalResponseLengthCode <- makeCalcTotalLengthBlock("ANS_",
    ##                                                         localVars,
    ##                                                         "totalDepLength_")
    setSizeLines[[iNextLine]] <- substitute(cppMemberFunction(resize(ADindependentVars, totalIndependentLength_)))
    
    setADresponseVarsSizeLine <- substitute(cppMemberFunction(resize(ADresponseVars, totalDepLength_)))

    ## line to finish taping
    finishTapingCall <- cppLiteral('RETURN_TAPE_->Dependent(ADindependentVars, ADresponseVars);')

    ADoptimizeCalls <- list(
        # cppLiteral(paste0("std::cout<<\"about to optimize for ", className,"\"<<std::endl;")),
        # cppLiteral("std::cout<<\"size before optimize = \"<< RETURN_TAPE_->size_var() <<\"\\n\";"),
                            cppLiteral("RETURN_TAPE_->optimize();"))
        #                     cppLiteral("std::cout<<\"size after optimize = \"<< RETURN_TAPE_->size_var() <<\"\\n\";"))

    returnCall <- cppLiteral("return(RETURN_TAPE_);")

  initADdynamicVarsCode <- if(usesModelCalculate) {
                               cppLiteral("init_dynamicVars(*ADinfo.updaterNV(), ADdynamicVars);")
##                                 substitute(init_dynamicVars(NV, ADdynamicVars),
                               ##                                             list(NV = as.name(nodeFxnVector_name[1])))
                               } else
                                   quote(blank())
    
    copyDynamicVarsToModelCode <- if(usesModelCalculate) {
                                      cppLiteral("copy_dynamicVars_to_model(*ADinfo.updaterNV(), ADdynamicVars);")
##                                      substitute(copy_dynamicVars_to_model(NV, ADdynamicVars),
##                                             list(NV = as.name(nodeFxnVector_name[1])))
                                  } else
                                      quote(blank())
    
    ## Finally put together all the code, parse it into the nimble exprClass system,
    ## and add it to the result (CFT)
    allRcode <- do.call('call', c(list('{'),
                                  list(modelInitCode),
                                  calcTotalLengthCode,
                                  setSizeLines,
                                  dummyIndexNodeInfoCode,
                                  initADdynamicVarsCode,
                                  initADindependentVarsCode,
                                  list(initADptrCode,
                                       setRecordingFalseLine,
                                       tapingCallRCode,
                                       CppADindependentCode,
                                       setRecordingTrueLine,
                                       copyDynamicVarsToModelCode),
                                  copyIntoIndepVarCode,
                                  list(tapingCallRCode,
                                       calcTotalResponseLengthCode,
                                       setADresponseVarsSizeLine),
                                  copyFromDepVarCode,
                                  list(finishTapingCall),
                                  ADoptimizeCalls,
                                  list(returnCall)),
                        quote=TRUE)
    allCode <- RparseTree2ExprClasses(allRcode)
  CFT$code <- cppCodeBlock(code = allCode, objectDefs = localVars)
  CFT$args$addSymbol(cppVar(baseType = "nimbleCppADinfoClass", ref = TRUE, name = "ADinfo"))
    CFT
}

addADinfoObjects <- function(cppDef) {
  firstStatic <- TRUE
  globals <- NULL
  for(i in seq_along(cppDef$nimCompProc$compileInfos)) {
    ADinfoNames <- cppDef$nimCompProc$compileInfos[[i]]$typeEnv[['ADinfoNames']]
    ADinfoNames_calculate <- cppDef$nimCompProc$compileInfos[[i]]$typeEnv[['ADinfoNames_calculate']]
    ADinfoNames <- c(ADinfoNames, ADinfoNames_calculate)
    if(!is.null(ADinfoNames)) {
      for(ADinfoName in ADinfoNames)
        cppDef$objectDefs$addSymbol(cppVar(name = ADinfoName, ptr = 0, baseType = "nimbleCppADinfoClass"))
    }
    ADstaticInfoNames <- cppDef$nimCompProc$compileInfos[[i]]$typeEnv[['ADstaticInfoNames']]
    if(!is.null(ADstaticInfoNames)) {
      for(ADstaticInfoName in ADstaticInfoNames) {
        cppDef$objectDefs$addSymbol(cppVarFull(name = ADstaticInfoName, ptr = 0, baseType = "nimbleCppADinfoClass", static = TRUE))
        if(firstStatic) {
          globals <- cppGlobalObjects(name = paste0('staticGlobals_', cppDef$name),
                                      staticMembers = TRUE)
          firstStatic <- FALSE
        }
        globals$objectDefs[[ADstaticInfoName]] <-
          cppVarFull(baseType = 'nimbleCppADinfoClass',
                     name = paste0(cppDef$name,'::', ADstaticInfoName))
      }
    }
    
  }
  if(!is.null(globals))
    cppDef$neededTypeDefs[['staticTapeInfos']] <- globals
}

## makeStaticInitClass <- function(cppDef, derivMethods) {
##     if(length(derivMethods) == 0) return(NULL)
##     cppClass <- cppClassDef(name = paste0('initTape_', cppDef$name), useGenerator = FALSE)
##     globalsDef <- cppGlobalObjects(name = paste0('initTapeGlobals_', cppDef$name))
##     globalsDef$objectDefs[['staticInitClassObject']] <- cppVarFull(baseType = paste0('initTape_', cppDef$name),
##                                                                    name = paste0('initTape_', cppDef$name, '_Object_'))
##     initializerCodeList <- list()
##     initializerDef <- cppFunctionDef(name = paste0('initTape_', cppDef$name), returnType = emptyTypeInfo())
##     for(derivFun in derivMethods){
##         ADinfoNames <- cppDef$nimCompProc$compileInfos[[derivFun]]$typeEnv[['ADinfoNames']]
        
##         ## use of parse instead of substitute is so R CMD check doesn't flag CLASSNAME:: as an unmentioned dependency on package named CLASSNAME
##         initializerCodeList <- c(initializerCodeList,
##                                  parse(text = paste0("push_back(", cppDef$name, "::allADtapePtrs_, ",
##                                                      cppDef$name, "::", paste0(derivFun, "_callForADtaping_"), "() )"))[[1]])
##         ## initializerCodeList <- c(initializerCodeList, substitute(push_back(CLASSNAME::allADtapePtrs_, CLASSNAME::ADTAPINGNAME() ), list(CLASSNAME = as.name(cppDef$name),ADTAPINGNAME = as.name(paste0(derivFun, "_callForADtaping_")))))
##     }
##     initializerCode <- do.call('call', c('{', initializerCodeList), quote = TRUE)
##     initializerECcode <- RparseTree2ExprClasses(initializerCode)
##     initializerDef$code <- cppCodeBlock(code = initializerECcode, objectDefs = symbolTable())
##     cppClass$functionDefs[['initializer']] <- initializerDef
##     cppClass$globalObjectsDefs[['globals']] <- globalsDef
##     cppClass
## }

## makeADargumentTransferFunction <- function(newFunName = 'arguments2cppad', targetFunDef, independentVarNames, funIndex = 1, parentsSizeAndDims,
##                                            ADconstantsInfo) {
##     ## modeled closely parts of /*  */
##     ## needs to set the ADtapePtr to one element of the ADtape
##     TF <- RCfunctionDef$new() ## should it be static?
##     TF$returnType <- cppVarFull(baseType = 'nimbleCppADinfoClass', ref = TRUE, name = 'RETURN_OBJ')
##     TF$name <- newFunName
##     localVars <- symbolTable() 
##     isNode <- !inherits(parentsSizeAndDims, 'uninitializedField')
##     if(!isNode)
##       TF$args <- targetFunDef$args
##     else{
##       TF$args <- symbolTable()
##       indexNodeInfoSym <- targetFunDef$args$getSymbolObject('ARG1_INDEXEDNODEINFO__')
##       # indexNodeInfoSym$name <-'ARG1_INDEXEDNODEINFO__' ## to conform with original R function indexing
##       TF$args$addSymbol(indexNodeInfoSym)   
##     }
    
##     ## set up index vars (up to 6)
##     indexVarNames <- paste0(letters[9:14],'_')
##     nimbleSymTab <- targetFunDef$RCfunProc$compileInfo$newLocalSymTab

##     ## assign tape ptr code
##     assignTapePtrCode <- substitute(memberData(ADtapeSetup, ADtape) <- allADtapePtrs_[FUNINDEX], list(FUNINDEX = funIndex)) ## This will have to become a unique index in general. -1 added during output
    
##     ## create code to copy from arguments into the independentVars
##     numIndependentVars <- length(independentVarNames)
##     copyIntoIndepVarCode <- vector('list', numIndependentVars+1)
##     ## create the netIncrement_ variable and code to initialize it to 1
##     localVars$addSymbol( cppVar(name = 'netIncrement_', baseType = 'int') )
##     copyIntoIndepVarCode[[1]] <- quote(netIncrement_ <- 1) 
##     totalIndependentLength <- 0
##     subArgIndexedInfo <- function(x){
##       if(deparse(x[[1]])== 'getNodeFunctionIndexedInfo'){
##         x[[2]] <- parse(text = "ARG1_INDEXEDNODEINFO__")[[1]]
##       }
##       return(deparse(x))
##     }
##     maxSize <- 0
##     for(ivn in seq_along(independentVarNames)) {
##         thisName <- independentVarNames[ivn]
##         thisSym <- nimbleSymTab$getSymbolObject(thisName)
##         if(isNode){
##           nameSubList <- targetFunDef$RCfunProc$nameSubList
##           thisName <- names(nameSubList)[sapply(nameSubList, function(x) return(as.character(x) == thisName))]
##           thisModelElementNum <- as.numeric(gsub(".*([0-9]+)$", "\\1", thisName)) ## Extract 1, 2, etc. from end of arg name.
##           thisName <- sub("_[0-9]+$", "", thisName)
##           thisModelName <- paste0('model_', Rname2CppName(thisName)) ## Add model_ at beginning and remove _1, _2, etc. at end of arg name.
##           thisSizeAndDims <- parentsSizeAndDims[[thisName]][[thisModelElementNum]]
##           if(is.null(thisSizeAndDims)){
##             thisConstInfo <- ADconstantsInfo[[thisName]][[thisModelElementNum]]
##             copyIntoIndepVarCode[[ivn+1]] <- substitute({memberData(ADtapeSetup, independentVars)[netIncrement_] <- ARG1_INDEXEDNODEINFO__.info[INT]; netIncrement_ <- netIncrement_ + 1}, list(INT = thisConstInfo$indexColumn)) 
##             totalIndependentLength <- totalIndependentLength + 1
##             next
##           }
##         }
##         if(thisSym$nDim > 0) {
##             thisSizes <- thisSym$size
##             if(isNode){
##               sizeList <- list()
##               for(i in 1:length(thisSizeAndDims$lengths)){
##                 if(thisSizeAndDims$lengths[i] == 1){
##                   if(deparse(thisSizeAndDims$indexExpr[[i]][[1]]) == 'getNodeFunctionIndexedInfo'){
##                     thisSizeAndDims$indexExpr[[i]][[1]] <- parse(text = paste0(
##                       'ARG1_INDEXEDNODEINFO__.info[', thisSizeAndDims$indexExpr[[i]][[3]], ']'))[[1]]
##                   }
##                   sizeList[[i]] <-  list(thisSizeAndDims$indexExpr[[i]][[1]], thisSizeAndDims$indexExpr[[i]][[1]])
##                 }
##                 else{
##                   sizeList[[i]] <-  list(thisSizeAndDims$indexExpr[[i]][[1]], thisSizeAndDims$indexExpr[[i]][[2]])
##                 }
##               }
##             }
##             else{
##               sizeList <- lapply(thisSizes, function(x) list(1, x))
##             }
##             names(sizeList) <- indexVarNames[1:length(sizeList)]
##             if(length(sizeList) > maxSize) maxSize <- length(sizeList)
##             newRcode <- makeCopyingCodeBlock(quote(memberData(ADtapeSetup, independentVars)), as.name(thisName), sizeList, indicesRHS = TRUE, incrementIndex = quote(netIncrement_), isNode)
##             copyIntoIndepVarCode[[ivn+1]] <- newRcode 
##             totalIndependentLength <- totalIndependentLength + prod(thisSizes)
##         } 
##         else {
##           if(isNode){
##             indexBracketInfo <- paste0('[', paste0(sapply(parentsSizeAndDims[[thisName]][[thisModelElementNum]]$indexExpr,
##                                                           subArgIndexedInfo), collapse = ', '),']')
##             indexName <- paste0("cppLiteral('(**", thisModelName, ")')", indexBracketInfo)
##             RHS <- parse(text = substitute(INDEXNAME, list(INDEXNAME = as.name(indexName))))[[1]]
##           }
##           else{
##             RHS <- as.name(thisName)
##           } 
##           copyIntoIndepVarCode[[ivn+1]] <- substitute({memberData(ADtapeSetup, independentVars)[netIncrement_] <- RHS; netIncrement_ <- netIncrement_ + 1}, list(RHS = RHS)) 
##           totalIndependentLength <- totalIndependentLength + 1
##         }
##     }
##     setSizeLine <- substitute(cppMemberFunction(resize(memberData(ADtapeSetup, independentVars), TIL)), list(TIL = totalIndependentLength))
##     returnCall <- cppLiteral("return(ADtapeSetup);")
    
##     if(maxSize > 0){
##       for(ivn in 1:maxSize)
##         localVars$addSymbol( cppVar(name = indexVarNames[ivn], baseType = 'int') )    
##     }
    
##     allRcode <- do.call('call', c(list('{'), list(setSizeLine), list(assignTapePtrCode), copyIntoIndepVarCode, list(returnCall)), quote=TRUE)
##     allCode <- RparseTree2ExprClasses(allRcode)
##     TF$code <- cppCodeBlock(code = allCode, objectDefs = localVars)
##     TF
## }

## 1. Use myADtapePtrs_
## 2. Don't assume declared known lengths.
makeADargumentTransferFunction2 <- function(newFunName = 'arguments2cppad',
                                            targetFunDef,
                                            callForTapingName,
                                            independentVarNames,
                                            funIndex = 1,
                                            parentsSizeAndDims,
                                            ADconstantsInfo,
                                            useModelInfo = list(),
                                            derivControl = list(),
                                            metaTape = FALSE) {
    if(!metaTape) {
      ADtape_independentVarsName <- as.name('independentVars')
      ## ADtape_dynamicVarsName <- as.name('dynamicVars')
      update_dynamicVars_funName <- as.name("update_dynamicVars")
    } else {
      ADtape_independentVarsName <- as.name('independentVars_meta')
      ## ADtape_dynamicVarsName <- as.name('dynamicVars_meta')
      update_dynamicVars_funName <- as.name("update_dynamicVars_meta")
    }
    nodeFxnVector_name <- useModelInfo[['nodeFxnVector_name']]
    usesModelCalculate <- length(nodeFxnVector_name) > 0    ## modeled closely parts of /*  */
    ## needs to set the ADtapePtr to one element of the ADtape
    TF <- RCfunctionDef$new() ## should it be static?
    TF$returnType <- cppVarFull(baseType = 'nimbleCppADinfoClass', ref = TRUE, name = 'RETURN_OBJ')
    if(metaTape)
      TF$template <- cppVarFull(name = character(), baseType = 'template', templateArgs = list('class TYPE_'))

    TF$name <- newFunName
    localVars <- symbolTable() 
    isNode <- !inherits(parentsSizeAndDims, 'uninitializedField')
    if(isNode) warning("makeADargumentTransferFunction_2 not yet updated for isNode = TRUE")
    if(!isNode) {
      TF$args <- targetFunDef$args$copy()
      if(metaTape) {
        symNames <- TF$args$getSymbolNames()
        newTempSymbols <- list()
        for(sn in symNames) {
          oldSym <- TF$args$getSymbolObject(sn)
          newSym <- cppVarSym2templateTypeCppVarSym(oldSym)
          newTempSymbols[[ sn ]] <- oldSym$copy()
          TF$args$addSymbol(newSym, allowReplace = TRUE)
        }
      }
    } else {
      TF$args <- symbolTable()
      indexNodeInfoSym <- targetFunDef$args$getSymbolObject('ARG1_INDEXEDNODEINFO__')
      # indexNodeInfoSym$name <-'ARG1_INDEXEDNODEINFO__' ## to conform with original R function indexing
      TF$args$addSymbol(indexNodeInfoSym)
    }
    
    ## set up index vars (up to 6)
    indexVarNames <- paste0(letters[9:14],'_')
    nimbleSymTab <- targetFunDef$RCfunProc$compileInfo$newLocalSymTab

    if(!metaTape) {
      copyToDoublesLines <- quote(blank())
      callForTapingNames <- TF$args$getSymbolNames()
    } else {
      copyToDoublesLines <- list()
      callForTapingNames <- character()
        for(ivn in seq_along(independentVarNames)) {
          thisName <- independentVarNames[ivn]
          oldSym <- newTempSymbols[[ thisName ]]
          newName <- paste0(thisName, "_double_temp_")
          oldSym$name <- newName
          oldSym$ref <- FALSE
          callForTapingNames <- c(callForTapingNames, newName)
          localVars$addSymbol(oldSym)
          newRline <- "NEWV__.setSize(ORIGV__.getSizeVec(), false, false); copy_CppADdouble_to_double(ORIGV__.getPtr(), ORIGV__.getPtr() + ORIGV__.size(), NEWV__.getPtr());"
          newRline <- gsub("NEWV__", newName, newRline)
          newRline <- gsub("ORIGV__", thisName, newRline)
          newRline <- substitute(cppLiteral(LINE), list(LINE = newRline))
          copyToDoublesLines[[ length(copyToDoublesLines)+1 ]] <- newRline
        }
        copyToDoublesLines <- do.call("call", c(list("{"), copyToDoublesLines))
    }

    if(!isNode) {
        TF$args$addSymbol(cppVar(baseType = "nimbleCppADinfoClass", ref = TRUE, name = "ADinfo"))
    }
    
    ## record tape if needed
    runCallForTapingCode <- do.call('call',
                                    c(list(callForTapingName),
                                      lapply(callForTapingNames, as.name),
                                      quote(ADinfo)),
                                    quote = TRUE)
    recordIfNeededCode <- substitute(
      if(!memberData(ADinfo, ADtape)) {
        COPYTODOUBLESLINES
        memberData(ADinfo, ADtape) <- RUNCALLFORTAPING
      },
      list(RUNCALLFORTAPING = runCallForTapingCode,
           COPYTODOUBLESLINES = copyToDoublesLines))
    ## recordIfNeededCode <- substitute(
    ##   if(!myADtapePtrs_[FUNINDEX]) {
    ##     COPYTODOUBLESLINES
    ##     myADtapePtrs_[FUNINDEX] <- RUNCALLFORTAPING
    ##   },
    ##   list(FUNINDEX = funIndex,
    ##        RUNCALLFORTAPING = runCallForTapingCode,
    ##        COPYTODOUBLESLINES = copyToDoublesLines))
    
    ## assign tape ptr code
    ## assignTapePtrCode <- substitute(memberData(ADtapeSetup, ADtape) <- myADtapePtrs_[FUNINDEX], list(FUNINDEX = funIndex)) ## This will have to become a unique index in general. -1 added during output
    
    ## create code to copy from arguments into the independentVars
    numIndependentVars <- length(independentVarNames)
  copyIntoIndepVarCode <- vector('list', numIndependentVars+1)

     ## create the netIncrement_ variable and code to initialize it to 1
    localVars$addSymbol( cppVar(name = 'netIncrement_', baseType = 'int') )
    copyIntoIndepVarCode[[1]] <- quote(netIncrement_ <- 1) 
 ##   totalIndependentLength <- 0
    subArgIndexedInfo <- function(x){
      if(deparse(x[[1]])== 'getNodeFunctionIndexedInfo'){
        x[[2]] <- parse(text = "ARG1_INDEXEDNODEINFO__")[[1]]
      }
      return(deparse(x))
    }
    maxSize <- 0
    for(ivn in seq_along(independentVarNames)) {
        thisName <- independentVarNames[ivn]
        thisSym <- nimbleSymTab$getSymbolObject(thisName)
        if(isNode){
          nameSubList <- targetFunDef$RCfunProc$nameSubList
          thisName <- names(nameSubList)[sapply(nameSubList, function(x) return(as.character(x) == thisName))]
          thisModelElementNum <- as.numeric(gsub(".*([0-9]+)$", "\\1", thisName)) ## Extract 1, 2, etc. from end of arg name.
          thisName <- sub("_[0-9]+$", "", thisName)
          thisModelName <- paste0('model_', Rname2CppName(thisName)) ## Add model_ at beginning and remove _1, _2, etc. at end of arg name.
          thisSizeAndDims <- parentsSizeAndDims[[thisName]][[thisModelElementNum]]
          if(is.null(thisSizeAndDims)){
            thisConstInfo <- ADconstantsInfo[[thisName]][[thisModelElementNum]]
            copyIntoIndepVarCode[[ivn+1]] <- substitute(
            {
              memberData(ADinfo, IVN)[netIncrement_] <- ARG1_INDEXEDNODEINFO__.info[INT]
##              memberData(ADtapeSetup, IVN)[netIncrement_] <- ARG1_INDEXEDNODEINFO__.info[INT]
              netIncrement_ <- netIncrement_ + 1
            },
            list(INT = thisConstInfo$indexColumn,
                 IVN = ADtape_independentVarsName)) 
            totalIndependentLength <- totalIndependentLength + 1
            next
          }
        }
        if(thisSym$nDim > 0) {
          ## thisSizes <- thisSym$size
            if(isNode){
              sizeList <- list()
              for(i in 1:length(thisSizeAndDims$lengths)){
                if(thisSizeAndDims$lengths[i] == 1){
                  if(deparse(thisSizeAndDims$indexExpr[[i]][[1]]) == 'getNodeFunctionIndexedInfo'){
                    thisSizeAndDims$indexExpr[[i]][[1]] <- parse(text = paste0(
                      'ARG1_INDEXEDNODEINFO__.info[', thisSizeAndDims$indexExpr[[i]][[3]], ']'))[[1]]
                  }
                  sizeList[[i]] <-  list(thisSizeAndDims$indexExpr[[i]][[1]], thisSizeAndDims$indexExpr[[i]][[1]])
                }
                else{
                  sizeList[[i]] <-  list(thisSizeAndDims$indexExpr[[i]][[1]], thisSizeAndDims$indexExpr[[i]][[2]])
                }
              }
            }
            else{
              sizeList <- lapply(1:thisSym$nDim,
                                 function(x) list(1,
                                                  substitute(dim(RHS)[INDEX],
                                                             list(RHS = as.name(thisName),
                                                                  INDEX = x))))
            }
            names(sizeList) <- indexVarNames[1:length(sizeList)]
            if(length(sizeList) > maxSize) maxSize <- length(sizeList)
          newRcode <- makeCopyingCodeBlock(
              substitute(
 ##                 memberData(ADtapeSetup, IVN),
                  memberData(ADinfo, IVN),
                  list(IVN = ADtape_independentVarsName)),
            as.name(thisName),
            sizeList,
            indicesRHS = TRUE,
            incrementIndex = quote(netIncrement_),
            isNode)
          copyIntoIndepVarCode[[ivn+1]] <- newRcode 
##            totalIndependentLength <- totalIndependentLength + prod(thisSizes)
        }
        
        else {
          if(isNode){
            indexBracketInfo <- paste0('[', paste0(sapply(parentsSizeAndDims[[thisName]][[thisModelElementNum]]$indexExpr,
                                                          subArgIndexedInfo), collapse = ', '),']')
            indexName <- paste0("cppLiteral('(**", thisModelName, ")')", indexBracketInfo)
            RHS <- parse(text = substitute(INDEXNAME, list(INDEXNAME = as.name(indexName))))[[1]]
          }
          else{
            RHS <- as.name(thisName)
          } 
          copyIntoIndepVarCode[[ivn+1]] <- substitute(
          {
              ## memberData(ADtapeSetup, IVN)[netIncrement_] <- RHS
              memberData(ADinfo, IVN)[netIncrement_] <- RHS
              netIncrement_ <- netIncrement_ + 1
          },
          list(RHS = RHS, IVN = ADtape_independentVarsName)) 
          ##          totalIndependentLength <- totalIndependentLength + 1
        }
    }
  ##   setSizeLine <- substitute(cppMemberFunction(resize(memberData(ADtapeSetup, independentVars), TIL)), list(TIL = totalIndependentLength))
  localVars$addSymbol( cppVar(name = "totalIndependentVarLength_", baseType = "int"))
  
  calcTotalLengthCode <- makeCalcTotalLengthBlock(independentVarNames,
                                                  nimbleSymTab,
                                                  "totalIndependentVarLength_",
                                                  0) ## from extraInput scheme: if(usesModelCalculate) 1 else 0)
    setSizeLine <- substitute(
      cppMemberFunction(resize(
        ## memberData(ADtapeSetup, IVN),
        memberData(ADinfo, IVN),
        totalIndependentVarLength_)),
      list(IVN = ADtape_independentVarsName))
##    returnCall <- cppLiteral("return(ADtapeSetup);")
    returnCall <- cppLiteral("return(ADinfo);")
    if(maxSize > 0){
      for(ivn in 1:maxSize)
        localVars$addSymbol( cppVar(name = indexVarNames[ivn], baseType = 'int') )    
    }

  dynamicVarsLine <- if(usesModelCalculate) {
    substitute(UDV(ADinfo),
               list(UDV = update_dynamicVars_funName))
  } else {
    quote(blank())
  }
    
    if(!metaTape) {
      setMetaFlagLine <- quote(blank())
    } else {
      ## setMetaFlagLine <- cppLiteral("ADtapeSetup.metaFlag = true;")
      setMetaFlagLine <- cppLiteral("ADinfo.metaFlag = true;")
    }
    
  allRcode <- do.call('call', c(list('{'),
                                list(recordIfNeededCode),
                                calcTotalLengthCode,
                                list(setSizeLine),
##                                list(assignTapePtrCode),
                                copyIntoIndepVarCode,
                                list(dynamicVarsLine, setMetaFlagLine),
                                list(returnCall)),
                      quote=TRUE)
    allCode <- RparseTree2ExprClasses(allRcode)
    TF$code <- cppCodeBlock(code = allCode, objectDefs = localVars)
    TF
}

## Generate a code block to determine the total length of a bunch of variables that may be
## scalar or non-scalar.
makeCalcTotalLengthBlock <- function(independentVarNames,
                                    symTab,
                                    totalLengthName,
                                    initLength = 0) {
  numIndependentVars <- length(independentVarNames)
  calcTotalLengthLines <- vector('list', numIndependentVars + 1)
  calcTotalLengthLines[[1]] <- substitute(TOTALLENGTH <- INITLENGTH,  ## Potential for extraInputDummy dimension here.
                                          list(TOTALLENGTH = as.name(totalLengthName),
                                               INITLENGTH = as.numeric(initLength))) 
  for(ivn in seq_along(independentVarNames)) {
      thisName <- independentVarNames[ivn]
      thisSym <- symTab$getSymbolObject(thisName)
    if(thisSym$nDim > 0) 
      calcTotalLengthLines[[ivn+1]] <- substitute(cppLiteral(LITERAL_CODE),
                                                  list(LITERAL_CODE =
                                                         paste0(totalLengthName, " += ", thisName, ".size();"))) 
    else
      calcTotalLengthLines[[ivn+1]] <- substitute(cppLiteral(LITERAL_CODE),
                                                  list(LITERAL_CODE =
                                                         paste0("++",totalLengthName,";")))
  }
  calcTotalLengthLines
}

## Generate a block of code for copying to or from CppAD objects, to or from original C++ objects
## On the CppAD side, we are always flattening to 1D.
##
## The code this generates is embedded in the ADtapingFunction made by makeADtapingFunction
##
## Note this does some work similar to BUGScontextClass::embedCodeInForLoop
makeCopyingCodeBlock <- function(LHSvar, 
                                 RHSvar, 
                                 indexList, 
                                 indicesRHS = TRUE,
                                 incrementIndex, 
                                 isNode = FALSE) {
  indexNames <- rev(names(indexList))
  indexedBracketExpr <- do.call('call', c(list('[', as.name('TO_BE_REPLACED')),
                                          lapply(rev(indexNames), as.name)), ## use rev() to force column-major order for results
                                quote = TRUE)
  if(indicesRHS) {
    if(isNode)
      RHS <- eval(
        substitute(
          substitute(
            indexedBracketExpr, 
            list(TO_BE_REPLACED = cppLiteral(paste0('(**model_', deparse(RHSvar), ')')))),
          list(indexedBracketExpr = indexedBracketExpr)
        ))
    else 
      RHS <- eval(
        substitute(
          substitute(
            indexedBracketExpr, 
            list(TO_BE_REPLACED = RHSvar)),
          list(indexedBracketExpr = indexedBracketExpr)
        ))
    LHS <- substitute(A[i], 
                      list(A = LHSvar,
                           i = incrementIndex))
  } else {
    LHS <- eval(
      substitute(
        substitute(
          indexedBracketExpr, 
          list(TO_BE_REPLACED = LHSvar)),
        list(indexedBracketExpr = indexedBracketExpr)
      ))
    RHS <- substitute(A[i],
                      list(A = RHSvar, 
                           i = incrementIndex))
  }
  innerCode <- substitute(
    {
      LHS <- RHS
      cppLiteral(incrementCode)
    },
    list(LHS = LHS,
         RHS = RHS,
         incrementCode = paste0("++", incrementIndex,";")))
##         incrementIndex = incrementIndex))
  for(i in length(indexList):1) {
    newForLoop <- 
      substitute(
        for(NEWINDEX_ in NEWSTART_:NEWEND_) INNERCODE, 
        list(NEWINDEX_ = as.name(indexNames[i]),
             NEWSTART_ = rev(indexList)[[i]][[1]],
             NEWEND_ = rev(indexList)[[i]][[2]],
             INNERCODE = innerCode))
    innerCode <- newForLoop
  }
  innerCode
}
