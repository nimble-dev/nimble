## This class represents a C++ BUGS model class
cppBUGSmodelClass <- setRefClass('cppBUGSmodelClass',
                                 contains = 'cppNamedObjectsClass',
                                 fields = list(
                                     modelDef = 'ANY',               ## modelDefClass object
                                     model = 'ANY',                  ## RmodelBaseClass object, needed for the nodeFunctions
                                     nodeFuns = 'ANY',		     ## list of cppNimbleFunctionClass objects
                                     CmodelValuesClassName = 'ANY'   ## character
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
                                         addObject('defaultModelValues_', cppVar(name = 'defaultModelValues_', baseType = CmodelValuesClassName))
                                     },
                                     ## buildConstructorFunctionDef adds a cppFunctionDef object for the constructor to the functionDefs field (inherited from cppNamespaceClass via cppNamedObjectsClass).  Most of this comprises the namedObjects assignments generated from  namedObjectsConstructorCodeBlock() inherited from cppNamedObjectsClass
                                     buildConstructorFunctionDef = function() {

                                         lines12 <- cppLiteral(c("defaultModelValues_.resize(1);",
                                                                 "pointAtAll(&defaultModelValues_, 0);",
                                                                 "modelValues_ = static_cast<Values *>(&defaultModelValues_);")) 
                                         
                                         code <- putCodeLinesInBrackets(list(lines12, namedObjectsConstructorCodeBlock())) 
                                         conFunDef <- cppFunctionDef(name = name,
                                                                     returnType = emptyTypeInfo(),
                                                                     code = cppCodeBlock(code = code, skipBrackets = TRUE))
                                         functionDefs[['constructor']] <<- conFunDef
                                     },
                                     ## buildPointAtAll adds a cppFunctionDef for pointAtAll to the functionDefs field.  This consists of lines like one__over_sqrt_tau__between = &(values->one__over_sqrt_tau__between_Vec[i]);  These point each model variable at the corresponding element of a modelValues object.
                                     buildPointAtAll = function() {
                                         template = quote( {NAME = cppReference( values_ %->% VECNAME[i_] ) } )[[2]]

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
                                                               
                                                                  args = list(cppVar(ptr = 1, baseType = CmodelValuesClassName, name = 'values_'),
                                                                      cppInt('i_')),
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
                                         
                                     },
                                    
                                     buildAll = function(buildNodeDefs = TRUE, where = globalenv(), ...) {
                                         makeCppNames() 
                                         buildVars()
                                         buildConstructorFunctionDef()
                                         buildSEXPgenerator(finalizer = 'namedObjects_Finalizer')
                                         buildPointAtAll()
                                         if(buildNodeDefs) buildNodes(where = where, debugCpp = debugCpp)
                                     }
                                     )
                                 )

