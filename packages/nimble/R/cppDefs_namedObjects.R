## C++ definition for a namedObjects class
## This is a permanaent C++ base class to provide one interface
## to get external pointers to class member data by character names

cppNamedObjectsClass <- setRefClass('cppNamedObjectsClass',
                               contains = 'cppClassDef',
                               fields = list(Rnames2CppNames = 'ANY'),
                               methods = list(
                                   initialize = function(...) {
                                   		Rnames2CppNames <<- list()
                                       inheritance <<- c(inheritance, list('NamedObjects'))
                                       Hincludes <<- c(Hincludes, nimbleIncludeFile("NamedObjects.h"))
                                       callSuper(...)
                                   },
                                   namedObjectsConstructorCodeBlock = function() {
                                       template <- quote( {namedObjects[RNAME] = cppReference(CNAME)} )[[2]]
                                       codeLines <- vector('list', length = length(Rnames2CppNames))
                                       for(i in seq_along(Rnames2CppNames))
                                           codeLines[[i]] <- codeSubstitute(template, list(RNAME = names(Rnames2CppNames)[i],
                                                                                           CNAME = as.name(Rnames2CppNames[[i]])))
                                       cppCodeBlock(code = putCodeLinesInBrackets(codeLines), skipBrackets = TRUE)
                                   }##,
                                   ## namedObjects has its own version of this so it can return both a base class (namedObjects) ptr and a derived class ptr
                                   ## buildSEXPgenerator = function(finalizer = NULL) { ## build a function that will provide a new object and return an external pointer
                                   ##     CBobjectDefs <- list(cppVar(name = 'newObj', baseType = name, ptr = 1),
                                   ##                          Sans = cppSEXP(name = 'Sans'),
                                   ##                          SextPtrBase = cppSEXP(name = 'SextPtrBase'),
                                   ##                          SextPtrDerived = cppSEXP(name = 'SextPtrDerived'));
                                   ##     newCodeLine <- cppLiteral(c(paste0('newObj = new ', name,';'),
                                   ##                                 'PROTECT(SextPtrDerived = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));',
                                   ##                                 'PROTECT(SextPtrBase = R_MakeExternalPtr(dynamic_cast<NamedObjects*>(newObj), R_NilValue, R_NilValue));',
                                   ##                                 'PROTECT(Sans = allocVector(VECSXP, 2));',
                                   ##                                 'SET_VECTOR_ELT(Sans, 0, SextPtrDerived);',
                                   ##                                 'SET_VECTOR_ELT(Sans, 1, SextPtrBase);'))
                                   ##     notificationLine <- if(nimbleOptions()$messagesWhenBuildingOrFinalizingCppObjects)
                                   ##                             paste0('std::cout<< \"In generator for ', name, '. Created at pointer \" << R_ExternalPtrAddr(Sans) << \"\\n\";')
                                   ##                         else character(0)
                                   ##     if(is.null(finalizer)) finalizer <- paste0(name,'_Finalizer')
                                   ##     codeLines <- substitute({
                                   ##         UNPROTECT(3)
                                   ##         return(Sans)
                                   ##     }, list(TYPE = as.name(name), FINALIZER = as.name(finalizer)))
                                   ##     allCode <- putCodeLinesInBrackets(list(newCodeLine, cppLiteral(notificationLine), codeLines))
                                   ##     SEXPgeneratorFun <<- cppFunctionDef(name = paste0('new_',name),
                                   ##                                         args = list(),
                                   ##                                         code = cppCodeBlock(code = allCode, objectDefs = CBobjectDefs, skipBrackets = TRUE),
                                   ##                                         returnType = cppSEXP(),
                                   ##                                         externC = TRUE)
                                   ## }
                                   )
                               )
