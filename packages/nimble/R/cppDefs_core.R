## This file contains reference classes used for C++ units such as
## cppDefinition (base class for all others)
## cppNamespace (base class for classes)
## cppClass
##

## Base class for C++ file and compilation information
cppDefinition <- setRefClass('cppDefinition', 
                             fields = list(
                                 filename = 'ANY',	#'character',  ## what filename (to which .h and .cpp will be appended) is this definition in
                                 CPPusings = 'ANY',	#''character',
                                 neededTypeDefs = 'ANY',	#''list',
                                 Hincludes = 'list',
                                 CPPincludes = 'list',
                                 nimbleProject = 'ANY'),  
                             methods = list(
                                 initialize = function(..., project) {
                                 	filename <<- character()
                                 	if(!is.character(CPPusings))	CPPusings <<- character()
                                 	if(!is.list(neededTypeDefs))	neededTypeDefs <<- list()
                                     nimbleProject <<- if(missing(project)) NULL else project
                                     callSuper(...)
                                 },
                                 getHincludes = function() {return(Hincludes)},
                                 getCPPincludes = function() {return(CPPincludes)},
                                 getCPPusings = function() {return(CPPusings)},
                                 getDefs = function() {return(list(.self))} ## return all objects to be included.  This allows adjunct objects like SEXPinterfaceFuns to be included
                                 )
                             )


cppGlobalObjects <- setRefClass('cppGlobalObject',
                                contains = 'cppDefinition',
                                field = list(
                                    name = 'ANY',
                                    objectDefs = 'ANY', ## a list or symbolTable
                                    staticMembers = 'ANY'
                                ),
                                methods = list(
                                    initialize = function(...) {name <<- character(); objectDefs <<-list(); staticMembers <<- FALSE; callSuper(...)},
                                    generate = function(declaration = FALSE) {
                                        ## For globals for class static members, we want no declaration
                                        ## Otherwise we want a declaration and put extern in front
                                        if(staticMembers & declaration) return(character())
                                        objectDefsToUse <- if(inherits(objectDefs, 'symbolTable')) objectDefs$symbols else objectDefs
                                        output <- paste0(generateAll(objectDefsToUse, declaration = declaration),';')
                                        if(declaration) output <- paste("extern ", output)
                                        output
                                    }
                                    ))
## class for C++ namespaces.
## A namespace includes, objects, classes, functions, typedefs, and other namespaces
## This is incomplete.  The typeDefs and nested namespaces are not used yet.
## The other components are used and so have been more developed.
cppNamespace <- setRefClass('cppNamespace',
                            contains = 'cppDefinition',
                            fields = list(
                                name = 'ANY',	    ## character,
                                objectDefs = 'ANY', ## list or symbolTable
                                functionDefs = 'ANY'),
                            methods = list(
                                initialize = function(...) {name <<- character();functionDefs <<- list(); objectDefs <<- list(); callSuper(...)}, ## By default a list, but can be a symbolTable
                                addObject = function(newName, newObj) objectDefs[[newName]] <<- newObj,
                                addFunction = function(newName, newFun) functionDefs[[newName]] <<- newFun,
                                generate = function() {
                                    objectDefsToUse <- if(inherits(objectDefs, 'symbolTable')) objectDefs$symbols else objectDefs
                                    output <- c(generateNameSpaceHeader(name$generate()),
                                                generateObjectDefs(objectDefsToUse),
                                                generateAll(functionDefs, declaration = TRUE),
                                                '};'
                                                )
                                }
                                )
                            )


nimbleCppInheritanceInfo <- list(
    Values = c('NamedObjects'),
    ModelBase = c('NamedObjects'),
    nodeFun = c('NamedObjects'),
    derivNodeFun = c('nodeFun')
)

## C++ class object.
## A class is like a namespace with inheritance
## At the moment everything is public. 
## This class can build cppFunction objects for a generator function and a finalizer function
## The generator can be called via .Call to return an external pointer to a new object of the class
## The finalizer is the finalizer assigned to the object when the external pointer is made
cppClassDef <- setRefClass('cppClassDef',
                           contains = 'cppNamespace',
                           fields = list(
                               inheritance = 'list',           ## classes to be declared as public 
                               ancestors = 'list',             ## classes inherited by inherited classes, needed to make all cast pointers
                               extPtrTypes = 'ANY',
                               private = 'list',		# 'list'. This field is a placeholder for future functionality.  Currently everything is generated as public
                               useGenerator = 'ANY',		#'logical',
                               SEXPgeneratorFun = 'ANY', ## These will be cppFunctionDefs.
                               SEXPfinalizerFun = 'ANY',
                               globalObjectsDefs = 'ANY'
                           ),
                           methods = list(
                               initialize = function(...) {
                                   useGenerator <<- TRUE
                                   globalObjectsDefs <<- list()
                                   Hincludes <<-	c(Hincludes, '<Rinternals.h>')	
                                   CPPincludes <<-	c(CPPincludes, '<iostream>') 
                                   callSuper(...)
                               },
                               getHincludes = function() {
                                   Hinc <- c(Hincludes,
                                             if(!inherits(SEXPgeneratorFun, 'uninitializedField')) SEXPgeneratorFun$getHincludes(),
                                             if(!inherits(SEXPfinalizerFun, 'uninitializedField')) SEXPfinalizerFun$getHincludes(),
                                             unlist(lapply(functionDefs, function(x) x$getHincludes()), recursive = FALSE))
                                   Hinc
                               },
                               getCPPincludes = function() {
                                   CPPinc <- c(CPPincludes,
                                               if(!inherits(SEXPgeneratorFun, 'uninitializedField')) SEXPgeneratorFun$getCPPincludes(),
                                               if(!inherits(SEXPfinalizerFun, 'uninitializedField')) SEXPfinalizerFun$getCPPincludes(),
                                               unlist(lapply(functionDefs, function(x) x$getCPPincludes()), recursive = FALSE))
                                   CPPinc
                               },
                               getCPPusings = function() {
                                   CPPuse <- unique(c(CPPusings,
                                                      if(!inherits(SEXPgeneratorFun, 'uninitializedField')) SEXPgeneratorFun$getCPPusings(),
                                                      if(!inherits(SEXPfinalizerFun, 'uninitializedField')) SEXPfinalizerFun$getCPPusings(),
                                                      unlist(lapply(functionDefs, function(x) x$getCPPusings()))))
                                   CPPuse
                               },
                               getDefs = function() {
                                   ans <- if(useGenerator) {
                                       if(inherits(SEXPgeneratorFun, 'uninitializedField')) stop('Trying to getDefs from a CppClassDef with useGenerator==TRUE but SEXPgeneratorFun not defined')
                                       if(inherits(SEXPfinalizerFun, 'uninitializedField'))
                                           list(.self, SEXPgeneratorFun)
                                       else
                                           list(.self, SEXPgeneratorFun, SEXPfinalizerFun)
                                   } else {
                                       list(.self)
                                   }
                                   if(length(globalObjectsDefs) > 0) ans <- c(ans, globalObjectsDefs)
                                   ans
                               },
                               addInheritance = function(newI) inheritance <<- c(inheritance, newI),
                               addAncestors = function(newI) ancestors <<- c(ancestors, newI), 
                               setPrivate = function(name) private[[name]] <<- TRUE,
                               generate = function(declaration = TRUE, definition = FALSE, ...) {
                                   if(declaration) {
                                       objectDefsToUse <- if(inherits(objectDefs, 'symbolTable')) objectDefs$symbols else objectDefs
                                       output <- c(generateClassHeader(name, inheritance),
                                                   list('public:'), ## In the future we can separate public and private
                                                   lapply(generateObjectDefs(objectDefsToUse), function(x) if(length(x)==0) '' else pasteSemicolon(x, indent = '  ')),
                                                   generateAll(functionDefs, declaration = TRUE),
                                                   '};'
                                               )
                                   } else {
                                       output <- generateAll(functionDefs, scopes = name)
                                   }
                                   output
                               },
                               getExtPtrTypeIndex = function() {
                                   structure(1:(length(extPtrTypes)+1), names = c(name, extPtrTypes))
                               },
                               buildSEXPgenerator = function(finalizer = NULL) { ## build a function that will provide a new object and return an external pointer
                                   addinheritance <- nimbleCppInheritanceInfo[c(unlist(inheritance), unlist(ancestors))]
                                   allinheritance <- unique(unlist(c(inheritance, ancestors, addinheritance)))
                                   extPtrTypes <<- allinheritance
                                   numBaseClasses <- length(allinheritance)
                                   SallinheritanceNames <- paste0('S', allinheritance, '_EXTPTR')
                                   baseClassPtrObjectDefs <- lapply(SallinheritanceNames, cppSEXP) 
                                   CBobjectDefs <- list(cppVar(name = 'newObj', baseType = name, ptr = 1),
                                                        Sans = cppSEXP(name = 'Sans'),
                                                        Sderived = cppSEXP(name = 'Sderived_EXTPTR'))
                                   CBobjectDefs <- c(CBobjectDefs, baseClassPtrObjectDefs)
                                   newCodeLine <- cppLiteral(c(paste0('newObj = new ', name,';'),
                                                               'Sderived_EXTPTR = PROTECT(R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));'))
                                   baseClassCastLines <- mapply(
                                       function(baseClassName, Sname)
                                           paste0('PROTECT(',Sname,'= R_MakeExternalPtr(dynamic_cast<',baseClassName,'*>(newObj), R_NilValue, R_NilValue));'),
                                       baseClassName = allinheritance,
                                       Sname = SallinheritanceNames)
                                   baseClassCastLines <- cppLiteral(baseClassCastLines)
                                   allocVectorLine <- cppLiteral(paste0('Sans = PROTECT(Rf_allocVector(VECSXP,',numBaseClasses + 1, '));'))

                                   packListLines <- mapply(function(Sname, i) paste0('SET_VECTOR_ELT(Sans,',i,',',Sname,');'),
                                                           Sname = c('Sderived_EXTPTR', SallinheritanceNames),
                                                           i = 0:(length(SallinheritanceNames)))
                                   packListLines <- cppLiteral(packListLines)
                                   
                                   notificationLine <- if(getNimbleOption('messagesWhenBuildingOrFinalizingCppObjects'))
                                                           paste0('std::cout<< \"In generator for ', name, '. Created at pointer \" << newObj << \"\\n\";')
                                                       else character(0)
                                   if(is.null(finalizer)) finalizer <- paste0(name,'_Finalizer')
                                   codeLines <- substitute({
                                       ## Finalizer registration now happens through nimble's finalizer mapping system.
                                       UNPROTECT(NUMPROT)
                                       return(Sans)
                                   }, list(TYPE = as.name(name), FINALIZER = as.name(finalizer), NUMPROT = length(allinheritance) + 2))
                                   allCodeList <- list(newCodeLine, cppLiteral(notificationLine),  baseClassCastLines, allocVectorLine, packListLines, codeLines)
                                   allCode <- putCodeLinesInBrackets(allCodeList)
                                   SEXPgeneratorFun <<- cppFunctionDef(name = paste0('new_',name),
                                                                       args = list(),
                                                                       code = cppCodeBlock(code = allCode, objectDefs = CBobjectDefs, skipBrackets = TRUE),
                                                                       returnType = cppSEXP(),
                                                                       externC = TRUE)
                               },
                               buildSEXPfinalizer = function() {
                                   CBobjectDefs <- list(cppVar(name = 'oldObj', baseType = name, ptr = 1))
                                   inputArgs <- list(cppSEXP(name = 'Sv'))
                                   notificationLine <- if(getNimbleOption('messagesWhenBuildingOrFinalizingCppObjects'))
                                       paste0('std::cout<< \"In finalizer for ', name, ' with pointer \" << R_ExternalPtrAddr(Sv) << \"\\n\";')
                                   else character(0)
                                   castLine <- paste0('oldObj = static_cast<',name,' *>(R_ExternalPtrAddr(Sv));')
                                   deleteLine <- c('if(oldObj) delete oldObj;', 'R_ClearExternalPtr(Sv);')
                                   code <- putCodeLinesInBrackets(list(cppLiteral(c(notificationLine, castLine, deleteLine))))
                                   SEXPfinalizerFun <<- cppFunctionDef(
                                       name = paste0(name,'_Finalizer'),
                                       args = inputArgs,
                                       returnType = cppVoid(),
                                       code = cppCodeBlock(code = code, objectDefs = CBobjectDefs, skipBrackets = TRUE)
                                       )
                               }
                               )
                           )

stripUnusedTypeDefs <- function(cppOutput) {
  cppOrigOutput <- cppOutput
  ans <- try({
    cppOutput <- unlist(cppOutput)
    typeDefLines <- grep("^typedef", cppOutput)
    lines_to_remove <- integer()
    for(i in typeDefLines) {
      pieces <- strsplit(cppOutput[typeDefLines[i]], ' ')[[1]]
      typedef_name <- pieces[length(pieces)]
      typedef_name <- gsub(";","",typedef_name)
      i_typedef_name <- grep(typedef_name, cppOutput)
      if(length(i_typedef_name)==1 && i_typedef_name[1]==i) {
        lines_to_remove <- c(lines_to_remove, i)
      }
    }
    if(length(lines_to_remove)>0) cppOutput<-cppOutput[-lines_to_remove]
    cppOutput
  })
  if(inherits(ans, 'try-error')) cppOrigOutput else ans
}

## A cppCodeBlock is an arbitrary collection of parse tree and other cppCodeBlocks (defined below)
## The parse tree can be either an R parse tree or one of our exprClass objects
cppCodeBlock <- setRefClass('cppCodeBlock',
                            fields = list(typeDefs = 'ANY', objectDefs = 'ANY', code = 'ANY', skipBrackets = 'ANY',
                                          cppADCode = 'ANY', generatorSymTab = 'ANY'),
                            methods = list(
                                generate = function(indent = '', ...) {
                                    ## TRUE was an original possible value.
                                    ## It was deprecated for the new value, 2L.
                                    ## This originally was version 2. Now it is the version.
                                    ## if(isTRUE(cppADCode)){
                                    ##     oldCppADCode <- nimbleUserNamespace$cppADCode
                                    ##     nimbleUserNamespace$cppADCode <- TRUE
                                    ##     on.exit(nimbleUserNamespace$cppADCode <- oldCppADCode)
                                    ## }
                                    if(identical(cppADCode, 2L)){
                                        oldCppADCode <- nimbleUserNamespace$cppADCode
                                        nimbleUserNamespace$cppADCode <- 2L
                                        on.exit(nimbleUserNamespace$cppADCode <- oldCppADCode)
                                    }
                                    if(inherits(typeDefs, 'uninitializedField')) typeDefs <<- list()
                                    typeDefsToUse <- if(inherits(typeDefs, 'symbolTable')) typeDefs$symbols else typeDefs
                                    if(length(typeDefsToUse) > 0) {
                                        outputCppCode <- paste0(indent, generateObjectDefs(typeDefsToUse),';')
                                    } else outputCppCode <- list()
                                    if(inherits(objectDefs, 'uninitializedField')) objectDefs <<- list()
                                    objectDefsToUse <- if(inherits(objectDefs, 'symbolTable')) objectDefs$symbols else objectDefs
                                    if(length(objectDefsToUse) > 0) {
                                        outputCppCode <- c(outputCppCode, paste0(indent, generateObjectDefs(objectDefsToUse),';'))
                                    } 

                                    if(inherits(code, 'exprClass')) {
                                        if(inherits(generatorSymTab, 'symbolTable')) useSymTab <- generatorSymTab
                                        else if(!inherits(objectDefs, 'symbolTable')) stop('Error, with exprClass code in the cppCodeBlock, must have objectDefs be a symbolTable')
                                        else useSymTab <- objectDefs
                                        outputCppCode <- c(outputCppCode, nimGenerateCpp(code, useSymTab, indent = ' ', showBracket = FALSE))
                                    } else {
                                        outputCppCode <- c(outputCppCode, outputCppParseTree2(code, indent))
                                    }
                                    if(isTRUE(getNimbleOption("stripUnusedTypeDefs")))
                                       outputCppCode <- stripUnusedTypeDefs(outputCppCode)
                                    outputCppCode
                                  }
                                )
                            )

## C++ function definitions
##
cppFunctionDef <- setRefClass('cppFunctionDef',
                              contains = 'cppDefinition',
                              fields = list(name = 'ANY',	#'character',
                                  returnType = 'ANY',	#'cppVar', 
                                  args = 'ANY',  # list
                                  code = 'ANY',	#	'cppCodeBlock',
                                  externC = 'ANY',  # logical
                                  virtual = 'ANY', ## logical, only relevant for class members
                                  static = 'ANY',  ## ditto
                                  abstract = 'ANY',  ## ditto
                                  template = 'ANY',
                                  const = 'ANY'
                                            ),
                              methods = list(
                                  initialize = function(...) {
                                      name <<- character()
                                      CPPincludes <<- as.list( c(CPPincludes, '<iostream>') )
                                      callSuper(...)
                                      if(inherits(virtual, 'uninitializedField')) virtual <<- FALSE
                                      if(inherits(abstract, 'uninitializedField')) abstract <<- FALSE
                                      if(inherits(template, 'uninitializedField')) template <<- NULL
                                      if(inherits(static, 'uninitializedField')) static <<- FALSE                                          
                                      if(inherits(const, 'uninitializedField')) const <<- FALSE
                                  },
                                  generate = function(declaration = FALSE, scopes = character(), ...) {
                                      if(inherits(args, 'uninitializedField')) args <<- list()
                                      argsToUse <- if(inherits(args, 'symbolTable')) args$symbols else args
                                      if(declaration) {
                                          outputCode <- paste0(if(virtual) 'virtual ' else character(0),
                                                               generateFunctionHeader(returnType, name, argsToUse, scopes, template, static, ...),
                                                               if(const) ' const ' else character(0),
                                                               if(abstract) '= 0' else character(0), ';')
                                          if(!inherits(externC, 'uninitializedField' ) ){
                                            if(externC == TRUE)
                                              outputCode <- paste0('extern "C" ', outputCode)
                                          }
                                           return(outputCode) 
                                      } else {
                                          code_is_empty <- inherits(code$code, 'uninitializedField')
                                          if(code_is_empty) {
                                              ## There is no code. This can occur for a nimbleFunctionVirtual, which is an abstract base class.
                                              if(abstract)
                                                  return(character(0))
                                          }
                                          c(paste(generateFunctionHeader(returnType, name, argsToUse, scopes, template, static = FALSE, ...),
                                                  if(const) ' const ' else character(0), '{'),
                                            if(!code_is_empty) code$generate(...) else '',
                                            list('}'))
                                      }
                                  }
                                  )
                              )

generateFunctionHeader <- function(returnType, name, args, scopes = character(), template = NULL, static = FALSE) {
    list(paste(paste0(if(is.null(template)) character() else paste(template$generate(),'\n'),
               paste0(if(static) 'static ' else character(), returnType$generate(printName = character()), collapse = '')),
               paste(c(scopes, name), collapse = '::'),
               '(',
               paste(unlist(lapply(args, function(x) x$generate())), collapse = ', '),
               ')'))
}

generateClassHeader <- function(ns, inheritance) {
    inheritancePart <- if(length(inheritance) > 0) {
        paste(':', paste('public', unlist(inheritance), collapse = ', '))
    } else NULL
    list(paste('class', ns, inheritancePart, '{'))
}

generateObjectDefs <- function(objectDefs, ...) {
    generateAll(objectDefs, ...)
}
