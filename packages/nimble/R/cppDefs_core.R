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


## class for C++ namespaces.
## A namespace includes, objects, classes, functions, typedefs, and other namespaces
## This is incomplete.  The typeDefs and nested namespaces are not used yet.
## The other components are used and so have been more developed.
cppNamespace <- setRefClass('cppNamespace',
                            contains = 'cppDefinition',
                            fields = list(
                                name = 'ANY',	#'character',
                                objectDefs = 'ANY', ## This one must be ANY because it could be a list or a symbolTable
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

## C++ class object.
## A class is like a namespace with inheritance
## At the moment everything is public. 
## This class can build cppFunction objects for a generator function and a finalizer function
## The generator can be called via .Call to return an external pointer to a new object of the class
## The finalizer is the finalizer assigned to the object when the external pointer is made
cppClassDef <- setRefClass('cppClassDef',
                           contains = 'cppNamespace',
                           fields = list(
                               inheritance = 'list',			#'list',
                               private = 'list',		#'list', ## a placeholder.  everything is public
                               useGenerator = 'ANY',		#'logical',
                               SEXPgeneratorFun = 'ANY', ## These will be cppFunctionDefs
                               SEXPfinalizerFun = 'ANY'),
                           methods = list(
                               initialize = function(...) {
                                   useGenerator <<- TRUE
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
                                   if(useGenerator) {
                                       if(inherits(SEXPfinalizerFun, 'uninitializedField'))
                                           list(.self, SEXPgeneratorFun)
                                       else
                                           list(.self, SEXPgeneratorFun, SEXPfinalizerFun)
                                   }
                                   else
                                       list(.self)
                               },
                               addInheritance = function(newI) inheritance <<- c(inheritance, newI),
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
                               buildSEXPgenerator = function(finalizer = NULL) { ## build a function that will provide a new object and return an external pointer
                                   CBobjectDefs <- list(cppVar(name = 'newObj', baseType = name, ptr = 1),
                                                      Sans = cppSEXP(name = 'Sans'));
                                   newCodeLine <- cppLiteral(c(paste0('newObj = new ', name,';'), 'PROTECT(Sans = R_MakeExternalPtr(newObj, R_NilValue, R_NilValue));'))
                                   notificationLine <- if(nimbleOptions()$messagesWhenBuildingOrFinalizingCppObjects)
                                                           paste0('std::cout<< \"In generator for ', name, '. Created at pointer \" << R_ExternalPtrAddr(Sans) << \"\\n\";')
                                                       else character(0)
                                   if(is.null(finalizer)) finalizer <- paste0(name,'_Finalizer')
                                   codeLines <- substitute({
                                       ## R_RegisterCFinalizerEx(Sans, cppReference(FINALIZER), FALSE) ## last argument is whether to call finalizer upon R exit.  If TRUE we can generate segfaults if the library has been dyn.unloaded, unfortunately
                                       UNPROTECT(1)
                                       return(Sans)
                                   }, list(TYPE = as.name(name), FINALIZER = as.name(finalizer)))
                                   allCode <- putCodeLinesInBrackets(list(newCodeLine, cppLiteral(notificationLine), codeLines))
                                   SEXPgeneratorFun <<- cppFunctionDef(name = paste0('new_',name),
                                                                       args = list(),
                                                                       code = cppCodeBlock(code = allCode, objectDefs = CBobjectDefs, skipBrackets = TRUE),
                                                                       returnType = cppSEXP(),
                                                                       externC = TRUE)
                               },
                               buildSEXPfinalizer = function() {
                                   CBobjectDefs <- list(cppVar(name = 'oldObj', baseType = name, ptr = 1))
                                   inputArgs <- list(cppSEXP(name = 'Sv'))
                                   notificationLine <- if(nimbleOptions()$messagesWhenBuildingOrFinalizingCppObjects)
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


## A cppCodeBlock is an arbitrary collection of parse tree and other cppCodeBlocks (defined below)
## The parse tree can be either an R parse tree or one of our exprClass objects
cppCodeBlock <- setRefClass('cppCodeBlock',
                            fields = list(objectDefs = 'ANY', code = 'ANY', skipBrackets = 'ANY'),			#'logical'),
                            methods = list(
                                generate = function(indent = '', ...) {
                                    if(inherits(objectDefs, 'uninitializedField')) objectDefs <<- list()
                                    objectDefsToUse <- if(inherits(objectDefs, 'symbolTable')) objectDefs$symbols else objectDefs
                                    if(length(objectDefsToUse) > 0) {
                                        outputCppCode <- paste0(indent, generateObjectDefs(objectDefsToUse),';')
                                    } else outputCppCode <- list()

                                    if(inherits(code, 'exprClass')) {
                                        if(!inherits(objectDefs, 'symbolTable')) stop('Error, with exprClass code in the cppCodeBlock, must have objectDefs be a symbolTable')
                                        outputCppCode <- c(outputCppCode, nimGenerateCpp(code, objectDefs, indent = ' ', showBracket = FALSE))
                                    } else {
                                        outputCppCode <- c(outputCppCode, outputCppParseTree2(code, indent))
                                    }
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
                                  args = 'ANY', 
                                  code = 'ANY',	#	'cppCodeBlock',
                                  externC = 'ANY',
                                  virtual = 'ANY',  ## only relevant for class members
                                  abstract = 'ANY', ## ditto
                                  const = 'ANY'     ## ditto
                                            ),
                              methods = list(
                                  initialize = function(...) {
                                      name <<- character()
                                      CPPincludes <<- as.list( c(CPPincludes, '<iostream>') )
                                      callSuper(...)
                                      if(inherits(virtual, 'uninitializedField')) virtual <<- FALSE
                                      if(inherits(abstract, 'uninitializedField')) abstract <<- FALSE
                                      if(inherits(const, 'uninitializedField')) const <<- FALSE
                                  },
                                  generate = function(declaration = FALSE, scopes = character(), ...) {
                                      if(inherits(args, 'uninitializedField')) args <<- list()
                                      argsToUse <- if(inherits(args, 'symbolTable')) args$symbols else args
                                      if(declaration) {
                                          outputCode <- paste0(if(virtual) 'virtual ' else character(0), generateFunctionHeader(returnType, name, argsToUse, scopes, ...), if(const) ' const ' else character(0), if(abstract) '= 0' else character(0), ';')
                                          if(!inherits(externC, 'uninitializedField' ) ){
                                            if(externC == TRUE)
                                              outputCode <- paste0('extern "C" ', outputCode)
                                          }
                                           return(outputCode) 
                                          
                                          
                                      } else {
                                          if(inherits(code$code, 'uninitializedField')) {
                                              ## There is no code. This can occur for a nimbleFunctionVirtual, which is an abstract base class.
                                              return(character(0))
                                          }
                                          c(paste(generateFunctionHeader(returnType, name, argsToUse, scopes, ...), if(const) ' const ' else character(0), '{'),
                                            code$generate(...),
                                            list('}'))
                                      }
                                  }
                                  )
                              )

generateFunctionHeader <- function(returnType, name, args, scopes = character()) {
    list(paste(returnType$generate(printName = character()), paste(c(scopes, name), collapse = '::'), '(',
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
