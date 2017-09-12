## Remove excess white spaces.  Keep one space.
cleanWhite <- function(s) gsub('[[:blank:]]+', ' ', s)

## This is a base class for c++ variables.
## To avoid having many unused fields in many cases, this can handle only ptrs and refs.
## If you need static, const, template arguments, or other adornments, use cppFullVar
## Below there are some wrappers for common cases like cppDouble, cppNimArrPtr, cppVoid, etc.
cppVar <- setRefClass('cppVar',
                             fields = list(
                                 baseType = 'ANY',	#'character',
                                 ptr = 'ANY',		#'numeric',
                                 ref = 'ANY',		#'logical',
                                 name = 'ANY'), 	#'character'),
                             methods = list(
                             	initialize = function(...){
               		           	  baseType <<- character()
               		           	  ptr <<- numeric()
                	          	  ref <<- logical()
                	          	  name <<- character()
					    		callSuper(...)
                             	},
                                 generate = function(printName = Rname2CppName(.self$name), ...) {
                                     ptrs <- if(length(ptr) > 0) paste(rep('*', ptr), collapse = '')
                                     if(length(printName) > 0) printName <- paste0(printName, collapse = ', ')
                                     cleanWhite(paste(baseType, ptrs, if(identical(ref, TRUE)) '&' else NULL, printName))
                                 },
                                 generateUse = function(...) {
                                     Rname2CppName(name)
                                 },
                                 generateUseDeref = function(...) { 
                                     paste0('(', paste(rep('*', max(0, ptr)), collapse = ''), Rname2CppName(name), ')') ## used to be ptr-asArg
                                 }
                                 )
                            )

## Here is the full version that can handle most c++ variable declarations.
## One thing that cannot be handled is function pointers or member pointers
cppVarFull <- setRefClass('cppVarFull',
                      contains = 'cppVar',
                      fields = list(
                          templateArgs = 'ANY', #'list',
                          baseScope = 'ANY',    #'list',
                          baseConst = 'ANY',    #'logical',
                          baseConstPtr = 'ANY', #'numeric',
                          const = 'ANY', 	#'logical',
                          static = 'ANY', 	#'logical',
                          arraySizes = 'ANY', 	#'integer',
                          constructor = 'ANY', 	#'character',
                          selfDereference = 'ANY',	#'logical'
                          silent = 'ANY'
                          ),
                      methods = list(
                          initialize = function(...) {
                              templateArgs <<- list()
                              baseScope <<- list()
                              baseConstPtr <<- numeric()
                              baseConst <<- logical()
                              const <<- logical()
                              static <<- logical()
                              arraySizes <<- integer()   
                              constructor <<- character()
                              silent <<- FALSE
                              selfDereference <<- FALSE
                              callSuper(...)
                          },
                          generateUse = function(deref, ...) {
                              if(missing(deref)) {
                                  if(selfDereference) generateUseDeref(...)
                                  else callSuper(...)
                              } else {
                                  if(deref) generateUseDeref(...)
                                  else callSuper(...)
                              }
                          },
                          generate = function(printName = Rname2CppName(.self$name), ...) {
                              if(silent) return(character())
                              bCP <- if(length(baseConst) > 0) { 
                                  if(length(baseConstPtr) > 0) paste(paste(rep('*', baseConstPtr), collapse = ''), 'const')
                                  else 'const'
                              }
                              baseTypePlusTemplate <- if(length(templateArgs)==0) baseType
                              else {
                                  expandedTemplateArgs <- unlist(lapply(templateArgs,
                                                                        function(x) {
                                                                            if(inherits(x, 'cppVar')) return(x$generate())
                                                                            return(as.character(x))
                                                                        }))
                                  paste0(baseType,'<', paste(expandedTemplateArgs, collapse = ', '), '>')
                              }
                              ptrs <- if(length(ptr) > 0) paste(rep('*', ptr), collapse = '')
                              if(length(printName) > 0) printName <- paste0(printName, collapse = ', ')
                              ans <- cleanWhite(paste(baseTypePlusTemplate, bCP, ptrs,  if(length(const) > 0) 'const', if(identical(ref, TRUE)) '&' else NULL, printName))
                              if(length(arraySizes) > 0) ans <- paste0(ans, '[', paste0(arraySizes, collapse =']['), ']')
                              ans <- paste0(ans, constructor)
                              if(length(static) > 0) if(static[1]) ans <- paste('static', ans)
                              ans
                          }
                          )
                      )

## Here are some wrappers for simple types

cppStrideType <- function(name = character(0), type = "Stride", strides,...) {
    if(length(strides) != 2) stop('Error in cppStrideType: expecting two strides')
    tA <- lapply(strides, function(x) if(is.na(x)) 'Dynamic' else x)
    cppVarFull(name = name, baseType = 'Stride', templateArgs = tA, ptr = 0, static = FALSE,...)
}

cppEigenMap <- function(name = character(0), type = 'double', eigMatrix = TRUE, strides = numeric(), constructor = '(0,0,0)', ...) {
    templateArgs <- list( paste0(if(eigMatrix) 'MatrixX' else 'ArrayXX',
                                 if(type == 'double') 'd' else if(type == 'integer') 'i' else 'b' ))
    if(length(strides) > 0) templateArgs[[2]] <- cppStrideType(strides = strides)
    cppVarFull(name = name,
               baseType = 'Map',
               templateArgs = templateArgs,
               constructor = constructor,
               ptr = 0,
               static = FALSE,
               ...)
}

emptyTypeInfo <- function() cppVar(baseType = character()) ## for return type of constructors and destructors

cppDouble <- function(name = character(0), ...) cppVar(name = name, baseType = 'double', ...)
cppInt <-  function(name = character(0), ...) cppVar(name = name, baseType = 'int', ...)
cppVoid <- function(name = character(0), ...) cppVar(name = name, baseType = 'void', ...)
cppNimArr <- function(name = character(0), nDim = 1, type = 'double', ptr = 0, ...) cppVarFull(name = name,
                                                                          baseType = 'NimArr',
                                                                          templateArgs = list(nDim, type),
                                                                          ptr = ptr, static = FALSE, ...)
cppNimArrPtr <- function(name = character(0), nDim = 1, type = 'double', ptr = 1, ...){
    cppVarFull(name = name, selfDereference = TRUE,
                                baseType = 'NimArr',
                                templateArgs = list(nDim, type),
                                ptr = ptr, static = FALSE, ...)
                                                                                }
cppVecNimArr <- function(name = character(0), nDim = 1, type = 'double',...) {
    cppVarFull(name = name,
               baseType = 'VecNimArr',
               templateArgs = list(nDim, type),
               ptr = 0,
               static = FALSE, ...)
}

cppVecNimArrPtr <- function(name = character(0), nDim = 1, type = 'double', ptr = 1,...) {
    cppVarFull(name = name,
               baseType = 'VecNimArr',
               templateArgs = list(nDim, type),
               ptr = ptr,
               static = FALSE, ...)
}

cppSEXP <- function(name = character(0), ...) cppVar(name = name, baseType = 'SEXP', ptr = 0, ...)

cppNodeFunctionVector <- function(name = character(0), ...) cppVar(name = name, baseType = 'NodeVectorClassNew', ptr = 0, ...) 

## to be defunct
cppModelVariableAccessorVector <- function(name = character(0), ...) cppVar(name = name, baseType = 'ManyVariablesAccessor', ptr = 0, ...) 
## to be defunct
cppModelValuesAccessorVector <- function(name = character(0), ...) cppVar(name = name, baseType = 'ManyModelValuesAccessor', ptr = 0, ...) 

cppModelVariableMapAccessorVector <- function(name = character(0), ...) cppVar(name = name, baseType = 'ManyVariablesMapAccessor', ptr = 0, ...) 

cppModelValuesMapAccessorVector <- function(name = character(0), ...) cppVar(name = name, baseType = 'ManyModelValuesMapAccessor', ptr = 0, ...) 

cppVecVoidPtr <- function(name = character(0), ...) cppVar(name = name, baseType = 'vector<void*>', ...)
