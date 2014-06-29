
###		buildModelInterface(refName, symbolTable, basePtrCall)
###					Builds a newreference class around a function which returns a
###					pointer to a Model.
###					refname will be the name of the new class
###
###					symbolTable will be a used to determine element types
###					basePtrCall is the C++ call which builds the C++ object
###					i.e. in R .Call(basePtrCall) should return the pointer
###					to our NimbleFunction Class
###					IMPORTANT NOTE: symbolTable MUST match form object pointed to by
###					basePtrCall! In particular, if the symbolTable says an element is a
###					scalar (i.e. nDims = 0) and it's really a NimArr (or vice versa), 
###					attempting to access this element will cause R to crash!

###					Once the new class is built, a new object is create via refname$new()
###					Access to the elements is provided through nimbleModel$Varname
###					Access to the modelValue elements is provided through nimbleModel$modelValues$Varname
###					In this interface, it is assumed that the pointers associated with nimbleModel$Varname
###					are double pointers, i.e. pointers to pointers to NimArr's
###					IMPORTANT NOTE: the function which builds the C++ object (i.e. basePtrCall)
###					MUST initialize the pointers correctly, or accessing will cause R to crash!
###					Default initialization is assumed that elements are pointing to row 1 of
###					the modelValues, but this is not necessary (as long as it is a double pointer to something

getMVptr <- function(rPtr)
  .Call("getModelValuesPtrFromModel", rPtr)
getMVName <- function(modelValuePtr)
  .Call("getMVBuildName", modelValuePtr)

CmodelBaseClass <- setRefClass('CmodelBaseClass',
                               contains = 'modelBaseClass',
                               fields = list(
                                   Rmodel = 'ANY',
                                   cppNames = 'ANY',
                                   cppCopyTypes = 'ANY', ## At the given moment these will all be 'numeric', but the system allows more flexibility
                                   ##CnodeFunClasses = 'list',
                                   compiledModel = 'ANY'
                                   ),
                               methods = list(
                                   show = function() {
                                       cat('CmodelBaseClass object\n')
                                   },
                                   setModel = function(model) {
                                       ## This is creating a circular reference, so be careful with show(), and with Rstudio
                                       Rmodel <<- model 
                                       model$CobjectInterface <- .self
                                   },
                                   copyFromModel = function(model) { ## Could potentially be generated customized for each model
                                       if(missing(model)) model <- Rmodel
                                       for(v in cppNames) {
                                           if(cppCopyTypes[[v]] == 'numeric') {
                                               .self[[v]] <<- model[[v]]
                                               next
                                           }
                                       }
                                   },
                                   setupNodes = function(where = classEnvironment, dll = NULL) {
                                       ## assumes we have the model set
                                       ## there is a field nodes in the modelBaseClass
                                       ##
                                       ## 1. generate CnodeFunClasses
                                       ##     - by iterating through the nodeGenerators in the Rmodel
                                       ##browser()
                                       for(i in names(Rmodel$nodeFunctions)) {
                                           nodes[[i]] <<- nimbleProject$instantiateNimbleFunction(Rmodel$nodeFunctions[[i]], dll = dll)
                                       }
                                       ## for(i in seq_along(Rmodel$nodeGenerators)) {
                                       ##     nodeGenName <- names(Rmodel$nodeGenerators)[i]
                                       ##     nfName <- environment(compiledModel$nodeFuns[[nodeGenName]]$Rgenerator)$refName
                                       ##     CnodeFunClasses[[nfName]] <<- compiledModel$nodeFuns[[nodeGenName]]$Rgenerator
                                           
                                       ##     ## using eval(substitute is purely to avoid a warning about <- instead of <<-
                                       ##     ## nfName is the name of the class created in buildNimbleFxnInterface, i.e. the interface class
                                       ##     ## I don't think this is used anywhere else.
                                       ##     eval(substitute(environment(RM_$nodeGenerators[[nodeGenName]])$CnimbleFunClassName <- nfName, list(RM_ = as.name('Rmodel'))))
                                       ## }
                                       ## ## 2. generate each node object
                                       ## ##     - by iterating through either the nodeFunctions or the instances of the nodeGenerators
                                       ## ##     - the former have names but we must look up the nodeGenerator
                                       ## ##     - the latter don't have names, but the nodeGenerator is there
                                       ## ##         - We will insert a CnimbleFunClassName into the nodeGenerator and iterator over nodeFunctions
                                       ## for(i in names(Rmodel$nodeFunctions)) {
                                       ##     ## We could have a better structure here
                                       ##     nfName <- environment(nf_getRefClassObject(Rmodel$nodeFunctions[[i]])$.generatorFunction)$CnimbleFunClassName
                                       ##     nodes[[i]] <<- CnodeFunClasses[[nfName]](Rmodel$nodeFunctions[[i]], dll = dll)
                                       ##     ## avoid <<- warning:
                                       ##     assign('.CobjectInterface', nodes[[i]], envir = environment(Rmodel$nodeFunctions[[i]]))
                                       ## }
                                   }
                                   )
                               )

makeModelCppCopyTypes <- function(symTab) {
    ans <- list()
    for(s in symTab$symbols) {
        ans[[s$name]] <- 'numeric' ## only one option right now, in a model
    }
    ans
}



makeModelBindingFields <- function(symTab) {
  fieldList = list(.basePtr = "ANY", .modelValues_Ptr = "ANY", .DUMMY = "ANY")  
  vNames = names(symTab$symbols)
  for(vn in vNames){
    ptrName = paste(".", vn, "_Ptr", sep = "")
    fieldList[[ptrName]] <- "ANY"
    eval(substitute( fieldList$VARNAME <- function(x){
      if(missing(x) ) 
        getNimValues(VPTR, 2)
      else
        setNimValues(VPTR, x, 2, allowResize = FALSE)
    }, list(VPTR = as.name(ptrName), VARNAME = vn) ) )
  }
  return(fieldList)
}

buildModelInterface <- function(refName, compiledModel, basePtrCall, project = NULL, dll = NULL, where = globalenv()){
    defaults <- list()
    if(inherits(compiledModel, 'symbolTable')) {
        symTab <- compiledModel
        defaults$cm <- NULL
        warning('No compiled model provided, so interface will be incomplete')
    } else {
        symTab <- compiledModel$model$modelDef$symTab
        defaults$cm <- compiledModel
    }
    defaults$cppCT <- makeModelCppCopyTypes(symTab)
    defaults$project <- project
    
    # resolve basePtrCall rather than looking it up later.
  if(!is.null(dll) && is.character(basePtrCall))
     basePtrCall = getNativeSymbolInfo(basePtrCall, dll)
  else
     warning("creating an initialization method that calls a C routine without any DLL information")

  eval(substitute( newClass <-  setRefClass(refName,
                                            fields = FIELDS,
                                            contains = 'CmodelBaseClass',
                                            methods = list(initialize = function(model, defaults, ..., dll = NULL) { ## model is an optional R model
                                                callSuper()
                                                .basePtr <<- .Call(BPTRCALL)
                                                .modelValues_Ptr <<- getMVptr(.basePtr)
                                                defaultModelValues <<- CmodelValues$new(existingPtr = .modelValues_Ptr, buildCall = getMVName(.modelValues_Ptr) )
                                                modelDef <<- model$modelDef
                                                graph <<- model$graph
                                                vars <<- model$vars
                                                isDataVars <<- model$isDataVars
                                                nimbleProject <<- defaults$project
                                                for(v in ls(model$isDataEnv)) isDataEnv[[v]] <<- model$isDataEnv[[v]]
                                                setData(modelDef$constantsList, warnAboutMissingNames = FALSE)
                                                cppNames <<- .Call("getAvailableNames", .basePtr) ## or could get this from R objects
                                                cppCopyTypes <<- defaults$cppCT
                                                compiledModel <<- defaults$cm
                                                for(vn in cppNames)
                                                    {
                                                        vPtrName <- paste(".", vn, "_Ptr", sep = "")
                                                        eval(substitute(.DUMMY <<- newObjElementPtr(.basePtr, vn), list(.DUMMY = as.name(vPtrName)) ) ) 
                                                    }      
                                                if(!missing(model)) {
                                                    setModel(model)
                                                    copyFromModel()
                                                    setupNodes(dll = dll)
                                                }
                                            },
                                                show = function() {
                                                    writeLines(paste0("Derived CmodelBaseClass created by buildModelInterface for model ", modelDef$name))
                                                }),
                                            where = where
  ), list(BPTRCALL = basePtrCall, FIELDS = makeModelBindingFields(symTab), where = where ) ) )

    ans <- function(model, where = globalenv(), dll = NULL, ...) {
        newClass$new(model, defaults, classEnvironment = where, dll = dll, ...) ## this means defaults will be in the closure of this function and hence found
    } 
    formals(ans)$where = where
    formals(ans)$dll = dll
    ans
}
