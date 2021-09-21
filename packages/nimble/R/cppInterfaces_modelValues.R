###     CmodelValues    Reference class to access the modelValues objects
###
###                     New CModelValues object is created with 
###						"CModelValues$new(type, exisitingptr, cppTypeName)"
###						If type provided, object is built around
###						genSymbolTablefromModel(type)	                        
###						Otherwise, both exisitingptr and cppTypeName is required
###						exisitingptr should actually be a pointer to the C++ object
###						we want (not sure if this is the best way to do it: see how we 
###						handle reference classes later). cppTypeName names the class,
###						but is not required to match anything (unlike "type")

###						To access the elements of a CmodelValues object, we can call
###						CModelObject[["VARNAME"]]. This will return a list, with each
###						element of the list corresponding with that row of the VecNimArr
###						Values can also be set with CModelObject[["VARNAME"]] <- list

###						Another way we would like to access the list is through 
###						CModelObject[i, "VARNAME"]. We would like this to return
###						a new CModelObject of the same type, which the i^{th} row is filled out
###						for variable VARNAME, and all other entries are blank. HOWEVER
###						this requires that CModelObject was built using type, not existingptr
###						Need to discuss what to do if built around existingptr



CmodelValues <- setRefClass(
    Class = 'CmodelValues',
    fields = list(
        dll = 'ANY',
        extptr = 'ANY',
        namedObjectsPtr = 'ANY',
        extptrCall = 'ANY',
        varNames = 'ANY',
        componentExtptrs = 'ANY',
        blankAns = 'ANY',
        GID_map = 'ANY',
        symTab = 'ANY',
        initialized = 'ANY',
    	sizes = function(){
    		if(length(varNames) == 0)
    			return(NULL)
    		sizeList <- list()
    		for(compName in varNames ){
    			sizeList[[compName]] <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getMVsize, componentExtptrs[[compName]]))
    		}
    	return(sizeList)
    	}),
    methods = list(
        show = function() {
            writeLines(paste0("CmodelValues object with variables: ", paste(varNames, collapse = ", "), "."))
        },
        getVarNames = function(includeLogProb = FALSE){
            if(includeLogProb)
                return(varNames)
            return(varNames[!grepl('logProb_', varNames)])
        },
        expandNodeNames = function (nodeNames, returnType = "names", flatIndices = TRUE)  {
            stop('There was a call to expandNodeNames for a compiled modelValues object.\n  This is deprecated and should not have occurred.\n  Please contact nimble developers at the nimble-users google group or at nimble.stats@gmail.com to let them know this happened.  \n Thank you.')
        },
        finalizeInternal = function() {
            finalize()
            extptr <<- NULL
            namedObjectsPtr <<- NULL
        },
        finalize = function() {
            nimbleInternalFunctions$nimbleFinalize(namedObjectsPtr) ##
        },
        initialize = function(buildCall, existingPtr, initialized = FALSE, dll) {
            if(missing(existingPtr) ) {
                if(is.character(buildCall)) {
                    warning("a call to getNativeSymbolInfo with only a name and no DLL")
                }
               
                # avoid R CMD check problem with registration
                ## notice that buildCall is the result of getNativeSymbolInfo using the dll from nimbleProject$instantiateCmodelValues
                ## only other calling point is from cppInterfaces_models.R, and in that case existingPtr is provided
                extptrlist <- eval(parse(text = ".Call(buildCall)"))
                extptr <<- extptrlist[[1]]
                namedObjectsPtr <<- extptrlist[[3]] ## order should come from the cppDef, but cheating here to get it right
                eval(call('.Call',nimbleUserNamespace$sessionSpecificDll$register_namedObjects_Finalizer, namedObjectsPtr, dll[['handle']], 'modelValues'))
            }
            else{
                extptr <<- existingPtr
            }
            dll <<- dll
            initialized <<- initialized
            if(missing(buildCall) ) 
                stop("Cannot build object without buildCall!")
            extptrCall <<- buildCall
            varNames <<- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getAvailableNames, extptr))      
            componentExtptrs <<- vector(mode = 'list', length = length(varNames))
            names(componentExtptrs) <<- varNames
            blankAns <<- componentExtptrs
            for(comp in varNames) 
                componentExtptrs[[comp]] <<- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getModelObjectPtr, extptr, comp))
        },
        resize = function(rows){	
        	for(ptr in componentExtptrs)
        	jnk <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$setVecNimArrRows, ptr, as.integer(rows), TRUE))
        	jnk <- eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$manualSetNRows, extptr, as.integer(rows)) )  	
        	},
        getSize = function() {
            if(length(componentExtptrs))
                eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getNRow, componentExtptrs[[1]])) else 0
        } ## formerly getCRows(componentExtptrs[[1]])		}
        )
    )


setMethod('[', 'CmodelValues',
          function(x, i, j) {
              if(missing(i))
                  i = x$varNames 
              if(missing(j) ) 
                  j = 1:cGetNRow(x)
              if(length(i) == 1){
                  ptr = x$componentExtptrs[[i]]
                  if(is.null(ptr) ) 
                      stop(paste('variable', i, ' not found in modelValues') ) 
                  if(length(j) == 1){
                      output = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getMVElement, ptr, as.integer(j)) )
                      return(output) 
                  }
                  output = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getMVElementAsList, ptr, as.integer(j)) )
                  return(output)
              }
              output <- list() 	
              for(cmp in i)
                  output[[cmp]] <- x[cmp, j]
              return(output)
          }
          )

setMethod('[<-', 'CmodelValues',
          function(x, i, j, value){
              if(missing(i) ) 
                  i = x$varNames
              if(missing(j) ) 
                  j = 1:getsize(x)
              if(length(i) == 1){
                  ptr = x$componentExtptrs[[i]]
                  if(is.null(ptr) ) 
                      stop(paste('variable', i, ' not found in modelValues') ) 
                  if(length(j) == 1){
                      storage.mode(value) <- 'numeric'
                      eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$setMVElement, ptr, as.integer(j), value ))
                      return(x)
                  }
                  for(jj in j)
                      storage.mode(value[[jj]]) <- 'numeric'
                  eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$setMVElementFromList, ptr, value, as.integer(j)) )
                  return(x)
              }
              cmpNames = names(value)
              if( !all(cmpNames %in% x$varNames) ) 
                  stop('Warning: names of modelValue elements do not match')
              for(n in cmpNames)
                  x[n, j] <- value[[n]]
              return(x)
          })


setMethod('[[', 'CmodelValues',
          function(x, i, ...) {
              if(missing(i) )
                  i <- x$varNames
              if(length(i) == 1){
                  k <- getsize(x)
                  if(k == 0) return(list())
                  ptr = x$componentExtptrs[[i]]
                  if(is.null(ptr) ) 
                      stop(paste('variable', i, ' not found in modelValues') )
                  if(k == 1){
                      output = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getMVElement, ptr, as.integer(1)) )
                      return(output) 
                  }
                  output = eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$getMVElementAsList, ptr, as.integer(1:k)) )
                  return(output)
              }
              output <- list() 	
              for(cmp in i)
                  output[[cmp]] <- x[[cmp]]
              return(output)
          })

setMethod('[[<-', 'CmodelValues',
          function(x, i,..., value){
  				if(missing(i) ) 
					i = x$varNames
				k <- getsize(x)
				if(length(i) == 1){
	              	ptr = x$componentExtptrs[[i]]
	              	if(is.null(ptr) ) 
	              		stop(paste('variable', i, ' not found in modelValues') ) 
					if(k == 1){
						eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$setMVElement, ptr, as.integer(1), as.numeric(value) ))
						return(x)
					}
				eval(call('.Call', nimbleUserNamespace$sessionSpecificDll$setMVElementFromList, ptr, value, as.integer(1:k) ))
				return(x)
				}
			cmpNames = names(value)
			if( !all(cmpNames %in% x$varNames) ) 
				stop('Warning: names of modelValue elements do not match')
			for(n in cmpNames)
				x[[n]] <- value[[n]]
			return(x)
          })
