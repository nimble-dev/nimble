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
        extptrCall = 'ANY',
        varNames = 'ANY',
        componentExtptrs = 'ANY',
        blankAns = 'ANY',
        ##.nodePtrs_byGID = 'ANY',
        GID_map = 'ANY',
        symTab = 'ANY',
        initialized = 'ANY',
    	sizes = function(){
    		if(length(varNames) == 0)
    			return(NULL)
    		sizeList <- list()
    		for(compName in varNames ){
    			sizeList[[compName]] <- .Call(getNativeSymbolInfo('getMVsize', dll), componentExtptrs[[compName]])
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
            return(GID_map$expandNodeNames(nodeNames = nodeNames, returnType = returnType, 
                                           flatIndices = flatIndices))
        },
        initialize = function(buildCall, existingPtr, initialized = FALSE, dll) {
            if(missing(existingPtr) ) {
                if(is.character(buildCall)) {
                    warning("a call to getNativeSymbolInfo with only a name and no DLL")
                }
                # avoid R CMD check problem with registration
                ## notice that buildCall is the result of getNativeSymbolInfo using the dll from nimbleProject$instantiateCmodelValues
                ## only other calling point is from cppInterfaces_models.R, and in that case existingPtr is provided
                extptr <<- eval(parse(text = ".Call(buildCall)"))
#                extptr <<- .Call(buildCall) 
            }
            else{
                extptr <<- existingPtr
            }
            dll <<- dll
            initialized <<- initialized
            if(missing(buildCall) ) 
                break("Cannot build object without buildCall!")
            extptrCall <<- buildCall
            varNames <<- .Call(getNativeSymbolInfo('getAvailableNames', dll), extptr)      
            componentExtptrs <<- vector(mode = 'list', length = length(varNames))
            names(componentExtptrs) <<- varNames
            blankAns <<- componentExtptrs
            for(comp in varNames) 
                componentExtptrs[[comp]] <<- .Call(getNativeSymbolInfo('getModelObjectPtr', dll), extptr, comp)
                
            ## .nodePtrs_byGID <<- new('numberedObjects')
            ## GID_map <<- makeMV_GID_Map(.self)
            ## if(length(sizes) > 0){
            ## 	varLengths <- sapply(sizes, prod)
            ## 	totLength <- sum(varLengths)
            ##     	.nodePtrs_byGID$resize(totLength)
            ## 	index = 1
            ## 	for(i in seq_along(varNames)){
            ##     	  	vName <- varNames[i]
            ##     	  	.Call('populateNumberedObject_withSingleModelValuesAccessors', extptr, vName, as.integer(expandNodeNames(vName, returnType = 'ids')), as.integer(1) , .nodePtrs_byGID$.ptr)
            ## 		index = index + varLengths[i]
            ## 	}
            ## }
        },
        resize = function(rows){	
        	for(ptr in componentExtptrs)
        	jnk <- .Call(getNativeSymbolInfo("setNumListRows", dll), ptr, as.integer(rows), TRUE)
        	jnk <- .Call(getNativeSymbolInfo('manualSetNRows', dll), extptr, as.integer(rows) )  	
        	},
        getSize = function() {	.Call(getNativeSymbolInfo("getNRow", dll), componentExtptrs[[1]]) } ## formerly getCRows(componentExtptrs[[1]])		}
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
              		output = .Call(getNativeSymbolInfo('getMVElement', x$dll), ptr, as.integer(j) )
              		return(output) 
              		}
              	output = .Call(getNativeSymbolInfo('getMVElementAsList', x$dll), ptr, as.integer(j) )
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
						.Call(getNativeSymbolInfo('setMVElement', x$dll), ptr, as.integer(j), value )
						return(x)
					}
				for(jj in j)
					storage.mode(value[[jj]]) <- 'numeric'
				.Call(getNativeSymbolInfo('setMVElementFromList', x$dll), ptr, value, as.integer(j) )
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
                      output = .Call(getNativeSymbolInfo('getMVElement', x$dll), ptr, as.integer(1) )
                      return(output) 
                  }
                  output = .Call(getNativeSymbolInfo('getMVElementAsList', x$dll), ptr, as.integer(1:k) )
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
						.Call(getNativeSymbolInfo('setMVElement', x$dll), ptr, as.integer(1), as.numeric(value) )
						return(x)
					}
				.Call(getNativeSymbolInfo('setMVElementFromList', x$dll), ptr, as.numeric(value), as.integer(1:k) )
				return(x)
				}
			cmpNames = names(value)
			if( !all(cmpNames %in% x$varNames) ) 
				stop('Warning: names of modelValue elements do not match')
			for(n in cmpNames)
				x[[n]] <- value[[n]]
			return(x)
          })
