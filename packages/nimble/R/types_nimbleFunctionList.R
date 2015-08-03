
nimPointerList <- setRefClass('nimPointerList', ## A base class for a list of objects that, in C++ will be a vector of pointers
                              fields = list(
                                  contentsList = 'ANY', 		#'list',
                                  baseClass = 'ANY', ## R object with information about the pointer-to type
                                  CobjectInterface = 'ANY'
                                  ),
                              methods = list(
                                  initialize = function(virtualBaseClass = NULL, length = 0) {
                                  	contentsList <<- list()
                                      if(!is.null(virtualBaseClass)) {
                                          baseClass <<- virtualBaseClass
                                          if(length > 0) contentsList <<- vector('list', length)
                                      }
                                  }
                                  )
                               )
                              

nimbleFunctionList <- setRefClass('nimbleFunctionList',
                                  contains = 'nimPointerList',
                                  methods = list(
                                      isBaseClassValid = function(x) {
                                          if(!is.nf(x)) return(FALSE)
                                          if(is.null(nfGetDefVar(x, 'contains'))) return(FALSE)
                                          if(!identical(nfGetDefVar(x, 'contains'), baseClass)) return(FALSE)
                                          TRUE
                                      },
                                      checkAllContents = function(x) {
                                          all( unlist( lapply( contentsList, isBaseClassValid ) ) ) 
                                      }
                                  ))

setMethod('[[', 'nimPointerList',
          function(x, i) {
              x$contentsList[[i]]
          })

setMethod('[[<-', 'nimPointerList',
          function(x, i, value) {
              x$contentsList[[i]] <- value
              x
          })

setMethod('[[<-', 'nimbleFunctionList',
          function(x, i, value) {
              if(!x$isBaseClassValid(value)) {stop('An element being put in a nimbleFunctionList is not valid. It may not have the right contains (base class)', call. = FALSE)}
              x$contentsList[[i]] <- value
              x
          })

checkNimbleFunctionListCpp <- function(nfl) {
    otherProblem <- try(
        {
            baseClassName <- environment(nfl$baseClass)$name
            containedClassNames <- lapply( nfl$contentsList, function(x) if(is.list(x)) x[[1]]$compiledNodeFun$inheritance else x$compiledNodeFun$inheritance)
            ok <- lapply(containedClassNames, function(x) any(unlist(x) == baseClassName))
            if(all(unlist(ok))) return(TRUE) else return(FALSE)
        })
    return(FALSE)
}

length.nimPointerList <- function(pointerList) {
	return(length(pointerList$contentsList))
}
