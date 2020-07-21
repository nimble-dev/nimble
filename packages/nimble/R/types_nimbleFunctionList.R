
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
                              
#' Create a list of nimbleFunctions
#' 
#' Create an empty list of nimbleFunctions that all will inherit from a base class.
#' 
#' @author NIMBLE development team
#' @export
#' @details
#' See the User Manual for information about creating and populating a \code{nimbleFunctionList}.
nimbleFunctionList <- setRefClass('nimbleFunctionList',
                                  contains = 'nimPointerList',
                                  methods = list(
                                      isBaseClassValid = function(x) {
                                          if(!is.nf(x)) return(FALSE)
                                          this_contains <- nfGetDefVar(x, 'contains')
                                          done <- FALSE
                                          ok <- FALSE
                                          baseClassName <- environment(baseClass)[['name']]
                                          while(!done) {
                                              if(is.null(this_contains)) {
                                                  done <- TRUE
                                              }
                                              else {
                                                  this_name <- environment(this_contains)[['name']]
                                                  if(identical(this_name, baseClassName)) {
                                                      ok <- TRUE
                                                      done <- TRUE
                                                  }
                                              }
                                              ## recurse up inheritance path
                                              if(!done)
                                                  this_contains <- nfGetDefVar(this_contains, 'contains')
                                          }
                                          ok
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
    return(TRUE)
    ## We have disabled this test because it does not look up an multiple inheritance tree.
    ## At this point, that would take some effort.
    ## However, inheritance validity is checked by the `[[<-` method for a nimbleFunctionList,
    ## so there "shouldn't" be a downstream problem here if everything checked ok when the nimbleFunctionList was populated. 
    otherProblem <- try(
        {
            baseClassName <- environment(nfl$baseClass)$name
            containedClassNames <- lapply( nfl$contentsList, function(x) if(is.list(x)) x[[1]]$compiledNodeFun$inheritance else x$compiledNodeFun$inheritance)
            ok <- lapply(containedClassNames, function(x) any(unlist(x) == baseClassName))
            if(all(unlist(ok))) return(TRUE) else return(FALSE)
        })
    return(FALSE)
}

length.nimPointerList <- function(x, ...) {
	return(length(x$contentsList))
}
