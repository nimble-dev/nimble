
#' EXPERIMENTAL: Turn a function into a model macro builder
#' A model macro expands one line of code in a nimbleModel into one or
#' more new lines.  This supports compact programming by defining
#' re-usable modules.  \code{model_macro_builder} takes as input a
#' function that constructs new lines of model code from the original
#' line of code.  It returns a function suitable for internal use by
#' \code{nimbleModel} that arranges arguments for input function.  Macros
#' are an experimental feature and are available only after setting
#' \code{nimbleOptions(enableModelMacros = TRUE)}.
#'
#' @param fun A function written to construct new lines of model code.
#'
#' @param use3pieces (TRUE or FALSE) Should the arguments from the input
#' line be split into pieces for the LHS (left-hand side), RHS
#' (right-hand side, possibly further split depending on
#' \code{unpackArgs}), and \code{stoch} (TRUE if the line uses a
#' \code{~}, FALSE otherwise)?  See details and examples.
#'
#' @param unpackArgs (TRUE or FALSE) Should arguments be passed as a list
#' (FALSE) or as separate arguments (TRUE)?  See details and examples.
#'
#' @details The arguments \code{use3pieces} and \code{unpackArgs}
#' indicate how \code{fun} expects to have arguments arranged from an
#' input line of code (processed by \code{nimbleModel}).
#'
#' Consider the defaults \code{use3pieces = TRUE} and \code{unpackArgs =
#' TRUE}, for a macro called \code{macro1}.  In this case, the line of
#' model code \code{x ~ macro1(arg1 = z[1:10], arg2 = "hello")} will be
#' passed to \code{fun} as \code{fun(stoch = TRUE, LHS = x, arg1 =
#' z[1:10], arg2 = "hello")}.
#'
#' If \code{use3pieces = TRUE} but \code{unpackArgs = FALSE}, then the
#' RHS will be passed as is, without unpacking its arguments into
#' separate arguments to \code{fun}.  In this case, \code{x ~ macro1(arg1
#' = z[1:10], arg2 = "hello")} will be passed to \code{fun} as
#' \code{fun(stoch = TRUE, LHS = x, RHS = macro1(arg1 = z[1:10], arg2 =
#' "hello"))}.
#'
#' If \code{use3pieces = FALSE} and \code{unpackArgs = FALSE}, the entire
#' line of code is passed as a single object.  In this case, \code{x ~
#' macro1(arg1 = z[1:10], arg2 = "hello")} will be passed to \code{fun}
#' as \code{fun(x ~ macro1(arg1 = z[1:10], arg2 = "hello"))}.  It is also
#' possible in this case to pass a macro without using a \code{~} or
#' \code{<-}.  For example, the line \code{macro1(arg1 = z[1:10], arg2 =
#' "hello")} will be passed to \code{fun} as \code{fun(macro1(arg1 =
#' z[1:10], arg2 = "hello"))}.
#'
#' If \code{use3pieces = FALSE} and \code{unpackArgs = TRUE}, it
#' won't make sense to anticipate a declaration using \code{~} or \code{<-}.  Ins#' tead, arguments from an arbitrary call will be passed as separate arguments.  #' For example, the line \code{macro1(arg1 = z[1:10], arg2 = "hello")} will be pa#' ssed to \code{fun} as \code{fun(arg1 = z[1:10], arg2 = "hello")}.
#'
#' It is extremely useful to be familiar with processing R code as an
#' object to write \code{fun} correctly.  Functions such as
#' \code{\link{substitute}} and \code{\link{as.name}}
#' (e.g. \code{as.name('~')}), \code{\link{quote}}, \code{\link{parse}}
#' and \code{\link{deparse}} are particularly handy.
#'
#' Multiple lines of new code should be contained in \code{ {} }. Extra
#' curly braces are not a problem. See example 2.
#'
#' Macro expansion is done recursively: One macro can return code that
#' invokes another macro.
#' 
#' @return A list with a named element \code{code} that contains the
#' replacement code.
#'
#' @export
#' 
#' @examples
#' nimbleOptions(enableModelMacros = TRUE)
#' nimbleOptions(enableMacroComments = FALSE)
#' nimbleOptions(verbose = FALSE)
#' 
#' ## Example 1: Say one is tired of writing "for" loops.
#' ## This macro will generate a "for" loop with dnorm declarations
#' all_dnorm <- model_macro_builder(
#'     function(stoch, LHS, RHSvar, start, end, sd = 1, modelInfo, .env) {
#'         newCode <- substitute(
#'             for(i in START:END) {
#'                 LHS[i] ~ dnorm(RHSvar[i], SD)
#'             },
#'             list(START = start,
#'                  END = end,
#'                  LHS = LHS,
#'                  RHSvar = RHSvar,
#'                  SD = sd))
#'         list(code = newCode)
#'     },
#'     use3pieces = TRUE,
#'     unpackArgs = TRUE 
#' )
#' 
#' model1 <- nimbleModel(
#'     nimbleCode(
#'     {
#'         ## Create a "for" loop of dnorm declarations by invoking the macro
#'         x ~ all_dnorm(mu, start = 1, end = 10)
#'     }
#'     ))
#' 
#' ## show code from expansion of macro
#' model1$getCode()
#' ## The result should be:
#' ## {
#' ##     for (i in 1:10) {
#' ##         x[i] ~ dnorm(mu[i], 1)
#' ##     }
#' ## }
#' 
#' ## Example 2: Say one is tired of writing priors.
#' ## This macro will generate a set of priors in one statement
#' flat_normal_priors <- model_macro_builder(
#'     function(..., modelInfo, .env) {
#'         allVars <- list(...)
#'         priorDeclarations <- lapply(allVars,
#'                                     function(x)
#'                                         substitute(VAR ~ dnorm(0, sd = 1000),
#'                                                    list(VAR = x)))
#'         newCode <- quote({})
#'         newCode[2:(length(allVars)+1)] <- priorDeclarations
#'         list(code = newCode)
#'     },
#'     use3pieces = FALSE,
#'     unpackArgs = TRUE
#' )
#' 
#' model2 <- nimbleModel(
#'     nimbleCode(
#'     {
#'         flat_normal_priors(mu, beta, gamma)
#'     }
#'     ))
#' 
#' ## show code from expansion of macro
#' model2$getCode()
#' ## The result should be:
#' ## {
#' ##    mu ~ dnorm(0, sd = 1000)
#' ##    beta ~ dnorm(0, sd = 1000)
#' ##    gamma ~ dnorm(0, sd = 1000)
#' ## }
model_macro_builder <- function(fun,
                                use3pieces = TRUE,
                                unpackArgs = TRUE ) {
    if(use3pieces) {
        wrapper <- function(code, modelInfo, .env) {
            args <- as.list(code)
            args[[1]] <- args[[1]] == '~'
            names(args) <- c('stoch','LHS','RHS')
            if(unpackArgs) {
                RHSargs <- as.list(args[[3]])[-1]
                args <- c(args[1:2], RHSargs)
            }
            ## Since code contains unevaluated expressions, they need
            ## to be wrapped in another layer of quote() that do.call
            ## will strip.  Otherwise there will be an attempted
            ## evaluation of variable names.
            args <-  lapply(args,
                            function(x)
                                substitute(quote(X), list(X = x)))
            do.call(fun, c(args, list(modelInfo = modelInfo, .env = .env)))
        }
    } else {
        if(unpackArgs) {
            wrapper <- function(code, modelInfo, .env) {
                args <- as.list(code)[-1]
                args <-  lapply(args,
                            function(x)
                                substitute(quote(X), list(X = x)))
                do.call(fun, c(args, list(modelInfo = modelInfo, .env = .env)))
            }
        } else {
            wrapper <- fun
        }
    }
    ans <- structure(list(process = wrapper), class = "model_macro")
    ans
}

## This function recurses through a block of code and expands any submodels
codeProcessModelMacros <- function(code,
                                   modelInfo,
                                   env = parent.frame(),
                                   recursionLabels = character()) {
    expandRecursionLabels <- function(possibleMacroName,
                                      labels = character()) {
        paste0(possibleMacroName,
               if(length(labels) > 0)
                   paste0('(expanded from ',
                          paste(labels, collapse = '-->'),
                          ')')
               else
                   character()
               )
    }
    codeLength <- length(code)
    ## First check if this is the start of a curly-bracketed block
    if(code[[1]] == '{') {
        if(codeLength > 1)
            ## Recurse on each line
            for(i in 2:codeLength){
              macroOutput <- codeProcessModelMacros(code = code[[i]],
                                                    modelInfo = modelInfo,
                                                    env = env,
                                                    recursionLabels
                                                    )
              code[[i]] <- macroOutput$code
              modelInfo <- macroOutput$modelInfo
            }
        return(list(code=code, modelInfo=modelInfo))
    }
    ## If this is a for loop, recurse on the body of the loop
    if(code[[1]] == 'for') {
      macroOutput <- codeProcessModelMacros(code[[4]],
                                            modelInfo = modelInfo,
                                            env = env,
                                            recursionLabels
                                           )
      code[[4]] <- macroOutput$code
      return(list(code=code, modelInfo=macroOutput$modelInfo))
    }
    ## Check if this line invokes a submodel.
    ## This can be done in two ways:
    ## (i) node1 [<- | ~] <macro name>(...)
    ## or
    ## (ii) <macro name>(args)
    ##
    ## The first version is more BUGS-like.
    ## The second version allows more full control.

    ## Initialize possibleMacroName assuming version (ii):
    possibleMacroName <- safeDeparse(code[[1]], warn = TRUE)
    ## If it is really version (i), possibleMacroName will be
    ## ~ or <- and should be updated to the call on the right-hand side:
    if(possibleMacroName %in% c('<-', '~') && !is.name(code[[3]])) {
        possibleMacroName <- safeDeparse(code[[3]][[1]], warn = TRUE)
    }
    if(exists(possibleMacroName)) { ## may need to provide an envir argument
        possibleMacro <- get(possibleMacroName) ## ditto
        if(inherits(possibleMacro, "model_macro")) {
            expandedInfo <- try(possibleMacro$process(code, 
                                                      modelInfo=modelInfo,
                                                      .env = env
                                                      ))
            if(inherits(expandedInfo, 'try-error'))
                stop(paste0("Model macro ",
                            expandRecursionLabels(
                                possibleMacroName,
                                recursionLabels
                            ),
                            " failed."),
                     call. = FALSE)
            if(!is.list(expandedInfo))
                stop(paste0("Model macro ",
                            expandRecursionLabels(
                                possibleMacroName,
                                recursionLabels
                            ),
                            " should return a list with an element named ",
                            "'code'.  It did not return a list."),
                     call. = FALSE)
            if(!is.call(expandedInfo[['code']]))
                stop(paste0("Model macro ",
                            expandRecursionLabels(
                                possibleMacroName,
                                recursionLabels
                            ),
                            " should return a list with an element named ",
                            "'code' that is a call."),
                     call. = FALSE)
            
            # Check for newly created parameters
            newPars <- list(newMacroPars(code, expandedInfo$code))
            names(newPars) <- possibleMacroName

            # Add automatic comments showing code added by macros
            if(getNimbleOption("enableMacroComments")){
              if(getNimbleOption("codeInMacroComments")){
                macroComment <- paste("#", deparse(code))
              } else {
                macroComment <- paste("#", possibleMacroName)
              }
              macroEnd <- "# ----"
              if(length(recursionLabels > 0)){
                spacer <- paste(rep("  ", length(recursionLabels)), collapse="")
                hashes <- paste(rep("#", length(recursionLabels)+1), collapse="")
                if(getNimbleOption("codeInMacroComments")){  
                  macroComment <- paste0(spacer, hashes, " ", deparse(code))
                } else {
                  macroComment <- paste0(spacer, hashes, " ", possibleMacroName)
                }
                macroEnd <- paste0(spacer, hashes, " ----")
              }
              macroStartLine <- substitute(MACRO, list(MACRO = macroComment))
              macroEndLine <- substitute(END, list(END = macroEnd))
              expandedInfo$code <- as.call(c(list(quote(`{`)), list(macroStartLine, expandedInfo$code, macroEndLine)))
            }

            curPars <- modelInfo$parameters
            expandedInfo$modelInfo$parameters <- c(curPars, newPars)

            ## Return object is a list so we can ossibly extract other
            ## content in the future.  We recurse on the returned code
            ## to expand macros that it might contain.
            macroOutput <- codeProcessModelMacros(expandedInfo$code,
                                           modelInfo = expandedInfo$modelInfo,
                                           env = env,
                                           recursionLabels = c(recursionLabels, possibleMacroName)
                                           )
            return(list(code=macroOutput$code, modelInfo=macroOutput$modelInfo))
        }
    }
    list(code=code, modelInfo=modelInfo)
}

# Index generator for all macros
macroIndexCreator <- labelFunctionCreator("i")

processModelMacros <- function(code, modelInfo, env){
  # No generated parameters before any macros run, so parameters = empty list
  macroIndexCreator(reset = TRUE)
  modelInfo$indexCreator <- macroIndexCreator
  modelInfo$parameters <- list()
  macroOutput <- codeProcessModelMacros(code=code, modelInfo=modelInfo, env=env)
  # Clean up extra brackets
  macroOutput$code <- removeExtraBrackets(macroOutput$code)
  # Convert factors to numeric
  macroOutput$modelInfo$constants <- convertFactorConstantsToNumeric(macroOutput$modelInfo$constants)  
  # Remove intermediate parameters from list and check for duplicates
  macroOutput$modelInfo$parameters <- checkMacroPars(macroOutput$modelInfo$parameters,
                                                     code, macroOutput$code)

  list(code = macroOutput$code, modelInfo = macroOutput$modelInfo)
}

convertFactorConstantsToNumeric <- function(constants){
  lapply(constants, function(x){
    if(is.factor(x)){
      x <- as.numeric(x)
    } else if(is.character(x)){
      if(is.null(dim(x))){
        x <- as.numeric(as.factor(x))
      } else {
        d <- dim(x)
        x <- as.numeric(as.factor(x))
        x <- array(x, dim=d)
      }
    }
    x
  })
}

# Get parameters in a line of code generated by a macro
getMacroParsInternal <- function(code){
  if(is.call(code)){
    code <- as.list(code)[2:length(code)]
    return(lapply(code, getMacroParsInternal))
  } else {
    return(sapply(code, function(x) if(is.numeric(x) || x == "") return(NULL) else return(deparse(x))))
  }
}

getMacroPars <- function(code){
  out <- getMacroParsInternal(code)
  out <- unlist(out)
  out <- out[out != ""]
  out <- out[!is.numeric(out)]
  unique(out)
}

processCodeLine <- function(code){
  #stopifnot(is.call(code))
  if(is.character(code)) return(list(LHS = NULL, RHS = NULL))
  if(is.name(code)){
    return(list(LHS = NULL, RHS = deparse(code)))
  }
  if(code[[1]] == "{"){
    return(lapply(as.list(code)[2:length(code)], processCodeLine))
  }
  if(code[[1]] == "for"){
    return(lapply(as.list(code)[2:length(code)], processCodeLine))
  } else {
    #if(isAssignment())
    if(as.character(code[[1]]) %in% c("~", "<-")){
      RHS <- getMacroPars(code[[3]]) # getRHS()
      LHS <- getMacroPars(code[[2]]) # getLHS()
    } else {
      RHS <- getMacroPars(code)
      LHS <- NULL
    }
  }

  list(LHS = LHS, RHS = RHS)
}

findAllListElementsByNameInternal <- function(inputList, name){
  if(is.list(inputList) && name %in% names(inputList)){
    return(inputList[[name]])
  } else if(is.list(inputList)){
    return(lapply(inputList, findAllListElementsByNameInternal, name = name))
  }
  return(NULL)
}


findAllListElementsByName <- function(inputList, name){
  stopifnot(is.list(inputList))
  out <- findAllListElementsByNameInternal(inputList, name)
  unlist(out)
}

processCode <- function(code){
  out <- processCodeLine(code)
  LHS = as.character(findAllListElementsByName(out, "LHS"))
  RHS = as.character(findAllListElementsByName(out, "RHS"))
  list(LHS = unique(LHS), RHS = unique(RHS))
}

# Figure out which parameters a macro added by comparing with the original
# line of code
newMacroPars <- function(startCode, endCode){
  startPar <- processCode(startCode)
  endPar <- processCode(endCode)
  LHS <- endPar$LHS[! endPar$LHS %in% unlist(startPar)]
  if(length(LHS) == 0) LHS <- NULL
  RHS <- endPar$RHS[! endPar$RHS %in% unlist(startPar)]
  if(length(RHS) == 0) RHS <- NULL
  list(LHS = LHS, RHS = RHS)
}

# Remove intermediate parameters and check for parameters generated
# more than once
checkMacroPars <- function(parameters, startCode, endCode){

  # If no macros return empty list
  if(length(parameters) == 0) return(parameters)
  
  # Find new parameters all at once
  all_pars <- newMacroPars(startCode, endCode)

  if(length(all_pars) == 0) return(NULL)

  # Remove intermediate parameters that don't end up final code
  final_pars <- lapply(parameters, function(x){
    list(LHS = x$LHS[x$LHS %in% all_pars$LHS],
         RHS = x$RHS[x$RHS %in% all_pars$RHS])
  })
  
  # Check for duplicate parameters from previous macros and warn
  if(length(final_pars) > 1){
  for (i in 2:length(final_pars)){
    this_macro <- final_pars[[i]]$LHS
    for (j in 1:(i-1)){
      dups <- this_macro %in% final_pars[[j]]$LHS
      if(any(dups)){

       msg <- paste0("  [Note] Macro '", names(final_pars[i]), "' declared LHS parameter(s): '",
                paste(this_macro[dups], collapse="', '"), "'\n         previously declared by macro '",
                names(final_pars[j]),"'")
        messageIfVerbose(msg)

        #warning(paste0("Macro '", names(final_pars[i]), "' generated parameter name(s) '",
        #        paste(this_macro[dups], collapse="', '"), "' previously generated by macro '",
        #        names(final_pars[j]),"'"), call.=FALSE)
      }
    }
  }
  }

  # Organize multiple instances of macros into a list
  unq_macros <- unique(names(final_pars))
  fp <- lapply(unq_macros, function(x){
    out <- final_pars[names(final_pars) == x]
    names(out) <- NULL
    out
  })
  names(fp) <- unq_macros
  final_pars <- fp

  final_pars
}


#' EXPERIMENTAL: Get list of parameter names generated by model macros
#'
#' Get a list of all parameter names (or certain categories of parameters)
#' generated by model macros in the model code.
#'
#' @param model A NIMBLE model object
#' @param includeLHS Include generated parameters on the left-hand side (LHS) of 
#'  assignments (\code{<-} or \code{~}) in the output
#' @param includeRHS Include generated parameters on the left-hand side (RHS) of 
#'  assignments (\code{<-} or \code{~}) in the output
#' @param includeDeterm Include deterministic generated parameters in the output
#' @param includeStoch Include stochastic generated parameters in the output
#' @param includeIndices Include index parameters generated for use in for loops
#'  in the output
#'
#' @details Some model macros will generate new parameters to be included in the
#'  output code. NIMBLE automatically detects these new parameters and 
#'  records them in the model object. This function allows easy access
#'  to this stored list, or subsets of it (see arguments).
#'
#' @return A named list of generated parameters, with the element names 
#' corresponding to the original source macro.
#'
#' @export
#'
#' @examples
#' nimbleOptions(enableModelMacros = TRUE)
#' nimbleOptions(enableMacroComments = FALSE)
#' nimbleOptions(verbose = FALSE)
#'
#' testMacro <- list(process = function(code, modelInfo, .env){
#'   code <- quote({
#'     for (i_ in 1:n){
#'       mu[i_] <- alpha + beta
#'       y[i_] ~ dnorm(0, sigma)
#'     }
#'     alpha ~ dnorm(0, 1)
#'   })
#'   list(code = code, modelInfo=modelInfo)
#' })
#'  class(testMacro) <- "model_macro"
#'
#' code <- nimbleCode({
#'   y[1:n] ~ testMacro()
#' })
#'
#' const <- list(y = rnorm(10), n = 10)
#'
#' mod <- nimbleModel(code, constants=const)
#'
#' mod$getMacroParameters()
#' # should be list(testMacro = list(c("mu", "alpha", "beta", "sigma")))
#'
#' mod$getMacroParameters(includeRHS = FALSE)
#' # should be list(testMacro = list(c("mu", "alpha")))
getMacroParameters <- function(model, includeLHS = TRUE, includeRHS = TRUE, includeDeterm = TRUE, 
                              includeStoch = TRUE, includeIndices = FALSE) {
  dc <- model$modelDef$declInfo
  out <- model$modelDef$macroParameters

  #RHS/LHS
  if(!includeRHS){
    out <- lapply(out, function(x){
      lapply(x, function(z){
        z$RHS <- NULL
        z
      })
    })
  }

  if(!includeLHS){
    out <- lapply(out, function(x){
      lapply(x, function(z){
        z$LHS <- NULL
        z
      })
    })
  }

  # Indices
  if(!includeIndices){
    idx <- sapply(unlist(lapply(dc, function(x) x$indexExpr)), deparse)
    out <- lapply(out, function(x){
      lapply(x, function(y){
        yout <- lapply(y, function(z){
          zout <- z[! z %in% idx]
          if(length(zout) == 0) zout <- NULL
          zout
        })
      })
    })
  }

  # Type
  varnames <- sapply(unlist(lapply(dc, function(x) x$targetVarExpr)), deparse)
    type <- unlist(lapply(dc, function(x) x$type))
    names(type) <- varnames
    determ_pars <- names(type)[type == "determ"]
    stoch_pars <- names(type)[type == "stoch"]
    
    if(!includeDeterm){
      out <- lapply(out, function(x){
        lapply(x, function(y){
          yout <- lapply(y, function(z){
            zout <- z[! z %in% determ_pars]
            if(length(zout) == 0) zout <- NULL
            zout
          })
        })
      })
    }

  if(!includeStoch){
    out <- lapply(out, function(x){
      lapply(x, function(y){
        yout <- lapply(y, function(z){
          zout <- z[! z %in% stoch_pars]
          if(length(zout) == 0) zout <- NULL
          zout
        })
      })
    })
  }

  # Cleanup
  lapply(out, function(x){
    lapply(x, function(y){
      y <- unlist(y)
      if(length(y) == 0) x <- NULL
      names(y) <- NULL
      unique(y)
    })
  })
}

# Add macro-generated inits to inits only if they don't already exist
addMacroInits <- function(origInits, macroInits){
  
  outInits <- origInits

  # If no macro-generated inits, do nothing
  if(is.null(macroInits)) return(outInits)
  if(length(macroInits) == 0) return(outInits)

  # If original inits are null, make them
  if(is.null(origInits)) outInits <- list()

  # Get inits to add
  # Only inits not already in original inits list
  addInits <- macroInits[!names(macroInits) %in% names(origInits)]
  
  # If there are some new inits, add them
  if(length(addInits) > 0){
    outInits[names(addInits)] <- addInits
  }
  
  outInits
}

# Remove extra brackets in BUGS code (typically resulting from macros)
removeExtraBrackets <- function(code){
  as.call(removeExtraBracketsInternal(code))
}

removeExtraBracketsInternal <- function(code){
  unlist(lapply(code, function(x){
    if(length(x) == 1) return(x)                       
    if(x[[1]] == "{") x <- as.list(x)[2:length(x)]
    if(is.list(x)){
      x <- removeExtraBracketsInternal(x)
    } else if(x[[1]] == "for"){
      x[[4]] <- removeExtraBrackets(x[[4]])
    }
    x
  }))
}
