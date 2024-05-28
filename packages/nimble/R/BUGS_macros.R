
#' EXPERIMENTAL: Turn a function into a model macro
#' A model macro expands one line of code in a nimbleModel into one or
#' more new lines.  This supports compact programming by defining
#' re-usable modules.  \code{model_macro_builder} takes as input a
#' function that constructs new lines of model code from the original
#' line of code.  It returns a function suitable for internal use by
#' \code{nimbleModel} that arranges arguments for input function.  Macros
#' are an experimental feature and are available only after setting
#' \code{nimbleOptions(enableModelMacros = TRUE)}.
#'
#' @param fun A function written to construct new lines of model code (see below).
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
#' won't make sense to anticipate a declaration using \code{~} or \code{<-}.
#' Instead, arguments from an arbitrary call will be passed as separate arguments.
#' For example, the line \code{macro1(arg1 = z[1:10], arg2 = "hello")} will be
#' passed to \code{fun} as \code{fun(arg1 = z[1:10], arg2 = "hello")}.
#'
#' In addition, the final two arguments of \code{fun} must be called \code{modelInfo}
#' and \code{.env} respectively. 
#' 
#' During macro processing, \code{nimbleModel} passes a named list to the \code{modelInfo} 
#' argument of \code{fun} containing, among other things, elements called
#' \code{constants} and \code{dimensions}. Macro developers can modify these
#' two elements (for example, to add a new constant needed for a macro) and
#' these changes will be reflected in the final model object. Note that currently
#' it is not possible for a macro to modify the data. Furthermore, if your macro add a new element to the
#' constants that \code{nimbleModel} then moves to the data, this new data will not be retained
#' in the final model object and thus will not be usable. 
#'
#' \code{nimbleModel} passes the R environment from which \code{nimbleModel} was
#' called to the \code{.env} argument.
#'
#' The \code{fun} function must return a named list with two elements:
#' \code{code}, the replacement code, and \code{modelInfo}, the \code{modelInfo} 
#' list described above. \code{modelInfo} must be in the output even if the macro
#' does not modify it.
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
#' @return A list of class \code{model_macro} with one element called \code{process},
#' which contains the macro function suitable for use by \code{nimbleModel}.
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
#'
#' ## Example 3: Macro that modifies constants
#' new_constant <- model_macro_builder(
#'    function(stoch, LHS, RHS, modelInfo, .env) {
#'      # number of elements
#'      n <- as.numeric(length(modelInfo$constants[[deparse(LHS)]]))
#'      code <- substitute({
#'        for (i in 1:N){
#'          L[i] ~ dnorm(mu[i], 1)
#'        }
#'      }, list(L = LHS, N = n))
#'
#'      # Add a new constant mu
#'      modelInfo$constants$mu <- rnorm(n, 0, 1)
#'
#'      list(code = code, modelInfo = modelInfo)
#'    },
#'    use3pieces = TRUE,
#'    unpackArgs = TRUE
#' )
#'
#' const <- list(y = rnorm(10))
#' code <- nimbleCode({
#'  y ~ new_constant()
#' })
#' 
#' mod <- nimbleModel(code = code, constants=const)
#' mod$getCode()
#' mod$getConstants() # new constant is here
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
## Called "internal" because there's a similarly-named wrapper function below that
## calls this function but also does other things not in the recursion
processMacrosInternal <- function(code,
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
              macroOutput <- processMacrosInternal(code = code[[i]],
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
      macroOutput <- processMacrosInternal(code[[4]],
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
    if(exists(possibleMacroName, envir = env)) { ## may need to provide an envir argument
        possibleMacro <- get(possibleMacroName, envir = env) ## ditto
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
            
            # Check for newly created parameters from macros and put them in a list element
            newPars <- list(newMacroPars(code, expandedInfo$code))
            # Name the list element after the source macro
            names(newPars) <- possibleMacroName

            # Add automatic comments showing code added by macros
            # Includes one line indicating the start of the code generated (macroComment)
            # by a particular macro, and one line showing the end (macroEnd)
            # Comments are simply character strings
            if(getNimbleOption("enableMacroComments")){
              # If this option is also set, then print the entire original
               # line of code with the macro in the comments instead of just the macro name
              if(getNimbleOption("codeInMacroComments")){
                codeOrName <- safeDeparse(code, warn = TRUE)
              } else {
                codeOrName <- possibleMacroName
              }
              spacer <- ""
              hashes <- "#"
              if(length(recursionLabels > 0)){
                spacer <- paste(rep("  ", length(recursionLabels)), collapse="")
                hashes <- paste(rep("#", length(recursionLabels)+1), collapse="")
              }
              macroComment <- paste0(spacer, hashes, " ", codeOrName)
              macroEnd <- paste0(spacer, hashes, " ----")

              # Add the starting and ending comments to the code
              macroStartLine <- substitute(MACRO, list(MACRO = macroComment))
              macroEndLine <- substitute(END, list(END = macroEnd))
              expandedInfo$code <- as.call(c(list(quote(`{`)), list(macroStartLine, expandedInfo$code, macroEndLine)))
            }
            
            # Add the new macro parameters to modelInfo
            curPars <- modelInfo$parameters
            expandedInfo$modelInfo$parameters <- c(curPars, newPars)

            ## Return object is a list so we can possibly extract other
            ## content in the future.  We recurse on the returned code
            ## to expand macros that it might contain.
            macroOutput <- processMacrosInternal(expandedInfo$code,
                                           modelInfo = expandedInfo$modelInfo,
                                           env = env,
                                           recursionLabels = c(recursionLabels, possibleMacroName)
                                           )
            return(list(code=macroOutput$code, modelInfo=macroOutput$modelInfo))
        }
    }
    list(code=code, modelInfo=modelInfo)
}

# Create index generator for all macros to use
# Used in e.g. a for loop macro that needs to create an index parameter
macroIndexCreator <- labelFunctionCreator("i")

# Top-level function that does all the macro processing
# Code: input code from user
# modelInfo: List containing key model info like constants, dimensions
# env: environment to operate in
codeProcessModelMacros <- function(code, modelInfo, env){
  # Reset the index generator and insert it into modelInfo
  macroIndexCreator(reset = TRUE)
  modelInfo$indexCreator <- macroIndexCreator
  # No macro generated parameters before any macros run, so parameters = empty list
  modelInfo$parameters <- list()
  # Recursively step through the code expanding macros
  # and adding information to modelInfo such as updated constants and
  # parameter names
  macroOutput <- processMacrosInternal(code=code, modelInfo=modelInfo, env=env)
  # Clean up extra brackets in output code
  macroOutput$code <- removeExtraBrackets(macroOutput$code)
  # Convert factors in constants to numeric
  macroOutput$modelInfo$constants <- convertFactorConstantsToNumeric(macroOutput$modelInfo$constants)  
  # Remove intermediate parameters from list of generated parameters
  # and check for duplicates
  macroOutput$modelInfo$parameters <- checkMacroPars(macroOutput$modelInfo$parameters,
                                                     code, macroOutput$code)
  list(code = macroOutput$code, modelInfo = macroOutput$modelInfo)
}

# Function that converts any factors or character vectors/matrices/arrays
# in the constants into numeric so they can be used by nimble
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

# Functions for finding parameters generated by macros-------------------------

# Figure out which new parameters a macro added
# This is done by finding the parameters in the starting (pre-macro) code line
# and the ending (post-expanded-macro) code, 
# then comparing the two
newMacroPars <- function(startCode, endCode){
  startPar <- getParametersFromCode(startCode) # starting parameters
  endPar <- getParametersFromCode(endCode)     # ending parameters
  # Find the differences, keeping the LHS and RHS of assignments separated
  LHS <- endPar$LHS[! endPar$LHS %in% unlist(startPar)]
  if(length(LHS) == 0) LHS <- NULL
  RHS <- endPar$RHS[! endPar$RHS %in% unlist(startPar)]
  if(length(RHS) == 0) RHS <- NULL
  list(LHS = LHS, RHS = RHS)
}

# Get all parameters from a code chunk
# Uses recursive function below internally
# Separates parameters by LHS and RHS of assingments
getParametersFromCode <- function(code){
  out <- getParametersFromCodeInternal(code)
  LHS = as.character(findAllListElementsByName(out, "LHS"))
  RHS = as.character(findAllListElementsByName(out, "RHS"))
  list(LHS = unique(LHS), RHS = unique(RHS))
}

# Recursive function to work through a code chunk and get all parameters
# Used above
getParametersFromCodeInternal <- function(code){
  # Check for length 0 code
  if(length(code) == 0) return(list(LHS = NULL, RHS = NULL))
  # If this code is a comment, return NULL
  if(is.character(code)) return(list(LHS = NULL, RHS = NULL))
  # If the code is just class(name) all by itself (no assignment), return it
  # Calling these "RHS" even though there is no sides
  if(is.name(code)){
    return(list(LHS = NULL, RHS = safeDeparse(code, warn = TRUE)))
  }
  # If we have several lines of code, iterate through them recursively
  if(code[[1]] == "{"){
    # If empty brackets return now
    if(length(code) == 1) return(list(LHS = NULL, RHS = NULL))
    return(lapply(as.list(code)[2:length(code)], getParametersFromCodeInternal))
  }
  # If a for loop, iterate through the lines recursively
  if(code[[1]] == "for"){
    return(lapply(as.list(code)[2:length(code)], getParametersFromCodeInternal))
  } else {
    # Otherwise  we should have a single line of code which is an assignment
    if(as.character(code[[1]]) %in% c("~", "<-")){
      # Separate the assignment into LHS and RHS side "pieces" and get the
      # parameters from each
      RHS <- getParsFromCodePiece(code[[3]])
      LHS <- getParsFromCodePiece(code[[2]])
    } else {
      RHS <- getParsFromCodePiece(code)
      LHS <- NULL
    }
  }

  list(LHS = LHS, RHS = RHS)
}

# Wrapper function for getting generated parameters that runs the
# recursive function below and does some cleanup on the final output
getParsFromCodePiece <- function(code){
  out <- all.vars(code)
  if(length(out) == 0) return(NULL)
  unique(out) # I don't think this is necessary, but just to be sure
}

# Recursively earch through a (possibly nested) list and find all elements
# with a particular name
findAllListElementsByNameInternal <- function(inputList, name){
  if(is.list(inputList) && name %in% names(inputList)){
    return(inputList[[name]])
  } else if(is.list(inputList)){
    return(lapply(inputList, findAllListElementsByNameInternal, name = name))
  }
  return(NULL)
}

# Calls recursive function above and does some cleanup on the output
findAllListElementsByName <- function(inputList, name){
  stopifnot(is.list(inputList))
  out <- findAllListElementsByNameInternal(inputList, name)
  unlist(out)
}

# Function to do some final cleanup on the list of generated macro parameters
# Remove "intermediate" parameters
# "Intermediate" parameters could occur if you have macros creating other macros.
# For example if you have a code-creation pattern like
# macro(parameter1) --> macro(parameter2) --> final code
# parameter2 would not show up in the final code and would not be something
# nimble would be saving/sampling. Thus my assumption is that it would not be of
# interest to users/developers and should not end up in the output of getMacroParameters().
# also check for parameters generated more than once on the LHS
# and a give a message
# Sometimes this could be a mistake, or it could be correct if e.g.
# you have two macros that generate the following two lines of code
# alpha[1:2] <- x
# alpha[3] <- y
checkMacroPars <- function(parameters, startCode, endCode){

  # If no macros return empty list
  if(length(parameters) == 0) return(parameters)
  
  # Find new parameters all at once: startCode is the entire original code
  # provided by the user, and endCode is the final code after all macros expanded
  all_pars <- newMacroPars(startCode, endCode)

  if(length(all_pars) == 0) return(NULL)

  # Parameters that exist in the input parameter list (generated by going through
  # each macro separately) but not in the list generated by comparing just the
  # start and final, complete end code are "intermediate" parameters that the
  # user shouldn't actually care about.
  # Remove intermediate parameters that don't end up final code
  final_pars <- lapply(parameters, function(x){
    list(LHS = x$LHS[x$LHS %in% all_pars$LHS],
         RHS = x$RHS[x$RHS %in% all_pars$RHS])
  })
  
  # Check for duplicate parameters on the LHS from previous macros and warn
  if(length(final_pars) > 1){
  for (i in 2:length(final_pars)){
    this_macro <- final_pars[[i]]$LHS
    for (j in 1:(i-1)){
      dups <- this_macro %in% final_pars[[j]]$LHS
      if(any(dups)){

       msg <- paste0("  [Note] Macro '", names(final_pars[i]), "' declared LHS parameter(s): '",
                paste(this_macro[dups], collapse="', '"), "'\n         previously declared by macro '",
                names(final_pars[j]),"'.",
                "\n         This is not a problem if parameter indices don't overlap")
        messageIfVerbose(msg)

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
#' @param includeIndexPars Include index parameters generated for use in for loops
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
                              includeStoch = TRUE, includeIndexPars = FALSE) {
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
  if(!includeIndexPars){
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

  # Type (deterministic/stochastic)
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
      y <- unlist(y, use.names = FALSE)
      if(length(y) == 0) x <- NULL
      unique(y)
    })
  })
}

# End functions dealing with getting parameters generated by macros------------

# Macros can auto-generate appropriate inits
# We want to add these to the existing inits, but also not overwrite
# any inits that already exist
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

# Remove extra curly brackets in BUGS code
# Typically macros generate a bunch of extra layers of {} that we don't want

# Call recursive function below and turn the result back into code
removeExtraBrackets <- function(code){
  as.call(removeExtraBracketsInternal(code))
}

# Turn code into list, remove extra brackets, and return the result
# Operates recursively
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
