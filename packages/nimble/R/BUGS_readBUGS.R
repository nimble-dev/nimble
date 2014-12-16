# code for creating BUGS model from a variety of input formats
# pieces written by Daniel Turek and Christopher Paciorek

BUGSmodel <- function(code, name, constants=list(), dimensions=list(), data=list(), inits=list(), returnModel=FALSE, where=globalenv(), debug=FALSE) {
    if(missing(name)) name <- deparse(substitute(code))
    if(length(constants) && sum(names(constants) == ""))
      stop("BUGSmodel: 'constants' must be a named list")
    if(length(dimensions) && sum(names(dimensions) == ""))
      stop("BUGSmodel: 'dimensions' must be a named list")
    if(length(data) && sum(names(data) == ""))
      stop("BUGSmodel: 'data' must be a named list")
    md <- modelDefClass$new(name = name)
    md$setupModel(code=code, constants=constants, dimensions=dimensions, debug=debug)
    if(!returnModel) return(md)
    model <- md$newModel(data=data, inits=inits, where=where)
}


#' Create a NIMBLE BUGS model
#' 
#' processes BUGS model code and optional constants, data, and initial values. Returns a NIMBLE model or model definition
#' 
#' @param code code for the model in the form returned by \code{modelCode} or (equivalently) \code{quote}
#' @param name optional character vector giving a name of the model for internal use.  If omitted, a name will be provided.
#' @param constants named list of constants in the model (not including data values).  Constants cannot be subsequently modified.
#' @param dimensions named list of dimensions for variables.  Only needed for variables used with empty indices in model code that are not provided in constants or data.
#' @param returnDef logical indicating whether the model should be returned (FALSE) or just the model definition (TRUE). 
#' @param debug logical indicating whether to put the user in a browser for debugging.  Intended for developer use.
#' @param where argument passed to \code{setRefClass}, indicating the environment in which the reference class definitions generated for the model and its modelValues should be created.  This is needed for managing package namespace issues during package loading and does not normally need to be provided by a user. 
#' @param data named list of values for the data nodes.  Data values can be subsequently modified.  Providing this argument also flags nodes as having data for purposes of algorithms that inspect model structure.
#' @author NIMBLE development team
#' @export
#' @details
#' See the User Manual or \code{help(modelBaseClass)} for information about manipulating NIMBLE models created by \code{nimbleModel}, including methods that operate on models, such as \code{getDependencies}.
#'
#' The user may need to provide dimensions for certain variables as in some cases NIMBLE cannot automatically determine the dimensions and sizes of variables. See the User Manual for more information.
#' @examples
#' modelCode <- BUGScode({
#'  x ~ dnorm(mu, sd = 1)
#'  mu ~ dnorm(0, sd = prior_sd)
#' })
#' constants = list(prior_sd = 1)
#' data = list(x = 4)
#' dimensions = list(
#' Rmodel <- nimbleModel(modelCode, constants = constants, data = data)
nimbleModel <- function(code, name, constants=list(), dimensions=list(), data=list(), inits=list(), returnDef = FALSE, where=globalenv(), debug=FALSE)
    BUGSmodel(code, name, constants, dimensions, data, inits, returnModel = !returnDef, where, debug)


#' Turn BUGS model code into an object for use in \code{nimbleModel} or \code{readBUGSmodel}
#'
#' Simply keeps model code as an R call object, the form needed by \code{nimbleModel} and optionally usable by \code{readBUGSmodel}
#' 
#' @param code expression providing the code for the model 
#' @author Daniel Turek
#' @export
#' @details It is equivalent to use the R function \code{quote}.  \code{modelCode} is simply provided as a more readable alternative for NIMBLE users not familiar with \code{quote}.
#' @examples
#' modelCode <- BUGScode({
#'  x ~ dnorm(mu, sd = 1)
#'  mu ~ dnorm(0, sd = prior_sd)
#' })
nimbleCode <- function(code) {
  code <- substitute(code)
  return(code)
}


BUGScode <- nimbleCode

processVarBlock <- function(lines) {
  # processes a var block from a BUGS file, determining variable names, dimensions, and sizes
  # at this point, sizes may have unevaluated variables in them

  # helper functions
  getDim <- function(vec) {
    if(length(vec) > 1)
      return(length(strsplit(vec[2], ";")[[1]])) else return(0)
  }
  
  getSize <- function(vec) {
    if(length(vec) > 1)
      return(strsplit(vec[2], ";"))  else return("0")
  }

  lines <- gsub("#.*", "", lines)
  lines <- gsub(";", "", lines)
  lines <- gsub("[[:space:]]", "", lines)
  lines <- paste(lines, collapse = "")
  # replace commas in brackets so can split variables
  chars <- strsplit(lines, integer(0))[[1]]
  nch <- length(chars)
  inBrackets <- FALSE
  for(i in seq_len(nch)) {
    if(inBrackets && chars[i] == ",")
      chars[i] = ";"
    if(chars[i] == "[") inBrackets <- TRUE
    if(chars[i] == "]") inBrackets <- FALSE
  }
  lines <- paste(chars, collapse = "")
  lines <- gsub("\\]", "", lines)
  
  pieces <- unlist(strsplit(lines, ","))
  pieces <- strsplit(pieces, "\\[")
  # variable names are in front of '[' (if there is an '[')
  varNames <- sapply(pieces, "[[", 1)
  dim <- sapply(pieces, getDim)
  size <- sapply(pieces, getSize)
  names(dim) <- varNames
  names(size) <- varNames
  return(list(varNames = varNames, dim = dim, size = size))
}

processModelFile <- function(fileName) {
  # processes a BUGS model file (.bug), splitting into var, data, and code blocks
  
  codeLines <- readLines(fileName)
  # extract lines corresponding to var, data, code blocks
  # var used for dimension info and data lines sourced in environment of the data input file objects
  codeLines <- paste(codeLines, collapse = "\n")
  codeLines <- gsub("/\\*.*?\\*/", "", codeLines)  # remove C-style comment blocks
  codeLines <- gsub("#.*?\n", "\n", codeLines) # remove R-style comments
  codeLines <- paste("\n", codeLines, collapse = "") # make sure first block occurs after a \n so regex below works ok; this allows me to not mistakenly find 'var', 'data', etc as names of nodes
  varBlockRegEx = "\n\\s*var\\s*(\n*.*?)(\n+\\s*(data|model|const).*)"
  dataBlockRegEx = "\n\\s*data\\s*\\{(.*?)\\}(\n+\\s*(var|model|const).*)"
  modelBlockRegEx = "\n\\s*model\\s*\\{(.*?)\\}\\s*\n+\\s*(var|data|const).*"

  if(length(grep(varBlockRegEx, codeLines))) {
    varLines <- gsub(varBlockRegEx, "\\1", codeLines)
    codeLines <- gsub(varBlockRegEx, "\\2", codeLines)
  } else varLines = NULL
  if(length(grep(dataBlockRegEx, codeLines))) {
    dataLines <- gsub(dataBlockRegEx, "\\1", codeLines)
    codeLines <- gsub(dataBlockRegEx, "\\2", codeLines)
  } else dataLines = NULL
  if(!length(grep(modelBlockRegEx, codeLines))) # model block is last block
    modelBlockRegEx = "model\\s*\\{(.*?)\\}\\s*\n*\\s*$"
  modelLines <- gsub(modelBlockRegEx, "\\1", codeLines)  # removes 'model' and whitespace at begin/end
  modelLines <- paste("{\n", modelLines, "\n}\n", collapse = "")

  return(list(modelLines = modelLines, varLines = varLines, dataLines = dataLines))
}


mergeMultiLineStatementsAndParse <- function(text) {
  # deals with BUGS syntax that allows multi-line statements where first line appears
  # to be valid full statement (e.g., where '+' starts the 2nd line)
  text <- unlist( strsplit(text, "\n") )  
  firstNonWhiteSpaceIndex <- regexpr("[^[:blank:]]", text)
  firstNonWhiteSpaceChar <- substr(text, firstNonWhiteSpaceIndex, firstNonWhiteSpaceIndex)
  mergeUpward <- firstNonWhiteSpaceChar %in% c('+', '-', '*', '/')
  if(length(text) > 1) {
    for(i in seq.int(length(text), 2, by = -1)) {
      if(mergeUpward[i]) {
        text[i-1] <- paste(text[i-1], substring(text[i], firstNonWhiteSpaceIndex[i]) )
      }
    }
  }
  text <- text[!mergeUpward]
  return(parse(text = text)[[1]])
}


#' Create a NIMBLE BUGS model from a variety of input formats, including BUGS model files
#' 
#' \code{readBUGSmodel} processes inputs providing the model and values for constants, data, initial values of the model in a variety of forms, returning a NIMBLE BUGS R model
#' 
#' @param model one of (1) a character string giving the file name containing the BUGS model code, with relative or absolute path, (2) an R function whose body is the BUGS model code, or (3) the output of \code{modelCode}. If a file name, the file can contain a 'var' block and 'data' block in the manner of the JAGS versions of the BUGS examples but should not contain references to other input data files nor a const block. The '.bug' or '.txt' extension can be excluded.
#' @param data (optional) (1) character string giving the file name for an R file providing the input constants and data as R code [assigning individual objects or as a named list], with relative or absolute path, or (2) a named list providing the input constants and data. If neither is provided, the function will look for a file named \{modelName\}-data including extensions .R, .r, or .txt.
#' @param inits (optional) (1) character string giving the file name for an R file providing the input constants and data as R code [assigning individual objects or as a named list], with relative or absolute path, or (2) a named list providing the input constants and data
#' @param dir (optional) character string giving the directory where the (optional) files are located
#' @param useInits boolean indicating whether to set the initial values, either based on \code{inits} or by finding the '-inits' file corresponding to the input model file
#' @return return returns a NIMBLE BUGS R model
#' @details Note that \code{readBUGSmodel} should handle most common ways of providing information on a model as used in BUGS and JAGS but does not handle input model files that refer to additional files containing data. Please see the BUGS examples provided with JAGS (\url{http://sourceforge.net/projects/mcmc-jags/files/Examples/}) for examples of supported formats. Also, \code{readBUGSmodel} takes both constants and data via the 'data' argument, unlike \code{nimbleModel}, in which these are distinguished. The reason for allowing both to be given via 'data' is for backwards compatibility with the BUGS examples, in which constants and data are not distinguished.
#' @author Christopher Paciorek
#' @export
#' @examples
#' modelCode <- BUGScode({
#'  x ~ dnorm(mu, sd = 1)
#'  mu ~ dnorm(0, sd = prior_sd)
#' })
#' data = list(prior_sd = 1, x = 4)
#' Rmodel <- readBUGSmodel(modelCode, data = data, inits = list(mu = 0))
#' Rmodel$setData(data['x'])
#' Rmodel[['mu']]
#' Rmodel$nodes[['x']]$calculate()
readBUGSmodel <- function(model, data = NULL, inits = NULL, dir = NULL, useInits = TRUE, useData = TRUE, debug = FALSE) {

  # helper function
  doEval <- function(vec, env) {
    out <- rep(0, length(vec))
    if(vec[1] == "0") return(numeric(0))
    for(i in seq_along(vec))
      out[i] <- eval(parse(text = vec[i]), env)
    return(out)
  }

  dimOrLengthAlt <- function(x) {
  # returns sizes of vectors/matrices/arrays with difference from dimOrLength that the value for a scalar is 1
    tmp <- dimOrLength(x)
    if(!length(tmp)) tmp <- 1
    return(tmp)
  }


  # process model information

  modelFileOutput <- modelName <- NULL
  if(is.function(model)) model <- mergeMultiLineStatementsAndParse(deparse(body(model)))
  if(is.character(model)) {
    if(!is.null(dir) && dir == "") modelFile <- model else modelFile <- file.path(dir, model)  # check for "" avoids having "/model.bug" when user provides ""
    modelName <- gsub("\\..*", "", basename(model))
    if(!file.exists(modelFile)) {
      possibleNames <- c(paste0(modelFile, '.bug'), paste0(modelFile, '.txt'))
      fileExistence <- file.exists(possibleNames)
      if(!sum(fileExistence)) {
        stop("readBUGSmodel: 'model' input does not reference an existing file.")
      } else {
        if(sum(fileExistence) > 1)
          warning("readBUGSmodel: multiple possible model files; using .bug file.")
        modelFile <- possibleNames[which(fileExistence)[1]]
      }
    }
    modelFileOutput <- processModelFile(modelFile)
    model <- mergeMultiLineStatementsAndParse(modelFileOutput$modelLines)
  }
  if(! class(model) == "{")
    stop("readBUGSmodel: cannot process 'model' input.")
    
  # process initial values

  if(useInits) {
    initsFile <-  NULL
    if(is.character(inits)) {
      initsFile <- file.path(dir, inits)
      if(!file.exists(initsFile)) 
        stop("readBUGSmodel: 'inits' input does not reference an existing file.")
    }
    if(is.null(inits)) {
      possibleNames <- c(
                         file.path(dir, paste0(modelName, "-init.R")),
                         file.path(dir, paste0(modelName, "-inits.R")),
                         file.path(dir, paste0(modelName, "-init.txt")),
                         file.path(dir, paste0(modelName, "-inits.txt")),
                         file.path(dir, paste0(modelName, "-init")),
                         file.path(dir, paste0(modelName, "-inits")))
      if(!Sys.info()['sysname'] %in% c("Darwin", "Windows")) # UNIX-like is case-sensitive
        possibleNames <- c(possibleNames,
                           file.path(dir, paste0(modelName, "-init.r")),
                           file.path(dir, paste0(modelName, "-inits.r")))
      fileExistence <- file.exists(possibleNames)
      if(sum(fileExistence) > 1)
        stop("readBUGSmodel: multiple possible initial value files; please pass as explicit 'inits' argument.")
      if(sum(fileExistence))
        initsFile <- possibleNames[which(fileExistence)[1]]
    }
    if(!is.null(initsFile)) {
      inits <- new.env()
      source(initsFile, inits)
      inits <- as.list(inits)
    }
  } else {
    inits <- NULL
  }
  if(!(is.null(inits) || is.list(inits)))
    stop("readBUGSmodel: invalid input for 'inits'.")

  # process var info
  varInfo <- NULL
  if(!is.null(modelFileOutput) && !is.null(modelFileOutput$varLines))
    varInfo = processVarBlock(strsplit(modelFileOutput$varLines, "\n")[[1]])

  # process data and constants input
  # since data and constants are mixed together in JAGS and BUGS, we take the same approach here (unfortunately)

  dataFile <-  NULL
  if(is.character(data)) {
    dataFile <- file.path(dir, data)
    if(!file.exists(dataFile)) 
      stop("readBUGSmodel: 'data' input does not reference an existing file.")
  }
  if(is.null(data)) {
    possibleNames <- c(
                       file.path(dir, paste0(modelName, "-data.R")),
                       file.path(dir, paste0(modelName, "-data.txt")),
                       file.path(dir, paste0(modelName, "-data")))
    if(!Sys.info()['sysname'] %in% c("Darwin", "Windows")) # UNIX-like is case-sensitive
      possibleNames <- c(possibleNames,
                         file.path(dir, paste0(modelName, "-data.r")))
    fileExistence <- file.exists(possibleNames)
    if(sum(fileExistence) > 1)
      stop("readBUGSmodel: multiple possible initial value files; please pass as explicit 'data' argument.")
    if(sum(fileExistence))
      dataFile <- possibleNames[which(fileExistence)[1]]
  }
  if(!is.null(dataFile)) {
    data <- new.env()
    source(dataFile, data)
  }
  if(is.list(data)) {
    if(length(data) && sum(names(data) == ""))
      stop("readBUGSmodel: 'data' must be a named list")
    data <- list2env(data)  # need as environment for later use 
  }
  
  if(!(is.null(data) || is.environment(data)))
    stop("readBUGSmodel: invalid input for 'data'.")

  if(!is.null(modelFileOutput) && !is.null(modelFileOutput$dataLines)) {
    # process data block in context of data objects
    if(is.null(data))
      data = new.env()
  
    # create vectors/matrices/arrays for all objects in var block in case data block tries to fill objects
    vars <- varInfo$varNames[varInfo$dim > 0]
    vars <- vars[!(vars %in% ls(data))]
    for(thisVar in vars) {
      dimInfo <- sapply(varInfo$size[[thisVar]], function(x) eval(parse(text = x), envir = data))
      tmp <- 0
      length(tmp) <- prod(dimInfo)
      dim(tmp) <- dimInfo
      assign(thisVar, tmp, envir = data)
    }

    origVars <- ls(data)
    eval(parse(text = modelFileOutput$dataLines), envir = data)
    newVars <- nf_assignmentLHSvars(parse(text = modelFileOutput$dataLines)[[1]])
    data <- as.list(data)[c(origVars, newVars)]
  } else {
    data <- as.list(data)
  }

  # determine dimensions from data list and varInfo
  dims <- lapply(data, dimOrLength)
  if(!is.null(varInfo)) {
    env <- data
    sizeInfo <- lapply(varInfo$size, doEval, env)
    # by default, use sizes based on actual data objects and ignore info
    # in the var block if it conflicts
    newNames <- names(sizeInfo)[!(names(sizeInfo) %in% names(data))]
    dims[newNames] <- sizeInfo[newNames]
  }

  if(length(dims) && sum(names(dims) == ""))
    stop("readBUGSmodel: something is wrong; 'dims' object is not a named list")

  # create R model
  # 'data' will have constants and data, but BUGSmodel is written to be ok with this
  # we can't separate them before building model as we don't know names of nodes in model
  Rmodel <- nimbleModel(model, ifelse(is.null(modelName), 'model', modelName), constants = data, dimensions = dims, debug = debug)

  # now provide values for data nodes from 'data' list
  dataNodes <- names(data)[(names(data) %in% Rmodel$getVarNames())]
  data <- data[dataNodes]
  names(data) <- dataNodes
  Rmodel$setData(data)

  if(!is.null(inits)) {
    varNames <- names(inits)[names(inits) %in% Rmodel$getVarNames()]
    for(varName in varNames) {
      # check for isData in case a node is a mix of data and non-data and inits are supplied such
      # that they would overwrite the data nodes without this check
      Rmodel[[varName]][!Rmodel$isData(varName)] <- inits[[varName]][!Rmodel$isData(varName)] 
    }
  }
  return(Rmodel)
}
