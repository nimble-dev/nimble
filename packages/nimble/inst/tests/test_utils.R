### Functions for testing math, called from test_math.R
require(testthat)

gen_runFun <- function(input) {
  runFun <- function() {}
  formalsList <- vector('list', length(input$inputDim))
  formalsList <- lapply(input$inputDim, function(x) parse(text = paste0("double(", x, ")"))[[1]])
  names(formalsList) <- paste0('arg', seq_along(input$inputDim))
  formals(runFun) <- formalsList
  tmp <- quote({})
  tmp[[2]] <- input$expr
  tmp[[3]] <- quote(return(out))
  tmp[[4]] <- parse(text = paste0("returnType(double(", input$outputDim, "))"))[[1]]
  body(runFun) <- tmp
  return(runFun)
}

make_input <- function(dim, size = 3, logicalArg) {
  if(!logicalArg) rfun <- rnorm else rfun <- function(n) { rbinom(n, 1, .5) }
  if(dim == 0) return(rfun(1))
  if(dim == 1) return(rfun(size))
  if(dim == 2) return(matrix(rfun(size^2), size))
  stop("not set for dimension greater than 2")
}

test_math <- function(input, verbose = TRUE, size = 3) {
  if(verbose) cat("### Testing", input$name, "###\n")
  runFun <- gen_runFun(input)
  nfR <- nimbleFunction(  # formerly nfGen
      #       setup = TRUE,
             run = runFun)
  #nfR <- nfGen()
  project <- nimbleProjectClass(NULL, name = 'foo')
  nfC <- compileNimble(nfR, project = project)

  nArgs <- length(input$inputDim)
  logicalArgs <- rep(FALSE, nArgs)
  if("logicalArgs" %in% names(input))
    logicalArgs <- input$logicalArgs
  
  arg1 <- make_input(input$inputDim[1], size = size, logicalArgs[1])
  if(nArgs == 2)
    arg2 <- make_input(input$inputDim[2], size = size, logicalArgs[2])
  if("Rcode" %in% names(input)) {
    eval(input$Rcode)
  } else {
    eval(input$expr)
  }
  if(nArgs == 2) {
    out_nfR = nfR(arg1, arg2)
    out_nfC = nfC(arg1, arg2)
  } else {
    out_nfR = nfR(arg1)
    out_nfC = nfC(arg1)
  }
  attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
  if(is.logical(out)) out <- as.numeric(out)
  if(is.logical(out_nfR)) out_nfR <- as.numeric(out_nfR)
  try(test_that(paste0("Test of math (direct R calc vs. R nimbleFunction): ", input$name), expect_that(out, equals(out_nfR))))
  try(test_that(paste0("Test of math (direct R calc vs. C nimbleFunction): ", input$name), expect_that(out, equals(out_nfC))))
  # unload DLL as R doesn't like to have too many loaded
  dyn.unload(project$cppProjects[[1]]$getSOName())
  invisible(NULL)
}


### Function for testing MCMC called from test_mcmc.R


test_mcmc <- function(example, model, data = NULL, inits = NULL,
                      verbose = TRUE, numItsR = 5, numItsC = 1000,
                      basic = TRUE, exactSample = NULL, results = NULL, resultsTolerance = NULL,
                      numItsC_results = numItsC,
                      resampleData = FALSE,
                      topLevelValues = NULL, seed = 0, mcmcControl = NULL, samplers = NULL, removeAllDefaultSamplers = FALSE, 
                      doR = TRUE, doCpp = TRUE, returnSamples = FALSE, name = NULL) {
  # There are three modes of testing:
  # 1) basic = TRUE: compares R and C MCMC values and, if requested by passing values in 'exactSample', will compare results to actual samples (you'll need to make sure the seed matches what was used to generate those samples)
  # 2) if you pass 'results', it will compare MCMC output to known posterior summaries within tolerance specified in resultsTolerance
  # 3) resampleData = TRUE: runs initial MCMC to get top level nodes then simulates from the rest of the model, including data, to get known parameter values, and fits to the new data, comparing parameter estimates from MCMC with the known parameter values

    # samplers can be given individually for each node of interest or as a vector of nodes for univariate samplers or list of vectors of nodes for multivariate samplers
    # e.g.,
    # multiple univar samplers: samplers(type = 'RW', target = c('mu', 'x'))
    # single multivar sampler: samplers(type = "RW_block", target = c('x[1]', 'x[2]'))
    # single multivar sampler: samplers(type = "RW_block", target = 'x')
    # multiple multivar samplers: samplers(type = "RW_block", target = list('x', c('theta', 'mu')))

       
  setSampler <- function(var, spec) {
      currentTargets <- sapply(spec$samplerSpecs, function(x) x$target)
    # remove already defined scalar samplers
      inds <- which(unlist(var$target) %in% currentTargets)
      spec$removeSamplers(inds, print = FALSE)
    # look for cases where one is adding a blocked sampler specified on a variable and should remove scalar samplers for constituent nodes
      inds <- which(sapply(unlist(var$target), function(x) Rmodel$expandNodeNames(x)) %in% currentTargets)
#      inds <- which(sapply(spec$samplerSpecs, function(x)
#          gsub("\\[[0-9]+\\]", "", x$target))
#                         %in% var$target)
      spec$removeSamplers(inds, print = FALSE)

      if(is.list(var$target) && length(var$target) == 1) var$target <- var$target[[1]]
      if(length(var$target) == 1 || (var$type == "RW_block" && !is.list(var$target))) 
          tmp <- spec$addSampler(type = var$type, target = var$target, control = var$control, print = FALSE) else tmp <- sapply(var$target, function(x) spec$addSampler(type = var$type, target = x, control = var$control, print = FALSE))
  }
  
  if(is.null(name)) {
      if(!missing(example)) {
          name <- example
      } else {
            if(is.character(model)) {
                name <- model
            } else {
                  name <- 'unnamed case'
              }
        }
  }

  cat("===== Starting MCMC test for ", name, ". =====\n", sep = "")
  
  if(!missing(example)) {
    # classic-bugs example specified by name
  	dir = getBUGSexampleDir(example)
    if(missing(model)) model <- example
    Rmodel <- readBUGSmodel(model, dir = dir, data = data, inits = inits, useInits = TRUE,
                            check = FALSE)
  } else {
    # code, data and inits specified directly where 'model' contains the code
    example = deparse(substitute(model))
    if(missing(model)) stop("Neither BUGS example nor model code supplied.")
    Rmodel <- readBUGSmodel(model, data = data, inits = inits, dir = "", useInits = TRUE,
                            check = FALSE)
  }


  if(doCpp) {
      Cmodel <- compileNimble(Rmodel)
  }
  if(!is.null(mcmcControl)) mcmcspec <- configureMCMC(Rmodel, control = mcmcControl) else mcmcspec <- configureMCMC(Rmodel)
  if(removeAllDefaultSamplers) mcmcspec$removeSamplers()
  
  if(!is.null(samplers)) {
      sapply(samplers, setSampler, mcmcspec)
      if(verbose) {
          cat("Setting samplers to:\n")
          mcmcspec$getSamplers()
      }
  }
  
  vars <- Rmodel$getDependencies(Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE), stochOnly = TRUE, includeData = FALSE, downstream = TRUE)
  vars <- unique(removeIndexing(vars))
  mcmcspec$addMonitors(vars, print = FALSE)
  
  Rmcmc <- buildMCMC(mcmcspec)
  if(doCpp) {
      Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  }
      
  if(basic) {
    ## do short runs and compare R and C MCMC output
      if(doR) {
          set.seed(seed);
          RmcmcOut <- try(Rmcmc$run(numItsR))
          if(!is(RmcmcOut, "try-error")) {
            RmvSample  <- nfVar(Rmcmc, 'mvSamples')
            R_samples <- as.matrix(RmvSample)
          } else R_samples <- NULL
      }
      if(doCpp) {
          set.seed(seed)
          Cmcmc$run(numItsC)
          CmvSample <- nfVar(Cmcmc, 'mvSamples')    
          C_samples <- as.matrix(CmvSample)
          ## for some reason columns in different order in CmvSample...
          C_subSamples <- C_samples[seq_len(numItsR), attributes(R_samples)$dimnames[[2]], drop = FALSE]
      }

      if(doR && doCpp && !is.null(R_samples)) {
          context(paste0("testing ", example, " MCMC"))
          try(
              test_that(paste0("test of equality of output from R and C versions of ", example, " MCMC"), {
                  expect_that(R_samples, equals(C_subSamples), info = paste("R and C posterior samples are not equal"))
              })
              )
      }
      if(is.null(R_samples)) {
          cat("R MCMC failed.\n")
      }

      if(doCpp) {
          if(!is.null(exactSample)) {
              for(varName in names(exactSample))
                  try(
                      test_that(paste("Test of MCMC result against known samples for", example, ":", varName), {
                          expect_that(round(C_samples[seq_along(exactSample[[varName]]), varName], 8), equals(round(exactSample[[varName]], 8))) })
                      )
          }
      }
      
    summarize_posterior <- function(vals) 
      return(c(mean = mean(vals), sd = sd(vals), quantile(vals, .025), quantile(vals, .975)))

      if(doCpp) {
          if(verbose) {
              start <- round(numItsC / 2) + 1
              try(print(apply(C_samples[start:numItsC, , drop = FALSE], 2, summarize_posterior)))
          }
      }
  }

  ## assume doR and doCpp from here down
  if(!is.null(results)) { 
     # do (potentially) longer run and compare results to inputs given
    set.seed(seed)
    Cmcmc$run(numItsC_results)
    CmvSample <- nfVar(Cmcmc, 'mvSamples')
    postBurnin <- (round(numItsC_results/2)+1):numItsC_results
    C_samples <- as.matrix(CmvSample)[postBurnin, , drop = FALSE]
     for(metric in names(results)) {
      if(!metric %in% c('mean', 'median', 'sd', 'var', 'cov'))
        stop("Results input should be named list with the names indicating the summary metrics to be assessed, from amongst 'mean', 'median', 'sd', 'var', and 'cov'.")
      if(metric != 'cov') {
        postResult <- apply(C_samples, 2, metric)       
        for(varName in names(results[[metric]])) {
          varName <- gsub("_([0-9]+)", "\\[\\1\\]", varName) # allow users to use theta_1 instead of "theta[1]" in defining their lists
          samplesNames <- dimnames(C_samples)[[2]]
          if(!grepl(varName, "\\[", fixed = TRUE)) 
             samplesNames <- gsub("\\[.*\\]", "", samplesNames)
          matched <- which(varName == samplesNames) 
          diff <- abs(postResult[matched] - results[[metric]][[varName]])
          for(ind in seq_along(diff)) {
            strInfo <- ifelse(length(diff) > 1, paste0("[", ind, "]"), "")
            try(
              test_that(paste("Test of MCMC result against known posterior for", example, ":",  metric, "(", varName, strInfo, ")"), {
                expect_that(diff[ind], is_less_than(resultsTolerance[[metric]][[varName]][ind]))
              })
              )
          }
        }
      } else  { # 'cov'
        for(varName in names(results[[metric]])) {
          matched <- grep(varName, dimnames(C_samples)[[2]], fixed = TRUE)
          postResult <- cov(C_samples[ , matched])
           # next bit is on vectorized form of matrix so a bit awkward
          diff <- c(abs(postResult - results[[metric]][[varName]]))
          for(ind in seq_along(diff)) {
            strInfo <- ifelse(length(diff) > 1, paste0("[", ind, "]"), "")
            try(
              test_that(paste("Test of MCMC result against known posterior for", example, ":",  metric, "(", varName, ")", strInfo), {
                expect_that(diff[ind], is_less_than(resultsTolerance[[metric]][[varName]][ind]))
              })
              )
          }
        }
      }
    }
  }
  if(returnSamples) {
    if(exists('CmvSample'))
      returnVal <- as.matrix(CmvSample)
  } else returnVal <- NULL

  if(resampleData) {
    topNodes <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
    topNodesElements <- Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE,
                                            returnScalarComponents = TRUE)
    if(is.null(topLevelValues)) {
      postBurnin <- (round(numItsC/2)):numItsC
      if(is.null(results) && !basic) {
      # need to generate top-level node values so do a basic run
        set.seed(seed)
        Cmcmc$run(numItsC)
        CmvSample <- nfVar(Cmcmc, 'mvSamples')
        C_samples <- as.matrix(CmvSample)[postBurnin, ]
      }
      topLevelValues <- as.list(apply(C_samples[ , topNodesElements, drop = FALSE], 2, mean))
    }
    if(!is.list(topLevelValues)) {
      topLevelValues <- as.list(topLevelValues)
      if(sort(names(topLevelValues)) != sort(topNodesElements))
        stop("Values not provided for all top level nodes; possible name mismatch")
    }
    sapply(topNodesElements, function(x) Cmodel[[x]] <- topLevelValues[[x]])
    # check this works as side effect
    nontopNodes <- Rmodel$getDependencies(topNodes, self = FALSE, includeData = TRUE, downstream = TRUE, stochOnly = FALSE)
    # nonDataNodes <- Rmodel$getDependencies(topNodes, self = TRUE, includeData = FALSE, downstream = TRUE, stochOnly = TRUE)
    nonDataNodesElements <- Rmodel$getDependencies(topNodes, self = TRUE, includeData = FALSE, downstream = TRUE, stochOnly = TRUE, returnScalarComponents = TRUE)
    dataVars <- unique(removeIndexing(Rmodel$getDependencies(topNodes, dataOnly = TRUE, downstream = TRUE)))
    set.seed(seed)
    Cmodel$resetData()
    simulate(Cmodel, nontopNodes)

    dataList <- list()
    for(var in dataVars) {
      dataList[[var]] <- values(Cmodel, var)
      if(Cmodel$modelDef$varInfo[[var]]$nDim > 1)
        dim(dataList[[var]]) <- Cmodel$modelDef$varInfo[[var]]$maxs
    }
    Cmodel$setData(dataList)

    trueVals <- values(Cmodel, nonDataNodesElements)
    names(trueVals) <- nonDataNodesElements
    set.seed(seed)
    Cmcmc$run(numItsC_results)
    CmvSample <- nfVar(Cmcmc, 'mvSamples')
    
    postBurnin <- (round(numItsC_results/2)):numItsC
    C_samples <- as.matrix(CmvSample)[postBurnin, nonDataNodesElements, drop = FALSE]
    interval <- apply(C_samples, 2, quantile, c(.025, .975))
    interval <- interval[ , names(trueVals)]
    covered <- trueVals <= interval[2, ] & trueVals >= interval[1, ]
    coverage <- sum(covered) / length(nonDataNodesElements)
    tolerance <- 0.15
    if(verbose) 
      cat("Coverage for model", example, "is", coverage*100, "%.\n")
    miscoverage <- abs(coverage - 0.95)
    try(
      test_that(paste("Test of MCMC coverage on known parameter values for:", example), {
                expect_that(miscoverage, is_less_than(tolerance))
              })
      )
    if(miscoverage > tolerance || verbose) {
      cat("True values with 95% posterior interval:\n")
      print(cbind(trueVals, t(interval), covered))
    }
  }

  cat("===== Finished MCMC test for ", name, ". =====\n", sep = "")

  if(doCpp) {
    # DTL is looking at issue with segFault when close R
  #  dyn.unload(getNimbleProject(Rmodel)$cppProjects[[1]]$getSOName())
  #  dyn.unload(getNimbleProject(Rmcmc)$cppProjects[[1]]$getSOName())
  }
  return(returnVal)

}
