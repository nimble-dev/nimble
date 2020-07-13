## Mark tests that are know to fail with `if(RUN_FAILING_TESTS)`.
## By default these tests will not be run, but we will occasionally clean up by running them with
## At the moment only used in test-optim.R; otherwise we have
## mechanisms in various test files for wrapping failing tests
## in expect_failure.
## $ RUN_FAILING_TESTS=1 Rscript test-my-stuff.R
RUN_FAILING_TESTS <- (nchar(Sys.getenv('RUN_FAILING_TESTS')) != 0)

## We can get problems on Windows with system(cmd) if cmd contains
## paths with spaces in directory names.  Solution is to make a cmd
## string with only local names (to avoid spaces) and indicate what
## directory it should be done in.  setwd(dir) should handle dir with
## paths with spaces ok.
system.in.dir <- function(cmd, dir = '.') {
  curDir <- getwd()
  on.exit(setwd(curDir), add = TRUE)
  setwd(dir)
  if(.Platform$OS.type == "windows")
    shell(shQuote(cmd))
  else
    system(cmd)
}

## This sets up sink to also capture messages (in particular warnings).
sink_with_messages <- function(file, ...) {
  sinkfile <- file(file, open = 'wt')
  sink(sinkfile)
  sink(sinkfile, type = 'message')
}

## This is useful for working around scoping issues with nimbleFunctions using other nimbleFunctions.
temporarilyAssignInGlobalEnv <- function(value) {
  name <- deparse(substitute(value))
  assign(name, value, envir = .GlobalEnv)
  rmCommand <- substitute(remove(name, envir = .GlobalEnv))
  do.call('on.exit', list(rmCommand, add = TRUE), envir = parent.frame())
}

withTempProject <- function(code) {
  code <- substitute(code)
  project <- nimble:::nimbleProjectClass()
  nimbleOptions(nimbleProjectForTesting = project)
  on.exit({
    ## messages are suppressed by test_that, so "assign('.check', 1, globalenv())" can be used as a way to verify this code was called
    .project <- nimbleOptions('nimbleProjectForTesting')
    nimbleOptions(nimbleProject = NULL) ## clear this before clearCompiled step, in case clearCompiled() itself fails
    .project$clearCompiled()
  }, add = TRUE)
  eval(code)
}

expect_compiles <- function(..., info = NULL, link = FALSE, force01 = TRUE) {
  oldSCBL <- nimbleOptions('stopCompilationBeforeLinking')
  nimbleOptions(stopCompilationBeforeLinking = !link)
  oldForce01 <- nimbleOptions('force01')
  nimbleOptions(forceO1 = force01)
  on.exit({
    assign('.check', 1, globalenv())
    nimbleOptions(stopCompilationBeforeLinking = oldSCBL)
    nimbleOptions(forceO1 = oldForce01)
  }, add = TRUE)
  if(!link) {
    ans <- try(compileNimble(...)) ## expecting a thrown error
    expect_identical(as.character(ans), 'Error : safely stopping before linking\n', info = info)
  } else {
    ans <- compileNimble(...)
    ans
  }
}

gen_runFunCore <- function(input) {
  runFun <- function() {}
  formalsList <- input$args
  if(is.null(formalsList)) formalsList <- list()
  if(is.null(names(formalsList)))
    if(length(formalsList) > 0)
      names(formalsList) <- paste0('arg', seq_along(input$args))
  formals(runFun) <- formalsList
  tmp <- quote({})
  tmp[[2]] <- input$expr
  tmp[[3]] <- if(is.null(input[['return']]))
    quote(return(out))
  else
    input[['return']]
  tmp[[4]] <- substitute(returnType(OUT), list(OUT = input$outputType))
  body(runFun) <- tmp
  return(runFun)
}

## Indexes the names of a list of input lists for test_coreRfeature
indexNames <- function(x) {
  i <- 1
  lapply(x, function(z) {z$name <- paste(i, z$name); i <<- i + 1; z})
}

make_input <- function(dim, size = 3, logicalArg) {
  if(!logicalArg) rfun <- rnorm else rfun <- function(n) { rbinom(n, 1, .5) }
  if(dim == 0) return(rfun(1))
  if(dim == 1) return(rfun(size))
  if(dim == 2) return(matrix(rfun(size^2), size))
  stop("not set for dimension greater than 2")
}

wrap_if_matches <- function(pattern, string, wrapper, expr) {
  if (!is.null(pattern) && any(grepl(paste0('^', pattern, '$'), string))) {
    wrapper(expr)
  } else {
    expr
  }
}

weightedMetricFunc <- function(index, samples, weights, metric, samplesToWeightsMatch){
  samples <- samples[,index]
  weights <- exp(weights[,samplesToWeightsMatch[index]])/sum(exp(weights[,samplesToWeightsMatch[index]]))
  if(metric == "median"){
    ord <- order(samples)
    weights <- weights[ord]
    samples <- samples[ord]
    sumWts <- 0
    while(sumWts < .5){
      sumWts <- sumWts + weights[1]
      weights <- weights[-1]
    }
    return(samples[length(samples)-length(weights)])
  }
  wtMean <- weighted.mean(samples, weights)
  if(metric == "mean"){
    return(wtMean)
  }
  wtVar <- sum(weights*(samples - wtMean)^2)
  if(metric == "var"){
    return(wtVar)
  }
  if(metric == "sd"){
    return(sqrt(wtVar))
  }
}

compareFilesByLine <- function(trialResults, correctResults, main = "") {
  trialResults <- stripTestPlacementWarning(trialResults)
  correctResults <- stripTestPlacementWarning(correctResults)
  test_that(paste0(main, ': same number of output lines'),
            expect_equal(length(trialResults), length(correctResults)))
  
  linesToTest <- min(length(trialResults), length(correctResults))
  mapply(function(lineno, trialLine, correctLine) {
    test_that(paste0(main, ": output line #", lineno),
              expect_identical(trialLine, correctLine))
  }, 1:linesToTest, trialResults, correctResults)
  invisible(NULL)
}

compareFilesUsingDiff <- function(trialFile, correctFile, main = "") {
  if(main == "") main <- paste0(trialFile, ' and ', correctFile, ' do not match\n')
  diffOutput <- system2('diff', c(trialFile, correctFile), stdout = TRUE)
  test_that(paste0(main, paste0(diffOutput, collapse = '\n')),
            expect_true(length(diffOutput) == 0)
  )
  invisible(NULL)
}

stripTestPlacementWarning <- function(lines) {
  ## deal with Placing tests in `inst/tests/` is deprecated warning
  ## as it doesn't seem entirely predictable when/where it appears 
  coreLines <- grep("^Placing tests in", lines)
  addedLines <- lines[coreLines-1] == "In addition: Warning message:"
  totalLines <- c(coreLines-addedLines, coreLines)
  if(length(totalLines))
    return(lines[-totalLines]) else return(lines)
}

test_mcmc <- function(example, model, data = NULL, inits = NULL, ..., name = NULL, knownFailures = list(), expectWarnings = list(), avoidNestedTest = FALSE) {
  ## imitate processing test_mcmc_internal just to get a name for the test_that description
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
  name <- basename(name) ## name could be a pathed directory including tempdir(), which would change every time and hence appear as errors in line-by-line comparison with the gold file. So for futher purposes we use only the file name
  ## `missing(example)` does not work inside the test_that
  if(!missing(example)) {
    ## classic-bugs example specified by name
    dir = nimble:::getBUGSexampleDir(example)
    if(missing(model)) model <- example
    modelKnown <- TRUE
  } else {
    dir = ""
    modelKnown <- !missing(model)
  }
  
  if(avoidNestedTest) {  ## sometimes test_mcmc is called from within a test_that; this avoids report of empty test as of testthat v2.0.0
    expect_true(modelKnown, 'Neither BUGS example nor model code supplied.')
    Rmodel <- readBUGSmodel(model, data = data, inits = inits, dir = dir, useInits = TRUE,
                            check = FALSE)
    test_mcmc_internal(Rmodel, ..., name = name, knownFailures = knownFailures, expectWarnings = expectWarnings)
  } else {
    test_that(name, {
      expect_true(modelKnown, 'Neither BUGS example nor model code supplied.')
      Rmodel <- readBUGSmodel(model, data = data, inits = inits, dir = dir, useInits = TRUE,
                              check = FALSE)
      test_mcmc_internal(Rmodel, ..., name = name, knownFailures = knownFailures, expectWarnings = expectWarnings)
    })
  }
}


test_mcmc_internal <- function(Rmodel, ##data = NULL, inits = NULL,
                               verbose = nimbleOptions('verbose'), numItsR = 5, numItsC = 1000,
                               basic = TRUE, exactSample = NULL, results = NULL, resultsTolerance = NULL,
                               numItsC_results = numItsC,
                               resampleData = FALSE,
                               topLevelValues = NULL, seed = 0, mcmcControl = NULL, samplers = NULL, removeAllDefaultSamplers = FALSE,
                               doR = TRUE, doCpp = TRUE, returnSamples = FALSE, name = NULL, knownFailures = list(), expectWarnings = list()) {
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
  
  setSampler <- function(var, conf) {
    currentTargets <- sapply(conf$samplerConfs, function(x) x$target)
    # remove already defined scalar samplers
    inds <- which(unlist(var$target) %in% currentTargets)
    conf$removeSamplers(inds, print = FALSE)
    # look for cases where one is adding a blocked sampler specified on a variable and should remove scalar samplers for constituent nodes
    currentTargets <- sapply(conf$samplerConfs, function(x) x$target)
    inds <- which(sapply(unlist(var$target), function(x) Rmodel$expandNodeNames(x)) %in% currentTargets)
    conf$removeSamplers(inds, print = FALSE)
    
    if(is.list(var$target) && length(var$target) == 1) var$target <- var$target[[1]]
    if(length(var$target) == 1 || (var$type %in% c("RW_block", "RW_PF_block", "RW_llFunction_block") && !is.list(var$target)))
      tmp <- conf$addSampler(type = var$type, target = var$target, control = var$control, print = FALSE) else tmp <- sapply(var$target, function(x) conf$addSampler(type = var$type, target = x, control = var$control, print = FALSE))
  }
  
  wrap_if_matches('nameOK', names(knownFailures), expect_failure, {
    expect_false(is.null(name), info = 'name argument NULL')
  })
  
  ## leaving this message permanently on for now
  cat("===== Starting MCMC test for ", name, ". =====\n", sep = "") ## for log file, for comparison to gold file
  system(paste0("echo \"===== Starting MCMC test for ", name, ". =====\n\"", sep = "")) ## for travis log file, so it knows the process is not dead after 10 minutes of silence (message() does not work)
  
  if(doCpp) {
    Cmodel <- compileNimble(Rmodel)
    if(verbose) cat('done compiling model\n')
  }
  if(!is.null(mcmcControl)) mcmcConf <- configureMCMC(Rmodel, control = mcmcControl) else mcmcConf <- configureMCMC(Rmodel)
  if(removeAllDefaultSamplers) mcmcConf$removeSamplers()
  
  if(!is.null(samplers)) {
    sapply(samplers, setSampler, mcmcConf)
    cat("Setting samplers to:\n")
    print(mcmcConf$getSamplers())
  }
  
  vars <- Rmodel$getDependencies(Rmodel$getNodeNames(topOnly = TRUE, stochOnly = TRUE), stochOnly = TRUE, includeData = FALSE, downstream = TRUE)
  vars <- unique(nimble:::removeIndexing(vars))
  mcmcConf$addMonitors(vars, print = FALSE)
  
  Rmcmc <- buildMCMC(mcmcConf)
  if(doCpp) {
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  }
  
  if(basic) {
    ## do short runs and compare R and C MCMC output
    if(doR) {
      set.seed(seed)
      R_samples <- NULL
      ## need expect_error not expect_failure(expect_something()) because otherwise
      ## R error will stop execution
      wrap_if_matches('R MCMC', names(knownFailures), expect_error, {
        RmcmcOut <- Rmcmc$run(numItsR)
        RmvSample  <- nfVar(Rmcmc, 'mvSamples')
        R_samples <- as.matrix(RmvSample)
      })
    }
    if(doCpp) {
      set.seed(seed)
      Cmcmc$run(numItsC)
      CmvSample <- nfVar(Cmcmc, 'mvSamples')
      C_samples <- as.matrix(CmvSample)
      ## for some reason columns in different order in CmvSample...
      if(doR)
        C_subSamples <- C_samples[seq_len(numItsR), attributes(R_samples)$dimnames[[2]], drop = FALSE]
    }
    if(doR && doCpp && !is.null(R_samples)) {
      wrap_if_matches('R C samples match', names(knownFailures), expect_failure, {
        expect_equal(R_samples, C_subSamples, info = paste("R and C posterior samples are not equal"))
      })
    }
    
    if(doCpp) {
      if(!is.null(exactSample)) {
        for(varName in names(exactSample))
          wrap_if_matches('C samples match known samples', names(knownFailures), expect_failure, {
            expect_equal(round(C_samples[seq_along(exactSample[[varName]]), varName], 8),
                         round(exactSample[[varName]], 8),
                         info = paste0("Equality of compiled MCMC samples and known exact samples for variable ", varName))})
      }
    }
    
    summarize_posterior <- function(vals)
      return(c(mean = mean(vals), sd = sd(vals), quantile(vals, .025), quantile(vals, .975)))
    
    if(doCpp) {
      start <- round(numItsC / 2) + 1
      try(print(apply(C_samples[start:numItsC, , drop = FALSE], 2, summarize_posterior)))
    }
  }
  
  ## assume doR and doCpp from here down
  if(!is.null(results)) {
    ## do (potentially) longer run and compare results to inputs given
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
          samplesNames <- dimnames(C_samples)[[2]]
          if(!grepl("[", varName, fixed = TRUE))
            samplesNames <- gsub("\\[.*\\]", "", samplesNames)
          matched <- which(varName == samplesNames)
          diff <- abs(postResult[matched] - results[[metric]][[varName]])
          for(ind in seq_along(diff)) {
            strInfo <- ifelse(length(diff) > 1, paste0("[", ind, "]"), "")
            wrap_if_matches(paste('MCMC match to known posterior:', varName, metric, ind), names(knownFailures), expect_failure, {
              expect_true(diff[ind] < resultsTolerance[[metric]][[varName]][ind],
                          info = paste("Test of MCMC result against known posterior for :",  metric, "(", varName, strInfo, ")"))
            })
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
            wrap_if_matches(paste('MCMC match to known posterior:', varName, 'cov', ind), names(knownFailures), expect_failure, {
              expect_true(diff[ind] < resultsTolerance[[metric]][[varName]][ind],
                          info = paste("Test of MCMC result against known posterior for:",  metric, "(", varName, ")", strInfo))
            })
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
    nonDataNodesElements <- Rmodel$getDependencies(topNodes, self = TRUE, includeData = FALSE, downstream = TRUE, stochOnly = TRUE, returnScalarComponents = TRUE)
    dataVars <- unique(nimble:::removeIndexing(Rmodel$getDependencies(topNodes, dataOnly = TRUE, downstream = TRUE)))
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
    cat("Coverage for ", name, " is", coverage*100, "%.\n")
    miscoverage <- abs(coverage - 0.95)
    ## always print for purpose of goldfile
    # if(miscoverage > tolerance || verbose) {
    cat("True values with 95% posterior interval:\n")
    print(cbind(trueVals, t(interval), covered))
    # }
    wrap_if_matches('coverage', names(knownFailures), expect_failure, {
      expect_true(miscoverage < tolerance,
                  info = paste("Test of MCMC coverage on known parameter values for:", name))
    })
    
  }
  
  cat("===== Finished MCMC test for ", name, ". =====\n", sep = "")
  
  if(doCpp) {
    if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(Rmodel)
    }
  }
  return(returnVal)
}

test_resampler <- function(samplerName, wtsList, reps = 500, creps = reps){
  n <- sapply(wtsList, function(x){return(length(x))})
  output <- lapply(n, function(x){return(numeric(x))})
  avgCounts <- output
  samplerFunction <- getFromNamespace(samplerName, 'nimble')()
  for(rep in 1:reps){
    counts <- list()
    for(i in 1:length(wtsList)){
      output[[i]] <-    samplerFunction$run(wtsList[[i]])
      counts[[i]] <- numeric(length(output[[i]]))
      for(j in 1:n[i]){
        counts[[i]][j] <- length(which(output[[i]] == j))
      }
      avgCounts[[i]] <- avgCounts[[i]] + counts[[i]]
    }
  }
  expectedValue <- list(length(wtsList))
  for(i in 1:length(wtsList)){
    avgCounts[[i]] <- avgCounts[[i]]/reps
    expectedValue[[i]] <- n[i]*(wtsList[[i]]/sum(wtsList[[i]]))
    diffVec <- abs(expectedValue[[i]] - avgCounts[[i]])
    for(j in 1:n[i]){
      test_that(paste0("Test of accurate samples for uncompiled resampling
                      method ", samplerName, ", weight set ", i,
                       ", weight number ", j), 
                expect_lt(diffVec[j], sqrt(expectedValue[[i]][j]) + .01))
    }
  }
  avgCounts <- lapply(n, function(x){return(numeric(x))})
  compiledSamplerFunction <-  compileNimble(samplerFunction)
  for(rep in 1:reps){
    counts <- list()
    for(i in 1:length(wtsList)){
      output[[i]] <-    compiledSamplerFunction$run(wtsList[[i]])
      counts[[i]] <- numeric(length(output[[i]]))
      for(j in 1:n[i]){
        counts[[i]][j] <- length(which(output[[i]] == j))
      }
      avgCounts[[i]] <- avgCounts[[i]] + counts[[i]]
    }
  }
  expectedValue <- list(length(wtsList))
  for(i in 1:length(wtsList)){
    avgCounts[[i]] <- avgCounts[[i]]/reps
    expectedValue[[i]] <- n[i]*(wtsList[[i]]/sum(wtsList[[i]]))
    diffVec <- abs(expectedValue[[i]] - avgCounts[[i]])
    for(j in 1:n[i]){
      test_that(paste0("Test of accurate samples for compiled resampling
                       method ", samplerName, ", weight set ", i,
                       ", weight number ", j), 
                expect_lt(diffVec[j], sqrt(expectedValue[[i]][j]) + .01))
    }
  }
}


##### SET PARAMETERS #####
context("Testing that SMC functionality throws errors pointing to nimbleSMC.")

### particle filter testing follows similar steps to MCMC testing.
### for each example, we comapare filter output between R and C.  We also (where applicable) compare:
### 1) estimated values for latent states to true values for both weighted and equally weighted samples
### 2) estimated top level parameter values to known values (Liu-West filter, PMCMC)
### 3) estimated log-likelihood values to known values (for normal transition - observation
###    models where LL can be calculated analytically via KF)

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

oldWidth <- getOption("width")
options(width = 1000)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

### basic scalar latent node example, no top-level params

##### BEGIN TESTING #####


### scalar latent node and data node example, two top-level params

code <- nimbleCode({
  x[1] ~ dnorm(mean = mu0, sd = sigma_x);
  y[1] ~ dnorm(x[1], sd=sigma_y);
  for(i in 2:N){
    x[i] ~ dnorm(mean = a + b*x[i-1], sd = sigma_x);
    y[i] ~ dnorm(mean = x[i], sd=sigma_y);
  }
  sigma_x ~ T(dnorm(1, sd = .5), 0,);
  sigma_y ~ T(dnorm(.1, sd = .5), 0,);
  mu0 <- 0
})

set.seed(0)

a <- 0
b <- 1
N <- 5
sigma_x <- 1
sigma_y <- .1
x <- rep(NA, N)
y <- x
x[1] <- rnorm(1,0,sigma_x)
y[1] <- rnorm(1,x[1], sigma_y)
for(i in 2:N){
  x[i] <- rnorm(1,a+b*x[i-1], sigma_x)
  y[i] <- rnorm(1,x[i], sigma_y)
}

consts <- list(N=N)

testdata <- list(y=y)

inits <- list(sigma_x=sigma_x, sigma_y=sigma_y, x = x, a = a, b = b)
expect_error(
test_mcmc(model = code, name = 'scalar pmcmc', inits = inits, data = c(testdata, consts),  samplers = list(
  list(type = 'RW_PF', target = 'sigma_x', control = list(latents = 'x', m = 1000)),
  list(type = 'RW_PF', target = 'sigma_y', control = list(latents = 'x', m = 1000))),
  removeAllDefaultSamplers = TRUE, numItsC = 1000, results = list(
    mean = list( sigma_x = sigma_x,
                 sigma_y = sigma_y)),
  resultsTolerance = list(mean = list(sigma_x = .5,
                                      sigma_y = .5)))
)

expect_error(
test_mcmc(model = code, name = 'block pmcmc', inits = inits, data = c(testdata, consts),  samplers = list(
  list(type = 'RW_PF_block', target = c('sigma_x', 'sigma_y'), control = list(latents = 'x', m = 1000, resample = FALSE))),
  removeAllDefaultSamplers = TRUE, numItsC = 1000, results = list(
    mean = list(sigma_x = sigma_x,
                sigma_y = sigma_y)),
  resultsTolerance = list(mean = list(sigma_x = .5,
                                      sigma_y = .5)))
# , expectWarnings = list("C MCMC" = "NaNs produced"))
)
## Let's stop here to save testing time
## # test MCMC with longer runs and lower tolerance
## set.seed(0)
## N <- 50
## sigma_x <- 1
## sigma_y <- .1
## x <- rep(NA, N)
## y <- x
## x[1] <- rnorm(1,0,sigma_x)
## y[1] <- rnorm(1,x[1], sigma_y)
## for(i in 2:N){
##   x[i] <- rnorm(1,x[i-1], sigma_x)
##   y[i] <- rnorm(1,x[i], sigma_y)
## }

## consts <- list(N=N)

## testdata <- list(y=y)

## inits <- list(sigma_x=1, sigma_y=.1, x = x)

## test_mcmc(model = code, name = 'scalar pmcmc, more data', inits = inits, data = c(testdata, consts),  basic = FALSE, samplers = list(
##   list(type = 'RW_PF', target = 'sigma_x', control = list(latents = 'x', m = 1000, resample = FALSE)),
##   list(type = 'RW_PF', target = 'sigma_y', control = list(latents = 'x', m = 1000, resample = FALSE))),
##   removeAllDefaultSamplers = TRUE, numItsC = 1000, numItsC_results = 5000, results = list(
##   mean = list( sigma_x = sigma_x,
##               sigma_y = sigma_y)),
##   resultsTolerance = list(mean = list(sigma_x = .1,
##                                       sigma_y = .1)))

## test_mcmc(model = code, name = 'block pmcmc, more data', inits = inits, data = c(testdata, consts), basic = FALSE, samplers = list(
##   list(type = 'RW_PF_block', target = c('sigma_x', 'sigma_y'), control = list(latents = 'x', m = 1000, resample = FALSE))),
##   removeAllDefaultSamplers = TRUE, numItsC = 1000, numItsC_results = 5000, results = list(
##     mean = list(sigma_x = sigma_x,
##                 sigma_y = sigma_y)),
##   resultsTolerance = list(mean = list(sigma_x = .1,
##                                       sigma_y = .1)))




options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
options(width = oldWidth)