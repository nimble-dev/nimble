require(testthat)
require(methods)
require(nimble)

## These make it clear that error messages are expected.
if(FALSE) {
expect_failure <- function(...) {
    cat('BEGIN expected failure:\n', file = stderr())
    testthat:::expect_failure(...)
    argList <- list(...)
    cat(argList$info)
    cat('\nEND expected failure.\n', file = stderr())
}
}


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

expect_compiles <- function(..., info = NULL, link = FALSE, forceO1 = TRUE) {
    oldSCBL <- nimbleOptions('stopCompilationBeforeLinking')
    nimbleOptions(stopCompilationBeforeLinking = !link)
    oldForceO1 <- nimbleOptions('forceO1')
    nimbleOptions(forceO1 = forceO1)
    on.exit({
        assign('.check', 1, globalenv())
        nimbleOptions(stopCompilationBeforeLinking = oldSCBL)
        nimbleOptions(forceO1 = oldForceO1)
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

test_coreRfeature_batch <- function(input_batch, info = '', verbose = nimbleOptions('verbose'), dirName = NULL) {
        test_coreRfeature_batch_internal(input_batch, verbose = verbose, dirName)
}

test_coreRfeature_batch_internal <- function(input_batch, verbose = nimbleOptions('verbose'), dirName = NULL) { ## a lot like test_math but a bit more flexible
    names(input_batch) <- paste0('batch_case_', seq_along(input_batch))
    runFuns <- lapply(input_batch, gen_runFunCore)
    nfR <- nimbleFunction(setup = TRUE, methods = runFuns)()

    nfC <- compileNimble(nfR, dirName = dirName)

    for(i in seq_along(input_batch)) {
        input <- input_batch[[i]]
        if(verbose) cat("### Testing", input$name, "###\n")
        test_that(input$name, {
            funName <- names(input_batch)[i]
            if(is.null(input$expectWarnings))
                input$expectWarnings <- list()
            nArgs <- length(input$args)
            evalEnv <- new.env()
            eval(input$setArgVals, envir = evalEnv)
            savedArgs <- as.list(evalEnv)
            seedToUse <- if(is.null(input[['seed']])) 31415927 else input[['seed']]
            set.seed(seedToUse)
            wrap_if_matches("R eval", names(input$expectWarnings), expect_warning, {
                eval(input$expr, envir = evalEnv)
            })
            savedOutputs <- as.list(evalEnv)
            list2env(savedArgs, envir = evalEnv)
            ## with R ref classes, lookup of methods via `[[` does not work until it has been done via `$`
            ## force it via `$` here to allow simpler syntax below
            forceFind <- eval(substitute(nfR$FUNNAME, list(FUNNAME = as.name(funName))))
            forceFind <- eval(substitute(nfC$FUNNAME, list(FUNNAME = as.name(funName))))
            if(nArgs == 5) {
                set.seed(seedToUse)
                wrap_if_matches("R run", names(input$expectWarnings), expect_warning, {
                    out_nfR = nfR[[funName]](evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4, evalEnv$arg5)})
                list2env(savedArgs, envir = evalEnv)
                set.seed(seedToUse)
                out_nfC = nfC[[funName]](evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4, evalEnv$arg5)
            }  
            if(nArgs == 4) {
                set.seed(seedToUse)
                wrap_if_matches("R run", names(input$expectWarnings), expect_warning, {
                    out_nfR = nfR[[funName]](evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4)
                })
                list2env(savedArgs, envir = evalEnv)
                set.seed(seedToUse)
                out_nfC = nfC[[funName]](evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4)
            }  
            if(nArgs == 3) {
                set.seed(seedToUse)
                wrap_if_matches("R run", names(input$expectWarnings), expect_warning, {
                    out_nfR = nfR[[funName]](evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
                })
                list2env(savedArgs, envir = evalEnv)
                set.seed(seedToUse)
                out_nfC = nfC[[funName]](evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
            }  
            if(nArgs == 2) {
                set.seed(seedToUse)
                wrap_if_matches("R run", names(input$expectWarnings), expect_warning, {
                    out_nfR = nfR[[funName]](evalEnv$arg1, evalEnv$arg2)
                })
                list2env(savedArgs, envir = evalEnv)
                set.seed(seedToUse)
                out_nfC = nfC[[funName]](evalEnv$arg1, evalEnv$arg2)
            }
            if(nArgs == 1) {
                set.seed(seedToUse)
                wrap_if_matches("R run", names(input$expectWarnings), expect_warning, {
                    out_nfR = nfR[[funName]](evalEnv$arg1)
                })
                list2env(savedArgs, envir = evalEnv)
                set.seed(seedToUse)
                out_nfC = nfC[[funName]](evalEnv$arg1)
            }
            if(nArgs == 0) {
                set.seed(seedToUse)
                wrap_if_matches("R run", names(input$expectWarnings), expect_warning, {
                    out_nfR = nfR[[funName]]()
                })
                list2env(savedArgs, envir = evalEnv)
                set.seed(seedToUse)
                out_nfC = nfC[[funName]]()
            }
            out <- savedOutputs$out
            ## clear any attributes except dim
            dimOut <- attr(out, 'dim')
            dimOutR <- attr(out_nfR, 'dim')
            dimOutC <- attr(out_nfC, 'dim')
            attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
            if(!is.null(input[['storage.mode']]))
                storage.mode(out) <- storage.mode(out_nfR) <- storage.mode(out_nfC) <- input[['storage.mode']]
            attr(out, 'dim') <- dimOut
            attr(out_nfR, 'dim') <- dimOutR
            attr(out_nfC, 'dim') <- dimOutC
            checkEqual <- input[['checkEqual']]
            if(is.null(checkEqual)) checkEqual <- FALSE
            if(is.null(input[['return']])) { ## use default 'out' object
                if(!checkEqual) {
                    expect_identical(class(out), class(out_nfC), info = paste('iden tmp test of class', class(out), class(out_nfC)))
                    expect_identical(dim(out), dim(out_nfC), info = 'iden test of dim')
                    expect_identical(round(out, 10), round(out_nfC, 10), info='iden test of round')
                    expect_identical(out, out_nfR, info = "Identical test of coreRfeature (direct R vs. R nimbleFunction)")
                    wh <- which.max(abs(out - out_nfC))
                    expect_identical(out, out_nfC, info = paste("Identical test of coreRfeature (direct R vs. C++ nimbleFunction)", system('gcc --version', intern = T)[1], ' ', sprintf("%0.20f", out[wh]), " ", sprintf("%0.20f", out_nfC[wh])))
                } else {
                    expect_equal(out, out_nfR, info = "Equal test of coreRfeature (direct R vs. R nimbleFunction)")
                    expect_equal(out, out_nfC, info = "Equal test of coreRfeature (direct R vs. C++ nimbleFunction)")
                }
            } else { ## not using default return(out), so only compare out_nfR to out_nfC
                if(!checkEqual) {
                    expect_identical(out_nfC, out_nfR, info = "Identical test of coreRfeature (compiled vs. uncompied nimbleFunction)")
                } else {
                    expect_identical(out_nfC, out_nfR, info = "Equal test of coreRfeature (compiled vs. uncompied nimbleFunction)")
                }
            }
        })
    }
    ## unload DLL as R doesn't like to have too many loaded
    if(.Platform$OS.type != 'windows') nimble:::clearCompiled(nfR) ##dyn.unload(project$cppProjects[[1]]$getSOName())
    invisible(NULL)
}
        
test_coreRfeature <- function(input, verbose = nimbleOptions('verbose'), dirName = NULL) {
    test_that(input$name, {
        test_coreRfeature_internal(input, verbose, dirName)
    })
}

test_coreRfeature_internal <- function(input, verbose = nimbleOptions('verbose'), dirName = NULL) { ## a lot like test_math but a bit more flexible
  if(verbose) cat("### Testing", input$name, "###\n")
  runFun <- gen_runFunCore(input)
  nfR <- nimbleFunction(run = runFun)
  ## This try is safe because failure is caught by expect_equal below

  expectedCompilerFail <- FALSE
  if(!is.null(input[['expectedCompilerError']]))
      if(isTRUE(input[['expectedCompilerError']]))
          expectedCompilerError <- TRUE
  if(!expectedCompilerError) {
      nfC <- compileNimble(nfR, dirName = dirName)
  } else {
      expect_error(nfC <- compileNimble(nfR, dirName = dirName))
      return()
  }
  
  nArgs <- length(input$args)
  evalEnv <- new.env()
  eval(input$setArgVals, envir = evalEnv)
  savedArgs <- as.list(evalEnv)
  seedToUse <- if(is.null(input[['seed']])) 31415927 else input[['seed']]
  set.seed(seedToUse)
  eval(input$expr, envir = evalEnv)
  savedOutputs <- as.list(evalEnv)
  list2env(savedArgs, envir = evalEnv)
  if(nArgs == 5) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4, evalEnv$arg5)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4, evalEnv$arg5)
  }  
  if(nArgs == 4) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3, evalEnv$arg4)
  }  

  if(nArgs == 3) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2, evalEnv$arg3)
  }  
  if(nArgs == 2) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1, evalEnv$arg2)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1, evalEnv$arg2)
  }
  if(nArgs == 1) {
      set.seed(seedToUse)
      out_nfR = nfR(evalEnv$arg1)
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC(evalEnv$arg1)
  }
  if(nArgs == 0) {
      set.seed(seedToUse)
      out_nfR = nfR()
      list2env(savedArgs, envir = evalEnv)
      set.seed(seedToUse)
      out_nfC = nfC()
  }
  out <- savedOutputs$out
  ## clearn any attributes except dim
  dimOut <- attr(out, 'dim')
  dimOutR <- attr(out_nfR, 'dim')
  dimOutC <- attr(out_nfC, 'dim')
  attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
  attr(out, 'dim') <- dimOut
  attr(out_nfR, 'dim') <- dimOutR
  attr(out_nfC, 'dim') <- dimOutC
  checkEqual <- input[['checkEqual']]
  if(is.null(checkEqual)) checkEqual <- FALSE
  if(is.null(input[['return']])) { ## use default 'out' object
      if(!checkEqual) {
          expect_identical(out, out_nfR, info = paste0("FOO Identical test of coreRfeature (direct R vs. R nimbleFunction): ", input$name))
          expect_identical(out, out_nfC, info = paste0("FOO Identical test of coreRfeature (direct R vs. C++ nimbleFunction): ", input$name))
      } else {
          expect_equal(out, out_nfR, info = paste0("Equal test of coreRfeature (direct R vs. R nimbleFunction): ", input$name) )
          expect_equal(out, out_nfC, info = paste0("Equal test of coreRfeature (direct R vs. C++ nimbleFunction): ", input$name))
      }
  } else { ## not using default return(out), so only compare out_nfR to out_nfC
      if(!checkEqual) {
          expect_identical(out_nfC, out_nfR, info = paste0("Identical test of coreRfeature (compiled vs. uncompied nimbleFunction): ", input$name))
      } else {
          expect_identical(out_nfC, out_nfR, info = paste0("Equal test of coreRfeature (compiled vs. uncompied nimbleFunction): ", input$name))
      }
  }
  # unload DLL as R doesn't like to have too many loaded
  if(.Platform$OS.type != 'windows') nimble:::clearCompiled(nfR) ##dyn.unload(project$cppProjects[[1]]$getSOName())
  invisible(NULL)
}

gen_runFun <- function(param, logicalArgs, returnType = "double") {
    runFun <- function() {}
    types <- rep('double', length(param$inputDim))
    types[logicalArgs] <- 'logical'
    types <- paste0(types, '(', param$inputDim, ')')
    formalsList <- lapply(types, function(x) parse(text = x)[[1]])
    names(formalsList) <- paste0('arg', seq_along(param$inputDim))
    formals(runFun) <- formalsList
    tmp <- quote({})
    tmp[[2]] <- param$expr
    tmp[[3]] <- quote(return(out))
    tmp[[4]] <- parse(text = paste0("returnType(", returnType, "(", param$outputDim, "))"))[[1]]
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

wrap_if_matches <- function(pattern, string, wrapper, expr) {
    if (!is.null(pattern) && any(grepl(paste0('^', pattern, '$'), string))) {
        wrapper(expr)
    } else {
        expr
    }
}

wrap_if_true <- function(test, wrapper, expr, wrap_in_try = FALSE) {
  wrap <- if (isTRUE(wrap_in_try))
            function(x) try(x, silent = TRUE)
          else identity
  if (isTRUE(test)) wrap(wrapper(expr)) else wrap(expr)
}

## This is a parametrized test, where `param` is a list with names:
##   param$name - A descriptive test name.
##   param$expr - A quoted expression `quote(out <- some_function_of(arg1, arg2, ...))`.
##   param$Rcode - Optional R version of expr.
##   param$inputDim - A vector of dimensions of the input `arg`s.
##   param$outputDim - The dimension of the output `out.
##   param$xfail - Optional regular expression of tests that are expected to fail.
test_math <- function(param, caseName, verbose = nimbleOptions('verbose'), size = 3, dirName = NULL) {
    info <- paste0(caseName, ': ', param$name)
    ## in some cases, expect_error does not suppress error messages (I believe this has
    ## to do with how we trap errors in compilation), so make sure user realizes expectation
    if('knownFailureReport' %in% names(param) && param$knownFailureReport)
        cat("\nBegin expected error message:\n")
    test_that(info, {
        ## wrap_if_matches(param$xfail, paste0(info, ': compiles and runs'), expect_error, {
        test_math_internal(param, info, verbose, size, dirName)
        ## })
    })
    if('knownFailureReport' %in% names(param) && param$knownFailureReport)
        cat("End expected error message.\n")
    invisible(NULL)
}

test_math_internal <- function(param, info, verbose = nimbleOptions('verbose'), size = 3, dirName = NULL) {
  if(verbose) cat("### Testing", param$name, "###\n")
  nArgs <- length(param$inputDim)
  logicalArgs <- rep(FALSE, nArgs)
  if("logicalArgs" %in% names(param))
      logicalArgs <- param$logicalArgs
  returnType <- "double"
  if("returnType" %in% names(param))
      returnType <- param$returnType

  runFun <- gen_runFun(param, logicalArgs, returnType)
  wrap_if_matches(param$expectWarnings, "builds", expect_warning, {
      nfR <- nimbleFunction(  
          run = runFun)
  })
  
  info <- paste0(info, ": compiles")
  ## need expect_error not expect_failure(expect_something()) because otherwise
  ## R error will stop execution
  wrap_if_matches(param$knownFailure, info, expect_error, {
      nfC <- compileNimble(nfR, dirName = dirName)

      arg1 <- make_input(param$inputDim[1], size = size, logicalArgs[1])
      if(nArgs > 1)
          arg2 <- make_input(param$inputDim[2], size = size, logicalArgs[2])
      if(nArgs > 2)
          arg3 <- make_input(param$inputDim[3], size = size, logicalArgs[3])
      if(nArgs > 3)
          stop("test_math not set up for >3 args yet")
      
      if("Rcode" %in% names(param)) {      
          eval(param$Rcode)
      } else {
          eval(param$expr)
      }
      info <- paste0(info, ": runs")
      wrap_if_matches(param$knownFailure, info, expect_failure, {
          if(nArgs == 3) {
              expect_silent(out_nfR <- nfR(arg1, arg2, arg3))
              expect_silent(out_nfC <- nfC(arg1, arg2, arg3))
          }  
          if(nArgs == 2) {
              expect_silent(out_nfR <- nfR(arg1, arg2))
              expect_silent(out_nfC <- nfC(arg1, arg2))
          }
          if(nArgs == 1) {
              expect_silent(out_nfR <- nfR(arg1))
              expect_silent(out_nfC <- nfC(arg1))
          }
      
          attributes(out) <- attributes(out_nfR) <- attributes(out_nfC) <- NULL
          
          infoR <- paste0(info, ": R vs Nimble DSL")
          wrap_if_matches(param$knownFailure, infoR, expect_failure, {
              expect_equal(out, out_nfR, info = infoR)
          })
          infoC <- paste0(info, ": R vs Nimble Cpp")
          wrap_if_matches(param$knownFailure, infoC, expect_failure, {
              expect_equal(out, out_nfC, info = infoC)
          })
      })
      # Unload DLL as R doesn't like to have too many loaded.
      if(.Platform$OS.type != 'windows') nimble:::clearCompiled(nfR)
  })
  invisible(NULL)
}


### Function for testing MCMC called from test_mcmc.R

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


test_filter <- function(example, model, data = list(), inits = list(),
                        verbose = nimbleOptions('verbose'), numItsR = 3, numItsC = 10000,
                        basic = TRUE, exactSample = NULL, results = NULL, resultsTolerance = NULL,
                        numItsC_results = numItsC,
                        seed = 0, filterType = NULL, latentNodes = NULL, filterControl = NULL,
                        doubleCompare = FALSE, filterType2 = NULL,
                        doR = TRUE, doCpp = TRUE, returnSamples = FALSE, name = NULL, dirName = NULL, 
                        knownFailures = list()) {
    ## There are two modes of testing:
    ## 1) basic = TRUE: compares R and C Particle Filter likelihoods and sampled states
    ## 2) if you pass 'results', it will compare Filter output to known latent state posterior summaries, top-level parameter posterior summaries,
    ##    and likelihoods within tolerance specified in resultsTolerance.  Results are compared for both weighted and unweighted samples.
    ## filterType determines which filter to use for the model.  Valid options are: "bootstrap", "auxiliary", "LiuWest", "ensembleKF"
    ## filterControl specifies options to filter function, such as saveAll = TRUE/FALSE.

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

    ## keep this outside of test_that as use of missing within test_that triggers error with "'missing' can
    ## only be used for arguments"
    if(!missing(example)) {
        ## classic-bugs example specified by name
        dir <- getBUGSexampleDir(example)
        if(missing(model)) model <- example
    } else {
        ## code, data and inits specified directly where 'model' contains the code
        example = deparse(substitute(model))
        if(missing(model)) stop("Neither BUGS example nor model code supplied.")
        dir <- ""
    }
    returnVal <- NULL
    
    cat("===== Starting Filter test for ", name, " using ", filterType, ". =====\n", sep = "")

    test_that(name, {
        Rmodel <- readBUGSmodel(model, dir = dir, data = data, inits = inits, useInits = TRUE, check = FALSE)
        if(doCpp) {
            Cmodel <- compileNimble(Rmodel, dirName = dirName)
            if(verbose) cat('done compiling model\n')
        }
        if(verbose) cat("Building filter\n")
        if(filterType == "bootstrap"){
            if(!is.null(filterControl))  Rfilter <- buildBootstrapFilter(Rmodel, nodes = latentNodes, control = filterControl)
            else Rfilter <- buildBootstrapFilter(Rmodel, nodes = latentNodes, control = list(saveAll = TRUE, thresh = 0))
        }
        if(filterType == "auxiliary"){
            if(!is.null(filterControl))  Rfilter <- buildAuxiliaryFilter(Rmodel, nodes = latentNodes, control = filterControl)
            else Rfilter <- buildAuxiliaryFilter(Rmodel, nodes = latentNodes, control = list(saveAll = TRUE))
        }
        if(filterType == "LiuWest"){
            if(!is.null(filterControl))  Rfilter <- buildLiuWestFilter(Rmodel, nodes = latentNodes, control = filterControl)
            else Rfilter <- buildLiuWestFilter(Rmodel, nodes = latentNodes, control = list(saveAll = TRUE))
        }
        if(filterType == "ensembleKF"){
            if(!is.null(filterControl))  Rfilter <- buildEnsembleKF(Rmodel, nodes = latentNodes, control = filterControl)
            else Rfilter <- buildEnsembleKF(Rmodel, nodes = latentNodes, control = list(saveAll = TRUE))
        }
        saveAll <- TRUE 
        if(!is.null(filterControl) && exists('saveAll', filterControl))
            saveAll <- filterControl$saveAll
        
        if(doCpp) {
            Cfilter <- compileNimble(Rfilter, project = Rmodel, dirName = dirName)
        }

        if(basic) {
            ## do short runs and compare R and C filter output
            if(doR) {
                set.seed(seed)
                RfilterOut <- Rfilter$run(numItsR)
                if(filterType == "ensembleKF"){
                    RmvSample  <- nfVar(Rfilter, 'mvSamples')
                    R_samples <- as.matrix(RmvSample)
                }
                else{
                    RmvSample  <- nfVar(Rfilter, 'mvWSamples')
                    RmvSample2 <- nfVar(Rfilter, 'mvEWSamples')
                    R_samples <- as.matrix(RmvSample)
                    R_samples2 <- as.matrix(RmvSample2)
                    if(filterType != 'LiuWest'){
                        R_ESS <- Rfilter$returnESS()
                    }
                }
            } 
            if(doCpp) {
                set.seed(seed)
                CfilterOut <- Cfilter$run(numItsR)
                if(filterType == "ensembleKF"){
                    CmvSample <- nfVar(Cfilter, 'mvSamples')
                    C_samples <- as.matrix(CmvSample)
                    C_subSamples <- C_samples[, attributes(R_samples)$dimnames[[2]], drop = FALSE]
                }
                else{
                    CmvSample <- nfVar(Cfilter, 'mvWSamples')
                    CmvSample2 <- nfVar(Cfilter, 'mvEWSamples')
                    C_samples <- as.matrix(CmvSample)
                    C_samples2 <- as.matrix(CmvSample2)
                    C_subSamples <- C_samples[, attributes(R_samples)$dimnames[[2]], drop = FALSE]
                    C_subSamples2 <- C_samples2[, attributes(R_samples2)$dimnames[[2]], drop = FALSE]
                    if(filterType != 'LiuWest'){
                        C_ESS <- Rfilter$returnESS()
                        for(i in seq_along(length(C_ESS))){
                            wrap_if_matches('C ESS >= 0', names(knownFailures), expect_failure, {
                                expect_gte(C_ESS[i], 0)
                            })
                            wrap_if_matches('C ESS <= numIts', names(knownFailures), expect_failure, {
                                expect_lte(C_ESS[i], numItsR)
                            })
                        }
                    }
                }
                ## for some reason columns in different order in CmvSample...
            }
            if(doR && doCpp && !is.null(R_samples)) {
                ## context(paste0("testing ", example," ", filterType, " filter"))
                if(filterType == "ensembleKF"){
                    expect_equal(R_samples, C_subSamples, info = paste("R and C posterior samples are not equal"))
                }
                else{
                    expect_equal(R_samples, C_subSamples, info = paste("R and C weighted posterior samples are not equal"))
                    expect_equal(R_samples2, C_subSamples2, info = paste("R and C equally weighted posterior samples are not equal"))
                    expect_equal(RfilterOut, CfilterOut, info = paste("R and C log likelihood estimates are not equal"))
                    if(filterType != 'LiuWest'){
                        wrap_if_matches('R C ESS match', names(knownFailures), expect_failure, {
                            expect_equal(R_ESS, C_ESS, info = paste("R and C ESS are not equal"))
                        })
                    }
                }
            }

            if(doCpp) {
                if(!is.null(exactSample)) {
                    for(varName in names(exactSample))
                        expect_equal(round(C_samples[seq_along(exactSample[[varName]]), varName], 8), round(exactSample[[varName]], 8), info = paste0("filter result does not match known samples for: ", varName))
                }
            }

            summarize_posterior <- function(vals)
                return(c(mean = mean(vals), sd = sd(vals), quantile(vals, .025), quantile(vals, .975)))

            if(doCpp) {
                ## if(verbose) {
                ##   try(print(apply(C_samples[, , drop = FALSE], 2, summarize_posterior)))
                ## }
            }
        }

        ## assume doR and doCpp from here down
        if(!is.null(results)) {
            ## do (potentially) longer run and compare results to inputs given
            set.seed(seed)
            Cll <- Cfilter$run(numItsC_results)
            for(wMetric in c(TRUE, FALSE)){
                weightedOutput <- 'unweighted'
                if(filterType == "ensembleKF")
                    CfilterSample <- nfVar(Cfilter, 'mvSamples')
                else{
                    if(wMetric){
                        CfilterSample <- nfVar(Cfilter, 'mvWSamples')
                        weightedOutput <-"weighted"
                    }
                    else
                        CfilterSample <- nfVar(Cfilter, 'mvEWSamples')
                }

                C_samples <- as.matrix(CfilterSample)[, , drop = FALSE]
                if(weightedOutput == "weighted"){
                    wtIndices <- grep("^wts\\[", dimnames(C_samples)[[2]])
                    C_weights <- as.matrix(C_samples[,wtIndices, drop = FALSE])
                    C_samples <- as.matrix(C_samples[,-wtIndices, drop = FALSE])
                }
                latentNames <- Rmodel$expandNodeNames(latentNodes, sort = TRUE, returnScalarComponents = TRUE)
                if(weightedOutput == "weighted"){
                    samplesToWeightsMatch <- rep(dim(C_weights)[2], dim(C_samples)[2])
                    latentIndices <- match(latentNames, dimnames(C_samples)[[2]])
                    latentSampLength <- length(latentNames)
                    if(!saveAll) {  ## added without careful checking; may not be robust
                        latentIndices <- latentIndices[!is.na(latentIndices)]
                        latentSampLength <- 1
                    }
                    latentDim <- latentSampLength/dim(C_weights)[2]
                    samplesToWeightsMatch[latentIndices] <- rep(1:dim(C_weights)[2], each = latentDim )
                }
                for(metric in names(results)) {
                    if(!metric %in% c('mean', 'median', 'sd', 'var', 'cov', 'll'))
                        stop("Results input should be named list with the names indicating the summary metrics to be assessed, from amongst 'mean', 'median', 'sd', 'var', 'cov', and 'll'.")
                    if(!(metric %in% c('cov', 'll'))) {
                        if(weightedOutput == "weighted"){
                            postResult <- sapply(1:dim(C_samples)[2], weightedMetricFunc, metric = metric, weights = C_weights, samples = C_samples, samplesToWeightsMatch)
                        }
                        else
                            postResult <- apply(C_samples, 2, metric)
                        for(varName in names(results[[metric]])) {
                            samplesNames <- dimnames(C_samples)[[2]]
                            if(!grepl(varName, "[", fixed = TRUE))
                                samplesNames <- gsub("\\[.*\\]", "", samplesNames)
                            matched <- which(varName == samplesNames)
                            if(!saveAll) {  ## added without careful checking; may not be robust
                                diff <- abs(postResult[matched] - results[[metric]][[varName]][length(results[[metric]][[varName]])])
                            } else {
                                diff <- abs(postResult[matched] - results[[metric]][[varName]])
                            }
                            for(ind in seq_along(diff)) {
                                strInfo <- ifelse(length(diff) > 1, paste0("[", ind, "]"), "")
                                expect_lt(diff[ind], resultsTolerance[[metric]][[varName]][ind],
                                          label = paste0("filter posterior result against known posterior for:", weightedOutput,  metric, "(", varName, strInfo, ")"))
                            }
                        }
                    } else if (metric == 'cov' ) {
                        for(varName in names(results[[metric]])) {
                            matched <- grep(varName, dimnames(C_samples)[[2]], fixed = TRUE)
                            ##             if(weightedOutput == "weighted")
                            ##               postResult <- cov.wt(C_samples[, matched], wt = )  #weighted covariance not currently implemented
                            ##             else
                            postResult <- cov(C_samples[ , matched])
                            ## next bit is on vectorized form of matrix so a bit awkward
                            diff <- c(abs(postResult - results[[metric]][[varName]]))
                            for(ind in seq_along(diff)) {
                                strInfo <- ifelse(length(diff) > 1, paste0("[", ind, "]"), "")
                                expect_lt(diff[ind], resultsTolerance[[metric]][[varName]][ind],
                                          label = paste0("filter posterior result against known posterior for", example, ":",  metric, "(", varName, ")", strInfo))
                            }
                        }
                    }
                    else {  # ll (log likelihood)
                        diff <- abs(Cll - results[[metric]][[1]][1])
                        expect_lt(diff, resultsTolerance[[metric]][[1]][1], label = paste0("filter log-likelihood result against known log-likelihood for", example, ":",  metric))
                    }
                }
            }
        }
        try(print(apply(as.matrix(C_samples), 2, summarize_posterior)))  ## print summaries of equally weighted samples
        if(returnSamples) {
            if(exists('CmvSample'))
                returnVal <- as.matrix(CmvSample)
        } 
        if(doCpp) {
            if(.Platform$OS.type != 'windows') 
                nimble:::clearCompiled(Rmodel)
        }

    })
    cat("===== Finished ", filterType, " filter test for ", name, ". =====\n", sep = "")
    
    return(returnVal)
}

## Testing for correct behavior of different resampling methods
## used within PFs.
##   samplerName - A string with the name of 
##                 the resampling function to be tested.
##   wtsList - A list, where each element is a vector of weights to use for 
##             testing, given as input to the resampler functions.
##   reps - An integer, the number of repetitions to conduct.
##
## For each provided set of weights in wtsList, the test will produce
## 'reps' number of samples according to those weights.  For each set of
## samples, the number of times each element was sampled is recorded.  These
## recorded counts are averaged over the 'reps' number of samples, and then
## the averages are compared to the expected count of each element 
## (i.e. wts*length(wts)).

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

test_size <- function(input, verbose = nimbleOptions('verbose')) {
    if(is.null(input$expectPassWithConst)) input$expectPassWithConst <- input$expectPass
    if(is.null(input$knownProblem)) input$knownProblem <- FALSE
    if(is.null(input$knownProblemWithConst)) input$knownProblemWithConst <- input$knownProblem
    if(is.null(input$expectWarn)) input$expectWarn <- FALSE
    if(is.null(input$expectWarnWithConst)) input$expectWarnWithConst <- input$expectWarn

    if(verbose) cat("### Testing", input$name, " with RHS variable ###\n")
    code <- quote({
        m <- nimbleModel(code = input$expr, data = input$data, inits = input$inits)
        calculate(m)  ## Calculates from scratch.
        calculate(m)  ## Uses cached value.
    })
    message = paste(input$name, 'with RHS variable', ifelse(input$expectPass, 'works', 'fails'), 'as expected')
    if (input$knownProblem) message = paste(message, 'marked as KNOWN ISSUE')
    if(xor(input$expectPass, input$knownProblem)) {
        if(input$expectWarn) {
            test_that(message, expect_warning(eval(code)))
        } else test_that(message, expect_silent(eval(code)))
    } else {
        test_that(message, expect_error(suppressWarnings(eval(code))))
    }

    if(verbose) cat("### Testing", input$name, "with RHS constant ###\n")
    code <- quote({
        m <- nimbleModel(code = input$expr, data = input$data, constants = input$inits)
        calculate(m)  ## Calculates from scratch.
        calculate(m)  ## Uses cached value.
    })
    message = paste(input$name, 'with RHS constant', ifelse(input$expectPassWithConst, 'works', 'fails'), 'as expected')
    if (input$knownProblemWithConst) message = paste(message, 'marked as KNOWN ISSUE')
    if(xor(input$expectPassWithConst, input$knownProblemWithConst)) {
        if(input$expectWarnWithConst) {
            test_that(message, expect_warning(eval(code)))
        } else test_that(message, expect_silent(eval(code)))
    } else {
        ## As of testthat 3.0, warnings bubble up so need to deal with them by suppression.
        test_that(message, expect_error(suppressWarnings(eval(code))))
    }

    invisible(NULL)
}

# could redo test_size to always expect specific error, but not taking time to do that now
test_size_specific_error <- function(input, verbose = nimbleOptions('verbose')) {
    if(verbose) cat("### Testing", input$name, "###\n")
    test_that(paste0("Test 1 of size/dimension check: ", input$name), {
        expect_error(nimbleModel(code = input$expr, data = input$data, inits = input$inits),
                     regexp = input$correctErrorMsg, info = paste("Result does not match", input$expectPass))
    })

    invisible(NULL)
}


test_getParam <- function(distCall, dist = NULL) {
    distCallText <- deparse(distCall)
    test_that(distCallText, {
        gpScalar <- nimbleFunction(
            setup = function(model, node, param) {},
            run = function() {
                ans1 <- model$getParam(node, param)
                ans2 <- getParam(model, node, param) ## to become model$getParam(node, param)
                if(ans1 != ans2) stop('oops, ans1 != ans2')
                return(ans1)
                returnType(double())
            })
        if(is.null(dist)) dist <- nimble:::getDistributionInfo(as.character(distCall[[1]])) else dist <- nimble:::getDistributionInfo(dist)
        code <- substitute({x ~ DISTCALL}, list(DISTCALL = distCall))
        m <- nimbleModel( code = code )
        cm <- compileNimble(m)
        gpFuns <- list()
        expectedResults <- list()
        altParams <- dist$altParams
        altParamNames <- names(altParams)
        distCallText <- deparse(distCall)
        
        reqdArgs <- dist$reqdArgs ## these are canonical
        exprs <- dist$exprs
        alts <- dist$alts
        providedArgs <- names(distCall)
        providedArgs <- providedArgs[providedArgs != ""]
        whichExpr <- NULL
        ## figure out which way arguments were provided in distCall
        for(i in seq_along(exprs)) {
            if(all(providedArgs %in% alts[[i]])) whichExpr <- i
        }
        if(is.null(whichExpr)) {
            if(all(providedArgs %in% reqdArgs)) whichExpr <- 0
        }
        expect_equal(is.null(whichExpr), FALSE, 'args not found')

        
        ## exprs give expressions for calculating reqdArgs from alts

        ## altParams give expressions for calculating individual alts from reqdArgs
        
        ## put reqd in evalEnv, which means using exprs for the alts as needed
        ## if testing on something provided, grab what was provided.
        ## if testing on something not provided, if it is reqd then use it directly
        ## otherwise calculate it from altParams

        evalEnv <- new.env()

        for(i in seq_along(distCall)) {
            if(names(distCall)[i] != "") assign(names(distCall)[i], distCall[[i]], envir = evalEnv)
        }
        if(whichExpr > 0) {  ## what was provided was not canonical
            for(i in seq_along(exprs[[whichExpr]])) {
                assign(names(exprs[[whichExpr]])[i], eval(exprs[[whichExpr]][[i]], envir = evalEnv), envir = evalEnv)
            }
        }

        ## check recovery of alternative param names from what was provided
        for(i in seq_along(altParamNames)) {
            gpFuns[[i]] <- gpScalar(m, 'x', altParamNames[i])
            if(altParamNames[i] %in% providedArgs) ## it was provided so simply eval the name
                expectedResults[[i]] <- eval(as.name(altParamNames[i]), envir = evalEnv)
            else  ## it wasn't provided so eval the expression to calculate it from reqdArgs
                expectedResults[[i]] <- eval(altParams[[i]], envir = evalEnv)
            expect_equal(gpFuns[[i]]$run(), expectedResults[[i]], info = paste('error in uncompiled use',
                                                                        altParamNames[i]))
            expect_equal(m$getParam('x', altParamNames[i]), expectedResults[[i]],
                         info = paste('error in R getParam', altParamNames[i]))
            expect_equal(cm$getParam('x', altParamNames[i]), expectedResults[[i]],
                         info = paste('error in C getParam', altParamNames[i]))
        }

        resultsNames <- altParamNames
        nextI <- length(expectedResults)+1
        for(i in seq_along(reqdArgs)) {
            gpFuns[[nextI]] <- gpScalar(m, 'x', reqdArgs[i])
            expectedResults[[nextI]] <- eval(as.name(reqdArgs[i]), envir = evalEnv) ## it was already calculated into evalEnv above
            expect_equal(gpFuns[[nextI]]$run(), expectedResults[[nextI]],
                         info = paste('error in uncompiled reqd', reqdArgs[i]))
            expect_equal(m$getParam('x', reqdArgs[i]), expectedResults[[nextI]],
                         info = paste('error in R model reqd', reqdArgs[i]))
            expect_equal(cm$getParam('x', reqdArgs[i]), expectedResults[[nextI]],
                         info = paste('error in C model reqd', reqdArgs[i]))
            resultsNames[nextI] <- reqdArgs[i]
            nextI <- nextI + 1
        }
        
        compiled <- do.call('compileNimble', c(list(m), gpFuns, list(resetFunctions = TRUE)))
        for(i in seq_along(expectedResults)) {
            expect_equal(compiled[[i+1]]$run(), expectedResults[[i]],
                         info = paste('error in compiled', resultsNames[i]))
        }
        if(.Platform$OS.type != 'windows') nimble:::clearCompiled(m)
    })
    invisible(NULL)
}


test_getBound <- function(model, cmodel, test, node, bnd, truth, info) {
    test_that(paste0("getBound test: ", info), {
        rtest <- test(model, node, bnd)
        project <- nimble:::nimbleProjectClass(NULL, name = 'foo')
        ctest <- compileNimble(rtest, project = project)
        
        out1 <- model$getBound(node, bnd)
        out2 <- getBound(model, node, bnd)
        out3 <- cmodel$getBound(node, bnd)
        out4 <- getBound(cmodel, node, bnd)
        nfOutput <- rtest$run()
        cnfOutput <- ctest$run()
        
        expect_equal(truth, out1, info = paste0("mismatch of true bound with getBound result: ", info))
        expect_equal(out1, out2, info = paste0("function vs. method getBound call mismatch for uncompiled model with: ", info))
        expect_equal(out1, out2, info = paste0("function vs. method getBound call mismatch for compiled model with: ", info))
        expect_equal(out1, out3, info = paste0("uncompiled vs. compiled getBound call mismatch with: ", info))
        expect_equal(out1, nfOutput[1], info = paste0("direct vs. nimbleFunction getBound call mismatch for uncompiled model with: ", info))
        expect_equal(nfOutput[1], nfOutput[2], info = paste0("function vs. method getBound call mismatch for uncompiled nimbleFunction with: ", info))
        expect_equal(out3, cnfOutput[1], info = paste0("direct vs. nimbleFunction getBound call mismatch for compiled model with: ", info))
        expect_equal(cnfOutput[1], cnfOutput[2], info = paste0("function vs. method getBound call mismatch for compiled nimbleFunction with: ", info))
                                        #   if(.Platform$OS.type != 'windows') dyn.unload(project$cppProjects[[1]]$getSOName())
        })
    invisible(NULL)
}

## Nick's version of a nf that embeds deriv of model calculate
testCompiledModelDerivsNimFxn <- nimbleFunction(
  setup = function(model, calcNodes, wrtNodes, order){
  },
  run = function(){
    ansList <- nimDerivs(model$calculate(calcNodes), wrt = wrtNodes, order = order)
    returnType(ADNimbleList())
    return(ansList)
  }
)

## Chris' version of a nf that embeds deriv of model calculate
derivsNimbleFunction <- nimbleFunction(
  setup = function(model, calcNodes, wrt) {},
  run = function(x = double(1),
                 order = double(1),
                 reset = logical(0, default = FALSE)) {
    values(model, wrt) <<- x  
    ans <- nimDerivs(model$calculate(calcNodes), wrt = wrt, order = order, reset = reset)
    return(ans)
    returnType(ADNimbleList())
  }  ## don't need enableDerivs if call nimDerivs directly, but would if just have model$calc in nf
)

## nf for double-taping
derivsNimbleFunctionMeta <- nimbleFunction(
  setup = function(model, calcNodes, wrt, reset = FALSE) {
    innerWrtVec <- seq_along(model$expandNodeNames(wrt, returnScalarComponents = TRUE))
    d <- length(innerWrtVec)  
    allUpdateNodes <- makeUpdateNodes(wrt, calcNodes, model)
    updateNodes <- allUpdateNodes$updateNodes
    constantNodes <- allUpdateNodes$constantNodes
  },
  run = function(x = double(1)) {
    values(model, wrt) <<- x
    ans <- model$calculate(calcNodes)
    return(ans)
    returnType(double())
  },
  methods = list(
    ## inner first-order deriv
    derivs1Run = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = innerWrtVec, order = 1, reset = reset,
                         updateNodes = updateNodes, constantNodes = constantNodes, model = model)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## inner second-order deriv
    derivs2Run = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = innerWrtVec, order = 2, reset = reset,
                         updateNodes = updateNodes, constantNodes = constantNodes, model = model)
        ## not clear why can't do ans$hessian[,,1] with double(2) returnType
      return(ans$hessian)
      returnType(double(3))
    },
    ## outer arbitrary-order deriv calling inner first order
    metaDerivs1Run = function(x = double(1),
                             order = double(1),
                             reset = logical(0, default = FALSE)) {
      wrtVec <- 1:length(x)
      if(length(wrtVec) != d) stop("inner and outer wrt mismatch")  
      ans <- nimDerivs(derivs1Run(x), wrt = wrtVec, order = order, reset = reset,
                       updateNodes = updateNodes, constantNodes = constantNodes, model = model)
      return(ans)
      returnType(ADNimbleList())
    }, 
    ## outer arbitrary-order deriv calling inner second order
    metaDerivs2Run = function(x = double(1),
                             order = double(1),
                             reset = logical(0, default = FALSE)) {
      wrtVec <- 1:length(x)
      if(length(wrtVec) != d) stop("inner and outer wrt mismatch")  
      ans <- nimDerivs(derivs2Run(x), wrt = wrtVec, order = order, reset = reset,
                       updateNodes = updateNodes, constantNodes = constantNodes, model = model)
      return(ans)
      returnType(ADNimbleList())
    }
  ),
  enableDerivs = list(run = list(),
                      derivs1Run = list(),
                      derivs2Run = list())
)


derivsNimbleFunctionParamTransform <- nimbleFunction(
    setup = function(model, calcNodes, wrt) {
        wrtNodesAsScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
        my_parameterTransform <- parameterTransform(model, wrtNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        nimDerivs_wrt <- 1:d
        allUpdateNodes <- makeUpdateNodes(wrt, calcNodes, model)
        updateNodes   <- allUpdateNodes$updateNodes
        constantNodes <- allUpdateNodes$constantNodes
    },
    run = function(x = double(1),
                   order = double(1),
                   reset = logical(0, default = FALSE)) {
        transformed_x <- my_parameterTransform$transform(x)  # transform(x)
        ans <- nimDerivs(inverseTransformStoreCalculate(transformed_x), order = order, wrt = nimDerivs_wrt,
                         model = model, updateNodes = updateNodes,
                         constantNodes = constantNodes, reset = reset)

        return(ans)
        returnType(ADNimbleList())
    },
    methods = list(
        inverseTransformStoreCalculate = function(transformed_x = double(1)) {
            values(model, wrt) <<- my_parameterTransform$inverseTransform(transformed_x)
            lp <- model$calculate(calcNodes)
            returnType(double())
            return(lp)
        }
    ), enableDerivs = 'inverseTransformStoreCalculate'
)


derivsNimbleFunctionParamTransformMeta <- nimbleFunction(
    setup = function(model, calcNodes, wrt, reset = FALSE) {
        wrtNodesAsScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
        my_parameterTransform <- parameterTransform(model, wrtNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        nimDerivs_wrt <- 1:d
        allUpdateNodes <- makeUpdateNodes(wrt, calcNodes, model)
        updateNodes   <- allUpdateNodes$updateNodes
        constantNodes <- allUpdateNodes$constantNodes
    },
    ## formerly inverseTransformStoreCalculate
    run = function(transformed_x = double(1)) {
            values(model, wrt) <<- my_parameterTransform$inverseTransform(transformed_x)
            lp <- model$calculate(calcNodes)
            returnType(double())
            return(lp)
    },
    methods = list(
        derivs1Run = function(transformed_x = double(1)) {
            if(length(transformed_x) != d) stop("mismatch of x and wrtVec")
            ans <- nimDerivs(run(transformed_x), wrt = nimDerivs_wrt, order = 1, reset = reset,
                             updateNodes = updateNodes, constantNodes = constantNodes, model = model)
            return(ans$jacobian[1,])
            returnType(double(1))
        },
        derivs2Run = function(transformed_x = double(1)) {
            if(length(transformed_x) != d) stop("mismatch of x and wrtVec")
            ans <- nimDerivs(run(transformed_x), wrt = nimDerivs_wrt, order = 2, reset = reset,
                         updateNodes = updateNodes, constantNodes = constantNodes, model = model)
            ## not clear why can't do ans$hessian[,,1] with double(2) returnType
            return(ans$hessian)
            returnType(double(3))
        },
        metaDerivs1Run = function(x = double(1),
                                  order = double(1),
                                  reset = logical(0, default = FALSE)) {
            transformed_x <- my_parameterTransform$transform(x)
            if(length(transformed_x) != d) stop("mismatch of x and wrtVec")
            ans <- nimDerivs(derivs1Run(transformed_x), wrt = nimDerivs_wrt, order = order, reset = reset,
                             updateNodes = updateNodes, constantNodes = constantNodes, model = model)
            return(ans)
            returnType(ADNimbleList())
        }, 
        ## outer arbitrary-order deriv calling inner second order
        metaDerivs2Run = function(x = double(1),
                                  order = double(1),
                                  reset = logical(0, default = FALSE)) {
            transformed_x <- my_parameterTransform$transform(x)
            if(length(transformed_x) != d) stop("mismatch of x and wrtVec")
            ans <- nimDerivs(derivs2Run(transformed_x), wrt = nimDerivs_wrt, order = order, reset = reset,
                             updateNodes = updateNodes, constantNodes = constantNodes, model = model)
            return(ans)
            returnType(ADNimbleList())
        }       
    ),
    enableDerivs = list(run = list(),
                        derivs1Run = list(),
                        derivs2Run = list())
)




## For use with R-based derivative of compiled model calculate when useFasterRderivs is TRUE.
calcNodesForDerivs <- nimbleFunction(
        setup = function(model, calcNodes, wrt){
        },
    run = function(x = double(1)){
            values(model, wrt) <<- x
            ans <- calculate(model, calcNodes)
            return(ans)
            returnType(double())
        },
    enableDerivs = 'run')


## Tests taking derivatives of calls to model$calculate(nodes) (or equivalently calculate(model, nodes))
## Arguments:
##   model:         The uncompiled nimbleModel object to use in the call to calculate(model, nodes).
##   name:          The name of the model being tested.
##   calcNodeNames: A list, each element of which should be a character vector.   List elements  
##                  will be iterated through, and each element will be used as the 'nodes' argument
##                  in the call to calculate(model, nodes).
##   wrt:           A list, each element of which should be a character vector.  List elements will be iterated
##                  through, and each element will be used as the 'wrt' argument in a call to nimDerivs(calculate(model, nodes), wrt)
##   testR:         A logical argument.  If TRUE, the R version of nimDerivs will be checked for correct derivative calculations.
##                  This is accomplished by comparing derivatives calculated using the chain rule to derivatives of a function that
##                  wraps a call to calculate(model, nodes).
##   testCompiled:  A logical argument.  Currently only checks whether the model can compile.
##   tolerance:     A numeric argument, the tolerance to use when comparing wrapperDerivs to chainRuleDerivs.
##   verbose:       A logical argument.  Currently serves no purpose.
test_ADModelCalculate_nick <- function(model, name = NULL, calcNodeNames = NULL, wrt = NULL, order = c(0,1,2), 
                                  testCompiled = TRUE, tolerance = .001,  verbose = TRUE, gc = FALSE){
  temporarilyAssignInGlobalEnv(model)  

  if(testCompiled){
    expect_message(cModel <- compileNimble(model))
  }
  for(i in seq_along(calcNodeNames)){
    for(j in seq_along(wrt)){
      test_that(paste('R derivs of calculate function work for model', name, ', for calcNodes ', i, 
                      'and wrt ', j), {
                        wrapperDerivs <- nimDerivs(model$calculate(calcNodeNames[[i]]), wrt = wrt[[j]], order = order)
                        if(testCompiled){
                          print(calcNodeNames[[i]])
                          print(wrt[[j]])
                          testFunctionInstance <- testCompiledModelDerivsNimFxn(model, calcNodeNames[[i]], wrt[[j]], order)
                            if(gc) gc()
                            expect_message(ctestFunctionInstance <- compileNimble(testFunctionInstance, project =  model, resetFunctions = TRUE))
                            if(gc) gc()
                            cDerivs <- ctestFunctionInstance$run()
                          if(0 %in% order) expect_equal(wrapperDerivs$value, cDerivs$value, tolerance = tolerance)
                          if(1 %in% order) expect_equal(wrapperDerivs$jacobian, cDerivs$jacobian, tolerance = tolerance)
                          if(2 %in% order) expect_equal(wrapperDerivs$hessian, cDerivs$hessian, tolerance = tolerance)
                        }
                      })
    }
  }
}

## Chris' version of test_ADModelCalculate
## By default test a standardized set of {wrt, calcNodes} pairs representing common use cases (MAP, max lik, EB),
## unless user provides 'wrt' and 'calcNodes'.
test_ADModelCalculate <- function(model, name = 'unknown', x = 'given', xNew = NULL, calcNodes = NULL, wrt = NULL,
                                  newUpdateNodes = NULL, newConstantNodes = NULL,
                                  relTol = c(1e-15, 1e-8, 1e-3, 1e-3), absTolThreshold = 0, useFasterRderivs = FALSE, useParamTransform = FALSE,
                                  checkDoubleTape = TRUE, checkCompiledValuesIdentical = TRUE, checkDoubleUncHessian = TRUE,
                                  doAllUncHessian = TRUE, seed = 1, verbose = FALSE, debug = FALSE){
    if(!is.null(seed))
        set.seed(seed)
    initsHandling <- x
    xNewIn <- xNew
    ## Save model state so can restore for later use cases below.
    mv <- modelValues(model)
    nodes <- model$getNodeNames()
    nimCopy(model, mv, nodes, nodes, rowTo = 1, logProb = TRUE)

    if(is.null(wrt) && is.null(calcNodes)) {
        
        ## HMC/MAP use case
        if(verbose) cat("============================================\ntesting HMC/MAP-based scenario\n--------------------------------------------\n")
        calcNodes <- model$getNodeNames()
        wrt <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        tmp <- values(model, wrt)
        ## Hopefully values in (0,1) will always be legitimate for our test models.
        if(initsHandling %in% c('given','prior')) {
            x <- tmp
        } else if(initsHandling == 'random') {
            x <- runif(length(tmp))
        } 
        if(initsHandling == 'prior') {
            model$simulate(wrt)
            xNew <- values(model, wrt)
            nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        } else xNew <- runif(length(tmp))
        if(!is.null(xNewIn)) {
            wrtScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
            for(nm in names(xNewIn)) {
                wh <- grep(paste0("(^",nm,"$)|(^",nm,"\\[)"), wrtScalars)
                xNew[wh] <- c(xNewIn[[nm]])
            }
        }
        try(test_ADModelCalculate_internal(model, name = name, x = x, xNew = xNew, calcNodes = calcNodes, wrt = wrt,
                                           newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes,
                                           savedMV = mv, relTol = relTol, absTolThreshold = absTolThreshold,
                                           useFasterRderivs =  useFasterRderivs, useParamTransform = useParamTransform,
                                           checkDoubleTape = checkDoubleTape,
                                           checkCompiledValuesIdentical = checkCompiledValuesIdentical,
                                           checkDoubleUncHessian = checkDoubleUncHessian, doAllUncHessian = doAllUncHessian,
                                       verbose = verbose, debug = debug))
        ## max. lik. use case
        if(!is.null(seed))   ## There has weirdness where whether a test fails or passes modifies the RNG state, affecting the next call to test_ADModelCalculate_internal.
            set.seed(seed+1)
        if(verbose) cat("============================================\ntesting ML-based scenario\n--------------------------------------------\n")
        nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        calcNodes <- model$getNodeNames()
        topNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
        latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE, includeData = FALSE)
        calcNodes <- calcNodes[!calcNodes %in% c(topNodes, latentNodes)]  # should be data + deterministic
        wrt <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE) # will include hyps if present, but derivs wrt those should be zero
        tmp <- values(model, wrt)
        if(initsHandling %in% c('given','prior')) {
            x <- tmp
        } else if(initsHandling == 'random') {
            x <- runif(length(tmp))
        } 
        if(initsHandling == 'prior') {
            model$simulate(wrt)
            xNew <- values(model, wrt)
            nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        } else xNew <- runif(length(tmp))
        if(!is.null(xNewIn)) {
            wrtScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
            for(nm in names(xNewIn)) {
                wh <- grep(paste0("(^",nm,"$)|(^",nm,"\\[)"), wrtScalars)
                xNew[wh] <- c(xNewIn[[nm]])
            }
        }
        try(test_ADModelCalculate_internal(model, name = name, x = x, xNew = xNew, calcNodes = calcNodes, wrt = wrt, 
                                           newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes,
                                           savedMV =mv, relTol = relTol, absTolThreshold = absTolThreshold,
                                           useFasterRderivs =  useFasterRderivs, useParamTransform = useParamTransform,
                                           checkDoubleTape = checkDoubleTape,
                                           checkCompiledValuesIdentical = checkCompiledValuesIdentical,
                                           checkDoubleUncHessian = checkDoubleUncHessian, doAllUncHessian = doAllUncHessian,
                                           verbose = verbose, debug = debug))

        if(!is.null(seed))   ## There has weirdness where whether a test fails or passes modifies the RNG state, affecting the next call to test_ADModelCalculate_internal.
            set.seed(seed+2)
        ## modular HMC/MAP use case
        if(verbose) cat("============================================\ntesting HMC/MAP partial-based scenario\n--------------------------------------------\n")
        nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        calcNodes <- model$getNodeNames()
        wrt <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        wrtIdx <- sample(seq_along(wrt), round(length(wrt)/2), replace = FALSE)
        ## sample full wrt in case there are constraints built in, then subset wrt
        if(!length(wrtIdx))
            wrtIdx <- seq_along(wrt)
        tmp <- values(model, wrt)
        if(initsHandling %in% c('given','prior')) {
            x <- tmp
        } else if(initsHandling == 'random') {
            x <- runif(length(tmp))
        } 
        if(initsHandling == 'prior') {
            model$simulate(wrt)
            xNew <- values(model, wrt)
        } else xNew <- runif(length(tmp))
        if(!is.null(xNewIn)) {
            wrtScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
            for(nm in names(xNewIn)) {
                wh <- grep(paste0("(^",nm,"$)|(^",nm,"\\[)"), wrtScalars)
                xNew[wh] <- c(xNewIn[[nm]])
            }
        }
        wrtSub <- wrt[wrtIdx]
        ## get correct subset of x, xNew
        values(model, wrt) <- x
        x <- values(model, wrtSub)
        values(model, wrt) <- xNew
        xNew <- values(model, wrtSub)
        nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        try(test_ADModelCalculate_internal(model, name = name, x = x, xNew = xNew, calcNodes = calcNodes, wrt = wrtSub, 
                                           newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes,
                                           savedMV = mv, relTol = relTol, absTolThreshold = absTolThreshold,
                                           useFasterRderivs =  useFasterRderivs, useParamTransform = useParamTransform,
                                           checkDoubleTape = checkDoubleTape,
                                           checkCompiledValuesIdentical = checkCompiledValuesIdentical,
                                           checkDoubleUncHessian = checkDoubleUncHessian, doAllUncHessian = doAllUncHessian,
                                           verbose = verbose, debug = debug))

        if(!is.null(seed))   ## There has weirdness where whether a test fails or passes modifies the RNG state, affecting the next call to test_ADModelCalculate_internal.
            set.seed(seed+3)
        ## conditional max. lik. use case
        if(verbose) cat("============================================\ntesting ML partial-based scenario\n--------------------------------------------\n")
        nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        calcNodes <- model$getNodeNames()
        topNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
        latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE, includeData = FALSE)
        calcNodes <- calcNodes[!calcNodes %in% c(topNodes, latentNodes)]  # should be data + deterministic
        wrt <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        wrtIdx <- sample(seq_along(wrt), round(length(wrt)/2), replace = FALSE)
        ## sample full wrt in case there are constraints built in, then subset wrt
        if(!length(wrtIdx))
            wrtIdx <- seq_along(wrt)
        tmp <- values(model, wrt)
        if(initsHandling %in% c('given','prior')) {
            x <- tmp
        } else if(initsHandling == 'random') {
            x <- runif(length(tmp))
        } 
        if(initsHandling == 'prior') {
            model$simulate(wrt)
            xNew <- values(model, wrt)
        } else xNew <- runif(length(tmp))
        if(!is.null(xNewIn)) {
            wrtScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
            for(nm in names(xNewIn)) {
                wh <- grep(paste0("(^",nm,"$)|(^",nm,"\\[)"), wrtScalars)
                xNew[wh] <- c(xNewIn[[nm]])
            }
        }
        wrtSub <- wrt[wrtIdx]
        ## get correct subset of x, xNew
        values(model, wrt) <- x
        x <- values(model, wrtSub)
        values(model, wrt) <- xNew
        xNew <- values(model, wrtSub)
        nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        try(test_ADModelCalculate_internal(model, name = name, x = x, xNew = xNew, calcNodes = calcNodes, wrt = wrtSub, 
                                           newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes,
                                           savedMV = mv, relTol = relTol, absTolThreshold = absTolThreshold,
                                           useFasterRderivs =  useFasterRderivs, useParamTransform = useParamTransform,
                                           checkDoubleTape = checkDoubleTape,
                                           checkCompiledValuesIdentical = checkCompiledValuesIdentical,
                                           checkDoubleUncHessian = checkDoubleUncHessian, doAllUncHessian = doAllUncHessian,
                                           verbose = verbose, debug = debug))

        if(!is.null(seed))   ## There has weirdness where whether a test fails or passes modifies the RNG state, affecting the next call to test_ADModelCalculate_internal.
            set.seed(seed+4)
        ## empirical Bayes use case (though not actually integrating over any latent nodes)
        if(verbose) cat('============================================\ntesting EB-based scenario\n--------------------------------------------\n')
        nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        calcNodes <- model$getNodeNames()
        topNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
        calcNodes <- calcNodes[!calcNodes %in% topNodes]  # EB doesn't use hyperpriors
        wrt <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        tmp <- values(model, wrt)
        if(initsHandling %in% c('given','prior')) {
            x <- tmp
        } else if(initsHandling == 'random') {
            x <- runif(length(tmp))
        } 
        if(initsHandling == 'prior') {
            model$simulate(wrt)
            xNew <- values(model, wrt)
            nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        } else xNew <- runif(length(tmp))
        if(!is.null(xNewIn)) {
            wrtScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
            for(nm in names(xNewIn)) {
                wh <- grep(paste0("(^",nm,"$)|(^",nm,"\\[)"), wrtScalars)
                xNew[wh] <- c(xNewIn[[nm]])
            }
        }
        try(test_ADModelCalculate_internal(model, name = name, x = x, xNew = xNew, calcNodes = calcNodes, wrt = wrt, 
                                           newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes,
                                           savedMV = mv, relTol = relTol, absTolThreshold = absTolThreshold,
                                           useFasterRderivs =  useFasterRderivs, useParamTransform = useParamTransform,
                                           checkDoubleTape = checkDoubleTape,
                                           checkCompiledValuesIdentical = checkCompiledValuesIdentical,
                                           checkDoubleUncHessian = checkDoubleUncHessian, doAllUncHessian = doAllUncHessian,
                                           verbose = verbose, debug = debug))
    } else {
        if(is.null(calcNodes)) calcNodes <- model$getNodeNames()
        if(is.null(wrt)) wrt <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        ## Apply test to user-provided sets of nodes
        tmp <- values(model, wrt)
        if(initsHandling %in% c('given','prior')) {
            x <- tmp
        } else if(initsHandling == 'random') {
            x <- runif(length(tmp))
        } 
        if(initsHandling == 'prior') {
            model$simulate(wrt)
            xNew <- values(model, wrt)
            nimCopy(mv, model, nodes, nodes, row = 1, logProb = TRUE)
        } else xNew <- runif(length(tmp))
        if(!is.null(xNewIn)) {
            wrtScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
            for(nm in names(xNewIn)) {
                wh <- grep(paste0("(^",nm,"$)|(^",nm,"\\[)"), wrtScalars)
                xNew[wh] <- c(xNewIn[[nm]])
            }
        }
        try(test_ADModelCalculate_internal(model, name = name, x = x, xNew = xNew, calcNodes = calcNodes, wrt = wrt, 
                                           newUpdateNodes = newUpdateNodes, newConstantNodes = newConstantNodes,
                                           relTol = relTol, absTolThreshold = absTolThreshold,
                                           useFasterRderivs =  useFasterRderivs, useParamTransform = useParamTransform,
                                           checkDoubleTape = checkDoubleTape,
                                           checkCompiledValuesIdentical = checkCompiledValuesIdentical,
                                           checkDoubleUncHessian = checkDoubleUncHessian, doAllUncHessian = doAllUncHessian,
                                           verbose = verbose, debug = debug))
    }
}


## This does the core assessment, by default running with various sets of order values to be able to assess
## forward and backward mode and to assess whether values in the model are updated.
test_ADModelCalculate_internal <- function(model, name = 'unknown', xOrig = NULL, xNew = NULL,
                                           calcNodes = NULL, wrt = NULL, savedMV = NULL, 
                                           newUpdateNodes = NULL, newConstantNodes = NULL, 
                                           relTol = c(1e-15, 1e-8, 1e-3, 1e-3), absTolThreshold = 0, useFasterRderivs = FALSE,
                                           useParamTransform = FALSE, checkDoubleTape = TRUE, 
                                           checkCompiledValuesIdentical = TRUE, checkDoubleUncHessian = TRUE,
                                           doAllUncHessian = TRUE,
                                           verbose = FALSE, debug = FALSE){

    saved_edition <- edition_get()
    local_edition(3)
    on.exit(local_edition(saved_edition))
    
    test_that(paste0("Derivatives of calculate for model ", name), {
        if(exists('paciorek') && paciorek == 0) browser()
        if(is.null(calcNodes))
            calcNodes <- model$getNodeNames()

        nodes  <- model$getNodeNames()
        if(!(exists('CobjectInterface', model) && !is(model$CobjectInterface, 'uninitializedField'))) {
            cModel <- compileNimble(model)
        } else {
            cModel <- eval(quote(model$CobjectInterface))
            nimCopy(savedMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
        }
        tmpMV <- modelValues(model)

        if(is.null(wrt)) wrt <- calcNodes
        wrt <- model$expandNodeNames(wrt, unique = FALSE)
        discrete <- sapply(wrt, function(x) model$isDiscrete(x))
        wrt <- wrt[is.na(discrete) | !discrete]  # NAs from deterministic nodes, which in most cases should be continuous

        otherNodes <- model$getNodeNames()
        otherNodes <- otherNodes[!otherNodes %in% wrt]

        checkSquareMatrix <- function(node) {
            val <- model[[node]]
            d <- dim(val)
            if(!is.null(d) && length(d) == 2 && d[1] == d[2]) return(TRUE) else return(FALSE)
        }        
        
        allUpdateNodes <- makeUpdateNodes(wrt, calcNodes, model)
        updateNodes <- allUpdateNodes$updateNodes
        constantNodes <- allUpdateNodes$constantNodes

        updateNodesDeps <- model$getDependencies(updateNodes)
        constantNodesDeps <- model$getDependencies(constantNodes)
        
        if(useParamTransform) {
            rDerivs <- derivsNimbleFunctionParamTransform(model, calcNodes = calcNodes, wrt = wrt)
        } else rDerivs <- derivsNimbleFunction(model, calcNodes = calcNodes, wrt = wrt)

        ## Seem to need resetFunctions because of use of my_parameterTransform twice in paramTransform case.
        cDerivs <- compileNimble(rDerivs, project = model, resetFunctions = TRUE)

        if(checkDoubleTape) {
            if(useParamTransform) {
                rDerivsMeta <- derivsNimbleFunctionParamTransformMeta(model, calcNodes = calcNodes, wrt = wrt)
                rDerivsMetaReset <- derivsNimbleFunctionParamTransformMeta(model, calcNodes = calcNodes, wrt = wrt, reset = TRUE)
            } else {
                rDerivsMeta <- derivsNimbleFunctionMeta(model, calcNodes = calcNodes, wrt = wrt)
                rDerivsMetaReset <- derivsNimbleFunctionMeta(model, calcNodes = calcNodes, wrt = wrt, reset = TRUE)
            }
            cDerivsMeta <- compileNimble(rDerivsMeta, project = model, resetFunctions = TRUE)
            cDerivsMetaReset <- compileNimble(rDerivsMetaReset, project = model, resetFunctions = TRUE)
        }
        
        if(useFasterRderivs) {
            ## Set up a nf so R derivs use a model calculate that is done fully in compiled code (cModel$calculate loops over nodes in R)
            if(useParamTransform) {
                ## Need wrapper so that we are calling nimDerivs on a function call and not a nf method
                wrapper <- function(x) {
                    cDerivs$inverseTransformStoreCalculate(x)
                }
            } else {
                rCalcNodes <- calcNodesForDerivs(model, calcNodes = calcNodes, wrt = wrt)  
                cCalcNodes <- compileNimble(rCalcNodes, project = model)
                
                ## temporarilyAssignInGlobalEnv(cCalcNodes)
                
                ## Need wrapper so that we are calling nimDerivs on a function call and not a nf method
                wrapper <- function(x) {
                    cCalcNodes$run(x)
                }
            }
            if(checkDoubleTape) {
                wrapperMeta1 <- function(x) {
                    ans <- nimDerivs(wrapper(x), order = 1, reset = FALSE)
                    return(ans$jacobian[1, ])
                }
                wrapperMeta1Reset <- function(x) {
                    ans <- nimDerivs(wrapper(x), order = 1, reset = TRUE)
                    return(ans$jacobian[1, ])
                }
                wrapperMeta2 <- function(x) {
                    ans <- nimDerivs(wrapper(x), order = 2, reset = FALSE)
                    return(ans$hessian)
                }
                wrapperMeta2Reset <- function(x) {
                    ans <- nimDerivs(wrapper(x), order = 2, reset = TRUE)
                    return(ans$hessian)
                }
            }
        }

        xList <- list(xOrig)
        if(!is.null(xNew)) {
            xList[[2]] <- xNew
            xList[[3]] <- xNew
        }

        for(case in 1:2) {
            for(idx in seq_along(xList)) {
                if(exists('paciorek') && paciorek == idx) browser()
                if(verbose) {
                    if(case == 1) {
                        if(idx == 1) {
                            cat("Testing initial wrt values with initial constantNodes\n")
                            cat("  Using wrt: ", wrt, "\n")
                        }
                        if(idx == 2) cat("Testing new wrt values with initial constantNodes\n")
                        if(idx == 3 && length(updateNodes)) {
                             cat("Testing new updateNode values with initial constantNodes\n")
                            cat("  Using updateNodes: ", updateNodes, "\n")
                        }
                    } else {
                        if(idx == 1 && length(constantNodes)) {
                            cat("Testing initial wrt values with new constantNodes\n")
                            cat("  Using constantNodes: ", constantNodes, "\n")
                        }
                        if(idx == 2 && length(constantNodes)) cat("Testing new wrt values with new constantNodes\n")
                        if(idx == 3 && length(constantNodes) && length(updateNodes)) cat("Testing new updateNode values with new constantNodes\n")
                    }
                }

                if(idx == 3) {
                    if(length(updateNodes)) {
                        values(model, updateNodes) <- runif(length(updateNodes))
                        values(cModel, updateNodes) <- values(model, updateNodes)
                        ## Overwrite for tricky cases (p.d. matrices, integers, etc.)
                        if(!is.null(newUpdateNodes)) {
                            for(nm in names(newUpdateNodes)) {
                                ## Expansion is needed so that fill by element, which retains dimensionality in uncompiled model.
                                values(model, model$expandNodeNames(nm, returnScalarComponents = TRUE)) <- newUpdateNodes[[nm]]
                                values(cModel, nm) <- values(model, nm)
                            }
                        }
                        model$calculate(updateNodesDeps)
                        cModel$calculate(updateNodesDeps)
                    } else next
                }

                reset <- FALSE
                if(case == 2) {
                    reset <- TRUE
                    if(length(constantNodes)) {
                        values(model, constantNodes) <- runif(length(constantNodes))
                        values(cModel, constantNodes) <- values(model, constantNodes)
                        
                        if(!is.null(newConstantNodes)) {
                            for(nm in names(newConstantNodes)) {
                                ## Expansion is needed so that fill by element, which retains dimensionality in uncompiled model.
                                values(model, model$expandNodeNames(nm, returnScalarComponents = TRUE)) <- newConstantNodes[[nm]]
                                values(cModel, nm) <- values(model, nm)
                            }
                        }
                        model$calculate(constantNodesDeps)
                        cModel$calculate(constantNodesDeps)
                    } else next 
                }
                
                x <- xList[[idx]]

                nimCopy(model, tmpMV, nodes, nodes, row = 1, logProb = TRUE)
                ## Ensure both models are consistent.
                nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                
                rWrt_orig <- cWrt_orig <- values(model, wrt)
                
                ## Store current logProb and non-wrt values to check that order=c(1,2) doesn't change them.
                ## Don't calculate model as want to assess possibility model is out-of-state.
                ## model$calculate()
                ## cModel$calculate()
                
                rLogProb_orig <- cLogProb_orig <- model$getLogProb(calcNodes)
                rVals_orig <- cVals_orig <- values(model, otherNodes)
                
                if(useFasterRderivs) {
                    inputx <- x
                    if(useParamTransform)
                        inputx <- rDerivs$my_parameterTransform$transform(x)

                    
                    rOutput012 <- nimDerivs(wrapper(inputx), order = 0:2, reset = reset)
                    rVals012 <- values(cModel, otherNodes)
                    rLogProb012 <- cModel$getLogProb(calcNodes)
                    rWrt012 <- values(cModel, wrt)
                    
                    nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                    rOutput01 <- nimDerivs(wrapper(inputx), order = 0:1, reset = reset)
                    rLogProb01 <- cModel$getLogProb(calcNodes)
                    rVals01 <- values(cModel, otherNodes)
                    rWrt01 <- values(cModel, wrt)

                    if(doAllUncHessian) { 
                        nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                        rOutput12 <- nimDerivs(wrapper(inputx), order = 1:2, model = cModel, reset = reset)
                        rVals12 <- values(cModel, otherNodes)
                        rLogProb12 <- cModel$getLogProb(calcNodes)
                        rWrt12 <- values(cModel, wrt)
                        
                        nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                        rOutput02 <- nimDerivs(wrapper(inputx), order = c(0,2), reset = reset)
                        rVals02 <- values(cModel, otherNodes)
                        rLogProb02 <- cModel$getLogProb(calcNodes)
                        rWrt02 <- values(cModel, wrt)
                    }
                    
                    if(checkDoubleTape) {
                        ## Note that because inner deriv is order 1 or 2, don't expect model to be updated,
                        ## so need to do this before 01, 012 cases below.
                        nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                        if(reset) {
                            rOutput1d <- nimDerivs(wrapperMeta1Reset(inputx), order = 0, model = cModel, reset = reset)
                        } else rOutput1d <- nimDerivs(wrapperMeta1(inputx), order = 0, model = cModel, reset = reset)
                        rVals1d <- values(cModel, otherNodes)
                        rLogProb1d <- cModel$getLogProb(calcNodes)
                        rWrt1d <- values(cModel, wrt)

                        if(doAllUncHessian) { 
                            nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                            if(reset) {
                                rOutput2d <- nimDerivs(wrapperMeta2Reset(inputx), order = 0, model = cModel, reset = reset)
                            } else rOutput2d <- nimDerivs(wrapperMeta2(inputx), order = 0, model = cModel, reset = reset)
                            rVals2d <- values(cModel, otherNodes)
                            rLogProb2d <- cModel$getLogProb(calcNodes)
                            rWrt2d <- values(cModel, wrt)
                        }
                        
                        if(checkDoubleUncHessian) {
                            nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                            if(reset) {
                                rOutput2d11 <- nimDerivs(wrapperMeta1Reset(inputx), order = 1, model = cModel, reset = reset)
                            } else rOutput2d11 <- nimDerivs(wrapperMeta1(inputx), order = 1, model = cModel, reset = reset)
                            
                            rVals2d11 <- values(cModel, otherNodes)
                            rLogProb2d11 <- cModel$getLogProb(calcNodes)
                            rWrt2d11 <- values(cModel, wrt)
                        }
                    }

                    rLogProb_new <- wrapper(inputx)
                    rVals_new <- values(cModel, otherNodes)

                    ## now reset cModel for use in compiled derivs
                    nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)

                } else {
                         
                    rOutput012 <- rDerivs$run(x, 0:2, reset = reset)
                    rVals012 <- values(model, otherNodes)
                    rLogProb012 <- model$getLogProb(calcNodes)
                    rWrt012 <- values(model, wrt)

                    nimCopy(tmpMV, model, nodes, nodes, row = 1, logProb = TRUE)
                    rOutput01 <- rDerivs$run(x, 0:1, reset = reset)
                    rLogProb01 <- model$getLogProb(calcNodes)
                    rVals01 <- values(model, otherNodes)
                    rWrt01 <- values(model, wrt)

                    if(doAllUncHessian) {
                        nimCopy(tmpMV, model, nodes, nodes, row = 1, logProb = TRUE)
                        rOutput12 <- rDerivs$run(x, 1:2, reset = reset)
                        rVals12 <- values(model, otherNodes)
                        rLogProb12 <- model$getLogProb(calcNodes)
                        rWrt12 <- values(model, wrt)
                        
                        nimCopy(tmpMV, model, nodes, nodes, row = 1, logProb = TRUE)
                        rOutput02 <- rDerivs$run(x, c(0,2), reset = reset)
                        rVals02 <- values(model, otherNodes)
                        rLogProb02 <- model$getLogProb(calcNodes)
                        rWrt02 <- values(model, wrt)
                    }
                    
                    if(checkDoubleTape) {
                        ## Note that because inner deriv is order 1 or 2, don't expect model to be updated.
                        nimCopy(tmpMV, model, nodes, nodes, row = 1, logProb = TRUE)
                        if(reset) {
                            rOutput1d <- rDerivsMetaReset$metaDerivs1Run(x = x, order = 0, reset = reset)
                        } else rOutput1d <- rDerivsMeta$metaDerivs1Run(x = x, order = 0, reset = reset)
                        rVals1d <- values(model, otherNodes)
                        rLogProb1d <- model$getLogProb(calcNodes)
                        rWrt1d <- values(model, wrt)
                        
                        if(doAllUncHessian) {
                            nimCopy(tmpMV, model, nodes, nodes, row = 1, logProb = TRUE)
                            if(reset) {
                                rOutput2d <- rDerivsMetaReset$metaDerivs2Run(x = x, order = 0, reset = reset)
                            } else rOutput2d <- rDerivsMeta$metaDerivs2Run(x = x, order = 0, reset = reset)
                            rVals2d <- values(model, otherNodes)
                            rLogProb2d <- model$getLogProb(calcNodes)
                            rWrt2d <- values(model, wrt)
                        }
                            
                        if(checkDoubleUncHessian) {
                            nimCopy(tmpMV, model, nodes, nodes, row = 1, logProb = TRUE)
                            if(reset) {
                                rOutput2d11 <- rDerivsMetaReset$metaDerivs1Run(x = x, order = 1, reset = reset)
                            } else rOutput2d11 <- rDerivsMeta$metaDerivs1Run(x = x, order = 1, reset = reset)
                            rVals2d11 <- values(model, otherNodes)
                            rLogProb2d11 <- model$getLogProb(calcNodes)
                            rWrt2d11 <- values(model, wrt)
                        }
                    }

                    values(model, wrt) <- x
                    rLogProb_new <- model$calculate(calcNodes)
                    rVals_new <- values(model, otherNodes)
                }

                ## Without useParamTransform, wrt should be updated, because assignment into model is outside nimDerivs call,
                ## but with useParamTransform they should not,
                cOutput12 <- cDerivs$run(x, 1:2, reset = reset)
                cVals12 <- values(cModel, otherNodes)
                cLogProb12 <- cModel$getLogProb(calcNodes)
                cWrt12 <- values(cModel, wrt)

                nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                cOutput01 <- cDerivs$run(x, 0:1, reset = reset)
                cVals01 <- values(cModel, otherNodes)
                cLogProb01 <- cModel$getLogProb(calcNodes)
                cWrt01 <- values(cModel, wrt)

                nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                cOutput012 <- cDerivs$run(x, 0:2, reset = reset)
                cVals012 <- values(cModel, otherNodes)
                cLogProb012 <- cModel$getLogProb(calcNodes)
                cWrt012 <- values(cModel, wrt)

                nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                cOutput02 <- cDerivs$run(x, c(0,2), reset = reset)
                cVals02 <- values(cModel, otherNodes)
                cLogProb02 <- cModel$getLogProb(calcNodes)
                cWrt02 <- values(cModel, wrt)

                if(checkDoubleTape) {
                    ## Note that because inner deriv is order 1 or 2, don't expect model to be updated.
                    nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                    if(reset) {
                        cOutput1d <- cDerivsMetaReset$metaDerivs1Run(x = x, order = 0, reset = reset)
                    } else cOutput1d <- cDerivsMeta$metaDerivs1Run(x = x, order = 0, reset = reset)
                    cVals1d <- values(cModel, otherNodes)
                    cLogProb1d <- cModel$getLogProb(calcNodes)
                    cWrt1d <- values(cModel, wrt)
                    
                    nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                    if(reset) {
                        cOutput2d <- cDerivsMetaReset$metaDerivs2Run(x = x, order = 0, reset = reset)
                    } else cOutput2d <- cDerivsMeta$metaDerivs2Run(x = x, order = 0, reset = reset)
                    cVals2d <- values(cModel, otherNodes)
                    cLogProb2d <- cModel$getLogProb(calcNodes)
                    cWrt2d <- values(cModel, wrt)

                    nimCopy(tmpMV, cModel, nodes, nodes, row = 1, logProb = TRUE)
                    if(reset) {
                        cOutput2d11 <- cDerivsMetaReset$metaDerivs1Run(x = x, order = 1, reset = reset)
                    } else cOutput2d11 <- cDerivsMeta$metaDerivs1Run(x = x, order = 1, reset = reset)
                    cVals2d11 <- values(cModel, otherNodes)
                    cLogProb2d11 <- cModel$getLogProb(calcNodes)
                    cWrt2d11 <- values(cModel, wrt)
                }

                values(cModel, wrt) <- x
                cLogProb_new <- cModel$calculate(calcNodes)
                cVals_new <- values(cModel, otherNodes)
                
                ## Check results ##

                if(checkCompiledValuesIdentical) {
                    expect_fun <- expect_identical
                } else expect_fun <- nim_expect_equal
                
                ## Check that only requested orders provided.

                expect_identical(length(rOutput012$value), 1L)
                expect_identical(length(rOutput01$value), 1L)
                if(doAllUncHessian) {
                    expect_identical(length(rOutput12$value), 0L)
                    expect_identical(length(rOutput02$value), 1L)
                }
                expect_identical(length(cOutput01$value), 1L)
                expect_identical(length(cOutput12$value), 0L)
                expect_identical(length(cOutput012$value), 1L)
                expect_identical(length(cOutput02$value), 1L)

                expect_gte(length(rOutput012$jacobian), 1)
                expect_gte(length(rOutput01$jacobian), 1)
                if(doAllUncHessian) {
                    expect_gte(length(rOutput12$jacobian), 1)
                    expect_identical(length(rOutput02$jacobian), 0L)
                }
                expect_gte(length(cOutput01$jacobian), 1)
                expect_gte(length(cOutput12$jacobian), 1)
                expect_gte(length(cOutput012$jacobian), 1)
                expect_identical(length(cOutput02$jacobian), 0L)

                expect_gte(length(rOutput012$hessian), 1)
                expect_identical(length(rOutput01$hessian), 0L)
                if(doAllUncHessian) {
                    expect_gte(length(rOutput12$hessian), 1)
                    expect_gte(length(rOutput02$hessian), 1)
                }
                expect_identical(length(cOutput01$hessian), 0L)
                expect_gte(length(cOutput12$hessian), 1)
                expect_gte(length(cOutput012$hessian), 1)
                expect_gte(length(cOutput02$hessian), 1)

                if(checkDoubleTape) {
                    expect_gte(length(rOutput1d$value), 1)
                    expect_identical(length(rOutput1d$jacobian), 0L)
                    expect_identical(length(rOutput1d$hessian), 0L)
                    expect_gte(length(cOutput1d$value), 1)
                    expect_identical(length(cOutput1d$jacobian), 0L)
                    expect_identical(length(cOutput1d$hessian), 0L)

                    if(doAllUncHessian) {
                        expect_gte(length(rOutput2d$value), 1)
                        expect_identical(length(rOutput2d$jacobian), 0L)
                        expect_identical(length(rOutput2d$hessian), 0L)
                    }
                    expect_gte(length(cOutput2d$value), 1)
                    expect_identical(length(cOutput2d$jacobian), 0L)
                    expect_identical(length(cOutput2d$hessian), 0L)

                    if(checkDoubleUncHessian) {
                        expect_identical(length(rOutput2d11$value), 0L)
                        expect_gte(length(rOutput2d11$jacobian), 1)
                        expect_identical(length(rOutput2d11$hessian), 0L)
                    }
                    expect_identical(length(cOutput2d11$value), 0L)
                    expect_gte(length(cOutput2d11$jacobian), 1)
                    expect_identical(length(cOutput2d11$hessian), 0L)
                }
                
                ## 0th order 'derivative'
                expect_identical(rOutput01$value, rLogProb_new)
                expect_identical(rOutput012$value, rLogProb_new)
                if(doAllUncHessian) 
                    expect_identical(rOutput02$value, rLogProb_new)
                expect_fun(cOutput01$value, cLogProb_new)
                expect_fun(cOutput012$value, cLogProb_new)
                expect_fun(cOutput02$value, cLogProb_new)

                nim_expect_equal(rOutput01$value, cOutput01$value, tolerance = relTol[1], abs_threshold = absTolThreshold)
                nim_expect_equal(rOutput012$value, cOutput012$value, tolerance = relTol[1], abs_threshold = absTolThreshold)
                if(doAllUncHessian) 
                    nim_expect_equal(rOutput02$value, cOutput02$value, tolerance = relTol[1], abs_threshold = absTolThreshold)

                if(FALSE) {  ## not relevant if using nim_expect_equal
                ## expect_equal (via waldo::compare and waldo::num_equal) uses absolute tolerance if 'y' value <= tolerance.
                ## Consider creating a nim_equal operator that checks if max(abs(x-y)/abs(y)) > tolerance using
                ## relative tolerance unless y == 0.
                    if(mean(abs(cOutput01$value)) <= relTol[1])
                        warning("Using absolute tolerance for 01$value comparison.")
                    if(mean(abs(cOutput012$value)) <= relTol[1])
                        warning("Using absolute tolerance for 012$value comparison.")
                    if(mean(abs(cOutput02$value)) <= relTol[1])
                        warning("Using absolute tolerance for 02$value comparison.")
                }
                
                expect_identical(sum(is.na(rOutput01$value)), 0L, info = "NAs found in uncompiled 0th derivative")
                expect_identical(sum(is.na(cOutput01$value)), 0L, info = "NAs found in compiled 0th derivative")
                expect_identical(sum(is.na(rOutput012$value)), 0L, info = "NAs found in uncompiled 0th derivative")
                expect_identical(sum(is.na(cOutput012$value)), 0L, info = "NAs found in compiled 0th derivative")
                if(doAllUncHessian)
                    expect_identical(sum(is.na(rOutput02$value)), 0L, info = "NAs found in uncompiled 0th derivative")
                expect_identical(sum(is.na(cOutput02$value)), 0L, info = "NAs found in compiled 0th derivative")
                
                ## 1st derivative
                nim_expect_equal(rOutput01$jacobian, cOutput01$jacobian, tolerance = relTol[2], abs_threshold = absTolThreshold)
                expect_identical(sum(is.na(rOutput01$jacobian)), 0L, info = "NAs found in uncompiled 1st derivative")
                expect_identical(sum(is.na(cOutput01$jacobian)), 0L, info = "NAs found in compiled 1st derivative")

                if(doAllUncHessian) {
                    nim_expect_equal(rOutput12$jacobian, cOutput12$jacobian, tolerance = relTol[2], abs_threshold = absTolThreshold)
                    expect_identical(sum(is.na(rOutput12$jacobian)), 0L, info = "NAs found in uncompiled 1st derivative")
                }
                expect_identical(sum(is.na(cOutput12$jacobian)), 0L, info = "NAs found in compiled 1st derivative")

                nim_expect_equal(rOutput012$jacobian, cOutput012$jacobian, tolerance = relTol[2], abs_threshold = absTolThreshold)
                expect_identical(sum(is.na(rOutput012$jacobian)), 0L, info = "NAs found in uncompiled 1st derivative")
                expect_identical(sum(is.na(cOutput012$jacobian)), 0L, info = "NAs found in compiled 1st derivative")

                if(FALSE) {  ## not relevant if using nim_expect_equal
                    if(mean(abs(cOutput01$jacobian)) <= relTol[2])
                        warning("Using absolute tolerance for 01$jacobian comparison.")
                    if(mean(abs(cOutput12$jacobian)) <= relTol[2])
                        warning("Using absolute tolerance for 12$jacobian comparison.")
                    if(mean(abs(cOutput012$jacobian)) <= relTol[2])
                        warning("Using absolute tolerance for 012$jacobian comparison.")
                }
                
                ## explicit comparison of first derivs;
                ## both of these are reverse mode because 2nd order reverse also invokes first order reverse
                expect_identical(cOutput01$jacobian, cOutput012$jacobian)
                
                expect_identical(cOutput12$jacobian, cOutput012$jacobian)

                if(checkDoubleTape) {
                    nim_expect_equal(rOutput1d$value, cOutput1d$value, tolerance = relTol[2], abs_threshold = absTolThreshold)
                    expect_identical(sum(is.na(rOutput1d$value)), 0L, info = "NAs found in uncompiled double-taped 1st derivative")
                    expect_identical(sum(is.na(cOutput1d$value)), 0L, info = "NAs found in compiled double-taped 1st derivative")

                    ## explicit comparison to single-taped result
                    expect_fun(cOutput1d$value, c(cOutput012$jacobian))

                    if(FALSE) {  ## not relevant if using nim_expect_equal
                        if(mean(abs(cOutput1d$value)) <= relTol[2])
                            warning("Using absolute tolerance for 1d$value comparison.")
                    }
                }

                ## 2nd derivative
                if(doAllUncHessian) {
                    nim_expect_equal(rOutput12$hessian, cOutput12$hessian, tolerance = relTol[3], abs_threshold = absTolThreshold)
                    expect_identical(sum(is.na(rOutput12$hessian)), 0L, info = "NAs found in uncompiled 2nd derivative")
                }
                expect_identical(sum(is.na(cOutput12$hessian)), 0L, info = "NAs found in compiled 2nd derivative")

                nim_expect_equal(rOutput012$hessian, cOutput012$hessian, tolerance = relTol[3], abs_threshold = absTolThreshold)
                expect_identical(sum(is.na(rOutput012$hessian)), 0L, info = "NAs found in uncompiled 2nd derivative")
                expect_identical(sum(is.na(cOutput012$hessian)), 0L, info = "NAs found in compiled 2nd derivative")

                if(doAllUncHessian) {
                    nim_expect_equal(rOutput02$hessian, cOutput02$hessian, tolerance = relTol[3], abs_threshold = absTolThreshold)
                    expect_identical(sum(is.na(rOutput02$hessian)), 0L, info = "NAs found in uncompiled 2nd derivative")
                }
                expect_identical(sum(is.na(cOutput02$hessian)), 0L, info = "NAs found in compiled 2nd derivative")

                if(FALSE) {  ## not relevant if using nim_expect_equal
                    if(mean(abs(cOutput12$hessian)) <= relTol[3])
                        warning("Using absolute tolerance for 12$hessian comparison.")
                    if(mean(abs(cOutput012$hessian)) <= relTol[3])
                        warning("Using absolute tolerance for 012$hessian comparison.")
                    if(mean(abs(cOutput02$hessian)) <= relTol[3])
                        warning("Using absolute tolerance for 02$hessian comparison.")
                }
                
                expect_identical(cOutput12$hessian, cOutput012$hessian)
                expect_identical(cOutput02$hessian, cOutput012$hessian)


                if(checkDoubleTape) {
                    if(doAllUncHessian) {
                        nim_expect_equal(rOutput2d$value, cOutput2d$value, tolerance = relTol[3], abs_threshold = absTolThreshold)
                        expect_identical(sum(is.na(rOutput2d$value)), 0L, info = "NAs found in uncompiled double-taped 2nd derivative")
                    }
                    expect_identical(sum(is.na(cOutput2d$value)), 0L, info = "NAs found in compiled double-taped 2nd derivative")

                    if(checkDoubleUncHessian) {
                        nim_expect_equal(rOutput2d11$jacobian, cOutput2d11$jacobian, tolerance = relTol[4], abs_threshold = absTolThreshold)
                        expect_identical(sum(is.na(rOutput2d11$jacobian)), 0L, info = "NAs found in uncompiled double-taped 2nd derivative")
                    }
                    expect_identical(sum(is.na(cOutput2d11$jacobian)), 0L, info = "NAs found in compiled double-taped 2nd derivative")

                    ## explicit comparison to single-taped result
                    ## Not clear why 2d$value not identical to 012$hessian
                    nim_expect_equal(cOutput2d$value, c(cOutput012$hessian), tolerance = 1e-15, abs_threshold = absTolThreshold)
                    if(length(cOutput2d11$jacobian) == 1) cOutput2d11$jacobian <- c(cOutput2d11$jacobian)
                    expect_fun(cOutput2d11$jacobian, cOutput012$hessian[,,1])

                    if(FALSE) {  ## not relevant if using nim_expect_equal
                        if(mean(abs(cOutput2d$value)) <= relTol[3])
                            warning("Using absolute tolerance for 2d$value comparison.")
                        if(mean(abs(cOutput2d11$jacobian)) <= relTol[4])
                            warning("Using absolute tolerance for 2d11$jacobian comparison.")
                    }
                }

                ## wrt values should equal original wrt values if order !=0 or doubleTape
                ## because setting of wrt in model is done within nimDerivs call, so should obey our rules about when model state is altered.

                expect_identical(rWrt01, x)
                ## We provide the model to nimDerivs, so expect restoration except for non-paramTransform
                ## and non-fasterRderivs, where assignment is outside nimDerivs.
                if(doAllUncHessian) {
                    if(!useParamTransform && !useFasterRderivs) {
                        expect_identical(rWrt12, x)
                    } else expect_identical(rWrt12, rWrt_orig)
                }
                expect_identical(rWrt012, x)
                expect_identical(cWrt01, x)
                if(!useParamTransform) {
                    ## Assignment to model is outside nimDerivs call, so expect change.
                    expect_identical(cWrt12, x)
                } else expect_identical(cWrt12, cWrt_orig)
                expect_identical(cWrt012, x)

                if(checkDoubleTape) {
                    expect_identical(rWrt1d, rWrt_orig)
                    if(doAllUncHessian) 
                        expect_identical(rWrt2d, rWrt_orig)
                    if(checkDoubleUncHessian)
                        expect_identical(rWrt2d11, rWrt_orig)
                    expect_identical(cWrt1d, cWrt_orig)
                    expect_identical(cWrt2d, cWrt_orig)
                    expect_identical(cWrt2d11, cWrt_orig)
                }

                ## Also, should we take otherNodes and break into those that are in calcNodes and those not?
                ## Those in not in calcNodes should never be changed. So maybe nothing to check.
                
                ## model state - when order 0 is included, logProb and determistic nodes should be updated; otherwise not
                expect_identical(rLogProb01, rLogProb_new)
                expect_identical(rVals01, rVals_new)
                expect_identical(rLogProb012, rLogProb_new)
                expect_identical(rVals012, rVals_new)
                if(doAllUncHessian) {
                    expect_identical(rLogProb02, rLogProb_new)
                    expect_identical(rVals02, rVals_new)
                    expect_identical(rLogProb12, rLogProb_orig)
                    expect_identical(rVals12, rVals_orig)
                }
                
                if(checkDoubleTape) {
                    ## Double tapes here don't have order = 0 in inner tape, so model should not be updated since I do pass model into nimDerivs.
                    expect_identical(rLogProb1d, rLogProb_orig)
                    if(doAllUncHessian) 
                        expect_identical(rLogProb2d, rLogProb_orig)
                    if(checkDoubleUncHessian)
                        expect_identical(rLogProb2d11, rLogProb_orig)
                    expect_identical(rVals1d, rVals_orig)
                    if(doAllUncHessian) 
                        expect_identical(rVals2d, rVals_orig)
                    if(checkDoubleUncHessian)
                        expect_identical(rVals2d11, rVals_orig)
                }
                
                ## Not clear if next check should be expect_identical (in many cases they are identical);
                ## Check with PdV whether values from taped model could get into the compiled model.        
                
                expect_fun(cLogProb01, cLogProb_new) 
                expect_fun(cVals01, cVals_new)
                expect_fun(cLogProb012, cLogProb_new) 
                expect_fun(cVals012, cVals_new)
                expect_fun(cLogProb02, cLogProb_new) 
                expect_fun(cVals02, cVals_new)

                expect_fun(cLogProb12, cLogProb_orig)
                expect_fun(cVals12, cVals_orig)

                if(checkDoubleTape) {
                    ## Double tapes here don't have order = 0 in inner tape, so model should not be updated.
                    expect_fun(cLogProb1d, cLogProb_orig)
                    expect_fun(cVals1d, cVals_orig)
                    expect_fun(cLogProb2d, cLogProb_orig)
                    expect_fun(cVals2d, cVals_orig)
                    expect_fun(cLogProb2d11, cLogProb_orig)
                    expect_fun(cVals2d11, cVals_orig)
                }
            }
        }
    })
}


## Makes random vectors of wrt elements, following James Duncan's code
make_wrt <- function(argTypes, n_random = 10, allCombinations = FALSE) {
    ## always include each arg on its own, and all combinations of the args
    ## Note that for models with a large number of variables this might turn out to be too much.
    wrts <- as.list(names(argTypes))
    if(allCombinations) {
        if (length(argTypes) > 1)
            for (m in 2:length(argTypes)) {
                this_combn <- combn(names(argTypes), m)
                wrts <- c(
                    wrts,
                    unlist(apply(this_combn, 2, list), recursive = FALSE)
                )
            }
    }
    
  while (n_random > 0) {
    n_random  <- n_random - 1
    n <- sample(1:length(argTypes), 1) # how many of the args to use?
    ## grab a random subset of the args of length n
    args <- sample(argTypes, n)
    ## may repeat an arg up to 2 times
    reps <- sample(1:2, length(args), replace = TRUE)
    argSymbols <- lapply(args, nimble:::argType2symbol)
    this_wrt <- c()
    for (i in 1:length(args)) {
      while (reps[i] > 0) {
        reps[i] <- reps[i] - 1
        ## coin flip determines whether to index vectors/matrices
        use_indexing <- sample(c(TRUE, FALSE), 1)
        if (use_indexing && argSymbols[[i]]$nDim > 0) {
          rand_row <- sample(1:argSymbols[[i]]$size[1], size = 1)
          ## another coin flip determines whether to use : in indexing or not
          use_colon <- sample(c(TRUE, FALSE), 1)
          if (use_colon && rand_row < argSymbols[[i]]$size[1]) {
            end_row <- rand_row +
              sample(1:(argSymbols[[i]]$size[1] - rand_row), size = 1)
            rand_row <- paste0(rand_row, ':', end_row)
          }
          index <- rand_row
          if (argSymbols[[i]]$nDim == 2) {
            rand_col <- sample(1:argSymbols[[i]]$size[2], size = 1)
            ## one more coin flip to subscript second dimension
            use_colon_again <- sample(c(TRUE, FALSE), 1)
            if (use_colon_again && rand_col < argSymbols[[i]]$size[2]) {
              end_col <- rand_col +
                sample(1:(argSymbols[[i]]$size[2] - rand_col), size = 1)
              rand_col <- paste0(rand_col, ':', end_col)
            }
            index <- paste0(index, ',', rand_col)
          }
          this_wrt <- c(this_wrt, paste0(names(args)[i], '[', index, ']'))
        }
        ## if first coin flip was FALSE, just
        ## use the arg name without indexing
        else this_wrt <- c(this_wrt, names(args)[i])
      }
    }
    if (!is.null(this_wrt)) wrts <- c(wrts, list(unique(this_wrt)))
  }
  wrts <- unique(wrts)
}


makeADDistributionTestList <- function(distnList){
  argsList <- lapply(distnList$args, function(x){
    return(x)
  })
  ansList <- list(args = argsList,
                  expr = substitute(out <- nimDerivs(METHODEXPR, wrt = WRT, order = c(0,1,2)),
                                    list(METHODEXPR = as.call(c(list(quote(method1)),
                                                               lapply(names(distnList$args),
                                                                      function(x){return(parse(text = x)[[1]])}))),
                                         WRT = if(is.null(distnList$WRT)) names(distnList$args) else distnList$WRT
                                    )),
                  outputType = quote(ADNimbleList())
  )
  return(ansList)
}

makeADDistributionMethodTestList <- function(distnList){
  argsList <- lapply(distnList$args, function(x){
    return(x)
  })
  argsValsList <- list()
  for(iArg in seq_along(distnList$args)){
    argsValsList[[names(distnList$args)[iArg]]] <- parse(text = names(distnList$args)[iArg])[[1]]
  }
  ansList <- list(args = argsList,
                  expr = substitute({out <- numeric(2);
                                     out[1] <- DISTNEXPR;
                                     out[2] <- LOGDISTNEXPR;},
                                    list(DISTNEXPR = as.call(c(list(parse(text = distnList$distnName)[[1]]),
                                                               argsValsList,
                                                               list(log = FALSE))),
                                         LOGDISTNEXPR = as.call(c(list(parse(text = distnList$distnName)[[1]]),
                                                               argsValsList,
                                                               list(log = TRUE)))
                                         
                                    )),
                  outputType = quote(double(1, 2))
  )
  return(ansList)
}

testADDistribution <- function(ADfunGen, argsList, name, debug = FALSE){
    ADfun <- ADfunGen()
    CADfun <- compileNimble(ADfun)
    for(iArg in seq_along(argsList)){
      iOrdersToCheck <- argsList[[iArg]][['ordersToCheck']]
      if(is.null(iOrdersToCheck)) iOrdersToCheck <- 0:2 ## check all orders if not specified
      else argsList[[iArg]][['ordersToCheck']] <- NULL
      RfunCallList <- c(list(quote(ADfun$run)), argsList[[iArg]])
      CfunCallList <- c(list(quote(CADfun$run)), argsList[[iArg]])
      RderivsList <- eval(as.call(RfunCallList))
      CderivsList <- eval(as.call(CfunCallList))
      argValsText <- paste(sapply(names(argsList[[iArg]]), 
                          function(x){return(paste(x, " = ",
                          paste(argsList[[iArg]][[x]], collapse = ', ')))}), collapse = ', ')
      if(is.logical(debug) && debug == TRUE) browser()
      else if(is.numeric(debug) && debug == iArg) browser()
      if(0 %in% iOrdersToCheck)
        expect_equal(RderivsList$value, CderivsList$value, tolerance = 1e-6, 
                     info = paste("Values of", name , "not equal for arguments: ",
                                  argValsText, '.'))
      else
        print(paste("Skipping check of R and C++ `value` equality for ",
                    name, " with arguments: ", argValsText ))
      if(1 %in% iOrdersToCheck)
        expect_equal(RderivsList$jacobian, CderivsList$jacobian, tolerance = 1e-2,
                     info = paste("Jacobians of", name , "not equal for arguments: ",
                                  argValsText, '.'))
      else
        print(paste("Skipping check of R and C++ `jacobian` equality for ",
                    name, " with arguments: ", argValsText ))
      if(2 %in% iOrdersToCheck)
        expect_equal(RderivsList$hessian, CderivsList$hessian, tolerance = 1e-2,
                     info = paste("Hessians of", name , "not equal for arguments: ",
                                  argValsText, '.'))
      else
        print(paste("Skipping check of R and C++ `hessian` equality for ",
                    name, " with arguments: ", argValsText ))
    }
    nimble:::clearCompiled(CADfun)
}

expandNames <- function(var, ...) {
    tmp <- as.matrix(expand.grid(...))
    indChars <- apply(tmp, 1, paste0, collapse=', ')
    ## Use of sort here only puts names in same order as NIMBLE if have no double-digit indexes.
    sort(paste0(var, "[", indChars, "]"))
}

test_dynamic_indexing_model <- function(param) {
    test_that(param$case, test_dynamic_indexing_model_internal(param))
    invisible(NULL)
}
                
test_dynamic_indexing_model_internal <- function(param) {
        if(!is.null(param$expectError) && param$expectError) {
            expect_error(m <- nimbleModel(param$code, dimensions = param$dims, inits = param$inits, data = param$data, constants = param$constants), param$expectErrorMsg, info = "expected error not generated")
        } else {
            m <- nimbleModel(param$code, dimensions = param$dims, inits = param$inits, data = param$data, constants = param$constants)
            expect_true(inherits(m, 'modelBaseClass'), info = "problem creating model")
            for(i in seq_along(param$expectedDeps)) 
                expect_identical(m$getDependencies(param$expectedDeps[[i]]$parent, stochOnly = TRUE),
                    param$expectedDeps[[i]]$result, info = paste0("dependencies don't match expected in dependency of ", param$expectedDeps[[i]]$parent))
            cm <- compileNimble(m)
            expect_true(is.Cmodel(cm), info = "compiled model object improperly formed")
            expect_identical(calculate(m), calculate(cm), info = "problem with R vs. C calculate with initial indexes")
            for(i in seq_along(param$validIndexes)) {
                for(j in seq_along(param$validIndexes[[i]]$var)) {
                    m[[param$validIndexes[[i]]$var[j]]] <- param$validIndexes[[i]]$value[j]
                    cm[[param$validIndexes[[i]]$var[j]]] <- param$validIndexes[[i]]$value[j]
                }
                expect_true(is.numeric(calculate(cm)), 1, info = paste0("problem with C calculate with valid indexes, case: ", i))
                expect_true(is.numeric(calculateDiff(cm)), info = paste0("problem with C calculateDiff with valid indexes, case: ", i))              
                expect_identical(calculate(m), calculate(cm), info = paste0("problem with R vs. C calculate with valid indexes, case: ", i))
                deps <- m$getDependencies(param$invalidIndexes[[i]]$var, self = FALSE)
                expect_true(is.null(simulate(cm, deps, includeData = TRUE)), info = paste0("problem with C simulate with valid indexes, case: ", i))
                ## Reset values so can models have same values for comparisons in later iterations.
                cm$setInits(param$inits)
                cm$setData(param$data)
                cm$calculate()
            }
            for(i in seq_along(param$invalidIndexes)) {
                for(j in seq_along(param$invalidIndexes[[i]]$var)) {
                    m[[param$invalidIndexes[[i]]$var[j]]] <- param$invalidIndexes[[i]]$value[j]
                    
                    cm[[param$invalidIndexes[[i]]$var[j]]] <- param$invalidIndexes[[i]]$value[j]
                }
                expect_output(out <- calculate(m), "Dynamic index out of bounds", info = paste0("problem with lack of warning in R calculate with non-NA invalid indexes, case: ", i))
                expect_equal(out, NaN, info = paste0("problem with lack of NaN in R calculate with non-NA invalid indexes, case: ", i))
                expect_output(out <- calculate(cm), "Dynamic index out of bounds", info = paste0("problem with lack of warning in C calculate with invalid indexes, case: ", i))
                expect_equal(out, NaN, info = paste0("problem with lack of NaN in C calculate with invalid indexes, case: ", i))
                expect_output(out <- calculateDiff(cm), "Dynamic index out of bounds", info = paste0("problem with lack of warning in C calculateDiff with invalid indexes, case: ", i))
                expect_equal(out, NaN, info = paste0("problem with lack of NaN in C calculateDiff with invalid indexes, case: ", i))
                deps <- m$getDependencies(param$invalidIndexes[[i]]$var, self = FALSE)
                expect_output(simulate(cm, deps, includeData = TRUE), "Dynamic index out of bounds", info = paste0("problem with lack of warning in C simulate with invalid indexes, case: ", i))
                expect_true(sum(is.nan(values(cm, deps))) >= 1, info = paste0("problem with lack of NaN in C simulate with invalid indexes, case: ", i))
            }
            if(.Platform$OS.type != "windows") {
                nimble:::clearCompiled(m)
            }
        }
    invisible(NULL)
}


## utilities for saving test output to a reference file
## and making the test a comparison of the file
clearOldOutput <- function(filename) {
    if(file.exists(filename)) file.remove(filename)
}

appendOutput <- function(filename, case, caseName, casePrefix = "") {
    outputConnection <- file(filename, open = 'at')
    writeLines(caseName, con = outputConnection)
    outputAns <- lapply(case, function(x) writeLines(paste0(casePrefix, paste(x, collapse = " ")), con = outputConnection))
    close(outputConnection)
}

writeOutput <- function(cases, filename) {
    clearOldOutput(filename)
    for(i in seq_along(cases)) appendOutput(filename, cases[[i]], names(cases)[i], casePrefix = paste0(i,": "))
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

stripTestsPassedMessage <- function(lines) {
    stripLines <- grep("^Test passed", lines)
    if(length(stripLines))
       return(lines[-stripLines]) else return(lines)
}
    

compareFilesByLine <- function(trialResults, correctResults, main = "") {
    trialResults <- stripTestsPassedMessage(stripTestPlacementWarning(trialResults))
    correctResults <- stripTestPlacementWarning(correctResults)
    expect_equal(length(trialResults), length(correctResults))
    
    linesToTest <- min(length(trialResults), length(correctResults))
    mapply(function(lineno, trialLine, correctLine) {
        expect_identical(trialLine, correctLine)
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

## Create a nimbleFunction parametrization to be passed to gen_runFunCore().
##
## op:        An operator string.
## argTypes:  A character vector of argTypes (e.g. "double(0)").
##            If this is a named vector, then the names will be
##            interpreted as the formals to the constructed op call.
## more_args: A named list, e.g. list(log = 1), that will be added
##            as a formal to the output expr but not part of args.
##
## output_code: quoted code such as quote(Y^2) to be used with Y
##              substituted as the operation result.  E.g., if the
##              operation is `arg1 + arg2`, this would replace it with
##             (arg1 + arg2)^2.
##
## inner_codes: list of quoted code such as quote(X^2) to be used with
##              X substituted as the corresponding argument. The list
##              must be as long as the number of arugments, with NULL
##              entries indicating no substitution.  E.g. with
##              list(NULL, quote(X^2)), `arg1 + arg2` would be changed
##              to `art1 + arg2^2`.
make_op_param <- function(op, argTypes, more_args = NULL,
                          outer_code = NULL, inner_codes = NULL) {
  arg_names <- names(argTypes)

  if (is.null(arg_names)) {
    arg_names <- paste0('arg', 1:length(argTypes))
    op_args <- lapply(arg_names, as.name)
  } else {
    op_args <- sapply(arg_names, as.name, simplify = FALSE)
  }

  args_string <- paste0(arg_names, ' = ', argTypes, collapse = ' ')
  name <- paste(op, args_string)

  # Add inner funs around arguments
  if(!is.null(inner_codes)) {
    for(i in seq_along(op_args)) {
      if(!is.null(inner_codes[[i]])) {
        op_args[[i]] <- eval(substitute(
          substitute(INNER_CODE,
                     list(X = op_args[[i]])),
          list(INNER_CODE = inner_codes[[i]]))
        )
      }
    }
  }
  
  this_call <- as.call(c(
      substitute(FOO, list(FOO = as.name(op))),
      op_args, more_args
  ))

  if(is.null(outer_code)) {
      expr <- substitute(
          out <- THIS_CALL,
          list(THIS_CALL = this_call)
      )
  } else {
    expr <- eval(substitute(
      substitute(
        out <- OUTER_CODE,
        list(Y = this_call)),
      list(OUTER_CODE = outer_code))
      )
  }

  argTypesList <- as.list(argTypes)
  names(argTypesList) <- arg_names
  argTypesList <- lapply(argTypesList, function(arg) {
    parse(text = arg)[[1]]
  })

  list(
    name = name,
    expr = expr,
    args = argTypesList,
    outputType = parse(text = return_type_string(op, argTypes))[[1]]
  )
}

## Takes an operator and its input types as a character vector and
## creates a string representing the returnType for the operation.
##
## op:       An operator string
## argTypes: A character vector of argTypes (e.g. "double(0)".
##
return_type_string <- function(op, argTypes) {

  ## multivariate distributions ops.  These do not support recycling rule behavior, so the return type is always double(0)
  mvdist_ops <- names(nimble:::sizeCalls)[nimble:::sizeCalls == 'sizeScalarRecurseAllowMaps'] 
  if(op %in% mvdist_ops)
    return("double(0)")
  
  ## see ops handled by eigenize_recyclingRuleFunction in genCpp_eigenization.R
  recycling_rule_ops <- c(
    nimble:::scalar_distribution_dFuns,
    nimble:::scalar_distribution_pFuns,
    nimble:::scalar_distribution_qFuns,
    nimble:::scalar_distribution_rFuns,
    paste0(c('d', 'q', 'p', 'r'), 't'),
    paste0(c('d', 'q', 'p', 'r'), 'exp'),
    'bessel_k'
  )

  returnTypeCode <- nimble:::returnTypeHandling[[op]]

  if (is.null(returnTypeCode))
    if (!op %in% recycling_rule_ops)
      return(argTypes[1])
    else returnTypeCode <- 1

  scalarTypeString <- switch(
    returnTypeCode,
    'double', ## 1
    'integer', ## 2
    'logical'  ## 3
  )

  args <- lapply(
    argTypes, function(argType)
      nimble:::argType2symbol(parse(text = argType)[[1]])
  )

  if (is.null(scalarTypeString)) ## returnTypeCode is 4 or 5
    scalarTypeString <-
      if (length(argTypes) == 1)
        if (returnTypeCode == 5 && args[[1]]$type == 'logical') 'integer'
        else args[[1]]$type
      else if (length(argTypes) == 2) {
        aot <- nimble:::arithmeticOutputType(args[[1]]$type, args[[2]]$type)
        if (returnTypeCode == 5 && aot == 'logical') 'integer'
        else aot
      }

  reductionOperators <- c(
    nimble:::reductionUnaryOperators,
    nimble:::matrixSquareReductionOperators,
    nimble:::reductionBinaryOperatorsEither,
    'dmulti'
  )

  nDim <- if (op %in% reductionOperators) 0
          else max(sapply(args, `[[`, 'nDim'))

  if (nDim > 2)
    stop(
      'Testing does not currently support args with nDim > 2',
      call. = FALSE
    )

  sizes <- if (nDim == 0) 1
           else if (length(argTypes) == 1) args[[1]]$size
           else if (op %in% nimble:::matrixMultOperators)  {
             if (!length(argTypes) == 2)
               stop(
                 paste0(
                   'matrixMultOperators should only have 2 args but got ',
                   length(argTypes)
                 ), call. = FALSE
               )
             c(args[[1]]$size[1], args[[2]]$size[2])
           } else if (nDim == 2) {
             ## one arg is a matrix but this is not matrix multiplication
             ## so assume that the first arg with nDim > 1
             has_right_nDim <- sapply(args, function(arg) arg$nDim == nDim)
             args[has_right_nDim][[1]]$size
           } else {
             ## nDim is 1 so either recycling rule or simple vector operator
             max((sapply(args, `[[`, 'size')))
           }

    size_string <- if(all(is.na(sizes)))
                       ''
                   else if (nDim > 0)
                       paste0(', c(', paste(sizes, collapse = ', '), ')')
                   else ''

  return(paste0(scalarTypeString, '(', nDim, size_string, ')'))
}

## Takes an argSymbol and if argSymbol$size is NA adds default sizes.
add_missing_size <- function(argSymbol, vector_size = 3, matrix_size = c(3, 4)) {
  if (any(is.na(argSymbol$size))) {
    if (argSymbol$nDim == 1)
      argSymbol$size <- vector_size
    else if (argSymbol$nDim == 2)
      argSymbol$size <- matrix_size
  }
  invisible(argSymbol)
}

arg_type_2_input <- function(argType, input_gen_fun = NULL, size = NULL, return_function = FALSE) {
  argSymbol <- add_missing_size(
    nimble:::argType2symbol(argType)
  )
  type <- argSymbol$type
  nDim <- argSymbol$nDim
  if(is.null(size))
    size <- argSymbol$size
  
  if (is.null(input_gen_fun))
    input_gen_fun <- switch(
      type,
      "double"  = function(arg_size) rnorm(prod(arg_size)),
      "integer" = function(arg_size) rgeom(prod(arg_size), 0.5),
      "logical" = function(arg_size)
        sample(c(TRUE, FALSE), prod(arg_size), replace = TRUE)
    )
  ans <- function() {
    arg <- switch(
      nDim + 1,
      input_gen_fun(1), ## nDim is 0
      input_gen_fun(prod(size)), ## nDim is 1
      matrix(input_gen_fun(prod(size)), nrow = size[1], ncol = size[2]), ## nDim is 2
      array(input_gen_fun(prod(size)), dim = size) ## nDim is 3
    )
    if(is.null(arg))
      stop('Something went wrong while making test input.', call.=FALSE)
    arg
  }
  if(return_function)
    return(ans)
  ans()
}

modify_on_match <- function(x, pattern, key, value, env = parent.frame(), ...) {
  ## Modify any elements of a named list that match pattern.
  ##
  ## @param x A named list of lists.
  ## @param pattern A regex pattern to compare with `names(x)`.
  ## @param key The key to modify in any lists whose names match `pattern`.
  ## @param value The new value for `key`.
  ## @param env The environment in which to modify `x`.
  ## @param ... Additional arguments for `grepl`.
  for (name in names(x)) {
    if (grepl(pattern, name, ...)) {
      eval(substitute(x[[name]][[key]] <- value), env)
    }
  }
}

nim_all_equal <- function(x, y, tolerance = .Machine$double.eps^0.5, abs_threshold = 0, verbose = FALSE, info = "") {
  xlab <- deparse1(substitute(x))
  ylab <- deparse1(substitute(y))
  
  denom <- abs(y)
  ## Use absolute tolerance for sufficiently small values.
  ## This is necessary for y values exactly zero.
  ## In some cases (such as derivatives very near zero and affected by floating point errors)
  ## we also need to use absolute tolerance.
  denom[denom <= abs_threshold] <- 1
  rel_diff <- abs((x-y)/denom)
  result <- rel_diff < tolerance
  all_result <- all(result)
  if(verbose) {
    if(!all_result) {
      ord <- order(rel_diff, decreasing = TRUE)
      wh <- 1:min(length(rel_diff), 5)
      report <- cbind(x[ord[wh]], y[ord[wh]], rel_diff[ord[wh]])
      report <- report[report[,3] > tolerance, ]
      cat("\n******************\n")
      cat("Detected some values out of relative tolerance ", info, ": ", xlab, " ", ylab, ".\n")
      print(report)
      cat("******************\n")
    } 
  }
  all_result
}

nim_expect_equal <- function(x, y, tolerance = .Machine$double.eps^0.5)
  expect_true(nim_all_equal(x, y, tolerance, verbose = TRUE))
