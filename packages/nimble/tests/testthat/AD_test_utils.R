source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# An environment for storing some defaults for AD testing.
# Having this allows one to change these values globally.
# Each entry is tolerances for 0th, 1st, and 2nd order followed by an
# absolute threshold below which the tolerances will be used as absolute, not
# relative. See nim_all_equal.  It is applied to the second of the two valued
# to be compared.
ADtestEnv <- new.env()
resetTols <- function() {
  ADtestEnv$RRrelTol <<- c(1e-8, 1e-3, 1e-2, 1e-10)
  ADtestEnv$RCrelTol <<- c(1e-8, 1e-3, 1e-2, 1e-10)
  ADtestEnv$CCrelTol <<- c(rep(sqrt(.Machine$double.eps), 3), 1e-10)
}
resetTols()

#####################
## knownFailure utils
#####################

## test_param_name: an AD test parameterization name as produced by
##                  make_op_param()
##
## returns: a length 2 character vector with the op and args
##
split_test_param_name <- function(test_param_name) {
  first_space <- regexpr(' ', test_param_name)
  op <- substr(test_param_name, 1, first_space - 1)
  args <- substr(test_param_name, first_space + 1, nchar(test_param_name))
  return(c(op, args))
}

is_failure <- function(test_param_name, failure_name, knownFailures = list()) {
  if (length(knownFailures) == 0) return(FALSE)
  op_and_args <- split_test_param_name(test_param_name)
  op <- op_and_args[1]
  args <- op_and_args[2]
  if (!is.null(knownFailures[[op]])) {
    return(
      isTRUE(knownFailures[[op]][[args]][[failure_name]]) ||
      isTRUE(knownFailures[[op]][['*']][[failure_name]])
    )
  }
  return(FALSE)
}

## test_param_name: an AD test parameterization name as produced by
##                  make_op_param()
## knownFailures:   a list of known test failures, e.g. in the format of
##                  AD_knownFailures
## 
## returns: TRUE if the test param is known to lead to a compilation failure and
##          FALSE otherwise
## 
is_compilation_failure <- function(test_param_name, knownFailures = list()) {
  return(is_failure(test_param_name, 'compilation', knownFailures))
}

is_segfault_failure <- function(test_param_name, knownFailures = list()) {
  return(is_failure(test_param_name, 'segfault', knownFailures))
}

is_method_failure <- function(test_param_name, method_name,
                              output_name = c('value', 'jacobian', 'hessian'),
                              knownFailures = list()) {
  if (length(knownFailures) == 0) return(FALSE)
  output_name <- match.arg(output_name)
  op_and_args <- split_test_param_name(test_param_name)
  op <- op_and_args[1]
  args <- op_and_args[2]
  if (!is.null(knownFailures[[op]])) {
    this_arg_fails <- all_args_fail <- FALSE
    if (!is.null(knownFailures[[op]][[args]]))
      this_arg_fails <- identical(
        knownFailures[[op]][[args]][[output_name]],
        method_name
      )
    if (!is.null(knownFailures[[op]][['*']]))
      all_args_fail <- identical(
        knownFailures[[op]][['*']][[output_name]],
        method_name
      )
    return(this_arg_fails || all_args_fail)
  }
  return(FALSE)
}

all_equal_ignore_zeros <- function(v1, v2, ...) {
  if(length(v1) != length(v2)) return(FALSE)
  zv1 <- as.logical(v1 == 0) ## as.logical will strip any names
  zv2 <- as.logical(v2 == 0)
  nzboth <- !(zv1 & zv2)
  if(sum(nzboth)==0) return(FALSE) ## This is used to validate expected discrepancies except when they are zero, so if all are zero, we default to that there were discrepancies, which allows the test to pass.  
  all.equal(v1[nzboth], v2[nzboth], ...)
}

########################
## Main AD testing utils
########################

test_AD2_oneCall <- function(Robj, Cobj,
                             recordArgs,
                             testArgs,
                             order,
                             wrt = NULL,
                             check.equality = TRUE,
                             # In normal usage, the tols are set from test_AD2 and
                             # these defaults are not used
                             RRrelTol = c(1e-8, 1e-3, 1e-2, 1e-14),
                             RCrelTol = c(1e-15, 1e-8, 1e-3, 1e-14),
                             CCrelTol = c(rep(sqrt(.Machine$double.eps), 3), 1e-14),
                             Rmodel = NULL, Cmodel = NULL,
                             recordInits = NULL, testInits = NULL,
                             nodesToChange = NULL) {
  resRecord <- list()
  resTest <- list()

  useRmodel <- !is.null(Rmodel)
  useCmodel <- !is.null(Cmodel)

  RRabsThresh <- 0
  RCabsThresh <- 0
  CCabsThresh <- 0
  if(length(RRrelTol) > 3) RRabsThresh <- RRrelTol[4]
  if(length(RCrelTol) > 3) RCabsThresh <- RCrelTol[4]
  if(length(CCrelTol) > 3) CCabsThresh <- CCrelTol[4]

  if(is.null(names(RRrelTol))) names(RRrelTol)[1:3] <- c("0", "1", "2")
  if(is.null(names(RCrelTol))) names(RCrelTol)[1:3] <- c("0", "1", "2")
  if(is.null(names(CCrelTol))) names(CCrelTol)[1:3] <- c("0", "1", "2")

  setModelInits <- function(m, inits) {
    for(v in names(inits)) {
      m[[v]] <- inits[[v]]
    }
    m$calculate()
  }

  changeModelInits <- function(m, inits, nodesToChange) {
    for(node in nodesToChange) { # This will not include constant nodes and can include like 'x[2:3]'
      if(all.vars(parse(text = node))[1] %in% names(inits))
        eval(parse(text = paste0("m$", node, " <- inits$", node)))
    }
    m$calculate()
  }
  
  do_one_call <- function(fun, argList) {
    eval( as.call( c(fun, argList) ) )
  }
  
  do_one_set <- function(method, order, metaLevel = 0, name,
                         doR, doC, tapingLevel = c(1, 0, 2), fixedOrder = FALSE) {
    # tapingLevel: 0 for the (run) function itself, 1 for a tape of run, 2 for a tape of a tape of run
    # metaLevel: 0 if order 0 from the tape is order 0 result.  1 if order 0 gives order 1 result due to double taping.  2 if order 0 gives order 2 result due to double taping.
    # fixedOrder: FALSE is return object is result of nimDerivs call.
    #             0, 1, or 2 if return object is the value, jacobian or hessian directly.
    if(tapingLevel == 0) order <- 0
    if(isFALSE(fixedOrder))
      metaOrder <- order[order >= metaLevel] - metaLevel # e.g. for metaLevel 1, count 1 as 0
    else
      metaOrder <- fixedOrder
    # Follow construction is because Robj[[method]] doesn't work until
    # Robj$run (if method == "run") has been done.  It's a flaw in reference classes
    Rfun <- eval( parse( text = paste0("Robj$", method), keep.source = FALSE)[[1]])
    Cfun <- eval( parse( text = paste0("Cobj$", method), keep.source = FALSE)[[1]])
    extraArgs <- list()
    if(tapingLevel >= 1) {
      extraArgs$wrt <- wrt_all
      # If there is a tape and no fixed order, the order arg is needed
      if(isFALSE(fixedOrder))
        extraArgs$order <- metaOrder
      extraArgs$reset <- TRUE
    }
    if(tapingLevel == 2) {
      extraArgs$innerWrt = wrt_all
    }
    argList <- c(recordArgs, extraArgs)
    if(doR) {
      if(useRmodel) setModelInits(Rmodel, recordInits)
      RansRecord <- do_one_call(Rfun, argList)
    }
    if(doC) {
      if(useCmodel) setModelInits(Cmodel, recordInits)
      CansRecord <- do_one_call(Cfun, argList)
    }
    if(tapingLevel >= 1)
      extraArgs$reset <- FALSE

    argList <- c(testArgs, extraArgs)
    if(doR) {
      if(useRmodel) changeModelInits(Rmodel, testInits, nodesToChange)
      RansTest <- do_one_call(Rfun, argList)
    }
    if(doC) {
      if(useCmodel) changeModelInits(Cmodel, testInits, nodesToChange)
      CansTest <- do_one_call(Cfun, argList)
    }

    pieces <- c("0"="value", "1"="jacobian", "2"="hessian")
    rName <- paste0("R", name)
    cName <- paste0("C", name)

    if(!isFALSE(fixedOrder)) {
      oChar <- as.character(fixedOrder) # Assume fixedOrder given directly (no need to adjust for metaOrder)
      if(doR) {
        resRecord[[oChar]][[rName]] <<- RansRecord
        resTest[[oChar]][[rName]] <<- RansTest
      }
      if(doC) {
        resRecord[[oChar]][[cName]] <<- CansRecord
        resTest[[oChar]][[cName]] <<- CansTest
      }
    } else {      
      for(o in as.character(metaOrder)) {
        oChar <- as.character(as.numeric(o) + metaLevel)
        if((metaLevel == 1) & o == '1') {
          argsLength <- length(wrt_all)
          reorder_jac_jac <- function(x) {
            outLength <- length(x) / (argsLength*argsLength)
            as.numeric(x) |> array(c(outLength, argsLength, argsLength)) |> aperm(c(2, 3, 1))
          }
        }
        if(doR) {
          if((metaLevel == 1) & o == '1') {            
            resRecord[[oChar]][[rName]] <<- reorder_jac_jac(RansRecord[[ pieces[[o]] ]])
            resTest[[oChar]][[rName]] <<- reorder_jac_jac(RansRecord[[ pieces[[o]] ]])
          } else {
            resRecord[[oChar]][[rName]] <<- RansRecord[[ pieces[[o]] ]]
            resTest[[oChar]][[rName]] <<- RansTest[[ pieces[[o]] ]]
          }
        }
        if(doC) {
          if((metaLevel == 1) & o == '1') {            
            resRecord[[oChar]][[cName]] <<- reorder_jac_jac(CansRecord[[ pieces[[o]] ]])
            resTest[[oChar]][[cName]] <<- reorder_jac_jac(CansTest[[ pieces[[o]] ]])
          } else {
            resRecord[[oChar]][[cName]] <<- CansRecord[[ pieces[[o]] ]]
            resTest[[oChar]][[cName]] <<- CansTest[[ pieces[[o]] ]]
          }
        }
      }
    }
  }
  
  len_record <- sum(unlist(lapply(recordArgs, length)))
  len_test <- sum(unlist(lapply(testArgs, length)))
  if(len_record != len_test) {
    stop('lengths of recordArgs and testArgs must match.')
  }

  if(is.null(wrt))
    wrt_all <- 1:len_record
  else
    wrt_all <- wrt
  
  do_one_set("run", name = "run",
             doR = TRUE, doC = TRUE, tapingLevel = 0, fixedOrder = 0)
  do_one_set("derivsRun", order = order, metaLevel = 0, name = "derivsRun",
             doR = TRUE, doC = TRUE, tapingLevel = 1, fixedOrder = FALSE)
  do_one_set("value", name = "value",
             doR = FALSE, doC = TRUE, tapingLevel = 1, fixedOrder = 0)
  do_one_set("jac", name = "jacR1",
             doR = FALSE, doC = TRUE, tapingLevel = 1, fixedOrder = 1)
  do_one_set("jac2", name = "jacF1",
             doR = FALSE, doC = TRUE, tapingLevel = 1, fixedOrder = 1)
  do_one_set("hess", name = "hess",
             doR = FALSE, doC = TRUE, tapingLevel = 1, fixedOrder = 2)
  do_one_set("derivsValue", order = 0:2, name = "derivsValue",
             doR = FALSE, doC = TRUE, tapingLevel = 2, fixedOrder = FALSE)
  do_one_set("derivsJac", order = 1:2, name = "derivsJacR1", metaLevel = 1,
             doR = FALSE, doC = TRUE, tapingLevel = 2, fixedOrder = FALSE)
  do_one_set("derivsJac2", order = 1:2, name = "derivsJacF1", metaLevel = 1,
             doR = FALSE, doC = TRUE, tapingLevel = 2, fixedOrder = FALSE)
  do_one_set("derivsHess", order = 2, name = "derivsHess", metaLevel = 2,
             doR = FALSE, doC = TRUE, tapingLevel = 2, fixedOrder = FALSE)

  all_equal_list <- function(first, others, tol,
                             abs_threshold, info = "") {
    pass <- TRUE
    for(i in seq_along(others)) {
      pass <- pass && nim_all_equal(as.numeric(first), as.numeric(others[[i]]),
                                    tol, abs_threshold = abs_threshold,
                                    verbose = TRUE, info = info)
    }
    pass
  }

  if(check.equality) {
    pass <- TRUE
    for(res in list(resRecord, resTest))
      for(o in as.character(order)) {
        ansSet <- res[[o]]
        splitAnsSet <- split(ansSet,
                             c("C", "R")[as.integer(grepl("^R", names(ansSet)))+1])
        RansSet <- splitAnsSet[["R"]]
        CansSet <- splitAnsSet[["C"]]
        if(length(RansSet) > 1) {
          pass <- pass && all_equal_list(RansSet[[1]], RansSet[-1], tol = RRrelTol[[o]],
                                         abs_threshold = RRabsThresh,
                                         info = paste0("(RR order ", o,")"))
          if(!pass) {
            cat(paste('Some R-to-R derivatives do not match for order', o, '.\n'))
            browser()
          }
        }
        if(length(RansSet) > 0 && length(CansSet) > 0) {
          pass <- pass && all_equal_list(RansSet[[1]], CansSet, tol = RCrelTol[[o]],
                                         abs_threshold = RCabsThresh,
                                         info = paste0("(RC order ", o,")"))
          if(!pass) {
            cat(paste('Some C-to-R derivatives to not match for order', o, '.\n'))
            browser()
          }
        }
        if(length(CansSet) > 1) {
          pass <- pass && all_equal_list(CansSet[[1]], CansSet[-1], CCrelTol[[o]],
                                         abs_threshold = CCabsThresh,
                                         info = paste0("(CC order ", o, ")"))
          if(!pass) {
            cat(paste('Some C-to-C derivatives to not match for order', o, '.\n'))
            browser()
          }
        }
      }
    expect_true(pass)
  }
  list(resRecord = resRecord, resTest = resTest)
}

test_AD2 <- function(param, dir = file.path(tempdir(), "nimble_generatedCode"),
                     control = list(), verbose = nimbleOptions('verbose'),
                     catch_failures = FALSE, seed = NULL,
                     nimbleProject_name = '', return_compiled_nf = FALSE,
                     knownFailures = list()) {
  if (!is.null(param$debug) && param$debug) browser()
  if (verbose) cat(paste0("### Testing ", param$name, "\n"))

  ## by default, reset the seed for every test
  if(is.numeric(param[['seed']])) set.seed(param[['seed']])
     else if(is.numeric(seed)) set.seed(seed)

  nf <- nimbleFunction(
    setup = function() {},
    run = param$run,
    methods = param$methods,
    buildDerivs = param$buildDerivs
  )
  Robj <- nf()
  temporarilyAssignInGlobalEnv(Robj)

  if(!is.null(param$inputs)) {
    inputsRecord <- param$inputs[[1]]
    inputsTest <- param$inputs[[2]]
  } else {
    ## input_gen_funs might be generalized.
    opParam <- param$opParam
    if (is.null(param$input_gen_funs) || is.null(names(param$input_gen_funs)))
      if (length(param$input_gen_funs) <= 1) {
        inputsRecord <- lapply(opParam$args, arg_type_2_input,
                               input_gen_fun = param$input_gen_funs, size = param$size)
        inputsTest <- lapply(opParam$args, arg_type_2_input,
                             input_gen_fun = param$input_gen_funs, size = param$size)
      }
    else
      stop(
        'input_gen_funs of length greater than 1 must have names',
        call. = FALSE
      )
    else {
      inputsRecord <- sapply(
        names(opParam$args),
        function(name)
          arg_type_2_input(opParam$args[[name]],  input_gen_fun = param$input_gen_funs[[name]],
                           size = param$size[[name]]),
        simplify = FALSE
      )
      inputsTest <- sapply(
        names(opParam$args),
        function(name)
          arg_type_2_input(opParam$args[[name]],  input_gen_fun = param$input_gen_funs[[name]],
                           size = param$size[[name]]),
        simplify = FALSE
      )
    }
  }
  ##
  ## generate inputs that depend on the other inputs
  ##
  is_fun <- sapply(inputsRecord, is.function)
  inputsRecord[is_fun] <- lapply(
    inputsRecord[is_fun], function(fun) {
      eval(as.call(c(fun, inputsRecord[names(formals(fun))])))
    }
  )
  is_fun <- sapply(inputsTest, is.function)
  inputsTest[is_fun] <- lapply(
    inputsTest[is_fun], function(fun) {
      eval(as.call(c(fun, inputsTest[names(formals(fun))])))
    }
  )

  # Currently only whole arguments supported.
  # Hence all this does is create indices that skip over non-wrt args
  lens <- lapply(inputsRecord, length)
  if(!is.null(param$wrt_args)) {
    nextInd <- 1
    lenInds <- lapply(lens,
                      function(x) {
                        ans <- nextInd:(nextInd + x - 1); nextInd <<- nextInd + x; ans})
    wrt <- unlist(lenInds[param$wrt_args])
  } else {
    wrt <- 1:(sum(unlist(lens)))
  }
  
  if (is_segfault_failure(param$name, knownFailures)) {
    if (verbose) cat("## Skipping the rest of test before compilation",
                     "due to known segmentation fault\n")
    return(invisible(NULL))
  }

  if (!is.null(param$dir)) dir <- param$dir
  compilation_fails <- is_compilation_failure(param$name, knownFailures)
  Cobj <- param$Cobj ## user provided compiled nimbleFunction?
  if (is.null(Cobj)) {
    if (verbose) cat("## Compiling nimbleFunction \n")
    Cobj <- wrap_if_true(compilation_fails, expect_error, {
      compileNimble(
        Robj, dirName = dir, projectName = nimbleProject_name,
        control = control
      )
    }, wrap_in_try = isTRUE(catch_failures))
  }
  if (isTRUE(catch_failures) && inherits(Cobj, 'try-error')) {
    warning(
      paste0(
        'The test of ', opParam$name,
        ' failed to compile.\n'
      ),
      call. = FALSE,
      immediate. = TRUE
    )
    ## stop the test here because it didn't compile
    return(invisible(NULL))
  } else if (compilation_fails) {
      if (verbose) cat("## Compilation failed, as expected \n")
      return(invisible(NULL))
  } else {
    RRrelTol <- param$RRrelTol
    RCrelTol <- param$RCrelTol
    CCrelTol <- param$CCrelTol
    if(is.null(RRrelTol)) RRrelTol <- ADtestEnv$RRrelTol
    if(is.null(RCrelTol)) RCrelTol <- ADtestEnv$RCrelTol
    if(is.null(CCrelTol)) CCrelTol <- ADtestEnv$CCrelTol
        
    res <- try(test_AD2_oneCall(Robj = Robj, Cobj = Cobj,
                                recordArgs = inputsRecord, testArgs = inputsTest,
                                order = 0:2,
                                wrt = wrt,
                                RRrelTol = RRrelTol,
                                RCrelTol = RCrelTol,
                                CCrelTol = CCrelTol))
    if (inherits(res, 'try-error')) {
      msg <- paste(
        'Something failed in test', opParam$name, '.\n'
      )
      if (isTRUE(catch_failures)) ## continue to compilation
        warning(msg, call. = FALSE, immediate. = TRUE)
      else
        stop(msg, call. = FALSE) ## throw an error here
    }

  }
  if(!compilation_fails && getNimbleOption('useClearCompiledInADTesting')) {
    nimble:::clearCompiled(Cobj)
  }
  list(res = res)
}
## Take a test parameterization created by make_AD_test() or
## make_distribution_fun_AD_test(), generate a random input, and test for
## matching nimDerivs outputs from uncompiled and compiled versions of a
## nimbleFunction.
##
## param:              an test parameterization generated by make_AD_test() or
##                     make_distribution_fun_AD_test()
## dir:                passed to compileNimble() as the dirName argument
## control:            passed to compileNimble() as the control argument
## verbose:            if TRUE, print messages while testing
## catch_failures:     if TRUE, don't stop testing when a testthat expect_*
##                     fails
## seed:               seed to use in set.seed() before generating random inputs
## nimbleProject_name: passed to compileNimble() as the projectName argument
## return_compiled_nf: if TRUE, don't call clearCompiled() and include the
##                     compiled nimbleFunction instance in the output
##
## returns: a list with the randomly generated input and possibly the compiled
##          nimbleFunction instance, or NULL if the test has a known compilation
##          failure
test_AD <- function(param, dir = file.path(tempdir(), "nimble_generatedCode"),
                    control = list(), verbose = nimbleOptions('verbose'),
                    catch_failures = FALSE, seed = 0,
                    nimbleProject_name = '', return_compiled_nf = FALSE,
                    knownFailures = list()) {
  if (!is.null(param$debug) && param$debug) browser()
  if (verbose) cat(paste0("### Testing ", param$name, "\n"))

  ## by default, reset the seed for every test
  if (is.numeric(seed)) set.seed(seed)

  ## TODO: use nClass instead?
  nf <- nimbleFunction(
    setup = function() {},
    run = param$run,
    methods = param$methods,
    buildDerivs = param$buildDerivs
  )
  nfInst <- nf()
  temporarilyAssignInGlobalEnv(nfInst)

  ##
  ## generate inputs for the nimbleFunction methods
  ##
  opParam <- param$opParam
  if (is.null(param$input_gen_funs) || is.null(names(param$input_gen_funs)))
    if (length(param$input_gen_funs) <= 1)
      input <- lapply(opParam$args, arg_type_2_input, param$input_gen_funs)
    else
      stop(
        'input_gen_funs of length greater than 1 must have names',
        call. = FALSE
      )
  else {
    input <- sapply(
      names(opParam$args),
      function(name)
        arg_type_2_input(opParam$args[[name]], param$input_gen_funs[[name]]),
      simplify = FALSE
    )
  }
  ##
  ## generate inputs that depend on the other inputs
  ##
  is_fun <- sapply(input, is.function)
  input[is_fun] <- lapply(
    input[is_fun], function(fun) {
      eval(as.call(c(fun, input[names(formals(fun))])))
    }
  )

  if (is_segfault_failure(param$name, knownFailures)) {
    if (verbose) cat("## Skipping the rest of test before compilation",
                     "due to known segmentation fault\n")
    return(invisible(NULL))
  }

  ##
  ## compile the nimbleFunction
  ##
  if (!is.null(param$dir)) dir <- param$dir
  compilation_fails <- is_compilation_failure(param$name, knownFailures)
  CnfInst <- param$CnfInst ## user provided compiled nimbleFunction?
  if (is.null(CnfInst)) {
    if (verbose) cat("## Compiling nimbleFunction \n")
    CnfInst <- wrap_if_true(compilation_fails, expect_error, {
      compileNimble(
        nfInst, dirName = dir, projectName = nimbleProject_name,
        control = control
      )
    }, wrap_in_try = isTRUE(catch_failures))
  }
  if (isTRUE(catch_failures) && inherits(CnfInst, 'try-error')) {
    warning(
      paste0(
        'The test of ', opParam$name,
        ' failed to compile.\n', CnfInst[1]
      ),
      call. = FALSE,
      immediate. = TRUE
    )
    ## stop the test here because it didn't compile
    return(invisible(NULL))
  } else if (compilation_fails) {
    if (verbose) cat("## Compilation failed, as expected \n")
  } else {
    ##
    ## call R versions of nimbleFunction methods with generated input
    ##
    if (verbose) cat("## Calling R versions of nimbleFunction methods\n")
    Rderivs <- try(
      sapply(names(param$methods), function(method) {
        do.call(
          ## can't access nfInst[[method]] until $ has been used :(
          eval(substitute(nfInst$method, list(method = as.name(method)))), input
        )
      }, USE.NAMES = TRUE), silent = TRUE
    )
    if (inherits(Rderivs, 'try-error')) {
      msg <- paste(
        'Calling R version of test', opParam$name,
        'resulted in an error:\n', Rderivs[1]
      )
      if (isTRUE(catch_failures)) ## continue to compilation
        warning(msg, call. = FALSE, immediate. = TRUE)
      else
        stop(msg, call. = FALSE) ## throw an error here
    }

    ##
    ## call compiled nimbleFunction methods with generated input
    ##
    Cderivs <- sapply(names(param$methods), function(method) {
      do.call(
        ## same issue as with Rderivs
        eval(substitute(CnfInst$method, list(method = as.name(method)))),
        input
      )
    }, USE.NAMES = TRUE)
    if ('log' %in% names(opParam$args)) {
      input2 <- input
      input2$log <- as.numeric(!input$log)
      Rderivs2 <- try(
        sapply(names(param$methods), function(method) {
          do.call(nfInst[[method]], input2)
        }, USE.NAMES = TRUE), silent = TRUE
      )
      Cderivs2 <- sapply(names(param$methods), function(method) {
        do.call(CnfInst[[method]], input2)
      }, USE.NAMES = TRUE)
    }
    
    ##
    ## loop over test methods (each with a different wrt arg)
    ##
    ## set expect_equal tolerances
    tol1 <- if (is.null(param$tol1)) 1e-8 else param$tol1
    tol2 <- if (is.null(param$tol2)) 1e-6 else param$tol2
    tol3 <- if (is.null(param$tol3)) 1e-4 else param$tol3
    for (method_name in names(param$methods)) {
      info <- paste0(
        param$name, ', wrt = ',
        paste0(param$wrts[[method_name]], collapse = ', ')
      )
      if (verbose) {
        cat(paste0(
          "## Testing ", method_name, ': ',
          paste0(param$wrts[[method_name]], collapse = ', '),
          '\n'
        ))
      }
      ##
      ## test values
      ##
      value_test_fails <- is_method_failure(
        param$name, method_name, 'value', knownFailures
      )
      value_test <- wrap_if_true(value_test_fails, expect_failure, {
        if (verbose) cat("## Checking values\n")
        expect_equal(
          Cderivs[[method_name]]$value,
          Rderivs[[method_name]]$value,
          tolerance = tol1, info = info
        )
        if ('log' %in% names(opParam$args)) {
          if (verbose) cat("## Checking log behavior for values\n")
          expect_equal(
            Cderivs2[[method_name]]$value,
            Rderivs2[[method_name]]$value,
            tolerance = tol1, info = info
          )
          expect_false(isTRUE(all.equal(
            Rderivs[[method_name]]$value,
            Rderivs2[[method_name]]$value,
            tolerance = tol1)), info = info
          )
          expect_false(isTRUE(all.equal(
            Cderivs[[method_name]]$value,
            Cderivs2[[method_name]]$value,
            tolerance = tol1)), info = info
          )
        }
      }, wrap_in_try = isTRUE(catch_failures))
      if (isTRUE(catch_failures) && inherits(value_test, 'try-error')) {
        warning(
          paste0(
            'There was something wrong with the values of ',
            opParam$name, ' with wrt = c(',
            paste0(param$wrts[[method_name]], collapse = ', '), ').\n',
            value_test[1]
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (value_test_fails) {
        if (verbose) {
          cat(paste0(
            "## As expected, test of values failed for ", method_name, ' with wrt: ',
            paste0(param$wrts[[method_name]], collapse = ', '),
            '\n'
          ))
        }
        ## stop testing after an expected failure
        break
      }
      ##
      ## test jacobians
      ##
      jacobian_test_fails <- is_method_failure(
        param$name, method_name, 'jacobian', knownFailures
      )
      jacobian_test <- wrap_if_true(jacobian_test_fails, expect_failure, {
        if (verbose) cat("## Checking jacobians\n")
        expect_equal(
          Cderivs[[method_name]]$jacobian,
          Rderivs[[method_name]]$jacobian,
          tolerance = tol2, info = info
        )
        if ('log' %in% names(opParam$args)) {
          if (verbose) cat("## Checking log behavior for jacobians\n")
          expect_equal(
            Cderivs2[[method_name]]$jacobian,
            Rderivs2[[method_name]]$jacobian,
            tolerance = tol2, info = info
          )
          expect_false(isTRUE(all_equal_ignore_zeros(
            Rderivs[[method_name]]$jacobian,
            Rderivs2[[method_name]]$jacobian,
            tolerance = tol2)), info = info
          )
          expect_false(isTRUE(all_equal_ignore_zeros(
            Cderivs[[method_name]]$jacobian,
            Cderivs2[[method_name]]$jacobian,
            tolerance = tol2)), info = info
          )
        }
      }, wrap_in_try = isTRUE(catch_failures))
      if (isTRUE(catch_failures) && inherits(jacobian_test, 'try-error')) {
        warning(
          paste0(
            'There was something wrong with the jacobian of ',
            opParam$name, ' with wrt = c(',
            paste0(param$wrts[[method_name]], collapse = ', '), ').\n',
            jacobian_test[1]
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (jacobian_test_fails) {
        if (verbose) {
          cat(paste0(
            "## As expected, test of jacobian failed for ", method_name, ' with wrt: ',
            paste0(param$wrts[[method_name]], collapse = ', '),
            '\n'
          ))
        }
        ## stop testing after an expected failure
        break
      }
      ##
      ## test hessians
      ##
      hessian_test_fails <- is_method_failure(
        param$name, method_name, 'hessian', knownFailures
      )
      hessian_test <- wrap_if_true(hessian_test_fails, expect_failure, {
        if (verbose) cat("## Checking hessians\n")
        expect_equal(
          Cderivs[[method_name]]$hessian,
          Rderivs[[method_name]]$hessian,
          tolerance = tol3, info = info
        )
        if ('log' %in% names(opParam$args)) {
          if (verbose) cat("## Checking log behavior for hessians\n")
          expect_equal(
            Cderivs2[[method_name]]$hessian,
            Rderivs2[[method_name]]$hessian,
            tolerance = tol3, info = info
          )
          expect_false(isTRUE(all_equal_ignore_zeros(
            Rderivs[[method_name]]$hessian,
            Rderivs2[[method_name]]$hessian,
            tolerance = tol3)), info = info
          )
          expect_false(isTRUE(all_equal_ignore_zeros(
            Cderivs[[method_name]]$hessian,
            Cderivs2[[method_name]]$hessian,
            tolerance = tol3)), info = info
          )
        }
      }, wrap_in_try = isTRUE(catch_failures))
      if (isTRUE(catch_failures) && inherits(hessian_test, 'try-error')) {
        warning(
          paste0(
            'There was something wrong with the hessian of ',
            opParam$name, ' with wrt = c(',
            paste0(param$wrts[[method_name]], collapse = ', '), ').\n',
            hessian_test[1]
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      } else if (hessian_test_fails) {
        if (verbose) {
          cat(paste0(
            "## As expected, test of hessian failed for ", method_name, ' with wrt: ',
            paste0(param$wrts[[method_name]], collapse = ', '),
            '\n'
          ))
        }
        ## stop testing after an expected failure
        break
      }
    }
  }
  if (verbose) cat("### Test successful \n\n")
  if (return_compiled_nf)
    invisible(list(CnfInst = CnfInst, input = input))
  else if(!compilation_fails && getNimbleOption('useClearCompiledInADTesting')) {
    nimble:::clearCompiled(CnfInst)
    invisible(list(input = input))
  }
  invisible(NULL)
}

test_AD_batch <- function(batch, dir = file.path(tempdir(), "nimble_generatedCode"),
                          control = list(), verbose = nimbleOptions('verbose'),
                          catch_failures = FALSE, seed = 0,
                          nimbleProject_name = '', knownFailures = list(),
                          testFun = test_AD) {
  ## could try to do something more clever here, like putting the entire batch
  ## into one giant nimbleFunction generator so we only have to compile once
  ## (perhaps conditional on having no knownFailures in the entire batch)
  lapply(
    batch, testFun,
    dir, control, verbose, catch_failures, seed,
    nimbleProject_name, FALSE, knownFailures
  )
  invisible(NULL)
}

#########################################
## AD test parameterization builder utils
#########################################

## Takes a named list of `argTypes` and returns a list of character
## vectors, each of which is valid as the `wrt` argument of `nimDerivs()`.
## Each argument on its own and all combinations of the arguments will
## always be included, and then make_wrt will try to create up to `n_random`
## additional character vectors with random combinations of the arguments
## and indexing of those arguments when possible (i.e. for non-scalar args).
## n_arg_reps determines how many times an argument can be used in a given
## wrt character vector. By default, any argument will appear only 1 time.
make_wrt <- function(argTypes, n_random = 10, n_arg_reps = 1, sizes) {

  ## always include each arg on its own, and all combinations of the args
  wrts <- as.list(names(argTypes))
  if (length(argTypes) > 1)
    for (m in 2:length(argTypes)) {
      this_combn <- combn(names(argTypes), m)
      wrts <- c(
        wrts,
        unlist(apply(this_combn, 2, list), recursive = FALSE)
      )
    }

  argSymbols <- list()
  for(i in names(argTypes)) {
    argType <- argTypes[[i]]
    if(missing(sizes))
      argSymbols[[i]] <- add_missing_size(nimble:::argType2symbol(argType))
    else
      argSymbols[[i]] <- add_missing_size(nimble:::argType2symbol(argType), sizes[[i]])
  }

  ## argSymbols <- lapply(
  ##   argTypes, function(argType)
  ##     add_missing_size(nimble:::argType2symbol(argType))

  while (n_random > 0) {
    n_random  <- n_random - 1
    n <- sample(1:length(argTypes), 1) # how many of the args to use?
    ## grab a random subset of the args of length n
    args <- sample(argSymbols, n)
    ## may repeat an arg up to n_arg_reps times
    reps <- sample(1:n_arg_reps, length(args), replace = TRUE)
    this_wrt <- c()
    for (i in 1:length(args)) {
      while (reps[i] > 0) {
        reps[i] <- reps[i] - 1
        ## coin flip determines whether to index vectors/matrices
        use_indexing <- sample(c(TRUE, FALSE), 1)
        if (use_indexing && args[[i]]$nDim > 0) {
          rand_row <- sample(1:args[[i]]$size[1], size = 1)
          ## another coin flip determines whether to use : in indexing or not
          use_colon <- sample(c(TRUE, FALSE), 1)
          if (use_colon && rand_row < args[[i]]$size[1]) {
            end_row <- rand_row +
              sample(1:(args[[i]]$size[1] - rand_row), size = 1)
            rand_row <- paste0(rand_row, ':', end_row)
          }
          index <- rand_row
          if (args[[i]]$nDim == 2) {
            rand_col <- sample(1:args[[i]]$size[2], size = 1)
            ## one more coin flip to subscript second dimension
            use_colon_again <- sample(c(TRUE, FALSE), 1)
            if (use_colon_again && rand_col < args[[i]]$size[2]) {
              end_col <- rand_col +
                sample(1:(args[[i]]$size[2] - rand_col), size = 1)
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
  unique(wrts)
}

## A version 2 of make_AD_test.
#
# Construct nimbleFunction with following methods:
# value - simply calculate the function, the core of everything.
# derivsValue - get any derivatives of value.
# value2 - get the value from derivsValue. can be used for double (meta) taping.
# jac - get the jacobian from derivsValue via reverse order 1.  can double tape
# jac2 - get the jacobian from derivsValue via forward order 1 (same as en route to hessian). can double tape
# hess - get the hessian from derivsValue. can double tape.
# derivsJac - get any derivatives of jac
# derivsJac2 - get any derivatives of jac2
# derivsHess - get any derivatives of hess.
#
# For later testing (notation is "requested orders" --> "actual orders")
#   e.g. requesting order 0 from derivsJac gives order 1 of the actual function (value). 
# run: gives 0
# derivsRun: 0, 1, 2 --> 0, 1, 2
# value: gives 0 (tests F0)
# jac: gives 1    (tests F0R1)
# jac2: gives 1   (tests F0F1)
# hess: gives 2   (tests F0F1R2)
# derivsValue: 0, 1, 2 --> 0, 1, 2 (tests meta-taped F0 for later derivs)
# derivsJac: 0, 1 --> 1, 2         (tests meta-taped F0R1 for later derivs)
# derivsJac2: 0, 1 --> 1, 2        (tests meta-taped F0F1 for later derivs)
# derivsHess: 0 --> 2              (tests meta-taped F0F0R2 for later derivs)
make_AD_test2 <- function(op, argTypes, wrt_args = NULL,
                          input_gen_funs = NULL, more_args = NULL, seed = 0,
                          outer_code = NULL, inner_codes = NULL,
                          size = NULL, inputs = NULL,
                          includeModelArgs = FALSE) {
  if(!is.list(op)) {
    opParam <- make_op_param(op, argTypes, more_args = more_args,
                             outer_code = outer_code,  inner_codes = inner_codes)
  } else {
    opParam <- op
  }
  run <- gen_runFunCore(opParam)

  if(seed == 0) seed <- round(runif(1, 1, 10000))

  modelArg <- NA
  updateNodesArg <- NA
  constantNodesArg <- NA
  if(isTRUE(includeModelArgs)) {
    # If this is used, then some custom setup code
    # will need to be inserted to create model, updateNodes, and constantNodes
    modelArg <- as.name("model")
    updateNodesArg <- as.name("updateNodes")
    constantNodesArg <- as.name("constantNodes")
  }
  
  
  ## Make a set of methods. These need to be constructed
  ## so that they can have different argument names, numbers, and types.
  args_formals <- lapply(argTypes, function(argType) {
    parse(text = argType)[[1]]
  })
  if (is.null(names(args_formals)))
    names(args_formals) <- paste0('arg', 1:length(args_formals))

  runCall <- parse(text = paste0("run(", paste(names(args_formals), collapse=","), ")"),
                   keep.source = FALSE)[[1]]
  
  derivsRun <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             order = integer(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = order, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      return(ans)
      returnType(ADNimbleList())
    },
    list(RUNCALL = runCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(derivsRun)$srcref <- NULL
  formals(derivsRun) <- c(args_formals, formals(derivsRun))

  value <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 0, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      d1 <- dim(ans$value)[1]
      res <- numeric(length = d1)
      for(i in 1:d1)
        res[i] <- ans$value[i]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(value)$srcref <- NULL
  formals(value) <- c(args_formals, formals(value))
  
  jac <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 1, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      d1 <- dim(ans$jacobian)[1]
      d2  <- dim(ans$jacobian)[2]
      res <- numeric(length = d1*d2)
      for(i in 1:d1)
        for(j in 1:d2)
          res[i + (j-1)*d1] <- ans$jacobian[i, j]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(jac)$srcref <- NULL
  formals(jac) <- c(args_formals, formals(jac))
  
  jac2 <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      # Because this gets order 2, the order 1 comes from forward 1
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 1:2, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      d1 <- dim(ans$jacobian)[1]
      d2  <- dim(ans$jacobian)[2]
      res <- numeric(length = d1*d2)
      for(i in 1:d1)
        for(j in 1:d2)
          res[i + (j-1)*d1] <- ans$jacobian[i, j]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(jac2)$srcref <- NULL
  formals(jac2) <- c(args_formals, formals(jac2))

  hess <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 2, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      d1 <- dim(ans$hessian)[1]
      d2 <- dim(ans$hessian)[2]
      d3 <- dim(ans$hessian)[3]
      res <- numeric(length = d1*d2*d3)
      # There was a bug in pulling out results in a single line.
      # It is unrelated to AD so considered separate.
      for(i in 1:d1)
        for(j in 1:d2)
          for(k in 1:d3)
            res[i + (j-1)*d1 + (k-1)*d1*d2] <- ans$hessian[i, j, k]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(hess)$srcref <- NULL
  formals(hess) <- c(args_formals, formals(hess))

  valueCall <- parse(text = paste0("value(",
                                  paste(c(names(args_formals),
                                          "wrt = innerWrt",
                                          "reset=reset"), collapse=","),
                                  ")"),
                    keep.source = FALSE)[[1]]
  
  derivsValue <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             innerWrt = double(1),
             wrt = double(1),
             order = integer(1),
             reset = logical(0, default = FALSE) ) {
      ans <- nimDerivs(VALUECALL, wrt = wrt, order = order, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      return(ans)
      returnType(ADNimbleList())
    },
    list(VALUECALL = valueCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(derivsValue)$srcref <- NULL
  formals(derivsValue) <- c(args_formals, formals(derivsValue))

  
  jacCall <- parse(text = paste0("jac(",
                                  paste(c(names(args_formals),
                                          "wrt = innerWrt",
                                          "reset=reset"), collapse=","),
                                  ")"),
                    keep.source = FALSE)[[1]]
  
  derivsJac <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             innerWrt = double(1),
             wrt = double(1),
             order = integer(1),
             reset = logical(0, default = FALSE) ) {
      ans <- nimDerivs(JACCALL, wrt = wrt, order = order, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      return(ans)
      returnType(ADNimbleList())
    },
    list(JACCALL = jacCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(derivsJac)$srcref <- NULL
  formals(derivsJac) <- c(args_formals, formals(derivsJac))

  jac2Call <- parse(text = paste0("jac2(",
                                  paste(c(names(args_formals),
                                          "wrt = innerWrt",
                                          "reset=reset"), collapse=","),
                                  ")"),
                    keep.source = FALSE)[[1]]
  
  derivsJac2 <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             innerWrt = double(1),
             wrt = double(1),
             order = integer(1),
             reset = logical(0, default = FALSE) ) {
      ans <- nimDerivs(JAC2CALL, wrt = wrt, order = order, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      return(ans)
      returnType(ADNimbleList())
    },
    list(JAC2CALL = jac2Call, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
  ))
  attributes(derivsJac2)$srcref <- NULL
  formals(derivsJac2) <- c(args_formals, formals(derivsJac2))

  
  hessCall <- parse(text = paste0("hess(",
                                  paste(c(names(args_formals),
                                          "wrt = innerWrt",
                                          "reset=reset"), collapse=","),
                                  ")"),
                    keep.source = FALSE)[[1]]
  derivsHess  <-  eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             innerWrt = double(1),
             wrt = double(1),
             order = integer(1),
             reset = logical(0, default = FALSE) ) {
      ans <- nimDerivs(HESSCALL, wrt = wrt, order = order, reset = reset,
                       model = MODELARG, updateNodes = UPDATENODESARG,
                       constantNodes = CONSTANTNODESARG)
      return(ans)
      returnType(ADNimbleList())
    },
    list(HESSCALL = hessCall, MODELARG = modelArg,
         UPDATENODESARG = updateNodesArg,
         CONSTANTNODESARG = constantNodesArg)
    ))
  attributes(derivsHess)$srcref <- NULL
  formals(derivsHess) <- c(args_formals, formals(derivsHess))

  methods <- list(derivsRun = derivsRun,
                  value = value,
                  jac = jac,
                  jac2 = jac2,
                  hess = hess,
                  derivsValue = derivsValue,
                  derivsJac = derivsJac,
                  derivsJac2 = derivsJac2,
                  derivsHess = derivsHess
                  )
  
  list(
    name = opParam$name,
    opParam = opParam,
    run = run,
    methods = methods,
    buildDerivs = list(run = list()
                       ,
                       value = list(ignore = c('wrt', 'i')),
                       jac = list(ignore = c('wrt', 'i', 'j')),
                       jac2 = list(ignore = c('wrt', 'i', 'j')),
                       hess = list(ignore = c('wrt', 'i', 'j', 'k'))
                        ),
    wrt_args = wrt_args,
    input_gen_funs = input_gen_funs,
    size = size,
    seed = seed,
    inputs = inputs
  )
}


model_calculate_test_pieces <- make_AD_test2(
  op = list(
    expr = quote( {
      values(model, derivNodes) <<- arg1
      out <- model$calculate(calcNodes)
      return(out)
    }),
    args = list(arg1 = quote(double(1))),
    outputType = quote(double())
  ),
  argTypes = list(arg1 = "double(1)"),
  includeModelArgs = TRUE)

model_calculate_test <- nimbleFunction(
  setup = function(model, nodesList) {
    derivNodes <- nodesList$derivNodes
    updateNodes <- nodesList$updateNodes
    constantNodes <- nodesList$constantNodes
    calcNodes <- nodesList$calcNodes
    nNodes <- length(derivNodes)
  },
  run = model_calculate_test_pieces$run,
  methods = model_calculate_test_pieces$methods,
  buildDerivs = model_calculate_test_pieces$buildDerivs
)

setup_update_and_constant_nodes_for_tests <- function(model,
                                            derivNodes,
                                            forceConstantNodes = character(),
                                            forceUpdateNodes = character()) {
  ## "update" means "CppAD dynamic"
  ##  derivNodes <- model$expandNodeNames(derivNodes) # do not do this because do not want vector node names
  nNodes <- length(derivNodes)
  calcNodes <- model$getDependencies(derivNodes)
  ucNodes <- makeModelDerivsInfo(model, derivNodes, calcNodes, dataAsConstantNodes = TRUE)
  updateNodes <- ucNodes$updateNodes
  constantNodes <- ucNodes$constantNodes
  updateNodes <- setdiff(updateNodes, forceConstantNodes) # remove forceConstants from updates
  constantNodes <- setdiff(constantNodes, forceUpdateNodes) # remove forceUpdates from constants
  constantNodes <- union(constantNodes, forceConstantNodes) # add forceConstants to constants
  updateNodes <- union(updateNodes, forceUpdateNodes)     # add forceUpdates to updates
  list(derivNodes = derivNodes, updateNodes = updateNodes,
       constantNodes = constantNodes, calcNodes = calcNodes)
}

model_calculate_test_case <- function(Rmodel, Cmodel,
                                      deriv_nf, nodesList,
                                      v1, v2  ,
                                      order,
                                      varValues = list(),
                                      varValues2 = list(), ...) {
                                        # This sets up an instance of a checker fxn
  Rfxn <- deriv_nf( Rmodel, nodesList)
  Cfxn <- compileNimble(Rfxn, project = Rmodel)

  test_AD2_oneCall(Rfxn, Cfxn,
                   recordArgs = v1, testArgs = v2,
                   order = order,
                   Rmodel = Rmodel, Cmodel = Cmodel,
                   recordInits = varValues, testInits = varValues2,
                   nodesToChange = c(nodesList$updateNodes),
                   ...)
}

## Make a test parameterization to be used by test_AD. This method is primarily
## used by make_AD_test_batch() and make_distribution_fun_AD_test().
##
## op:             Character string, the operator that will be the focus of the
##                 test.
## argTypes:       Character vector of argType strings that, when parsed, can be
##                 passed to argType2symbol. If named, the names will as the formals
##                 of the nimbleFunction generator's run method and other methods.
##                 If not, formals are generated as arg1, arg2, etc.
## wrt_args:       Optional character vector of args to use in make_wrt(). If NULL,
##                 assumes that all the arguments should be used.
## input_gen_funs: A list of input generation functions which is simply passed to
##                 the output list. Should be NULL (use defaults found in
##                 arg_type_2_input()), length 1 (use same input gen mechanism for each
##                 argType, or a named list with names from among the argType names
##                 (possibly the sequentially generated names). This will be NULL
##                 when bulk generating the test params using make_AD_test_batch and
##                 added later via modify_on_match(). Used in the call to
##                 make_AD_test() in make_distribution_fun_AD_test(). 
## more_args:      A named list of additional fixed arguments to use in the
##                 generated operator call. E.g., if op = 'dnorm',
##                 argTypes = c('double(1, 5)', 'double(0)'), and
##                 more_args = list(log = 1), the call to make_op_param will include
##                 the expression dnorm(arg1, arg2, log = 1).
## seed:           A seed to use in set.seed().
##
## returns: A list with the following elements:
##          name:           character string, the operator and args
##          opParam:        result from make_op_param
##          run:            result from calling gen_runFunCore with opParam
##          methods:        a list of nimbleFunction method expressions that are
##                          the calls to nimDerivs() of the run method, each with
##                          a different wrt argument
##          buildDerivs:   list('run')
##          wrts:           a list of character vectors, each of which is the wrt
##                          argument for the corresponding method in methods
##          input_gen_funs: A list of random input generation functions to be
##                          used by arg_type_2_input(). 
make_AD_test <- function(op, argTypes, wrt_args = NULL,
                         input_gen_funs = NULL, more_args = NULL, seed = 0, ...) {
  ## set the seed for make_wrt
  if (is.numeric(seed)) set.seed(seed)
  opParam <- make_op_param(op, argTypes, more_args)

  run <- gen_runFunCore(opParam)
  method <- function() {
    {}
    returnType(ADNimbleList())
    return(outList)
  }
  formals_list <- lapply(argTypes, function(argType) {
    parse(text = argType)[[1]]
  })
  if (is.null(names(formals_list)))
    names(formals_list) <- paste0('arg', 1:length(formals_list))
  formals(method) <- formals_list

  body(method)[[2]] <-
    substitute(
      outList <- derivs(this_call),
      list(
        this_call = as.call(
          c(quote(run), lapply(names(formals_list), as.name))
        )
      )
    )

  if (is.null(wrt_args)) wrt_args_filter <- rep(TRUE, length(argTypes))
  else wrt_args_filter <- wrt_args

  dotArgs <- list(...)
  sizes = dotArgs$sizes
  if(is.null(sizes))
    wrts <- make_wrt(opParam$args[wrt_args_filter])
  else
    wrts <- make_wrt(opParam$args[wrt_args_filter], sizes = sizes)

  methods <- list()
  for (i in seq_along(wrts)) {
    method_i <- paste0('method', i)
    methods[[method_i]] <- method
    body(methods[[method_i]])[[2]][[3]][['wrt']] <- wrts[[i]]
  }

  if (is.null(wrt_args) || all(names(opParam$args) %in% wrt_args)) {
    method_no_wrt <- paste0('method', length(methods) + 1)
    methods[[method_no_wrt]] <- method
    wrts <- c(wrts, 'no wrt')
  }

  names(wrts) <- names(methods)

  ## method_all_wrt <- paste0('method', length(wrt) + 2)
  ## wrts[[method_all_wrt]] <- paste0('wrt ', paste0(wrt, collapse = ', '))
  ## methods[[method_all_wrt]] <- methods[[method_no_wrt]] <-  method
  ## body(methods[[method_all_wrt]])[[2]][[3]][['wrt']] <- wrt

  list(
    name = opParam$name,
    opParam = opParam,
    run = run,
    methods = methods,
    buildDerivs = 'run',
    wrts = wrts,
    input_gen_funs = input_gen_funs
  )
}

## ops: character vector of operator names
## argTypes: list of character vectors of argTypes
##           e.g. for a binary operator:
##             list(
##               c('double(1, 4)', 'double(0)'),
##               c('double(1, 4)', 'double(1, 4)')
##             )
make_AD_test_batch <- function(ops, argTypes, seed = 0,
                               maker = make_AD_test,
                               outer_code = NULL, inner_codes = NULL) {
  opTests <- vector(mode = 'list', length = length(ops) * length(argTypes))
  for(i in seq_along(ops)) {
    for(j in seq_along(argTypes)) {
      iOut <- (i-1) * length(argTypes) + j
      opTests[[iOut]] <- maker(op = ops[[i]], argTypes = argTypes[[j]],
                              seed = seed, outer_code = outer_code, inner_codes = inner_codes)
    }
  }
  
  ## opTests <- unlist(
  ##   recursive = FALSE,
  ##   x = lapply(
  ##     ops,
  ##     function(x) {
  ##       mapply(
  ##         maker,
  ##         argTypes = argTypes,
  ##         MoreArgs = list(op = x, seed = seed,
  ##                         outer_code = outer_code, inner_codes = inner_codes),
  ##         SIMPLIFY = FALSE
  ##       )
  ##     })
  ## )
  names(opTests) <- sapply(opTests, `[[`, 'name')
  invisible(opTests)
}

## Takes an element of distn_params list and returns a list of AD test
## parameterizations, each of which test_AD can use.
##
## distn_param: Probability distribution parameterization, which
##              must have the following fields/subfields:
##              - name
##              - variants
##              - args
##                - rand_variate
##                  - type
##              Additional args must also have type field.
## more_args:   Passed to make_op_param().
##
make_distribution_fun_AD_test <- function(distn_param, maker = make_AD_test) {
  distn_name <- distn_param$name
  ops <- sapply(distn_param$variants, paste0, distn_name, simplify = FALSE)

  rand_variate_idx <- which(names(distn_param$args) == 'rand_variate')

  argTypes <- sapply(distn_param$variants, function(variant) {
    op <- paste0(variant, distn_name)
    rand_variate_type <- distn_param$args$rand_variate$type
    first_argType <- switch(
      variant,
      d = rand_variate_type,
      p = rand_variate_type,
      q = c('double(0)')##, 'double(1, 4)')
    )
    first_arg_name <- switch(variant, d = 'x', p = 'q', q = 'p')
    ## need this complicated expand.grid call here because the argTypes might be
    ## character vectors (e.g. if rand_variate_type is c("double(0)", "double(1, 4)")
    ## then we create a test where we sample a scalar from the support and
    ## another where we sample a vector of length 4 from the support
    grid <- eval(as.call(c(
      expand.grid, list(first_argType),
      lapply(distn_param$args, `[[`, 'type')[-rand_variate_idx]
    )))
    argTypes <- as.list(data.frame(t(grid), stringsAsFactors=FALSE))
    lapply(argTypes, function(v) {
      names(v) <- c(
        first_arg_name,
        names(distn_param$args[-rand_variate_idx])
      )
      v
    })
  }, simplify = FALSE)
  for(i in seq_along(distn_param$args)) {
    if(is.null(distn_param$args[[i]]$size))
      distn_param$args[[i]]$size <- lapply(distn_param$args[[i]]$type,
                                           function(type) {
                                             if(is.character(type))
                                               type <- parse(text = type, keep.source=FALSE)[[1]]
                                             add_missing_size(nimble:::argType2symbol(type))$size
                                           })
  }
  argSizes <- sapply(distn_param$variants, function(variant) {
    rand_variate_size <- distn_param$args$rand_variate$size
    first_argSize <- switch(
      variant,
      d = rand_variate_size,
      p = rand_variate_size,
      q = 1
    )
    first_arg_name <- switch(variant, d = 'x', p = 'q', q = 'p')
    grid <- eval(as.call(c(
      expand.grid, list(first_argSize),
      lapply(distn_param$args, `[[`, 'size')[-rand_variate_idx]
    )))
    argSizes <- as.list(data.frame(t(grid), stringsAsFactors=FALSE))
    lapply(argSizes, function(v) {
      names(v) <- c(
        first_arg_name,
        names(distn_param$args[-rand_variate_idx])
      )
      v
    })
  }, simplify = FALSE)
  input_gen_funs <- lapply(distn_param$args[-rand_variate_idx], `[[`, 'input_gen_fun')
  input_gen_funs_list <- sapply(distn_param$variants, function(variant) {
    arg1_input_gen_fun <- switch(
      variant,
      d = list(x = distn_param$args$rand_variate$input_gen_fun),
      p = list(q = distn_param$args$rand_variate$input_gen_fun),
      q = list(p = runif)
    )
    c(arg1_input_gen_fun, input_gen_funs)
  }, simplify = FALSE)

  op_params <- unlist(lapply(distn_param$variants,
                             function(variant) {
                               ans <- list()
                               for(i in seq_along(argTypes[[variant]])) {
                                 these_argTypes <- argTypes[[variant]][[i]]
                                 these_argSizes  <- argSizes[[variant]][[i]]
                                 wrt_args = intersect(
                                   distn_param$wrt, names(these_argTypes)
                                 )
                                 new_ans <-  maker(
                                   ops[[variant]], these_argTypes,
                                   wrt_args = wrt_args,
                                   input_gen_funs = input_gen_funs_list[[variant]],
                                   more_args = distn_param$more_args[[variant]],
                                   size = these_argSizes
                                 )
                                 ans[[length(ans)+1]]<-new_ans
                               }
                               ans
                             }), recursive = FALSE)
    
  ## op_params <- unlist(
  ##   lapply(
  ##     distn_param$variants, function(variant) {
  ##       lapply(
  ##         argTypes[[variant]],
  ##         function(these_argTypes) {
  ##           wrt_args = intersect(
  ##             distn_param$wrt, names(these_argTypes)
  ##           )
  ##           maker(
  ##             ops[[variant]], these_argTypes,
  ##             wrt_args = wrt_args,
  ##             input_gen_funs = input_gen_funs_list[[variant]],
  ##             more_args = distn_param$more_args[[variant]]
  ##           )
  ##         }
  ##       )
  ##     }
  ##   ), recursive = FALSE
  ## )
  names(op_params) <- sapply(op_params, `[[`, 'name')
  return(op_params)
}

#############################
## input generation functions
#############################

## arg_size comes from arg$size where arg is a symbolBasic object
gen_pos_def_matrix <- function(arg_size, rfun = runif) {
  m <- arg_size[1] ## assumes matrix argType is square
  if(length(arg_size) == 1) # single length defined as arg_size = m^2
    m <- sqrt(m)
  mat <- diag(m)
  mat[lower.tri(mat, diag = TRUE)] <- rfun(m*(m + 1)/2)
  mat %*% t(mat)
}


## nf that embeds deriv of model calculate
derivsNimbleFunction <- nimbleFunction(
  setup = function(model, calcNodes, wrt) {},
  run = function(x = double(1),
                 order = double(1),
                 reset = logical(0, default = FALSE)) {
    values(model, wrt) <<- x  
    ans <- nimDerivs(model$calculate(calcNodes), wrt = wrt, order = order, reset = reset)
    return(ans)
    returnType(ADNimbleList())
  }  ## don't need buildDerivs if call nimDerivs directly, but would if just have model$calc in nf
)

## nf for double-taping
derivsNimbleFunctionMeta <- nimbleFunction(
  setup = function(model, calcNodes, wrt, reset = FALSE) {
    innerWrtVec <- seq_along(model$expandNodeNames(wrt, returnScalarComponents = TRUE))
    d <- length(innerWrtVec)  
    derivsInfo <- makeModelDerivsInfo(model, wrt, calcNodes)
    updateNodes <- derivsInfo$updateNodes
    constantNodes <- derivsInfo$constantNodes
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
  buildDerivs = list(run = list(),
                      derivs1Run = list(),
                      derivs2Run = list())
)


derivsNimbleFunctionParamTransform <- nimbleFunction(
    setup = function(model, calcNodes, wrt) {
        wrtNodesAsScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
        my_parameterTransform <- parameterTransform(model, wrtNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        nimDerivs_wrt <- 1:d
        derivsInfo <- makeModelDerivsInfo(model, wrt, calcNodes)
        updateNodes   <- derivsInfo$updateNodes
        constantNodes <- derivsInfo$constantNodes
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
    ), buildDerivs = 'inverseTransformStoreCalculate'
)


derivsNimbleFunctionParamTransformMeta <- nimbleFunction(
    setup = function(model, calcNodes, wrt, reset = FALSE) {
        wrtNodesAsScalars <- model$expandNodeNames(wrt, returnScalarComponents = TRUE)
        my_parameterTransform <- parameterTransform(model, wrtNodesAsScalars)
        d <- my_parameterTransform$getTransformedLength()
        nimDerivs_wrt <- 1:d
        derivsInfo <- makeModelDerivsInfo(model, wrt, calcNodes)
        updateNodes   <- derivsInfo$updateNodes
        constantNodes <- derivsInfo$constantNodes
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
    buildDerivs = list(run = list(),
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
    buildDerivs = 'run')


## By default test a standardized set of {wrt, calcNodes} pairs representing common use cases (MAP, max lik, EB),
## unless user provides 'wrt' and 'calcNodes'.
test_ADModelCalculate <- function(model, name = 'unknown', x = 'given', xNew = NULL, calcNodes = NULL, wrt = NULL,
                                  newUpdateNodes = NULL, newConstantNodes = NULL,
                                  relTol = c(1e-15, 1e-8, 1e-3, 1e-3, 1e-14), absTolThreshold = 0, useFasterRderivs = FALSE, useParamTransform = FALSE,
                                  checkDoubleTape = TRUE, checkCompiledValuesIdentical = TRUE, checkDoubleUncHessian = TRUE,
                                  doAllUncHessian = TRUE, seed = 1, verbose = FALSE, debug = FALSE){
    if(!is.null(seed))
        set.seed(seed)
    initsHandling <- x
    xNewIn <- xNew
    ## Save model state so can restore for later use cases below.
    mv <- modelValues(model)
    nodes <- model$getNodeNames(includeRHSonly = TRUE)
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
                                           relTol = c(1e-15, 1e-8, 1e-3, 1e-3, 1e-14), absTolThreshold = 0, useFasterRderivs = FALSE,
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

        nodes  <- model$getNodeNames(includeRHSonly = TRUE)
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
        
        derivsInfo <- makeModelDerivsInfo(model, wrt, calcNodes)
        updateNodes <- derivsInfo$updateNodes
        constantNodes <- derivsInfo$constantNodes

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
                } else {
                    expect_fun <- function(x,y) {
                            xlab <- deparse1(substitute(x))
                            ylab <- deparse1(substitute(y))
                            nim_expect_equal(x, y, tolerance = relTol[5], abs_threshold = absTolThreshold,
                                             xlab = xlab, ylab = ylab)
                    }
                }

                ## If go through transform and back, may not have identical.
                expect_fun_wrt <- ifelse(useParamTransform, expect_equal, expect_identical)
                
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
                    nim_expect_equal(cOutput2d$value, c(cOutput012$hessian), tolerance = relTol[5], abs_threshold = absTolThreshold)
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

                expect_fun_wrt(rWrt01, x)
                ## We provide the model to nimDerivs, so expect restoration except for non-paramTransform
                ## and non-fasterRderivs, where assignment is outside nimDerivs.
                if(doAllUncHessian) {
                    if(!useParamTransform && !useFasterRderivs) {
                        expect_fun_wrt(rWrt12, x)
                    } else expect_fun_wrt(rWrt12, rWrt_orig)
                }
                expect_fun_wrt(rWrt012, x)
                expect_fun_wrt(cWrt01, x)
                if(!useParamTransform) {
                    ## Assignment to model is outside nimDerivs call, so expect change.
                    expect_fun_wrt(cWrt12, x)
                } else expect_fun_wrt(cWrt12, cWrt_orig)
                expect_fun_wrt(cWrt012, x)

                if(checkDoubleTape) {
                    expect_fun_wrt(rWrt1d, rWrt_orig)
                    if(doAllUncHessian) 
                        expect_fun_wrt(rWrt2d, rWrt_orig)
                    if(checkDoubleUncHessian)
                        expect_fun_wrt(rWrt2d11, rWrt_orig)
                    expect_fun_wrt(cWrt1d, cWrt_orig)
                    expect_fun_wrt(cWrt2d, cWrt_orig)
                    expect_fun_wrt(cWrt2d11, cWrt_orig)
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


## ## Makes random vectors of wrt elements, following James Duncan's code
## make_wrt <- function(argTypes, n_random = 10, allCombinations = FALSE) {
##     ## always include each arg on its own, and all combinations of the args
##     ## Note that for models with a large number of variables this might turn out to be too much.
##     wrts <- as.list(names(argTypes))
##     if(allCombinations) {
##         if (length(argTypes) > 1)
##             for (m in 2:length(argTypes)) {
##                 this_combn <- combn(names(argTypes), m)
##                 wrts <- c(
##                     wrts,
##                     unlist(apply(this_combn, 2, list), recursive = FALSE)
##                 )
##             }
##     }

##   while (n_random > 0) {
##     n_random  <- n_random - 1
##     n <- sample(1:length(argTypes), 1) # how many of the args to use?
##     ## grab a random subset of the args of length n
##     args <- sample(argTypes, n)
##     ## may repeat an arg up to 2 times
##     reps <- sample(1:2, length(args), replace = TRUE)
##     argSymbols <- lapply(args, nimble:::argType2symbol)
##     this_wrt <- c()
##     for (i in 1:length(args)) {
##       while (reps[i] > 0) {
##         reps[i] <- reps[i] - 1
##         ## coin flip determines whether to index vectors/matrices
##         use_indexing <- sample(c(TRUE, FALSE), 1)
##         if (use_indexing && argSymbols[[i]]$nDim > 0) {
##           rand_row <- sample(1:argSymbols[[i]]$size[1], size = 1)
##           ## another coin flip determines whether to use : in indexing or not
##           use_colon <- sample(c(TRUE, FALSE), 1)
##           if (use_colon && rand_row < argSymbols[[i]]$size[1]) {
##             end_row <- rand_row +
##               sample(1:(argSymbols[[i]]$size[1] - rand_row), size = 1)
##             rand_row <- paste0(rand_row, ':', end_row)
##           }
##           index <- rand_row
##           if (argSymbols[[i]]$nDim == 2) {
##             rand_col <- sample(1:argSymbols[[i]]$size[2], size = 1)
##             ## one more coin flip to subscript second dimension
##             use_colon_again <- sample(c(TRUE, FALSE), 1)
##             if (use_colon_again && rand_col < argSymbols[[i]]$size[2]) {
##               end_col <- rand_col +
##                 sample(1:(argSymbols[[i]]$size[2] - rand_col), size = 1)
##               rand_col <- paste0(rand_col, ':', end_col)
##             }
##             index <- paste0(index, ',', rand_col)
##           }
##           this_wrt <- c(this_wrt, paste0(names(args)[i], '[', index, ']'))
##         }
##         ## if first coin flip was FALSE, just
##         ## use the arg name without indexing
##         else this_wrt <- c(this_wrt, names(args)[i])
##       }
##     }
##     if (!is.null(this_wrt)) wrts <- c(wrts, list(unique(this_wrt)))
##   }
##   wrts <- unique(wrts)
## }


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

test_ADDistribution <- function(ADfunGen, argsList, name, debug = FALSE){
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
    if(getNimbleOption('useClearCompiledInADTesting'))
        nimble:::clearCompiled(CADfun)
}
