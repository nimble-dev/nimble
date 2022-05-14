source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# An environment for storing some defaults for AD testing.
# Having this allows one to change these values globally.
ADtestEnv <- new.env()
resetTols <- function() {
  ADtestEnv$RRrelTol <<- c(1e-8, 1e-3, 1e-2)
  ADtestEnv$RCrelTol <<- c(1e-15, 1e-8, 1e-3)
  ADtestEnv$CCrelTol <<- rep(sqrt(.Machine$double.eps), 3)
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
                             RRrelTol = c(1e-8, 1e-3, 1e-2), 
                             RCrelTol = c(1e-15, 1e-8, 1e-3),
                             CCrelTol = rep(sqrt(.Machine$double.eps), 3)) {
  resRecord <- list()
  resTest <- list()

  if(is.null(names(RRrelTol))) names(RRrelTol) <- c("0", "1", "2")
  if(is.null(names(RCrelTol))) names(RCrelTol) <- c("0", "1", "2")
  if(is.null(names(CCrelTol))) names(CCrelTol) <- c("0", "1", "2")
  
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
    if(doR)
      RansRecord <- do_one_call(Rfun, argList)
    if(doC)
      CansRecord <- do_one_call(Cfun, argList)
    if(tapingLevel >= 1)
      extraArgs$reset <- FALSE

    argList <- c(testArgs, extraArgs)
    if(doR)
      RansTest <- do_one_call(Rfun, argList)
    if(doC)
      CansTest <- do_one_call(Cfun, argList)

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

  all_equal_list <- function(first, others, tol, info = "") {
    pass <- TRUE
    for(i in seq_along(others)) {
      pass <- pass && nim_all_equal(as.numeric(first), as.numeric(others[[i]]),
                                    tol, verbose = TRUE, info = info)
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
                                         info = paste0("(RR order ", o,")"))
          if(!pass) {
            cat(paste('Some R-to-R derivatives do not match for order',o))
            browser()
          }
        }
        if(length(RansSet) > 0 && length(CansSet) > 0) {
          pass <- pass && all_equal_list(RansSet[[1]], CansSet, tol = RCrelTol[[o]],
                                         info = paste0("(RC order ", o,")"))
          if(!pass) {
            cat(paste('Some C-to-R derivatives to not match for order',o))
            browser()
          }
        }
        if(length(CansSet) > 1) {
          pass <- pass && all_equal_list(CansSet[[1]], CansSet[-1], CCrelTol[[o]],
                                         info = paste0("(CC order ", o, ")"))
          if(!pass) {
            cat(paste('Some C-to-C derivatives to not match for order',o))
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
  if(!compilation_fails) {
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
  else if(!compilation_fails) {
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
make_wrt <- function(argTypes, n_random = 10, n_arg_reps = 1) {

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

  argSymbols <- lapply(
    argTypes, function(argType)
      add_missing_size(nimble:::argType2symbol(argType))
  )

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
                          size = NULL, inputs = NULL) {
  if(!is.list(op)) {
    opParam <- make_op_param(op, argTypes, more_args = more_args,
                             outer_code = outer_code,  inner_codes = inner_codes)
  } else {
    opParam <- op
  }
  run <- gen_runFunCore(opParam)

  if(seed == 0) seed <- round(runif(1, 1, 10000))
  
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
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = order, reset = reset)
      return(ans)
      returnType(ADNimbleList())
    },
    list(RUNCALL = runCall)
  ))
  attributes(derivsRun)$srcref <- NULL
  formals(derivsRun) <- c(args_formals, formals(derivsRun))

  value <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 0, reset = reset)
      d1 <- dim(ans$value)[1]
      res <- numeric(length = d1)
      for(i in 1:d1)
        res[i] <- ans$value[i]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall)
  ))
  attributes(value)$srcref <- NULL
  formals(value) <- c(args_formals, formals(value))
  
  jac <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 1, reset = reset)
      d1 <- dim(ans$jacobian)[1]
      d2  <- dim(ans$jacobian)[2]
      res <- numeric(length = d1*d2)
      for(i in 1:d1)
        for(j in 1:d2)
          res[i + (j-1)*d1] <- ans$jacobian[i, j]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall)
  ))
  attributes(jac)$srcref <- NULL
  formals(jac) <- c(args_formals, formals(jac))
  
  jac2 <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      # Because this gets order 2, the order 1 comes from forward 1
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 1:2, reset = reset)
      d1 <- dim(ans$jacobian)[1]
      d2  <- dim(ans$jacobian)[2]
      res <- numeric(length = d1*d2)
      for(i in 1:d1)
        for(j in 1:d2)
          res[i + (j-1)*d1] <- ans$jacobian[i, j]
      return(res)
      returnType(double(1))
    },
    list(RUNCALL = runCall)
  ))
  attributes(jac2)$srcref <- NULL
  formals(jac2) <- c(args_formals, formals(jac2))

  hess <- eval(substitute(
    function(#arg1 = double(1) etc to be inserted
             wrt = double(1),
             reset = logical(0, default = FALSE)) {
      ans <- nimDerivs(RUNCALL, wrt = wrt, order = 2, reset = reset)
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
    list(RUNCALL = runCall)
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
      ans <- nimDerivs(VALUECALL, wrt = wrt, order = order, reset = reset)
      return(ans)
      returnType(ADNimbleList())
    },
    list(VALUECALL = valueCall)
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
      ans <- nimDerivs(JACCALL, wrt = wrt, order = order, reset = reset)
      return(ans)
      returnType(ADNimbleList())
    },
    list(JACCALL = jacCall)
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
      ans <- nimDerivs(JAC2CALL, wrt = wrt, order = order, reset = reset)
      return(ans)
      returnType(ADNimbleList())
    },
    list(JAC2CALL = jac2Call)
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
      ans <- nimDerivs(HESSCALL, wrt = wrt, order = order, reset = reset)
      return(ans)
      returnType(ADNimbleList())
    },
    list(HESSCALL = hessCall)
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

  methods <- list()

  if (is.null(wrt_args)) wrt_args_filter <- rep(TRUE, length(argTypes))
  else wrt_args_filter <- wrt_args
  wrts <- make_wrt(opParam$args[wrt_args_filter])
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
gen_pos_def_matrix <- function(arg_size) {
  m <- arg_size[1] ## assumes matrix argType is square
  if(length(arg_size) == 1) # single length defined as arg_size = m^2
    m <- sqrt(m)
  mat <- diag(m)
  mat[lower.tri(mat, diag = TRUE)] <- runif(m*(m + 1)/2)
  mat %*% t(mat)
}
