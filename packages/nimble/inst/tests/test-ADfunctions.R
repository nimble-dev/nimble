source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

context("Testing of derivatives for nimbleFunctions.")

test_that('Derivatives of dnorm function correctly.',
  {
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(y = double(1)) {
        outList <- derivs(testMethod(y), wrt = c('x'))
        returnType(ADNimbleList())
        return(outList)
      },
      methods = list(
        testMethod = function(x = double(1, 2)) {
          out <- dnorm(x[1],0,1)
          returnType(double())
          return(out)
        }
      ), enableDerivs = list('testMethod')
    )
    ADfunInst <- ADfun1()
    x <- matrix(c(2, -2))
    Rderivs <- ADfunInst$run(x)
    temporarilyAssignInGlobalEnv(ADfunInst)  
    cADfunInst <- compileNimble(ADfunInst)
    cderivs <- cADfunInst$run(x)
    expect_equal(cderivs$value, Rderivs$value)
    expect_equal(cderivs$jacobian, Rderivs$jacobian)
    expect_equal(cderivs$hessian, Rderivs$hessian)
  }
)



test_that('Derivatives of x^2 function correctly.',
          {
            listOfADLists <- nimbleList(list1 = ADNimbleList(),
                                        list2 = ADNimbleList())
            temporarilyAssignInGlobalEnv(listOfADLists)  
            ADfun2 <- nimbleFunction(
              setup = function(){},
              run = function(y = double(1, c(2))) {
                ADlistList <- listOfADLists$new()
                ADlistList$list1 <- derivs(testMethod(y, y))
                ADlistList$list2 <- derivs(testMethod2(y, y))
                returnType(listOfADLists())
                return(ADlistList)
              },
              methods = list(
                testMethod = function(x = double(1, c(2)), y = double(1, c(2))) {
                  returnType(double(1, c(2)))
                  return(x[]^2)
                },
                testMethod2 = function(x = double(1, c(2)), y = double(1, c(2))) {
                  returnType(double(0))
                  return(2^x[1])
                }
                
              ), enableDerivs = list('testMethod', 'testMethod2')
            )
            ADfunInst <- ADfun2()
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs$list1$value, Rderivs$list1$value)
            expect_equal(cderivs$list1$jacobian, Rderivs$list1$jacobian, tolerance = 0.01)
            expect_equal(cderivs$list1$hessian, Rderivs$list1$hessian, tolerance = 0.01)
            expect_equal(cderivs$list2$value, Rderivs$list2$value)
            expect_equal(cderivs$list2$jacobian, Rderivs$list2$jacobian, tolerance = 0.01)
            expect_equal(cderivs$list2$hessian, Rderivs$list2$hessian, tolerance = 0.01)
          }
)

test_that('Derivatives of sum(log(x)) function correctly.',
          {
            ADfun3 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(1, c(2))) {
                outVals <- 2*derivs(testMethod(x))$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(1, c(2))) {
                  returnType(double(0))
                  return(sum(log(x[])))
                }
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun3()
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)


test_that('Derivatives of model$calculate() work in expressions.',
          {
            testModelCode <- nimbleCode({
              x[1:2] ~ dmnorm(zeroVec[1:2], diagMat[1:2, 1:2])
            })
            testModel <- nimbleModel(code = testModelCode,
                                     inits = list(x = rep(0, 2)),
                                     constants = list(zeroVec = rep(0, 2),
                                                      diagMat = diag(2)))
            ADfun4 <- nimbleFunction(
              setup = function(model){},
              run = function() {
                outVals <- 2*sum(derivs(model$calculate('x'), 
                                        wrt = 'x')$jacobian)
                returnType(double())
                return(outVals)
              }
            )
            ADfunInst <- ADfun4(testModel)
            Rderivs <- ADfunInst$run()
            temporarilyAssignInGlobalEnv(testModel)  
            temporarilyAssignInGlobalEnv(ADfunInst)  
            compileNimble(testModel)
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run()
            print(cderivs)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)

test_that('Derivatives of matrix multiplication function correctly.',
          {
            ADfun5 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(2, c(2, 2)),
                             y = double(2, c(2, 4))) {
                outVals <- derivs(testMethod(x, y), wrt = 'x', order = 1)$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(2, c(2, 2)),
                                      y = double(2, c(2, 4))) {
                  returnType(double(2, c(2, 4)))
                  return(x%*%y)
                }
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun5()
            x <- diag(2)
            set.seed(1)
            y <- matrix(rnorm(8), ncol = 4, nrow = 2)
            Rderivs <- ADfunInst$run(x, y)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x, y)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)



test_that('Derivatives of matrix exponentiation function correctly.',
          {
            ADfun6 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(2, c(2, 2))) {
                outVals <- derivs(testMethod(x), wrt = 'x', order = 1)$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(2, c(2, 2))) {
                  returnType(double(2, c(2, 2)))
                  return(exp(x))
                }
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun6()
            x <- diag(2)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)


## one for exp(mat)
## get order to owrk without c()!

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
                    nimbleProject_name = '', return_compiled_nf = FALSE) {
  if (!is.null(param$debug) && param$debug) browser()
  if (verbose) cat(paste0("### Testing ", param$name, "\n"))

  ## by default, reset the seed for every test
  if (is.numeric(seed)) set.seed(seed)

  nf <- nimbleFunction(
    setup = function() {},
    run = param$run,
    methods = param$methods,
    enableDerivs = param$enableDerivs
  )
  nfInst <- nf()

  opParam <- param$opParam

  if (is.null(param$arg_distns) || is.null(names(param$arg_distns)))
    if (length(param$arg_distns) <= 1)
      input <- lapply(opParam$args, arg_type_2_input, param$arg_distns)
    else
      stop(
        'arg_distns of length greater than 1 must have names',
        call. = FALSE
      )
  else {
    input <- sapply(
      names(opParam$args),
      function(name)
        arg_type_2_input(opParam$args[[name]], param$arg_distns[[name]]),
      simplify = FALSE
    )
  }

  is_fun <- sapply(input, is.function)
  input[is_fun] <- lapply(
    input[is_fun], function(fun) {
      eval(as.call(c(fun, input[names(formals(fun))])))
    }
  )

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

  temporarilyAssignInGlobalEnv(nfInst)

  if (!is.null(param$dir)) dir <- param$dir
  CnfInst <- param$CnfInst ## user provided compiled nimbleFunction?
  if (is.null(CnfInst)) {
    if (verbose) cat("## Compiling nimbleFunction \n")
    CnfInst <- wrap_if_true(isTRUE(param$knownFailures$compilation), expect_error, {
      compileNimble(nfInst, dirName = dir, projectName = nimbleProject_name)
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
  } else if (is.null(CnfInst)) { ## compilation failure was in knownFailures
    if (verbose) cat("## Compilation failed, as expected \n")
  } else {
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
    tol1 <- if (is.null(param$tol1)) 0.00001 else param$tol1
    tol2 <- if (is.null(param$tol2)) 0.0001 else param$tol2
    for (method_name in names(param$methods)) {
      if (verbose) {
        cat(paste0(
          "## Testing ", method_name, ': ',
          paste0(param$wrts[[method_name]], collapse = ', '),
          '\n'
        ))
      }
      value_test <- wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$value),
        param$knownFailures[[method_name]]$value, {
          if (verbose) cat("## Checking values\n")
          expect_equal(
            Cderivs[[method_name]]$value,
            Rderivs[[method_name]]$value
          )
          if ('log' %in% names(opParam$args)) {
            if (verbose) cat("## Checking log behavior for values\n")
            expect_equal(
              Cderivs2[[method_name]]$value,
              Rderivs2[[method_name]]$value
            )
            expect_false(isTRUE(all.equal(
              Rderivs[[method_name]]$value, Rderivs2[[method_name]]$value
            )))
            expect_false(isTRUE(all.equal(
              Cderivs[[method_name]]$value, Cderivs2[[method_name]]$value
            )))
          }
        }, wrap_in_try = isTRUE(catch_failures)
      )
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
      } else if (!is.null(param$knownFailures[[method_name]]$value)) {
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
      jacobian_test <- wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$jacobian),
        param$knownFailures[[method_name]]$jacobian, {
          if (verbose) cat("## Checking jacobians\n")
          expect_equal(
            Cderivs[[method_name]]$jacobian,
            Rderivs[[method_name]]$jacobian,
            tolerance = tol1
          )
          if ('log' %in% names(opParam$args)) {
            if (verbose) cat("## Checking log behavior for jacobians\n")
            expect_equal(
              Cderivs2[[method_name]]$jacobian,
              Rderivs2[[method_name]]$jacobian,
              tolerance = tol1
            )
            expect_false(isTRUE(all.equal(
              Rderivs[[method_name]]$jacobian, Rderivs2[[method_name]]$jacobian
            )))
            expect_false(isTRUE(all.equal(
              Cderivs[[method_name]]$jacobian, Cderivs2[[method_name]]$jacobian
            )))
          }
        }, wrap_in_try = isTRUE(catch_failures)
      )
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
      } else if (!is.null(param$knownFailures[[method_name]]$jacobian)) {
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
      hessian_test <- wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$hessian),
        param$knownFailures[[method_name]]$hessian, {
          if (verbose) cat("## Checking hessians\n")
          expect_equal(
            Cderivs[[method_name]]$hessian,
            Rderivs[[method_name]]$hessian,
            tolerance = tol2
          )
          if ('log' %in% names(opParam$args)) {
            if (verbose) cat("## Checking log behavior for hessians\n")
            expect_equal(
              Cderivs2[[method_name]]$hessian,
              Rderivs2[[method_name]]$hessian,
              tolerance = tol2
            )
            expect_false(isTRUE(all.equal(
              Rderivs[[method_name]]$hessian, Rderivs2[[method_name]]$hessian
            )))
            expect_false(isTRUE(all.equal(
              Cderivs[[method_name]]$hessian, Cderivs2[[method_name]]$hessian
            )))
          }
        }, wrap_in_try = isTRUE(catch_failures)
      )
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
      } else if (!is.null(param$knownFailures[[method_name]]$hessian)) {
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
  else if(!isTRUE(param$knownFailures$compilation)) {
    nimble:::clearCompiled(CnfInst)
    invisible(list(input = input))
  }
  invisible(NULL)
}

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

## Make a test parameterization to be used by test_AD. This method is primarily
## used by make_AD_test_batch() and make_distribution_fun_AD_test().
##
## op: Character string, the operator that will be the focus of the test.
## argTypes:   Character vector of argType strings that, when parsed, can be
##             passed to argType2symbol. If named, the names will as the formals
##             of the nimbleFunction generator's run method and other methods.
##             If not, formals are generated as arg1, arg2, etc.
## wrt_args:   Optional character vector of args to use in make_wrt(). If NULL,
##             assumes that all the arguments should be used.
## arg_distns: A list of input generation functions which is simply passed to
##             the output list. Should be NULL (use defaults found in
##             arg_type_2_input()), length 1 (use same input gen mechanism for each
##             argType, or a named list with names from among the argType names
##             (possibly the sequentially generated names). This will be NULL
##             when bulk generating the test params using make_AD_test_batch and
##             added later via modify_on_match(). Used in the call to
##             make_AD_test() in make_distribution_fun_AD_test(). 
## more_args:  A named list of additional fixed arguments to use in the
##             generated operator call. E.g., if op = 'dnorm',
##             argTypes = c('double(1, 5)', 'double(0)'), and
##             more_args = list(log = 1), the call to make_op_param will include
##             the expression dnorm(arg1, arg2, log = 1).
## seed:       A seed to use in set.seed().
##
## returns: A list with the following elements:
##          name         : character string, the operator
##          opParam      : result from make_op_param
##          run          : result from calling gen_runFunCore with opParam
##          methods      : a list of nimbleFunction method expressions that are
##                         the calls to nimDerivs() of the run method, each with
##                         a different wrt argument
##          enableDerivs : list('run')
##          wrts         : a list of character vectors, each of which is the wrt
##                         argument for the corresponding method in methods
##          arg_distns   : A list of random input generation functions to be
##                         used by arg_type_2_input(). 
make_AD_test <- function(op, argTypes, wrt_args = NULL,
                         arg_distns = NULL, more_args = NULL, seed = 0) {
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
    enableDerivs = list('run'),
    wrts = wrts,
    arg_distns = arg_distns
  )
}

## ops: character vector of operator names
## argTypes: list of character vectors of argTypes
##           e.g. for a binary operator:
##             list(
##               c('double(1, 4)', 'double(0)'),
##               c('double(1, 4)', 'double(1, 4)')
##             )
make_AD_test_batch <- function(ops, argTypes, seed = 0) {
  opTests <- unlist(
    recursive = FALSE,
    x = lapply(
      ops,
      function(x) {
        mapply(
          make_AD_test,
          argTypes = argTypes,
          MoreArgs = list(op = x, seed = seed),
          SIMPLIFY = FALSE
        )
      })
  )
  names(opTests) <- sapply(opTests, `[[`, 'name')
  invisible(opTests)
}

##################
## unary cwise ops
##################

## could instead use an inverted version of
## nimble:::specificCallReplacements in test_AD
gammafn <- gamma
lgammafn <- lgamma
ceil <- ceiling
ftrunc <- trunc

unaryArgs <- c('double(0)', 'double(1, 4)', 'double(2, c(3, 4))')
unaryOps <- c(
  '-', nimble:::unaryDoubleOperators,
  nimble:::unaryPromoteNoLogicalOperators
)

unaryOpTests <- make_AD_test_batch(
  unaryOps, unaryArgs
)

## tranform args
modify_on_match(unaryOpTests, '(log|sqrt) .+', 'arg_distns', function(x) abs(rnorm(x)))
modify_on_match(unaryOpTests, 'log1p .+', 'arg_distns', function(x) abs(rnorm(x)) - 1)
modify_on_match(unaryOpTests, '(logit|probit|cloglog) .+', 'arg_distns', runif)
modify_on_match(unaryOpTests, '(acos|asin|atanh) .+', 'arg_distns', function(x) runif(x, -1, 1))
modify_on_match(unaryOpTests, 'acosh .+', 'arg_distns', function(x) abs(rnorm(x)) + 1)

## see knownFailure below
modify_on_match(unaryOpTests, '^(gammafn|factorial) .+', 'arg_distns', function(x) abs(rnorm(x)))

## R implements lgamma as "the natural logarithm of _the absolute value of_ the
## gamma function", and since we use exp(lgamma(x)) for gammafn(x) we get the
## sign of the output flipped when gamma(x) < 0.
modify_on_match(
  unaryOpTests, '^(gammafn|factorial) .*', 'knownFailures',
  list(
    compilation = FALSE, ## no need to set this, but showing here by way of example
    input = list( # currently doesn't do anything, just a reminder
      cpp = -0.245, # any input x where gamma(x) < 0
      R = NULL ## if the R version had errors too, we could say so here
    )
  )
)

## the following unaryExpr operators are not currently supported
modify_on_match(
  unaryOpTests, '(nimRound|ftrunc|ceil|floor) .*',
  'knownFailures',
  list(
    compilation = TRUE
  )
)

## abs fails on negative inputs
## modify_on_match(unaryOpTests, 'abs .+', 'arg_distns', function(x) -abs(rnorm(x)))

## e.g. of how to set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
## modify_on_match(unaryOpTests, 'ilogit .+', 'tol2', 0.1)

######################
## unary reduction ops
######################

squaredNorm <- function(x) sum(x^2)

unaryReductionArgs <- c('double(1, 4)', 'double(2, c(3, 4))')

unaryReductionOps <- c(
  'sum', nimble:::reductionUnaryDoubleOperatorsEither,
  nimble:::reductionUnaryOperatorsArray
)

unaryReductionOpTests <- make_AD_test_batch(
  unaryReductionOps, unaryReductionArgs
)

## nimble doesn't support var of matrices, see sizeUnaryReduction
modify_on_match(
  unaryReductionOpTests, 'var double\\(2, c\\(3, 4\\)\\)',
  'knownFailures',
  list(
    compilation = TRUE
  )
)

#############
## binary ops
#############

## does not include combinations of vector and matrix
binaryArgs <- as.list(
  cbind(
    data.frame(t(expand.grid(unaryArgs[1], unaryArgs)), stringsAsFactors=FALSE),
    data.frame(t(expand.grid(unaryArgs[2:3], unaryArgs[1])), stringsAsFactors=FALSE)
  )
)
names(binaryArgs) <- NULL
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[2], 2)
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[3], 2)

binaryOps <- c(
  nimble:::binaryOrUnaryOperators,
  '/', '*', '%%'
)

binaryOpTests <- make_AD_test_batch(
  binaryOps, binaryArgs
)

modify_on_match(
  binaryOpTests, '%% .*',
  'knownFailures',
  list(
    compilation = TRUE
  )
)

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
##modify_on_match(binaryOpTests, '(\\+|-) double\\(0\\) double\\(1, 4\\)', 'tol2', 0.001)
##modify_on_match(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(0\\)', 'tol2', 0.001)
##modify_on_match(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(1, 4\\)', 'tol2', 0.001)

## runtime failures

## example of specifying when a particular method fails:
## modify_on_match(
##   binaryOpTests,
##   '\\+ double\\(0\\) double\\(0\\)',
##   'knownFailures',
##   list(
##     method3 = list( ## arg1, arg2
##       jacobian = expect_failure,
##       hessian = expect_failure
##     ),
##     method4 = list( ## no wrt
##       jacobian = expect_failure,
##       hessian = expect_failure
##     )
##   )
## )

###############
## pow-like ops
###############

powArgs <- list(
  c('double(0)', 'double(0)'),
  c('double(1, 4)', 'double(0)')
)

powOps <- c(
  'pow', '^'
)

powOpTests <- make_AD_test_batch(
  powOps, powArgs
)

#######################
## binary reduction ops
#######################

binaryReductionArgs <- list(
  c('double(1, 4)', 'double(1, 4)')
)

binaryReductionOps <- nimble:::reductionBinaryOperatorsEither

binaryReductionOpTests <- make_AD_test_batch(
  binaryReductionOps, binaryReductionArgs
)

##########################
## unary square matrix ops
##########################

squareMatrixArgs <- list('double(2, c(2, 2))', 'double(2, c(5, 5))')

squareMatrixOps <- c(nimble:::matrixSquareOperators,
                     nimble:::matrixSquareReductionOperators)

squareMatrixOpTests <- make_AD_test_batch(
  squareMatrixOps, squareMatrixArgs
)

gen_pos_def_matrix <- function(arg_size) {
  m <- arg_size[1] ## assumes matrix argType is square
  mat <- diag(1:m)
  mat[lower.tri(mat)] <- runif(m*(m - 1)/2)
  mat %*% t(mat)
}

modify_on_match(
  squareMatrixOpTests, 'chol .+', 'arg_distns',
  gen_pos_def_matrix
)

## ## compilation failures
## modify_on_match(
##   squareMatrixOpTests, '(logdet|det) .+',
##   'knownFailures', list(compilation = TRUE)
## )

####################
## binary matrix ops
####################

binaryMatrixOps <- nimble:::matrixMultOperators

binaryMatrixArgs <- as.list(
  data.frame(t(expand.grid(
    c('double(2, c(3, 4))', 'double(2, c(4, 4))', 'double(2, c(1, 4))'),
    c('double(2, c(4, 3))', 'double(2, c(4, 4))', 'double(2, c(4, 1))'))),
    stringsAsFactors=FALSE)
)

binaryMatrixOpTests <- make_AD_test_batch(
  binaryMatrixOps, binaryMatrixArgs
)

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
## modify_on_match(binaryMatrixOpTests, '%*% .+', 'tol2', 0.001)

#########################
## distribution functions
#########################

## create a base distribution specification which will
## be extended to specify the actual distributions
distn_base = list(
  ##
  ## name of distn_param entry should be abbreviated distribution name
  ## as listed in Rdist field of distributions_inputList.R (e.g. dnorm -> norm)
  ##
  name = '',
  ##
  ## variants should create a valid distribution function name
  ## when prepended to name
  ##
  variants = c('d'),

  args = list(
    ##
    ## the rand_variate arg must be included
    ##
    rand_variate = list( ## TODO: better naming convention
      ##
      ## distn is a function that should generate a number on the support,
      ## or return a function that uses the other args to generate random
      ## variates from this distribution. In the latter case, the formals of the
      ## returned function should be from among the other arg names. See the
      ## definition of binom for an example of how to use the parameters of the
      ## distribution.
      distn = NULL,
      ##
      ## The type field is required for all args.
      ##
      type = NULL
    )
    ##
    ## log is a special arg that will trigger a test to check that
    ## the _logFixed variant is not incorrectly used when included.
    ## See distributions_R.hpp. First we test with a fixed log arg
    ## in distn_tests but later we'll add a log argument to each
    ## distn_params entry in distn_with_log_tests, so don't add it here.
    ##
    ## log = list(
    ##   distn = function(n) sample(0:1, size = n, replace = TRUE),
    ##   type = 'double(0)'
    ## )
    ##
  ),
  ##
  ## You may want to pass in more_args to distribution function expr,
  ## but not as one of the nimbleFunction args, in which case they
  ## can be included as follows:
  ##
  ## more_args = list(
  ##   d = list(log = 1),
  ##   p = list(log.p = 1)
  ## ),
  ##
  ## wrt should be a character vector and can include arg names
  ## and any distribution function first argument names among the
  ## choices 'x', 'q', and 'p'.
  ##
  wrt = character(0)
)

## all distribution specifications will be added to this list and
## converted into an AD test parameterization by make_distribution_fun_AD_test()
distn_params <- list()

##
## Binomial distribution
##

distn_params[['binom_base']] <- distn_base
distn_params[['binom_base']]$name <- 'binom'
distn_params[['binom_base']]$args <- list(
  rand_variate = list(
    ##
    ## `size` here refers to the size parameter of the binomial distribution.
    ## Since the support depends on a parameter, return a function which
    ## test_AD will use with the realized value of the parameter `size`.
    ##
    distn = function(n) function(size, prob) rbinom(n, size, prob, replace = TRUE),
    type = c('double(0)', 'double(1, 5)')
  ),
  size = list(
    distn = function(n) sample(1:100, size = n, replace = TRUE),
    type = c('double(0)', 'double(1, 3)')
  ),
  prob = list(
    distn = function(n) runif(n),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['binom_base']]$wrt <- c('prob')

##
## Categorical distribution
##

distn_params[['cat_base']] <- distn_base
distn_params[['cat_base']]$name <- 'cat'
distn_params[['cat_base']]$args <- list(
  rand_variate = list(
    ## support depends on k = length(prob)
    distn = function(n) function(prob) sample(1:length(prob), n, replace = TRUE),
    type = c('double(0)')
  ),
  prob = list(
    distn = function(n) {prob <- runif(n); prob/sum(prob)},
    type = c('double(1, 5)', 'double(1, 8)')
  )
)
distn_params[['cat_base']]$wrt <- c('prob')

##
## Multinomial distribution
##

## no size arg
distn_params[['multi_no_size']] <- distn_base
distn_params[['multi_no_size']]$name <- 'multi'
distn_params[['multi_no_size']]$args <- list(
  rand_variate = list(
    distn = function(n) sample(1:10, n, replace = TRUE),
    type = c('double(1, 5)')
  ),
  prob = list(
    distn = function(n) {prob <- runif(n); prob/sum(prob)},
    type = c('double(1, 5)')
  )
)
distn_params[['multi_no_size']]$wrt = c('prob')

## including size
distn_params[['multi_with_size']] <- distn_params[['multi_no_size']]
distn_params[['multi_with_size']]$args[['rand_variate']]$distn <-
  function(n) function(size, prob) nimble::rmulti(n = 1, size, prob)
distn_params[['multi_with_size']]$args[['size']]  <- list(
  distn = function(n) sample(1:1000, n, replace = TRUE),
  type = c('double(1, 5)')
)

##
## Negative Binomial distribution
##

distn_params[['nbinom_base']] <- distn_base
distn_params[['nbinom_base']]$name <- 'nbinom'
distn_params[['nbinom_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(size, prob) rnbinom(n, size, prob),
    type = c('double(0)', 'double(1, 5)')
  ),
  size = list(
    distn = function(n) runif(n, max = 100),
    type = c('double(0)', 'double(1, 3)')
  ),
  prob = list(
    distn = function(n) runif(n),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['nbinom_base']]$wrt <- c('size', 'prob')

##
## Poisson distribution
##

distn_params[['pois_base']] <- distn_base
distn_params[['pois_base']]$name <- 'pois'
distn_params[['pois_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(lambda) rpois(n, lambda),
    type = c('double(0)', 'double(1, 5)')
  ),
  lambda = list(
    distn = function(n) runif(n, max = 100),
    type = c('double(0)', 'double(1, 3)')
  )
)
distn_params[['pois_base']]$wrt <- c('lambda')

##
## Beta distribution
##

distn_params[['beta_base']] <- distn_base
distn_params[['beta_base']]$name <- 'beta'
distn_params[['beta_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(shape1, shape2) rbeta(n, shape1, shape2),
    type = c('double(0)', 'double(1, 5)')
  ),
  shape1 = list(
    distn = function(n) runif(n, max = 10),
    type = c('double(0)', 'double(1, 3)')
  ),
  shape2 = list(
    distn = function(n) runif(n, max = 10),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['beta_base']]$wrt <- c('shape1', 'shape2', 'x')

##
## Chi-squared distribution
##

distn_params[['chisq_base']] <- distn_base
distn_params[['chisq_base']]$name <- 'chisq'
distn_params[['chisq_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(df) rchisq(n, df),
    type = c('double(0)', 'double(1, 5)')
  ),
  df = list(
    distn = function(n) runif(n, max = 100),
    type = c('double(0)', 'double(1, 3)')
  )
)
distn_params[['chisq_base']]$wrt <- c('x')

##
## Double Exponential (Laplace) distribution
##

## Here's an example of creating a distribution base
## specification and adding alternative parameterizations.

dexp_base <- distn_base
dexp_base$name <- 'dexp'
dexp_base$args <- list(
  rand_variate = list(
    distn = function(n) runif(n, min = -10, max = 10),
    type = c('double(0)', 'double(1, 5)')
  ),
  location = list(
    distn = function(n) runif(n, min = -10, max = 10),
    type = c('double(0)', 'double(1, 3)')
  )
)
dexp_base$wrt <- c('location', 'x')

## scale parameterization
distn_params[['dexp_scale']] <- dexp_base
distn_params[['dexp_scale']]$args[['scale']] <- list(
  distn = function(n) runif(n, max = 10),
  type = c('double(0)', 'double(1, 7)')
)
distn_params[['dexp_scale']]$wrt <- c(dexp_base$wrt, 'scale')

## rate parameterization
distn_params[['dexp_rate']] <- dexp_base
distn_params[['dexp_rate']]$args[['rate']] <-
  distn_params[['dexp_scale']]$args[['scale']]
distn_params[['dexp_rate']]$wrt <- c(dexp_base$wrt, 'rate')

##
## Exponential distribution
##

exp_base <- distn_base
exp_base$args = list(
  rand_variate = list(
    type = c('double(0)')
  )
)
exp_base$wrt <- c('x')

## base R version
distn_params[['exp_R']] <- exp_base
distn_params[['exp_R']]$name <- 'exp'
distn_params[['exp_R']]$args[['rand_variate']]$distn <- function(n) {
  function(rate) rexp(n, rate)
}
distn_params[['exp_R']]$args[['rate']] = list(
  distn = function(n) runif(n, max = 100),
  type = c('double(0)', 'double(1, 3)')
)
distn_params[['exp_R']]$wrt <- c(exp_base$wrt, 'rate')

## NIMBLE version, rate parameterization
distn_params[['exp_nimble_rate']] <- distn_params[['exp_R']]
distn_params[['exp_nimble_rate']]$name <- 'exp_nimble'

## NIMBLE version, scale parameterization
distn_params[['exp_nimble_scale']] <- exp_base
distn_params[['exp_nimble_scale']]$name <- 'exp_nimble'
distn_params[['exp_nimble_scale']]$args[['rand_variate']]$distn <- function(n) {
  function(rate) rexp(n, 1/scale)
}
distn_params[['exp_nimble_scale']]$args[['scale']]  <- list(
  distn = function(n) runif(n, max = 10),
  type = c('double(0)', 'double(1, 3)')
)
distn_params[['exp_nimble_scale']]$wrt <- c(exp_base$wrt, 'scale')

##
## Gamma distribution
##

distn_params[['gamma_scale']] <- distn_params[['exp_nimble_scale']]
distn_params[['gamma_rate']] <- distn_params[['exp_nimble_rate']]

## change the name
distn_params[['gamma_scale']]$name <-
  distn_params[['gamma_rate']]$name <- 'gamma'

## add the shape parameter
distn_params[['gamma_scale']]$args[['shape']] <-
  distn_params[['gamma_rate']]$args[['shape']] <- list(
    distn = function(n) runif(n, max = 10),
    type = c('double(0)', 'double(1, 7)')
  )

## add shape arg to wrt
distn_params[['gamma_scale']]$wrt <- c(
  distn_params[['gamma_scale']]$wrt, 'shape'
)
distn_params[['gamma_rate']]$wrt <- c(
  distn_params[['gamma_rate']]$wrt, 'shape'
)

## change how random variates are generated
distn_params[['gamma_rate']]$args[['rand_variate']]$distn <- function(n) {
  function(shape, rate) rgamma(n, shape, rate)
}
distn_params[['gamma_scale']]$args[['rand_variate']]$distn <- function(n) {
  function(shape, scale) rgamma(n, shape, scale=scale)
}

##
## Inverse Gamma distribution
##

distn_params[['invgamma_scale']] <- distn_params[['gamma_scale']]
distn_params[['invgamma_rate']] <- distn_params[['gamma_rate']]

distn_params[['invgamma_scale']]$name <-
  distn_params[['invgamma_rate']]$name <- 'invgamma'

##
## Sqrt Inverse Gamma distribution
##

distn_params[['sqrtinvgamma_scale']] <- distn_params[['gamma_scale']]
distn_params[['sqrtinvgamma_rate']] <- distn_params[['gamma_rate']]

distn_params[['sqrtinvgamma_scale']]$name <-
  distn_params[['sqrtinvgamma_rate']]$name <- 'sqrtinvgamma'

##
## Log-normal distribution
##

distn_params[['lnorm_base']] <- distn_base
distn_params[['lnorm_base']]$name <- 'lnorm'
distn_params[['lnorm_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(meanlog, sdlog) rlnorm(n, meanlog, sdlog),
    type = c('double(0)', 'double(1, 5)')
  ),
  meanlog = list(
    distn = function(n) runif(n, min = -100, max = 100),
    type = c('double(0)', 'double(1, 3)')
  ),
  sdlog = list(
    distn = function(n) runif(n, max = 100),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['lnorm_base']]$wrt <- c('meanlog', 'sd', 'x')

##
## Logistic distribution
##

distn_params[['logis_base']] <- distn_params[['dexp_scale']]
distn_params[['logis_base']]$name <- 'logis'
distn_params[['logis_base']]$args[['rand_variate']]$distn <- function(n) {
  function(location, scale) rlogis(n, location, scale)
}


##
## Normal distribution
##

distn_params[['norm_base']] <- distn_base
distn_params[['norm_base']]$name <- 'norm'
distn_params[['norm_base']]$variants = c('d')
distn_params[['norm_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(mean, sd) rnorm(n, mean, sd),
    type = c('double(0)', 'double(1, 5)')
  ),
  mean = list(
    distn = function(n) runif(n, min = -100, max = 100),
    type = c('double(0)', 'double(1, 3)')
  ),
  sd = list(
    distn = function(n) runif(n, max = 10),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['norm_base']]$wrt <- c('mean', 'sd', 'x')

##
## t-distribution
##

## usual R version
distn_params[['t_R']] <- distn_params[['chisq_base']]
distn_params[['t_R']]$name <- 't'
distn_params[['t_R']]$args[['rand_variate']]$distn <-
  function(n) function(df) rt(n, df)

## non-standard version
distn_params[['t_nonstandard']] <- distn_params[['t_R']]
distn_params[['t_nonstandard']]$name <- 't_nonstandard'

distn_params[['t_nonstandard']]$args[['mu']] <-
  distn_params[['norm_base']]$args[['mean']]

distn_params[['t_nonstandard']]$args[['sigma']] <-
  distn_params[['norm_base']]$args[['sd']]

distn_params[['t_nonstandard']]$wrt <- c(
  distn_params[['t_nonstandard']]$wrt, 'mu', 'sigma'
)

distn_params[['t_nonstandard']]$args[['rand_variate']]$distn <- function(n) {
  function(mu, sigma) rnorm(n, mu, sigma)
}

##
## Uniform distribution
##

distn_params[['unif_base']] <- distn_base
distn_params[['unif_base']]$name <- 'unif'

distn_params[['unif_base']]$args <- list(
  rand_variate = list(
    distn = function(n) function(min, max) runif(n, min, max),
    type = c('double(0)', 'double(1, 5)')
  ),
  min = list(
    distn = function(n) runif(n, min = -1000, max = 100),
    type = c('double(0)', 'double(1, 3)')
  ),
  max = list(
    distn = function(n) runif(n, min = 100, max = 1000),
    type = c('double(0)', 'double(1, 7)')
  )
)
distn_params[['unif_base']]$wrt <- c('min', 'max', 'x')

##
## Weibull distribution
##

distn_params[['weibull_base']] <- distn_params[['gamma_scale']]
distn_params[['weibull_base']]$name <- 'weibull'
distn_params[['weibull_base']]$args[['rand_variate']]$distn <- function(n) {
  function(shape, scale) rweibull(n, shape, scale)
}

##
## Dirichlet distribution
##

distn_params[['dirch_base']] <- distn_base
distn_params[['dirch_base']]$name <- 'dirch'
distn_params[['dirch_base']]$args <- list(
  rand_variate = list(
    distn = function(n) {prob <- runif(n); prob/sum(prob)},
    type = c('double(1, 5)')
  ),
  alpha = list(
    distn = function(n) runif(n, max = 1000),
    type = c('double(1, 5)')
  )
)
distn_params[['dirch_base']]$wrt = c('alpha', 'x')

##
## Multivariate Normal distribution
##

chol_base <- distn_base
chol_base$args <- list(
  cholesky = list(
    distn = gen_pos_def_matrix,
    type = c('double(2, c(5, 5))')
  )
)

distn_params[['mnorm_chol_base']] <- chol_base
distn_params[['mnorm_chol_base']]$name <- 'mnorm_chol'
distn_params[['mnorm_chol_base']]$args[['rand_variate']] <- list(
  distn = function(n) function(mean, cholesky)
    rmnorm_chol(n = 1, mean, cholesky),
  type = c('double(1, 5)')
)
distn_params[['mnorm_chol_base']]$args[['mean']] <- list(
  distn = function(n) runif(n, min = -10, max = 10),
  type = c('double(1, 5)')
)
distn_params[['mnorm_chol_base']]$wrt <- c('mean', 'cholesky', 'x')

##
## Multivariate t-distribution
##

distn_params[['mvt_chol_base']] <- chol_base
distn_params[['mvt_chol_base']]$name <- 'mvt_chol'
distn_params[['mvt_chol_base']]$args[['rand_variate']] <- list(
  distn = function(n) function(mu, cholesky, df)
    rmvt_chol(n = 1, mu, cholesky, df),
  type = c('double(1, 5)')
)
distn_params[['mvt_chol_base']]$args[['mu']] <- list(
  distn = function(n) runif(n, min = -10, max = 10),
  type = c('double(1, 5)')
)
distn_params[['mvt_chol_base']]$args[['df']] <- list(
  distn = function(n) runif(n, max = 100),
  type = c('double(0)')
)
distn_params[['mvt_chol_base']]$wrt <- c('mu', 'cholesky', 'x')

##
## Wishart distribution
##

distn_params[['wish_chol_base']] <- chol_base
distn_params[['wish_chol_base']]$name <- 'wish_chol'
distn_params[['wish_chol_base']]$args[['rand_variate']] <- list(
  distn = function(n) function(cholesky, df)
    rwish_chol(n = 1, cholesky, df),
  type = c('double(1, 5)')
)
distn_params[['wish_chol_base']]$args[['df']] <- list(
  distn = function(n) runif(n, max = 100),
  type = c('double(0)')
)
distn_params[['wish_chol_base']]$wrt <- c('cholesky', 'x')

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
##              Additional args must also have distn and type fields.
## more_args:   Passed to make_op_param().
##
make_distribution_fun_AD_test <- function(distn_param) {
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
  param_distns <- lapply(distn_param$args[-rand_variate_idx], `[[`, 'distn')
  arg_distns <- sapply(distn_param$variants, function(variant) {
    arg1_distn <- switch(
      variant,
      d = list(x = distn_param$args$rand_variate$distn),
      p = list(q = distn_param$args$rand_variate$distn),
      q = list(p = runif)
    )
    c(arg1_distn, param_distns)
  }, simplify = FALSE)
  op_params <- unlist(
    lapply(
      distn_param$variants, function(variant) {
        lapply(
          argTypes[[variant]],
          function(these_argTypes) {
            wrt_args = intersect(
              distn_param$wrt, names(these_argTypes)
            )
            make_AD_test(
              ops[[variant]], these_argTypes,
              wrt_args = wrt_args,
              arg_distns = arg_distns[[variant]],
              more_args = distn_param$more_args[[variant]]
            )
          }
        )
      }
    ), recursive = FALSE
  )
  names(op_params) <- sapply(op_params, `[[`, 'name')
  return(op_params)
}

## add more_args to include as part of distribution function expr
distn_params_log_1 <- lapply(distn_params, function(param) {
  more_args = lapply(
    param$variants, switch,
    d = list(log = 1),
    p = list(log.p = 1),
    q = list(log.p = 1)
  )
  names(more_args) <- param$variants
  param$more_args = more_args
  param
})

########################################
## distribution functions, fixed log arg
########################################

distn_tests <- unlist(
  lapply(distn_params_log_1, make_distribution_fun_AD_test),
  recursive = FALSE
)

## dexp is not processed by eigenize_recyclingRuleFunction()
distn_tests[['exp_R.dexp double(0) double(1, 3)']]$knownFailures <- list(
  compilation = TRUE
)

## scalar case: nimDerivs_dsqrtinvgamma is not defined
## vector recycling: log arg is not added to AST during compilation and as a
##                   result code generated by MAKE_RECYCLING_RULE_CLASS3_1scalar
##                   fails to compile
modify_on_match(
  distn_tests, 'dsqrtinvgamma',
  'knownFailures',
  list(
    compilation = TRUE
  )
)

## dt fails because of missing ncp argument
## dt_nonstandard fails because it has 3 parameters, but there is no
## MAKE_RECYCLING_RULE_CLASS4_1scalar macro yet
modify_on_match(
  distn_tests,
  '(dt|dt_nonstandard) (?!double\\(0\\) double\\(0\\))',
  'knownFailures',
  list(
    compilation = TRUE
  ),
  parent.frame(),
  perl = TRUE
)

modify_on_match(
  distn_tests,
  'dunif double\\(1, 5\\) double\\(1, 3\\) double\\(1, 7\\)',
  'knownFailures',
  list(
    method1 = list(
      value = expect_failure
    )
  )
)

###########################################
## distribution functions, variable log arg
###########################################

## add the arg 'log' to all the distn_params
distn_params_with_log <- lapply(distn_params, function(param) {
  param$args$log <- list(
    distn = function(n) sample(0:1, size = n, replace = TRUE),
    type = 'double(0)'
  )
  param
})

distn_with_log_tests <- unlist(
  lapply(distn_params_with_log, make_distribution_fun_AD_test),
  recursive = FALSE
)

distn_with_log_tests[['exp_R.dexp double(0) double(1, 3)']]$knownFailures <- list(
  compilation = TRUE
)

modify_on_match(
  distn_with_log_tests, 'dsqrtinvgamma',
  'knownFailures',
  list(
    compilation = TRUE
  )
)

#######################
## run all of the tests
#######################

invisible({
  sapply(unaryOpTests, test_AD)
  sapply(unaryReductionOpTests, test_AD)
  sapply(binaryOpTests, test_AD)
  sapply(powOpTests, test_AD)
  sapply(binaryReductionOpTests, test_AD)
  sapply(squareMatrixOpTests, test_AD)
  sapply(binaryMatrixOpTests, test_AD)
  sapply(distn_tests, test_AD)
  sapply(distn_with_log_tests, test_AD)
})
