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

test_AD <- function(param, dir = file.path(tempdir(), "nimble_generatedCode"),
                    control = list(), verbose = nimbleOptions('verbose'),
                    catch_failures = FALSE) {
  if (!is.null(param$debug) && param$debug) browser()
  if (verbose) cat(paste0("### Testing ", param$name, "\n"))

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
      input <- lapply(opParam$args, argType2input, param$arg_distns)
    else
      stop(
        'arg_distns of length greater than 1 must have names',
        call. = FALSE
      )
  else {
    input <- sapply(
      names(opParam$args),
      function(name)
        argType2input(opParam$args[[name]], param$arg_distns[[name]])
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
    warning(
      paste(
        'Calling R version of test', opParam$name,
        'resulted in an error:', Rderivs[1]
      ),
      immediate. = TRUE
    )
    return(invisible(NULL))
  }

  temporarilyAssignInGlobalEnv(nfInst)

  if (verbose) cat("## Compiling nimbleFunction \n")
  if (!is.null(param$dir)) dir <- param$dir
  CnfInst <- wrap_if_true(param$knownFailures$compilation, expect_error, {
    compileNimble(nfInst, dirName = dir)
  }, wrap_in_try = isTRUE(catch_failures))
  if (isTRUE(catch_failures) && inherits(CnfInst, 'try-error')) {
    warning(
      paste0(
        'The test of ', opParam$name,
        ' failed to compile.\n', CnfInst[1]
      ),
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
          immediate. = TRUE
        )
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
          immediate. = TRUE
        )
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
          immediate. = TRUE
        )
      }
    }
    nimble:::clearCompiled(CnfInst)
  }
  if (verbose) cat("### Test successful \n\n")
  invisible(NULL)
}

## Takes a named list of `argTypes` and returns a list of character
## vectors, each of which is valid as the `wrt` argument of `nimDerivs()`.
## Each argument on its own and all combinations of the arguments will
## always be included, and then make_wrt will try to create up to `n_random`
## additional character vectors with random combinations of the arguments
## and indexing of those arguments when possible (i.e. for non-scalar args).
make_wrt <- function(argTypes, n_random = 10) {

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

makeADtest <- function(op, argTypes, wrt_args = NULL,
                       arg_distns = NULL, more_args = NULL) {
  opParam <- makeOperatorParam(op, argTypes, more_args)

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

############
## unary ops
############

gammafn <- gamma
lgammafn <- lgamma
ceil <- ceiling

unaryArgs <- c('double(0)', 'double(1, 4)', 'double(2, c(3, 4))')
unaryOps <- c('-', 'sum',
              nimble:::unaryDoubleOperators,
              nimble:::unaryPromoteNoLogicalOperators,
              nimble:::reductionUnaryDoubleOperatorsEither,
              nimble:::reductionUnaryOperatorsArray)

unaryOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    unaryOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = unaryArgs,
        MoreArgs = list(op = x),
        SIMPLIFY = FALSE
      )
    })
)
names(unaryOpTests) <- sapply(unaryOpTests, `[[`, 'name')

## tranform args
modifyOnMatch(unaryOpTests, '(log|sqrt) .+', 'arg_distns', function(x) abs(rnorm(x)))
modifyOnMatch(unaryOpTests, 'log1p .+', 'arg_distns', function(x) abs(rnorm(x)) - 1)
modifyOnMatch(unaryOpTests, '(logit|probit|cloglog) .+', 'arg_distns', runif)
modifyOnMatch(unaryOpTests, '(acos|asin|atanh) .+', 'arg_distns', function(x) runif(x, -1, 1))
modifyOnMatch(unaryOpTests, 'acosh .+', 'arg_distns', function(x) abs(rnorm(x)) + 1)

## abs fails on negative inputs
modifyOnMatch(unaryOpTests, 'abs .+', 'arg_distns', function(x) -abs(rnorm(x)))

## e.g. of how to set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
## modifyOnMatch(unaryOpTests, 'ilogit .+', 'tol2', 0.1)

## compilation failures

modifyOnMatch(
  unaryOpTests,
  paste(
    '(logit|log1p|lfactorial|factorial|cloglog|atan|cosh|sinh|tanh|acosh|asinh|atanh)',
    '(double\\(1, 4\\)|double\\(2, c\\(3, 4\\)\\))'
  ),
  'knownFailures',
  list(compilation = TRUE)
)

modifyOnMatch(
  unaryOpTests, '(nimRound|ceil|floor|sd|var|probit|lgammafn|gammafn) .+',
  'knownFailures', list(compilation = TRUE)
)

ops_regex <- paste0(
  c('sum', nimble:::reductionUnaryDoubleOperatorsEither,
    nimble:::reductionUnaryOperatorsArray),
  collapse = '|'
)
modifyOnMatch(
  unaryOpTests,
  paste0('(', ops_regex, ') double\\(0\\)'),
  'knownFailures', list(compilation = TRUE)
)

set.seed(0)
sapply(unaryOpTests, test_AD)

#############
## binary ops
#############

binaryArgs <- as.list(
  cbind(
    data.frame(t(expand.grid(unaryArgs[1], unaryArgs)), stringsAsFactors=FALSE),
    data.frame(t(expand.grid(unaryArgs[2:3], unaryArgs[1])), stringsAsFactors=FALSE)
  )
)
names(binaryArgs) <- NULL
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[2], 2)
binaryArgs[[length(binaryArgs) + 1]] <- rep(unaryArgs[3], 2)

binaryOps <- c(nimble:::binaryOrUnaryOperators,
               nimble:::binaryMidOperators,
               nimble:::binaryLeftDoubleOperators,
               nimble:::reductionBinaryOperatorsEither)

binaryOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    binaryOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = binaryArgs,
        MoreArgs = list(op = x),
        SIMPLIFY = FALSE
      )
    })
)

names(binaryOpTests) <- sapply(binaryOpTests, `[[`, 'name')

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
##modifyOnMatch(binaryOpTests, '(\\+|-) double\\(0\\) double\\(1, 4\\)', 'tol2', 0.001)
##modifyOnMatch(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(0\\)', 'tol2', 0.001)
##modifyOnMatch(binaryOpTests, '(\\+|-) double\\(1, 4\\) double\\(1, 4\\)', 'tol2', 0.001)

## runtime failures

## example of specifying when a particular method fails:
modifyOnMatch(
  binaryOpTests,
  '\\+ double\\(0\\) double\\(0\\)',
  'knownFailures',
  list(
    method3 = list( ## arg1, arg2
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method4 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

sapply(binaryOpTests, test_AD)

##########################
## unary square matrix ops
##########################

squareMatrixArgs <- c('double(2, c(2, 2))', 'double(2, c(5, 5))')

squareMatrixOps <- c(nimble:::matrixSquareOperators,
                     nimble:::matrixSquareReductionOperators)

squareMatrixOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    squareMatrixOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = squareMatrixArgs,
        MoreArgs = list(op = x),
        SIMPLIFY = FALSE
      )
    })
)
names(squareMatrixOpTests) <- sapply(squareMatrixOpTests, `[[`, 'name')

modifyOnMatch(
  squareMatrixOpTests, 'chol .+', 'arg_distns',
  function(arg_size) {
    mat <- diag(1:arg_size[1]) ## assumes matrix is square
    mat[upper.tri(mat) | lower.tri(mat)] <- runif(1) ## chol only used upper tri
    mat
  }
)

## compilation failures
modifyOnMatch(
  squareMatrixOpTests, '(logdet|det) .+',
  'knownFailures', list(compilation = TRUE)
)

sapply(squareMatrixOpTests, test_AD)

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

binaryMatrixOpTests <- unlist(
  recursive = FALSE,
  x = lapply(
    binaryMatrixOps,
    function(x) {
      mapply(
        makeADtest,
        argTypes = binaryMatrixArgs,
        MoreArgs = list(op = x),
        SIMPLIFY = FALSE
      )
    })
)
names(binaryMatrixOpTests) <- sapply(binaryMatrixOpTests, `[[`, 'name')

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
modifyOnMatch(binaryMatrixOpTests, '%*% .+', 'tol2', 0.001)

sapply(binaryMatrixOpTests, test_AD)

#########################
## distribution functions
#########################

distn_params = list(
  list(
    ##
    ## name of distn_param entry should be abbreviated distribution name
    ## as used in distribution functions (e.g. dnorm -> norm)
    ##
    name = 'binom',
    ##
    ## variants should create a valid distribution function name
    ## when prepended to name
    ##
    variants = c('d'),
    args = list(
      ##
      ## the support arg must be included
      ##
      support = list( ## TODO: better naming convention
        ##
        ## This arg is required and should generate a number on the support,
        ## uniformly at random when possible.
        ## `size` here refers to the size parameter of the binomial distribution.
        ## Since the support depends on a parameter, return a function which
        ## test_AD will use with the realized value of the parameter `size`.
        ##
        distn = function(n) function(size) sample(1:size, n, replace = TRUE),
        ##
        ## The type field is required for all args.
        ##
        type = c('double(0)', 'double(1, 5)')
      ),
      size = list(
        ##
        ## For any arg (including support) distn is optional.
        ## See argType2input() for default argument generation distns.
        ##
        distn = function(n) sample(1:100, n, replace = TRUE),
        type = c('double(0)')##, 'double(1, 3)')
      ),
      prob = list(
        distn = function(n) runif(n),
        type = c('double(0)')##, 'double(1, 7)')
      )
      ## log is a special arg that will trigger a test to check that
      ## the _logFixed variant is not incorrectly used when included.
      ## See distributions_R.hpp.
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
    wrt = c('prob', 'p')
  )
)

## Takes an element of distn_params list and returns a list of AD test
## parameterizations, each of which test_AD can use.
##
## distn_param: Probability distribution parameterization, which
##              must have the following fields/subfields:
##              - name
##              - variants
##              - args
##                - support
##                  - type
##              Additional args must also have distn and type fields.
## more_args:   Passed to makeOperatorParam().
##
makeADdistnTest <- function(distn_param) {
  distn_name <- distn_param$name
  ops <- sapply(distn_param$variants, paste0, distn_name, simplify = FALSE)
  argTypes <- sapply(distn_param$variants, function(variant) {
    op <- paste0(variant, distn_name)
    support_type <- distn_param$args$support$type
    first_argType <- switch(
      variant,
      d = support_type,
      p = support_type,
      q = c('double(0)')##, 'double(1, 4)')
    )
    first_arg_name <- switch(variant, d = 'x', p = 'q', q = 'p')
    grid <- eval(as.call(c(
      expand.grid, list(first_argType),
      lapply(distn_param$args, `[[`, 'type')[-1]
    )))
    argTypes <- as.list(data.frame(t(grid), stringsAsFactors=FALSE))
    lapply(argTypes, function(v) {
      names(v) <- c(first_arg_name, names(distn_param$args)[-1])
      v
    })
  }, simplify = FALSE)
  param_distns <- lapply(distn_param$args[-1], `[[`, 'distn')
  arg_distns <- sapply(distn_param$variants, function(variant) {
    arg1_distn <- switch(
      variant,
      d = list(x = distn_param$args$support$distn),
      p = list(q = distn_param$args$support$distn),
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
            makeADtest(
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

distn_tests <- unlist(
  lapply(distn_params_log_1, makeADdistnTest),
  recursive = FALSE
)

lapply(distn_tests, test_AD)

#################################
## distribution functions log arg
#################################

## add the arg 'log' to all the distn_params
distn_params_with_log <- lapply(distn_params, function(param) {
  param$args$log <- list(
    distn = function(n) sample(0:1, size = n, replace = TRUE),
    type = 'double(0)'
  )
  param
})

distn_with_log_tests <- unlist(
  lapply(distn_params_with_log, makeADdistnTest),
  recursive = FALSE
)

lapply(distn_with_log_tests, test_AD)
