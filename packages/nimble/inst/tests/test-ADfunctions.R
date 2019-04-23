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
                    control = list(), verbose = nimbleOptions('verbose')) {
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
  if (is.null(param$arg_dists) || length(param$arg_dists) == 1)
    input <- lapply(opParam$args, argType2input, param$arg_dists)
  else {
    if (length(param$arg_dists) != length(opParam$args))
      stop(
        'arg_dists must be NULL, length 1, or the same length as the arguments',
        call. = FALSE
      )
    input <- mapply(argType2input, opParam$args, param$arg_dists)
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

  wrap_if_true(param$knownFailures$compilation, expect_error, {
    if (verbose) cat("## Compiling nimbleFunction \n")
    if (!is.null(param$dir)) dir <- param$dir
    CnfInst <- compileNimble(nfInst, dirName = dir)
    Cderivs <- sapply(names(param$methods), function(method) {
      do.call(
        ## same issue as with Rderivs
        eval(substitute(CnfInst$method, list(method = as.name(method)))),
        input
      )
    }, USE.NAMES = TRUE)
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
      wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$value),
        param$knownFailures[[method_name]]$value, {
          if (verbose) cat("## Checking values\n")
          expect_equal(
            Cderivs[[method_name]]$value,
            Rderivs[[method_name]]$value
          )
        }
      )
      wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$jacobian),
        param$knownFailures[[method_name]]$jacobian, {
          if (verbose) cat("## Checking jacobians\n")
          expect_equal(
            Cderivs[[method_name]]$jacobian,
            Rderivs[[method_name]]$jacobian,
            tolerance = tol1
          )
        }
      )
      wrap_if_true(
        !is.null(param$knownFailures[[method_name]]$hessian),
        param$knownFailures[[method_name]]$hessian, {
          if (verbose) cat("## Checking hessians\n")
          expect_equal(
            Cderivs[[method_name]]$hessian,
            Rderivs[[method_name]]$hessian,
            tolerance = tol2
          )
        }
      )
    }
    nimble:::clearCompiled(CnfInst)
  })
  if (verbose) cat("### Test successful \n\n")
  invisible(NULL)
}

## Takes an argType call and returns up to two components
## of the arg with respect to which to take derivatives
## in test_AD. One of the components will always be the
## argument itself and the other will be a random element
## within the arg (if the arg is non-scalar).
make_wrt <- function(argTypes, n_methods = 10) {

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

  while (n_methods > 0) {
    n_methods  <- n_methods - 1
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
        ## coin flip determines whether to subscript or not
        if (sample(c(TRUE, FALSE), 1))
          if (argSymbols[[i]]$nDim == 2)
            this_wrt <- c(
              this_wrt,
              paste0(
                names(args)[i], '[',
                sample(1:argSymbols[[i]]$size[1], size = 1) , ',',
                sample(1:argSymbols[[i]]$size[2], size = 1),
                ']'
              )
            )
          else ## nDim is 0 or 1
            this_wrt <- c(
              this_wrt,
              paste0(
                names(args)[i],
                '[', sample(1:argSymbols[[i]]$size, size = 1) , ']'
              )
            )
        else this_wrt <- c(this_wrt, names(args)[i])
      }
    }
    wrts <- c(wrts, list(unique(this_wrt)))
  }
  wrts <- unique(wrts)
}

makeADtest <- function(op, argTypes, wrt_args = NULL, arg_dists = NULL, more_args = NULL) {
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

  if (is.null(wrt_args)) wrt_args <- rep(TRUE, length(argTypes))
  wrts <- make_wrt(opParam$args[wrt_args])
  for (i in seq_along(wrts)) {
    method_i <- paste0('method', i)
    methods[[method_i]] <- method
    body(methods[[method_i]])[[2]][[3]][['wrt']] <- wrts[[i]]
  }

  method_no_wrt <- paste0('method', length(methods) + 1)
  methods[[method_no_wrt]] <- method

  wrts <- c(wrts, 'no wrt')
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
    arg_dists = arg_dists
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
modifyOnMatch(unaryOpTests, '(log|sqrt) .+', 'arg_dists', function(x) abs(rnorm(x)))
modifyOnMatch(unaryOpTests, 'log1p .+', 'arg_dists', function(x) abs(rnorm(x)) - 1)
modifyOnMatch(unaryOpTests, '(logit|probit|cloglog) .+', 'arg_dists', runif)
modifyOnMatch(unaryOpTests, '(acos|asin|atanh) .+', 'arg_dists', function(x) runif(x, -1, 1))
modifyOnMatch(unaryOpTests, 'acosh .+', 'arg_dists', function(x) abs(rnorm(x)) + 1)

## set tolerances (default is tol1 = 0.00001 and tol2 = 0.0001)
modifyOnMatch(unaryOpTests, 'ilogit .+', 'tol2', 0.1)
modifyOnMatch(unaryOpTests, '(factorial|exp|log1p|tan|acos|cos) .+', 'tol2', 0.001)

## runtime failures
## abs fails on negative inputs
modifyOnMatch(unaryOpTests, 'abs .+', 'arg_dists', function(x) -abs(rnorm(x)))
modifyOnMatch(
  unaryOpTests,
  'abs (double\\(1, 4\\)|double\\(2, c\\(3, 4\\)\\))',
  'knownFailure', '.*runs'
)

## wrt failures
modifyOnMatch(
  unaryOpTests, '.+ double\\(0\\)', 'knownFailures',
  list(
    method4 = list(
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)
modifyOnMatch(
  unaryOpTests, '.+ (double\\(1, 4\\)|double\\(2, c\\(3, 4\\)\\))',
  'knownFailures',
  list(
    method2 = list(
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method4 = list(
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

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
modifyOnMatch(
  binaryOpTests,
  '\\+ double\\(0\\) double\\(0\\)',
  'knownFailures',
  list(
    method5 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method6 = list( ## wrt arg1, arg1[1], arg2, arg2[1]
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

modifyOnMatch(
  binaryOpTests,
  '- double\\(0\\) double\\(0\\)',
  'knownFailures',
  list(
    method3 = list( ## wrt arg2
      jacobian = expect_failure
    ),
    method4 = list( ## wrt arg2[1]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method5 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method6 = list( ## wrt arg1, arg1[1], arg2, arg2[1]
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

modifyOnMatch(
  binaryOpTests,
  '\\+ double\\(0\\) double\\(1, 4\\)',
  'knownFailures',
  list(
    method3 = list( ## wrt arg2
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method4 = list( ## wrt arg2[*]
      jacobian = expect_failure
    ),
    method5 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method6 = list( ## wrt arg1, arg1[1], arg2, arg2[*]
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

modifyOnMatch(
  binaryOpTests,
  '\\+ double\\(1, 4\\) double\\(0\\)',
  'knownFailures',
  list(
    method2 = list( ## wrt arg1[*]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method3 = list( ## wrt arg2
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method4 = list( ## wrt arg2[1]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method5 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method6 = list( ## wrt arg1, arg1[*], arg2, arg2[1]
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

modifyOnMatch(
  binaryOpTests,
  '\\+ double\\(1, 4\\) double\\(1, 4\\)',
  'knownFailures',
  list(
    method2 = list( ## wrt arg1[*]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method4 = list( ## wrt arg2[*]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method5 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method6 = list( ## wrt arg1, arg1[*], arg2, arg2[*]
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

modifyOnMatch(
  binaryOpTests,
  '\\+ double\\(2, c\\(3, 4\\)\\) double\\(2, c\\(3, 4\\)\\)',
  'knownFailures',
  list(
    method2 = list( ## wrt arg1[*]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method4 = list( ## wrt arg2[*]
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method5 = list( ## no wrt
      jacobian = expect_failure,
      hessian = expect_failure
    ),
    method6 = list( ## wrt arg1, arg1[*], arg2, arg2[*]
      jacobian = expect_failure,
      hessian = expect_failure
    )
  )
)

## compilation failures
ops_regex <- paste0(
  c('\\+', '-'),
  collapse = '|'
)

modifyOnMatch(
  binaryOpTests,
  paste0('\\+ (',
         'double\\(2, c\\(3, 4\\)\\) double\\(0\\)|',
         'double\\(0\\) double\\(2, c\\(3, 4\\)\\))'),
  'knownFailures', list(compilation = TRUE)
)

modifyOnMatch(
  binaryOpTests,
  '/ double\\(0\\) double\\(2, c\\(3, 4\\)\\)',
  'knownFailures', list(compilation = TRUE)
)

sapply(binaryOpTests, test_AD)

##########################
## unary square matrix ops
##########################

squareMatrixArgs <- c('double(2, c(2, 2))', 'double(2, c(5, 5))')

squareMatrixOps <- c(nimble:::matrixSquareReductionOperators,
                     'inverse')

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

## compilation failures
modifyOnMatch(squareMatrixOpTests, '(logdet|det) .+', 'knownFailure', '.*compiles')

## runtime failures
modifyOnMatch(
  squareMatrixOpTests, 'inverse double\\(2, c\\(5, 5\\)\\)',
  'knownFailure', '.*compiles'
)

sapply(squareMatrixOpTests, test_AD, info = 'unary square matrix')

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

sapply(binaryMatrixOpTests, test_AD, info = 'binary matrix')

#########################
## distribution functions
#########################

dist_params = list(
  list(
    ## name of dist_param entry should be abbreviated distribution name
    ## as used in distribution functions (e.g. dnorm -> norm)
    name = 'binom',
    ## variants should create a valid distribution function name
    ## when prepended to name
    variants = c('d', 'p', 'q'),
    args = list(
      ## the support arg must be included
      support = list( ## TODO: better naming convention
        ## This arg is required and should generate a number on the support,
        ## uniformly at random when possible.
        ## `size` here refers to the size parameter of the binomial distribution.
        ## Since the support depends on a parameter, return a function which
        ## test_AD will use with the realized value of the parameter `size`.
        dist = function(n) function(size) sample(1:size, n, replace = TRUE),
        ## The type field is required for all args.
        type = c('double(0)', 'double(1, 5)')
      ),
      size = list(
        ## For any arg (including support) dist is optional.
        ## See argType2input() for default argument generation dists.
        dist = function(n) sample(1:100, n, replace = TRUE),
        type = c('double(0)', 'double(1, 3)')
      ),
      prob = list(
        dist = function(n) runif(n),
        type = c('double(0)', 'double(1, 7)')
      )
    ),
    ## wrt should be a character vector and can include arg names
    ## and any distribution function first argument names among the
    ## choices 'x', 'q', and 'p'.
    wrt = c('prob', 'p')
  )##,
  ##list(
  ##  name = 'beta',
  ##  variants = c('d', 'p', 'q', 'r'), ## prepended to name of this list
  ##  support = c(-Inf, Inf),
  ##  parameters = list(
  ##    mean = c('double(0)', 'double(1, 4)'),
  ##    sd = c('double(0)', 'double(1, 4)')
  ##  )
  ##),
  ##list(
  ##  name = 'chisq',
  ##  variants = c('d', 'p', 'q', 'r'), ## prepended to name of this list
  ##  support = c(-Inf, Inf),
  ##  parameters = list(
  ##    mean = c('double(0)', 'double(1, 4)'),
  ##    sd = c('double(0)', 'double(1, 4)')
  ##  )
  ##),
  ##list(
  ##  name = 'dexp',
  ##  variants = c('d', 'p', 'q', 'r'), ## prepended to name of this list
  ##  support = c(-Inf, Inf),
  ##  parameters = list(
  ##    mean = c('double(0)', 'double(1, 4)'),
  ##    sd = c('double(0)', 'double(1, 4)')
  ##  )
  ##),
  ##list(
  ##  name = 'norm',
  ##  variants = c('d', 'p', 'q', 'r'), ## prepended to name of this list
  ##  support = c(-Inf, Inf),
  ##  parameters = list(
  ##    mean = c('double(0)', 'double(1, 4)'),
  ##    sd = c('double(0)', 'double(1, 4)')
  ##  )
  ##)
)

## Takes an element of dist_params list and returns an AD test
## parameterization which test_AD can use.
##
## dist_param: Probability distribution parameterization, which
##             must have the following fields/subfields:
##             - name
##             - variants
##             - args
##               - support
##                 - type
##             Additional args must also have dist and type fields.
## more_args:  Passed to makeOperatorParam().
##
makeADdistTest <- function(dist_param, more_args = NULL) {
  dist_name <- dist_param$name
  ops <- sapply(dist_param$variants, paste0, dist_name, simplify = FALSE)
  argTypes <- sapply(dist_param$variants, function(variant) {
    op <- paste0(variant, dist_name)
    support_type <- dist_param$args$support$type
    first_argType <- switch(
      variant,
      d = support_type,
      p = support_type,
      q = c('double(0)', 'double(1, 4)')
    )
    first_arg_name <- switch(variant, d = 'x', p = 'q', q = 'p')
    grid <- eval(as.call(c(
      expand.grid, list(first_argType),
      lapply(dist_param$args, `[[`, 'type')[-1]
    )))
    argTypes <- as.list(data.frame(t(grid), stringsAsFactors=FALSE))
    lapply(argTypes, function(v) {
      names(v) <- c(first_arg_name, names(dist_param$args)[-1])
      v
    })
  }, simplify = FALSE)
  param_dists <- lapply(dist_param$args[-1], `[[`, 'dist')
  arg_dists <- sapply(dist_param$variants, function(variant) {
    arg1_dist <- switch(
      variant,
      d = list(x = dist_param$args$support$dist),
      p = list(q = dist_param$args$support$dist),
      q = list(p = runif)
    )
    c(arg1_dist, param_dists)
  }, simplify = FALSE)
  op_params <- unlist(
    lapply(
      dist_param$variants, function(variant) {
        lapply(
          argTypes[[variant]],
          function(these_argTypes) {
            wrt_args = intersect(
              dist_param$wrt, names(these_argTypes)
            )
            makeADtest(
              ops[[variant]], these_argTypes,
              wrt_args = wrt_args,
              arg_dists = arg_dists[[variant]],
              more_args = more_args
            )
          }
        )
      }
    ), recursive = FALSE
  )
  names(op_params) <- sapply(op_params, `[[`, 'name')
  return(op_params)
}

dist_tests <- unlist(
  lapply(dist_params, makeADdistTest, more_args = list(log = 1)),
  recursive = FALSE
)

lapply(dist_tests, test_AD)
