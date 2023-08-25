source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleOptions(enableModelMacros = TRUE)
nimbleOptions(enableMacroComments = FALSE)

context("Testing model macros")

test_that('Macro expansion 1',
{
    ## This converts a ~ testMacro(b)
    ## to inputs as (stoch = TRUE, LHS = a, b)
    testMacro <- nimble:::model_macro_builder(
        function(stoch, LHS, RHSarg, modelInfo, .env) {
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, 1)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, modelInfo=modelInfo)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )
    temporarilyAssignInGlobalEnv(testMacro)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(y[1])),
        modelInfo = list())$code,
        quote(x[1] ~ dnorm(y[1], 1))
    )
    
    model1 <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro(y[1])
        })
    )
    expect_identical(
        model1$getCode()[[2]],
        quote(x[1] ~ dnorm(y[1], 1))
    )

    model2 <- nimbleModel(
        nimbleCode({
            for(i in 1:3)
                x[i] ~ testMacro(y[i])
            for(j in 1:3) {
                a[j] ~ testMacro(b[j])
            }
            for(i in 1:3)
                for(j in 1:3)
                    f[i,j] ~ testMacro(g[i])
        })
    )
    expandedCode <- model2$getCode()
    expect_identical(
        expandedCode[[2]][[4]],
        quote(x[i] ~ dnorm(y[i], 1))
    )
    expect_identical(
        expandedCode[[3]][[4]][[2]],
        quote(a[j] ~ dnorm(b[j], 1))
    )
    expect_identical(
        expandedCode[[4]][[4]][[4]],
        quote(f[i, j] ~ dnorm(g[i], 1))
    ) 
})


test_that('Macro expansion 2',
{
    ## Like testMacro above but with unpackArgs = FALSE.
    ## This takes a ~ testMacro(b) with arguments (stoch = TRUE, LHS = a, RHS = b)
    testMacro <- nimble:::model_macro_builder(
        function(stoch, LHS, RHS, modelInfo, .env) {
            if(RHS[[1]] != "testMacro")
                stop("Problem with how testMacro was called")
            RHS[[1]] <- as.name('step') ## something valid
            ans <- substitute(
                OP(LHS, dnorm(RHS, 1)),
                list(LHS = LHS, RHS = RHS,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, modelInfo = modelInfo)
        },
        use3pieces = TRUE, ## default
        unpackArgs = FALSE  ## NON-default
    )
    
    temporarilyAssignInGlobalEnv(testMacro)

    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(y[1])),
        modelInfo = list())$code,
        quote(x[1] ~ dnorm(step(y[1]), 1))
    )
    
    
    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro(y[1])
        })
    )
    expect_identical(
        model$getCode()[[2]],
        quote(x[1] ~ dnorm(nimStep(y[1]), 1))
    )
})

test_that('Macro expansion 3',
{
    ## Like expansion 1, but with multiple named arguments.
    ## This takes a ~ testMacro(newVar = b, newIndex = 3)
    ## as arguments (stoch = TRUE, LHS = a, newVar = b, newIndex = 3)
    testMacro <- nimble:::model_macro_builder(
        function(stoch, LHS, newIndex, newVar, modelInfo, .env) {
        RHS <- substitute(
            X[I],
            list(X = newVar,
                 I = newIndex)
        )
        ans <- substitute(
            OP(LHS, dnorm(RHS, 1)),
            list(LHS = LHS, RHS = RHS,
                 OP = if(stoch) as.name('~') else as.name('<-'))
        )
        list(code = ans, modelInfo = modelInfo)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )

    temporarilyAssignInGlobalEnv(testMacro)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(newVar = z, newIndex = 5)),
            modelInfo = list()
        )$code,
        quote(x[1] ~ dnorm(z[5], 1))
    )

    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro(newVar = z, newIndex = 5)
        })
    )
    expect_identical(
        model$getCode()[[2]],
        quote(x[1] ~ dnorm(z[5], 1))
    )
})

test_that('Macro expansion 4',
{
    ## set use3pieces FALSE, unpackArgs FALSE.
    ## This takes a ~ testMacro(b) as a single argument,
    ## (code = a ~ testMacro(b)).
    ## It does not even need to be in a line with '~' or '<-'.
    testMacro <- nimble:::model_macro_builder(
        function(code, modelInfo, .env) {
            code[[3]][[1]] <- as.name('dnorm')
            list(code = code, modelInfo=modelInfo)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )

    temporarilyAssignInGlobalEnv(testMacro)

    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(y[1], 1)),
            modelInfo = list()
        )$code,
        quote(x[1] ~ dnorm(y[1], 1))
    )
    
    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro(y[1], 1)
        })
    )
    expect_identical(
        model$getCode()[[2]],
        quote(x[1] ~ dnorm(y[1], 1))
    )
})


test_that('Macro expansion 5',
{
    ## This takes a line of code split into arguments.
    ## It is designed for a line without '~' or '<-'.
    ## It takes testMacro(arg1, arg2, arg3)
    testMacro <- nimble:::model_macro_builder(
        function(arg1, arg2, arg3, modelInfo, .env) {
            code <- substitute(A <- B + C, list(A = arg1, B = arg2, C = arg3))
            list(code = code, modelInfo = modelInfo)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = TRUE  ## default
    )

    temporarilyAssignInGlobalEnv(testMacro)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(testMacro(x[1], y[2], z[i, j])),
            modelInfo = list()
        )$code,
        quote(x[1] <- y[2] + z[i, j])
    )

    model <- nimbleModel(
        nimbleCode({
            testMacro(a[1], b[2], c[3])
        })
    )
    expect_identical(
        model$getCode()[[2]],
        quote(a[1] <- b[2] + c[3])
    )
})

test_that('Macro expansion 6 (recursive macro expansion)',
{
    testMacroInner <- nimble:::model_macro_builder(
        function(stoch, LHS, RHSarg, modelInfo, ...) {
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, 1)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, modelInfo = modelInfo)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )
    temporarilyAssignInGlobalEnv(testMacroInner)
    
    ## a ~ testMacroOuter(b)
    ## becomes a ~ testMacroInner(b)
    testMacroOuter <- nimble:::model_macro_builder(
        function(code, modelInfo, ...) {
            code[[3]][[1]] <- as.name('testMacroInner')
            list(code = code, modelInfo = modelInfo)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )
    temporarilyAssignInGlobalEnv(testMacroOuter)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacroOuter(y[1])),
            modelInfo = list()
        )$code,
        quote(x[1] ~ dnorm(y[1], 1))
    )


    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacroOuter(y[1])
        })
    )

    expect_identical(
        model$getCode()[[2]],
        quote(x[1] ~ dnorm(y[1], 1))
    )
})

test_that(paste0('Macro expansion 7 (correct trapping of ',
                 'failure in recursive macro expansion)'),
{
    ## test of failure in recursive expansion
    ## test of recursive expansion:
    ## This a ~ testMacroInner(b)
    ## with inputs as (stoch = TRUE, LHS = a, b)
    testMacroInner <- nimble:::model_macro_builder(
        function(stoch, LHS, RHSarg, modelInfo, ...) {
            stop()
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, 1)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, modelInfo = modelInfo)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )
    temporarilyAssignInGlobalEnv(testMacroInner)
    ## a ~ testMacroOuter(b)
    ## becomes a ~ testMacroInner(b)
    testMacroOuter <- nimble:::model_macro_builder(
        function(code, modelInfo, .env) {
            code[[3]][[1]] <- as.name('testMacroInner')
            list(code = code, modelInfo = modelInfo)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )
    temporarilyAssignInGlobalEnv(testMacroOuter)

    #cat("\nTwo 'unused argument' error messages known to occur here:\n") 
    #cat("\nNot sure why these are supposed to fail:\n")

    expect_error(nimble:::codeProcessModelMacros(
                              quote(x[1] ~ testMacroOuter(y[1], 1)), modelInfo = list()),
                 "Model macro testMacroInner\\(expanded from testMacroOuter\\) failed.")

    #expect_error(ans <- nimbleModel(
    #        nimbleCode({
    #            x[1] ~ testMacroOuter(y[1], 1)
    #        })
    #    ), "Model macro testMacroInner\\(expanded from testMacroOuter\\) failed.")

    ## Replaced below with use of expect_error that does check error message.
    
    ## The expect_failure mechanism does not allow to check the error
    ## message returned from successful error-trapping.  Hence we use
    ## this try(), which is potentially dangerous for masking errors
    ## but in this case is safe because we follow it with an
    ## expectation that an error did occur.
    
    ## ans <- try(nimble:::codeProcessModelMacros(
    ##     quote(x[1] ~ testMacroOuter(y[1], 1)
    ##           ## second arg triggers
    ##           ## failure for testing
    ##           )))
    ## expect_true(inherits(ans, 'try-error'),
    ##             'failure occurred when intended')
    ## expect_identical(as.character(ans),
    ##                  "Error : Model macro testMacroInner(expanded from testMacroOuter) failed.\n")

    
    ## ans <- try(
    ##     nimbleModel(
    ##         nimbleCode({
    ##             x[1] ~ testMacroOuter(y[1], 1)
    ##         })
    ##     )
    ## )
    ## expect_true(inherits(ans, 'try-error'),
    ##             'failure occurred when intended')
    ## expect_identical(as.character(ans),
    ##                  "Error : Model macro testMacroInner(expanded from testMacroOuter) failed.\n")
})

test_that('duplicate variables from macro expansion error-trapped correctly',
{
    ## from roxygen example
    flat_normal_priors <- nimble:::model_macro_builder(
        function(..., modelInfo) {
          allVars <- list(...)
          allVars <- allVars[ !(names(allVars) %in% (c('.constants','.env'))) ]
            priorDeclarations <- lapply(allVars,
                                        function(x)
                                            substitute(VAR ~ dnorm(0, sd = 1000),
                                                       list(VAR = x)))
            newCode <- quote({})
            newCode[2:(length(allVars)+1)] <- priorDeclarations
            list(code = newCode, modelInfo = modelInfo)
        },
        use3pieces = FALSE,
        unpackArgs = TRUE
    )
    temporarilyAssignInGlobalEnv(flat_normal_priors)
    ## Safe use of try due to immediate next test
    expect_error(model <- nimbleModel(
        nimbleCode(
        {
            flat_normal_priors(mu, beta, gamma)
            mu ~ dexp(4)
        }
        )), "There are multiple definitions for node\\(s\\): mu.")
})

test_that('duplicate nested indices from macro expansion error-trapped correctly',
{
    all_dnorm <- nimble:::model_macro_builder(
        function(stoch, LHS, RHSvar, start, end, sd = 1, modelInfo, ...) {
            newCode <- substitute(
                for(i in START:END) {
                    LHS[i] ~ dnorm(RHSvar[i], SD)
                },
                list(START = start,
                     END = end,
                     LHS = LHS,
                     RHSvar = RHSvar,
                     SD = sd))
            list(code = newCode, modelInfo = modelInfo)
        },
        use3pieces = TRUE,
        unpackArgs = TRUE 
    )
    temporarilyAssignInGlobalEnv(all_dnorm)
    ## Safe use of try due to immediately next test
    expect_error(model <- nimbleModel(
        nimbleCode(
        {
            for(i in 1:3)
                x ~ all_dnorm(mu, start = 1, end = 10)        
        }
        )), "Variable i used multiple times as for loop index")
})

test_that('constants and env are accessed correctly, even in recursion',
{
  testEnv <- new.env()
  testEnv$var_in_env <- 5
  temporarilyAssignInGlobalEnv(testEnv)

  testMacroInner <- nimble:::model_macro_builder(
    function(stoch, LHS, RHSarg, modelInfo, .env) {
      constant_of_interest <- modelInfo$constants$constant_of_interest
      var_in_env <- eval(quote(var_in_env), envir = .env)
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, VAR_IN_ENV + CONSTANT_OF_INTEREST)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'),
                     VAR_IN_ENV = var_in_env,
                     CONSTANT_OF_INTEREST = constant_of_interest)
            )
            list(code = ans, modelInfo = modelInfo)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
        )
  temporarilyAssignInGlobalEnv(testMacroInner)


    ## a ~ testMacroOuter(b)
    ## becomes a ~ testMacroInner(b)
  testMacroOuter <- nimble:::model_macro_builder(
        function(code, modelInfo, .env) {
            code[[3]][[1]] <- as.name('testMacroInner')
            list(code = code, modelInfo = modelInfo)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
        )
  temporarilyAssignInGlobalEnv(testMacroOuter)

    macroOut <- nimble:::codeProcessModelMacros(
                   quote(x[1] ~ testMacroOuter(y[1])),
                   env = testEnv,
                   modelInfo = list(constants=list(constant_of_interest = 10))
        )

    expect_identical(macroOut$code,
                     quote(x[1] ~ dnorm(y[1], 5 + 10)))
    expect_identical(macroOut$modelInfo$constants,
                     list(constant_of_interest = 10))
                     
    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacroOuter(y[1])
        }),
        constants = list(constant_of_interest = 10),
        userEnv = testEnv
    )

    expect_identical(
        model$getCode()[[2]],
        quote(x[1] ~ dnorm(y[1], 5 + 10))
    )
})

test_that("macros can add constants", {

  inp_constants <- list(a = 1)

  testMacro <- list(process = function(code, modelInfo, .env){
      code[[3]] <- quote(dnorm(0, sd = 1))
      modelInfo$constants$b <- 2
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"


  temporarilyAssignInGlobalEnv(testMacro)

  code <- nimbleCode({
    x ~ testMacro()
  })
  
  modelInfo <- list(constants = inp_constants)
  out <- nimble:::processModelMacros(code, modelInfo = modelInfo)
  expect_equal(out$code,
    quote({
           x ~ dnorm(0, sd = 1)
        })
  )
  expect_equal(out$modelInfo$constants,
    list(a = 1, b = 2)
  )

})

test_that("processModelMacros converts factors to numeric", {
  inp_constants <- list(y = rnorm(3), x = factor(c("a","b","c")))

  code <- nimbleCode({
    for (i in 1:3){
      y[i] ~ dnorm(x[i], sd = 1)
    }
  })
  
  mod <- nimble:::processModelMacros(code, modelInfo=list(constants=inp_constants))

  expect_equal(
    mod$modelInfo$constants$x,
    as.numeric(inp_constants$x)
  )

})

test_that("convertFactorConstantsToNumeric converts factors to numeric", {

  # Make sure character matrices are handled
  inp_constants <- list(y = rnorm(3), x = factor(c("a","b","c")),
                        x2 = matrix(c("a","b","c","d"), 2, 2),
                        x3 = c("a","b","c"))
  out_constants <- nimble:::convertFactorConstantsToNumeric(inp_constants)
  expect_identical(out_constants,
                   list(y=inp_constants$y, x = c(1,2,3),
                        x2=matrix(c(1,2,3,4), 2, 2),
                        x3=c(1,2,3)))
})

test_that("processModelMacros makes index generator available", {

  testMacro <- list(process = function(code, modelInfo, .env){
      code[[3]] <- str2lang(modelInfo$indexCreator())
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"

  temporarilyAssignInGlobalEnv(testMacro)

  code <- nimbleCode({
    x ~ testMacro()
  })

  out <- nimble:::processModelMacros(code, modelInfo=list())
  
  expect_equal(out$code,
               quote({
                 x ~ i_1
               })
  )

})

test_that("generated parameters are stored in model definition",{  
  testMacro <- list(process = function(code, modelInfo, .env){
      code[[3]] <- quote(dnorm(alpha, 1))
      #modelInfo$parameters$testMacro <- "alpha"
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"

  temporarilyAssignInGlobalEnv(testMacro)

  code <- nimbleCode({
    x ~ testMacro()
  })

  mod <- nimbleModel(code, constants=list(), returnDef=TRUE)
  
  expect_equal(
    mod$macroParameters,
    list(testMacro = list(list(LHS = NULL, RHS = "alpha")))
  )

  mod <- nimbleModel(code, constants=list())
  expect_equal(
    mod$getMacroParameters(),
    list(testMacro = list("alpha"))
  )

})

test_that("getMacroPars finds parameters in a chunk of code", {

  inp <- quote(1 + 1)
  expect_equal(nimble:::getMacroPars(inp), NULL)

  inp <- quote(x ~ dnorm(alpha, beta))
  expect_equal(nimble:::getMacroPars(inp), c("x", "alpha", "beta"))

  inp <- quote(x[alpha] <- 1)
  expect_equal(nimble:::getMacroPars(inp), c("x", "alpha"))

  inp <- quote({x ~ dnorm(alpha, beta)
                y[gamma] <- 1
              })
  expect_equal(nimble:::getMacroPars(inp), c("x", "alpha", "beta", "y", "gamma"))

})

test_that("newMacroPars finds new parameters generated by a macro", {
  old <- quote(x ~ testMacro())
  new <- quote(x ~ dnorm(alpha, beta))

  expect_equal(nimble:::newMacroPars(old, new),
               list(LHS = NULL, RHS = c("alpha", "beta")))

})

test_that("checkMacroPars removes intermediate parameters and renames parameter list", {
  old <- quote({x ~ testMacro()
                y ~ testMacro()
               })
  new <- quote({x ~ dnorm(alpha, beta)
               y ~ dnorm(gamma, eps)
               })

  pars <- list(testMacro = list(LHS = NULL, RHS = c("alpha", "beta", "INT")),
               testMacro = list(LHS = NULL, RHS = c("gamma", "eps")))

  out <- nimble:::checkMacroPars(pars, old, new)

  expect_equal(names(out), c("testMacro"))
  expect_equal(length(out$testMacro), 2)
  expect_equal(out$testMacro, 
               list(list(LHS = NULL, RHS = c("alpha", "beta")),
                    list(LHS = NULL, RHS = c("gamma", "eps"))))

  old <- quote(x ~ testMacro())
  new <- quote(x ~ dnorm(alpha, beta))

  pars <- list(testMacro = list(LHS = NULL, RHS = c("alpha", "beta")))

  out <- nimble:::checkMacroPars(pars, old, new)

  expect_equal(out, list(testMacro = list(list(LHS = NULL, RHS = c("alpha", "beta")))))

  # when there are no macros
  old <- quote(x ~ dnorm(alpha, beta))
  new <- quote(x ~ dnorm(alpha, beta))
  pars <- list()

  out <- nimble:::checkMacroPars(pars, old, new)
  expect_equal(out, list())
})

test_that("nimbleModel creates and stores list of macro-generated parameters", {
  testMacro <- list(process = function(code, modelInfo, .env){
      code[[3]] <- quote(dnorm(alpha, beta))
      modelInfo$parameters$testMacro <- "alpha"
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"

  temporarilyAssignInGlobalEnv(testMacro)

  code <- nimbleCode({
    x ~ testMacro()
  })

  mod <- nimbleModel(code, constants=list())
  expect_equal(mod$getMacroParameters(),
               list(testMacro = list(c("alpha", "beta"))))

})

test_that("getMacroParameters() can pull out different subsets of parameters", {

  testMacro <- list(process = function(code, modelInfo, .env){
      code <- quote({
        for (i_ in 1:n){
          mu[i_] <- alpha + beta
          y[i_] ~ dnorm(0, sigma)
        }
        alpha ~ dnorm(0, 1)
      })
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"

  temporarilyAssignInGlobalEnv(testMacro)

  code <- nimbleCode({
    y[1:n] ~ testMacro()
  })

  const <- list(y = rnorm(10), n = 10)

  mod <- nimbleModel(code, constants=const)

  expect_equal(mod$getCode(),
                   quote({
                     for (i_ in 1:n){
                      mu[i_] <- alpha + beta
                      y[i_] ~ dnorm(0, sigma)
                     }
                     alpha ~ dnorm(0, 1)
                    })
                   )
  
  expect_equal(
    mod$getMacroParameters(),
    list(testMacro = list(c("mu", "alpha", "beta", "sigma")))
  )
  expect_equal(
    mod$getMacroParameters(includeRHS = FALSE),
    list(testMacro = list(c("mu", "alpha")))
  )
  expect_equal(
    mod$getMacroParameters(includeLHS = FALSE),
    list(testMacro = list(c("alpha", "beta", "sigma")))
  )
  expect_equal(
    mod$getMacroParameters(includeIndices = TRUE),
    list(testMacro = list(c("mu", "i_", "alpha", "beta", "sigma")))
  )
  expect_equal(
    mod$getMacroParameters(includeDeterm = FALSE),
    list(testMacro = list(c("alpha", "beta", "sigma")))
  )
  expect_equal(
    mod$getMacroParameters(includeStoch = FALSE),
    list(testMacro = list(c("mu", "beta", "sigma")))
  )
})

test_that("removeExtraBrackets cleans up output from codeProcessModelMacros",{

  code <- nimbleCode({
  {
    alpha <- 1
    for (i in 1:n){
      {
        for (j in 1:k){
          z <- 1
          {
          y <- 1
          }
        }
      }
    }
  }
  })
  
  expect_equal(
    nimble:::removeExtraBrackets(code),
    quote({
      alpha <- 1
      for (i in 1:n){
        for (j in 1:k){
          z <- 1
          y <- 1
        }
      }
    })
  )

})

test_that("comments are added to macro output if enableMacroComments = TRUE",{
  nimbleOptions(enableMacroComments = FALSE)
  
  testMacro <- list(process = function(code, modelInfo, .env){
      code <- quote({
        for (i_ in 1:n){
          mu[i_] <- testMacroInner()
          y[i_] ~ dnorm(0, sigma)
        }
        alpha ~ dnorm(0, 1)
        
      })
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"
  temporarilyAssignInGlobalEnv(testMacro)
  
  testMacroInner <- list(process = function(code, modelInfo, .env){
    code[[3]] <- quote(alpha + beta)
    list(code = code, modelInfo = modelInfo)
  })
  class(testMacroInner) <- "model_macro"
  temporarilyAssignInGlobalEnv(testMacroInner)

  code <- nimbleCode({
    y[1:n] ~ testMacro()
  })

  const <- list(y = rnorm(10), n = 10)

  mod <- nimbleModel(code, constants=const)

  expect_equal(mod$getCode(),
                   quote({
                     for (i_ in 1:n){
                      mu[i_] <- alpha + beta
                      y[i_] ~ dnorm(0, sigma)
                     }
                     alpha ~ dnorm(0, 1)
                    })
                   )

  nimbleOptions(enableMacroComments = TRUE)

  mod <- nimbleModel(code, constants=const)
  expect_equal(mod$getCode(),
                   quote({
                     "# testMacro"
                     for (i_ in 1:n){
                      "  ## testMacroInner"
                      mu[i_] <- alpha + beta
                      "  ## ----"
                      y[i_] ~ dnorm(0, sigma)
                     }
                     alpha ~ dnorm(0, 1)
                     "# ----"
                    })
                   )

  nimbleOptions(enableMacroComments = FALSE)

})

test_that("message given if same parameter is generated by macros multiple times on LHS", {
  testMacro <- list(process = function(code, modelInfo, .env){
      code <- quote(alpha ~ dnorm(0, 1))
      list(code = code, modelInfo=modelInfo)
    }
  )
  class(testMacro) <- "model_macro"
  temporarilyAssignInGlobalEnv(testMacro)
  
  code <- nimbleCode({
    testMacro()
    testMacro()
  })
  nimbleOptions(verbose = TRUE)
  expect_message(tryCatch({
    mod <- nimbleModel(code, constants=list())
  }, error = function(e) NULL), regexp = "LHS parameter")
  nimbleOptions(verbose = FALSE)
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
