source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleOptions(enableModelMacros = TRUE)

context("Testing model macros")

test_that('Macro expansion 1',
{
    ## This converts a ~ testMacro(b)
    ## to inputs as (stoch = TRUE, LHS = a, b)
    testMacro <- nimble:::model_macro_builder(
        function(stoch, LHS, RHSarg, .constants, .env) {
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, 1)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, constants=.constants)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )
    temporarilyAssignInGlobalEnv(testMacro)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(y[1]))
        ),
        list(code=quote(x[1] ~ dnorm(y[1], 1)), constants=list())
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
        function(stoch, LHS, RHS, .constants, .env) {
            if(RHS[[1]] != "testMacro")
                stop("Problem with how testMacro was called")
            RHS[[1]] <- as.name('step') ## something valid
            ans <- substitute(
                OP(LHS, dnorm(RHS, 1)),
                list(LHS = LHS, RHS = RHS,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, constants=.constants)
        },
        use3pieces = TRUE, ## default
        unpackArgs = FALSE  ## NON-default
    )
    
    temporarilyAssignInGlobalEnv(testMacro)

    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(y[1]))
        ),
        list(code=quote(x[1] ~ dnorm(step(y[1]), 1)), constants=list())
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
        function(stoch, LHS, newIndex, newVar, .constants, .env) {
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
        list(code = ans, constants=.constants)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )

    temporarilyAssignInGlobalEnv(testMacro)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(newVar = z, newIndex = 5))
        ),
        list(code=quote(x[1] ~ dnorm(z[5], 1)), constants=list())
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
        function(code, .constants, .env) {
            code[[3]][[1]] <- as.name('dnorm')
            list(code = code, constants=.constants)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )

    temporarilyAssignInGlobalEnv(testMacro)

    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro(y[1], 1))
        ),
        list(code=quote(x[1] ~ dnorm(y[1], 1)), constants=list())
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
        function(arg1, arg2, arg3, .constants, .env) {
            code <- substitute(A <- B + C, list(A = arg1, B = arg2, C = arg3))
            list(code = code, constants=.constants)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = TRUE  ## default
    )

    temporarilyAssignInGlobalEnv(testMacro)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(testMacro(x[1], y[2], z[i, j]))
        ),
        list(code=quote(x[1] <- y[2] + z[i, j]), constants=list())
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
        function(stoch, LHS, RHSarg, .constants, ...) {
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, 1)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans, constants=.constants)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )
    temporarilyAssignInGlobalEnv(testMacroInner)
    
    ## a ~ testMacroOuter(b)
    ## becomes a ~ testMacroInner(b)
    testMacroOuter <- nimble:::model_macro_builder(
        function(code, .constants, ...) {
            code[[3]][[1]] <- as.name('testMacroInner')
            list(code = code, constants=.constants)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )
    temporarilyAssignInGlobalEnv(testMacroOuter)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacroOuter(y[1]))
        ),
        list(code=quote(x[1] ~ dnorm(y[1], 1)), constants=list())
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
        function(stoch, LHS, RHSarg, ...) {
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, 1)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )
    temporarilyAssignInGlobalEnv(testMacroInner)
    ## a ~ testMacroOuter(b)
    ## becomes a ~ testMacroInner(b)
    testMacroOuter <- nimble:::model_macro_builder(
        function(code, .constants, .env) {
            code[[3]][[1]] <- as.name('testMacroInner')
            list(code = code, constants=.constants)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )
    temporarilyAssignInGlobalEnv(testMacroOuter)

    cat("\nTwo 'unused argument' error messages known to occur here:\n") 
    cat("\nNot sure why these are supposed to fail:\n")

    expect_error(nimble:::codeProcessModelMacros(
                              quote(x[1] ~ testMacroOuter(y[1], 1))),
                 "Model macro testMacroInner\\(expanded from testMacroOuter\\) failed.")

    expect_error(ans <- nimbleModel(
            nimbleCode({
                x[1] ~ testMacroOuter(y[1], 1)
            })
        ), "Model macro testMacroInner\\(expanded from testMacroOuter\\) failed.")

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
        function(...) {
          allVars <- list(...)
          allVars <- allVars[ !(names(allVars) %in% (c('.constants','.env'))) ]
            priorDeclarations <- lapply(allVars,
                                        function(x)
                                            substitute(VAR ~ dnorm(0, sd = 1000),
                                                       list(VAR = x)))
            newCode <- quote({})
            newCode[2:(length(allVars)+1)] <- priorDeclarations
            list(code = newCode)
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
        function(stoch, LHS, RHSvar, start, end, sd = 1, ...) {
            newCode <- substitute(
                for(i in START:END) {
                    LHS[i] ~ dnorm(RHSvar[i], SD)
                },
                list(START = start,
                     END = end,
                     LHS = LHS,
                     RHSvar = RHSvar,
                     SD = sd))
            list(code = newCode)
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
    function(stoch, LHS, RHSarg, .constants, .env) {
      constant_of_interest <- .constants$constant_of_interest
      var_in_env <- eval(quote(var_in_env), envir = .env)
            ans <- substitute(
                OP(LHS, dnorm(RHSarg, VAR_IN_ENV + CONSTANT_OF_INTEREST)),
                list(LHS = LHS, RHSarg = RHSarg,
                     OP = if(stoch) as.name('~') else as.name('<-'),
                     VAR_IN_ENV = var_in_env,
                     CONSTANT_OF_INTEREST = constant_of_interest)
            )
            list(code = ans, constants = .constants)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
        )
  temporarilyAssignInGlobalEnv(testMacroInner)


    ## a ~ testMacroOuter(b)
    ## becomes a ~ testMacroInner(b)
  testMacroOuter <- nimble:::model_macro_builder(
        function(code, .constants, .env) {
            code[[3]][[1]] <- as.name('testMacroInner')
            list(code = code, constants = .constants)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
        )
  temporarilyAssignInGlobalEnv(testMacroOuter)

    expect_identical(
        nimble:::codeProcessModelMacros(
                   quote(x[1] ~ testMacroOuter(y[1])),
                   env = testEnv,
                   constants = list(constant_of_interest = 10)
        ),
        list(code=quote(x[1] ~ dnorm(y[1], 5 + 10)),
             constants=list(constant_of_interest=10))
    )


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

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
