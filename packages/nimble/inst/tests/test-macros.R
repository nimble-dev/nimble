source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

nimbleOptions(enableModelMacros = TRUE)

context("Testing model macros")

test_that('Macro expansion 1',
{
    ## This converts a ~ testMacro1(b)
    ## to inputs as (stoch = TRUE, LHS = a, b)
    testMacro1 <- nimble:::model_macro_builder(
        function(stoch, LHS, RHSarg) {
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
    temporarilyAssignInGlobalEnv(testMacro1)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro1(y[1]))
        ),
        quote(x[1] ~ dnorm(y[1], 1))
    )
    
    model1 <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro1(y[1])
        })
    )
    expect_identical(
        model1$getCode()[[2]],
        quote(x[1] ~ dnorm(y[1], 1))
    )

    model2 <- nimbleModel(
        nimbleCode({
            for(i in 1:3)
                x[i] ~ testMacro1(y[i])
            for(j in 1:3) {
                a[j] ~ testMacro1(b[j])
            }
            for(i in 1:3)
                for(j in 1:3)
                    f[i,j] ~ testMacro1(g[i])
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
    ## Like testMacro1 but with unpackArgs = FALSE.
    ## This takes a ~ testMacro2(b) with arguments (stoch = TRUE, LHS = a, RHS = b)
    testMacro2 <- nimble:::model_macro_builder(
        function(stoch, LHS, RHS) {
            if(RHS[[1]] != "testMacro2") stop("Problem with how testMacro2 was called")
            RHS[[1]] <- as.name('step') ## something valid
            ans <- substitute(
                OP(LHS, dnorm(RHS, 1)),
                list(LHS = LHS, RHS = RHS,
                     OP = if(stoch) as.name('~') else as.name('<-'))
            )
            list(code = ans)
        },
        use3pieces = TRUE, ## default
        unpackArgs = FALSE  ## NON-default
    )
    
    temporarilyAssignInGlobalEnv(testMacro2)

    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro2(y[1]))
        ),
        quote(x[1] ~ dnorm(step(y[1]), 1))
    )
    
    
    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro2(y[1])
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
    ## This takes a ~ testMacro3(newVar = b, newIndex = 3)
    ## as arguments (stoch = TRUE, LHS = a, newVar = b, newIndex = 3)
    testMacro3 <- nimble:::model_macro_builder(
        function(stoch, LHS, newIndex, newVar) {
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
        list(code = ans)
        },
        use3pieces = TRUE, ## default
        unpackArgs = TRUE  ## default
    )

    temporarilyAssignInGlobalEnv(testMacro3)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro3(newVar = z, newIndex = 5))
        ),
        quote(x[1] ~ dnorm(z[5], 1))
    )

    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro3(newVar = z, newIndex = 5)
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
    ## This takes a ~ testMacro4(b) as a single argument,
    ## (code = a ~ testMacro4(b)).
    ## It does not even need to be in a line with '~' or '<-'.
    testMacro4 <- nimble:::model_macro_builder(
        function(code) {
            code[[3]][[1]] <- as.name('dnorm')
            list(code = code)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = FALSE  ## NON-default
    )

    temporarilyAssignInGlobalEnv(testMacro4)

    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(x[1] ~ testMacro4(y[1], 1))
        ),
        quote(x[1] ~ dnorm(y[1], 1))
    )
    
    model <- nimbleModel(
        nimbleCode({
            x[1] ~ testMacro4(y[1], 1)
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
    ## It takes testMacro5(arg1, arg2, arg3)
    testMacro5 <- nimble:::model_macro_builder(
        function(arg1, arg2, arg3) {
            code <- substitute(A <- B + C, list(A = arg1, B = arg2, C = arg3))
            list(code = code)
        },
        use3pieces = FALSE, ## NON-default
        unpackArgs = TRUE  ## default
    )

    temporarilyAssignInGlobalEnv(testMacro5)
    
    expect_identical(
        nimble:::codeProcessModelMacros(
            quote(testMacro5(x[1], y[2], z[i, j]))
        ),
        quote(x[1] <- y[2] + z[i, j])
    )

    model <- nimbleModel(
        nimbleCode({
            testMacro5(a[1], b[2], c[3])
        })
    )
    expect_identical(
        model$getCode()[[2]],
        quote(a[1] <- b[2] + c[3])
    )
})

