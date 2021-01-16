source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of model initializaion using initializeModel")


test_that('initializeModel works', {
    code <- nimbleCode({
        x[1] <- 0.5
        x[2] <- x[1] + 1
        x[3] <- x[2]^2
        x[4] ~ dnorm(x[3], 1)
        x[5] <- x[4] + 1
        y ~ dnorm(x[5], 1)
        x[7] <- y + 1
        x[8] <- x[7] + 10
        x[9] ~ dnorm(x[8], 1)
        x[10] ~ dnorm(x[9], 1)
    })

    Rmodel <- nimbleModel(code, data=list(y = 5), calculate = FALSE)

    Rinit <- initializeModel(Rmodel)
    Cmodel <- compileNimble(Rmodel)
    Cinit <- compileNimble(Rinit, project = Rmodel)

    correctValues <-  c(0.5, 1.5, 2.25, 3.512954, 4.512954, NA, 6, 16, 15.673767, 17.003566)

    set.seed(0)
    Rinit$run()
    expect_identical(as.numeric(round(Rmodel$x, 6)), correctValues)
    expect_true(is.na(Rmodel$x[6]))

    set.seed(0)
    Cinit$run()
    expect_identical(as.numeric(round(Cmodel$x, 6)), correctValues)
    expect_true(is.na(Cmodel$x[6]))
})


test_that('initializeModel recalculates all deterministic nodes in topological order', {
    set.seed(0)
    code <- nimbleCode({
        p <- ilogit(a)
        x ~ dbern(p)
    })
    constants <- list()
    data <- list()
    inits <- list(a = 10)
    Rmodel <- nimbleModel(code, constants, data, inits)
    Cmodel <- compileNimble(Rmodel)
    ##
    expect_true(is.na(Rmodel$calculate()))
    expect_true(is.na(Cmodel$calculate()))
    ##
    my_initializeModel <- initializeModel(Rmodel)
    Cmy_initializeModel <- compileNimble(my_initializeModel, project = Rmodel)
    ##
    my_initializeModel$run()
    Cmy_initializeModel$run()
    ##
    expect_equal(Rmodel$x, 1)
    expect_equal(Cmodel$x, 1)
    ##
    expect_equal(Rmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
    expect_equal(Cmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
    ##
    Rmodel$x <- NA
    Cmodel$x <- NA
    Rmodel$a <- -10
    Cmodel$a <- -10
    ##
    my_initializeModel$run()
    Cmy_initializeModel$run()
    ##
    expect_equal(Rmodel$x, 0)
    expect_equal(Cmodel$x, 0)
    expect_equal(Rmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
    expect_equal(Cmodel$getLogProb(), -0.0000453989, tol = 0.0000001)
})



options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)


