source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of modelValues interfaces")

test_that("copying of argument passing in nimbleFunctionInterface", {
    mvc <- modelValuesConf(vars = c('d1'),
                           types = c('double'),
                           sizes = list( d1 = 5 ))
    
    mv <- mvc(2)
    
    nfG <- nimbleFunction(
        setup = function(mv) {},
        run = function() {
            mv['d1', 1] <<- rnorm(5)
            ans <- mv['d1', 1]
            return(ans)
            returnType(double(1))
        }
    )
    nf <- nfG(mv)
    nf$run()
    
    cnf <- compileNimble(nf)
    cnf$run()
    
    cnf$mv['d1',1] <- as.numeric(1:5)
    ans <- cnf$mv['d1',1]
    
    expect_identical(ans, as.numeric(1:5))
})

test_that("as.matrix and matrix2mv work for modelValues",
{
    mvc <- modelValuesConf(vars = c('d1', 'd2', 'd3'),
                           types = c('double', 'double', 'double'),
                           sizes = list( d1 = 5, d2 = c(3, 5), d3 = c(3, 5, 7)))
    
    mv <- mvc(3)

    nfG <- nimbleFunction(
        setup = function(mv) {
            setupOutputs(mv)
        },
        run = function() {}
    )
    nf <- nfG(mv)
    cnf <- compileNimble(nf)

    for(i in 1:3) {
        mv$d1[[i]] <- rnorm(5)
        mv$d2[[i]] <- matrix(rnorm(15), nrow = 3)
        mv$d3[[i]] <- array(rnorm(3*5*7), dim = c(3, 5, 7))
    }
    ## check assignment by variable
    cnf$mv['d1'] <- mv['d1']
    expect_equal(cnf$mv['d1'], mv['d1'])
    cnf$mv['d2'] <- mv['d2']
    expect_equal(cnf$mv['d2'], mv['d2'])
    cnf$mv['d3'] <- mv['d3']
    expect_equal(cnf$mv['d3'], mv['d3'])

    ## get new values and check that as.matrix and matrix2mv are inverses
    for(i in 1:3) {
        mv$d1[[i]] <- rnorm(5)
        mv$d2[[i]] <- matrix(rnorm(15), nrow = 3)
        mv$d3[[i]] <- array(rnorm(3*5*7), dim = c(3, 5, 7))
    }
    nimble:::matrix2mv(as.matrix(mv), cnf$mv)
    expect_equal(cnf$mv['d1'], mv['d1'])
    expect_equal(cnf$mv['d2'], mv['d2'])
    expect_equal(cnf$mv['d3'], mv['d3'])
    newMat <- as.matrix(cnf$mv)
    expect_equal(as.matrix(mv), newMat)
}
)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
