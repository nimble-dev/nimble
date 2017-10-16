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

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
