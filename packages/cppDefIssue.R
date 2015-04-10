## see Perry's code below
library(nimble, lib.loc='/tmp/nim41')

dbl <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double(0))
        return(2*x)
    }
    )

cdbl <- compileNimble(dbl)

## FAILS

code1 <- nimbleCode({
    x ~ dnorm(0,1)
    y <- dbl(x)  # inclusion/exclusion reveals the issue
    v ~ dnorm(dbl(mu), 1)
    mu ~ dnorm(0,1)
})
m1 <- nimbleModel(code1, inits = list(x = 0, v=1, mu = 0.5))
cm1 <- compileNimble(m1)
#Error in envRefInferField(x, what, getClass(class(x)), selfEnv) : 
#  ‘filename’ is not a valid field or method name for reference class “RCfunInfoClass”

## WORKS FINE

code1 <- nimbleCode({
    x ~ dnorm(0,1)
    # y <- dbl(x)  # inclusion/exclusion reveals the issue
    v ~ dnorm(dbl(mu), 1)
    mu ~ dnorm(0,1)
})
m1 <- nimbleModel(code1, inits = list(x = 0, v=1, mu = 0.5))
cm1 <- compileNimble(m1)

#################
## Perry's code

library(nimble)
curdir <- getwd()

dbl <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double(0))
        return(2*x)
    }
    )

## not sure this should be compiled previously, but that's not the core issue
## cdbl <- compileNimble(dbl)

## I think the core issues should be a case like the following

useDbl <- nimbleFunction( ## proxy for a nodeFunction
    setup = function(){},
    run = function(y = double(0)) {
        returnType(double(0))
        return(2*dbl(y)) ## should be 4*y
    })

myUseDbl <- useDbl()
myUseDbl$run(4)
##CmyUseDbl <- compileNimble(myUseDbl, dirName = curdir, control = list(debug = TRUE))
##CmyUseDbl$run(4)
## That worked.  Why isn't it working when embedded in a model?

## Have a closer look in the compilation
CmyUseDbl <- compileNimble(myUseDbl, dirName = curdir, control = list(writeFiles = FALSE, compileCpp = FALSE, loadSO = FALSE))
names(CmyUseDbl$cppDefs)
lapply(CmyUseDbl$cppDefs, `[[`, 'filename')
CmyUseDbl$writeFiles('P1_nfRefClass60', stdout()) ## use whatever name came from previous line
CmyUseDbl$trace(writeFiles, browser)
CmyUseDbl$writeFiles('P1_nfRefClass60', stdout())




#############look at the model:
## FAILS
dbl <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double(0))
        return(2*x)
    }
    )

code1 <- nimbleCode({
    y <- dbl(x)  # inclusion/exclusion reveals the issue ## putting it first so it gets compiled first so I can get to debugging it more quickly
    x ~ dnorm(0,1)
    v ~ dnorm(dbl(mu), 1)
    mu ~ dnorm(0,1)
})
## Oh, now I get it the issue is when it is used twice
m1 <- nimbleModel(code1, inits = list(x = 0, v=1, mu = 0.5))
m1$nodes$lifted_dbl_mu$simulate
cm1 <- compileNimble(m1, dirName = curdir)
##cm1 <- compileNimble(m1, dirName = curdir, control = list(debug = TRUE)) ## triggers many browser()s along the way
##options(error = recover)
cm1 <- compileNimble(m1, dirName = curdir, control = list(debug = FALSE, writeFiles = FALSE, compileCpp = FALSE, loadSO = FALSE)) ## only do the R processing

names(cm1$cppDefs)
lapply(cm1$cppDefs, `[[`, 'filename')
cm1$writeFiles('P1_code', stdout()) ## previous line should have generated two unique names.  use each of them as input here.
cm1$writeFiles('P1_code_nfCode', stdout())
cm1$trace(writeFiles, browser)
cm1$writeFiles('P1_code_nfCode', stdout())
## n for a few lines until defs is defined
test <- defs[[6]]$getHincludes()
lapply(test, class)
test[[3]]$filename
test[[4]]$filename ## this is the bug, where did it get populated?

## ok, tracked down the error in populating this cppDef and fixed it

