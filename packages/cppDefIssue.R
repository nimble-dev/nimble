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

