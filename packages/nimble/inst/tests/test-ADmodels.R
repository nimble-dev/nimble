source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
context("Testing of derivatives for calculate() for nimbleModels")

test_that('nimbleModel with derivatives compiles correctly.',
          {
            ADCode1 <- nimbleCode({
              a[1] ~ dnorm(0, 1)
              a[2] ~ dnorm(0, 1)
              y[1] ~ dnorm(a[1], 1)
              y[2] ~ dnorm(a[2], 1)
            })
            
            ADMod1 <- nimbleModel(code = ADCode1, data = list(y = numeric(2)), dimensions = list(y = c(2)))
            temporarilyAssignInGlobalEnv(ADMod1)  
            cADMod1 <- compileNimble(ADMod1)
            ## derivs(calculate(model, nodes)) not fully implemented yet, this tests only that compilation is successful.
            }
)