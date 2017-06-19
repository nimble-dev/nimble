source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
context("Testing of derivatives for calculate() for nimbleModels")



test_that('Derivatives of model$calculate work for nimbleModel with scalar nodes.',
          {
            ADCode1 <- nimbleCode({
              a[1] ~ dnorm(0, 1)
              a[2] ~ dnorm(0, 1)
              y[1] ~ dnorm(a[1], 1)
              y[2] ~ dnorm(a[2], 1)
            })
            
            ADMod1 <- nimbleModel(code = ADCode1, data = list(y = numeric(2)), dimensions = list(y = c(2)),
                                  inits = list(a = c(1,1)))
            temporarilyAssignInGlobalEnv(ADMod1)  
            cADMod1 <- compileNimble(ADMod1)
            
            ## R derivatives are evaluated below.
            testFxn <- function(a){
              origa <- ADMod1$a
              ADMod1$a <- a
              outVal <- calculate(ADMod1, ADMod1$getDependencies('a'))
              ADMod1$a <- origa
              return(outVal)
            }
            rDerivs <- nimDerivs(testFxn(a = c(1,1)))
            rDerivs_chainRule <- nimDerivs(calculate(ADMod1, ADMod1$getDependencies('a')), wrtPars = 'a')
            expect_equal(rDerivs$value, rDerivs_chainRule$value)
            expect_equal(rDerivs$gradient, rDerivs_chainRule$gradient)
            }
)

test_that('Derivatives of model$calculate work for nimbleModel with multivariate nodes.',
          {
            
            ADCode2 <- nimbleCode({
              a[1:2] ~ dmnorm(c[1:2], diagMat[,])
              c[1:2] <- b[1:2] + c(1,1)
              b[1] ~ dnorm(1, 1)
              b[2] ~ dnorm(1, 1)
            })
            
            ADMod2 <- nimbleModel(
              code = ADCode2, dimensions = list(a = 2, b = 2, c = 2), constants = list(diagMat = diag(2)),
              inits = list(a = c(2.1, 1.2), b  = c(1,2)))
            temporarilyAssignInGlobalEnv(ADMod2)  
            # cADMod2 <- compileNimble(ADMod2)
            
            ## R derivatives are evaluated below.
            testFxn <- function(a, b){
              origa <- ADMod2$a
              ADMod2$a <- a
              origb <- ADMod2$b
              ADMod2$b <- b
              ADMod2$calculate(ADMod2$getDependencies('b'))
              outVal <- calculate(ADMod2, ADMod2$getDependencies('b'))
              ADMod2$b <- origb
              ADMod2$a <- origa
              return(outVal)
            }
            rDerivs <- nimDerivs(testFxn(a = c(2.1, 1.2), b = c(1, 2)))
            rDerivs_chainRule <- nimDerivs(calculate(ADMod2, ADMod2$getDependencies('b')), wrtPars = c('a','b'))
            expect_equal(rDerivs$value, rDerivs_chainRule$value)
            expect_equal(rDerivs$gradient, rDerivs_chainRule$gradient)
          }
)

