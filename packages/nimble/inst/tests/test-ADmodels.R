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
            expect_equal(rDerivs$hessian, rDerivs_chainRule$hessian, tolerance = .0001)
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
            expect_equal(rDerivs$hessian, rDerivs_chainRule$hessian, tolerance = .0001)
          }
)


test_that('Derivatives of model$calculate work for nimbleModel with a for loop.',
          {
            
            ADCode3 <- nimbleCode({
              for(i in 1:2){
                y[i, 1:2] ~ dmnorm(mu[1:2], sigma[1:2, 1:2])
              }
              meanVec[1:2] <- c(0,0)
              mu[1:2] ~ dmnorm(meanVec[1:2], diagMat[1:2, 1:2])
              sigma[1:2, 1:2] ~ dwish(diagMat[1:2, 1:2], 3)
            })
            
            simData <- matrix(0, nrow = 2, ncol = 2)
            ADMod3 <- nimbleModel(
              code = ADCode3, dimensions = list(a = 2, b = 2, c = 2), constants = list(diagMat = diag(2)),
              data = list(y = simData), inits = list(mu = c(-1.5, 0.8), sigma = diag(2)))
            temporarilyAssignInGlobalEnv(ADMod3)  
            #cADMod3 <- compileNimble(ADMod3)
            
            ## R derivatives are evaluated below.
            testFxn <- function(mu, sigma){
              origmu <- ADMod3$mu
              ADMod3$mu <- mu
              origsigma <- ADMod3$sigma
              ADMod3$sigma <- sigma
              ADMod3$calculate(ADMod3$getDependencies(c('mu', 'sigma')))
              outVal <- calculate(ADMod3, ADMod3$getDependencies(c('mu', 'sigma')))
              ADMod3$mu <- origmu
              ADMod3$sigma <- origsigma
              return(outVal)
            }
            rDerivs <- nimDerivs(testFxn(mu = c(-1.5, 0.8), sigma = diag(2)))
            rDerivs_chainRule <- nimDerivs(calculate(ADMod3, ADMod3$getDependencies(c('mu', 'sigma'))), wrtPars = c('mu', 'sigma'))
            expect_equal(rDerivs$value, rDerivs_chainRule$value)
            expect_equal(rDerivs$gradient, rDerivs_chainRule$gradient)
            expect_equal(rDerivs$hessian, rDerivs_chainRule$hessian, tolerance = .0001)
          }
)

makeCalcWrapperFunction <- function(calcNodeName, wrtName){
  argSaveLines <- list()
  argSaveLines[[1]] <- substitute(origVals <- list())
  for(i in seq_along(wrtName)){
    argSaveLines[[i+1]] <- substitute({origVals[[I]] <- model[[WRTNAME]];
                                      model[[WRTNAME]] <- WRTNAMEEXPR},
                                     list(I = i,
                                          WRTNAME = wrtName[[i]][1],
                                          WRTNAMEEXPR = as.name(wrtName[[i]][1])))
  }
  calcLine <- list(substitute(model$calculate(model$getDependencies(DEPNAMES)),
                              list(DEPNAMES = wrtName[[i]][1])))
  callLine <-  list(substitute(outVal <- calculate(model, CALCNODENAME),
                               list(CALCNODENAME = calcNodeName)))
  argReplaceLines <- list()
  for(i in seq_along(wrtName)){
    argReplaceLines[[i]] <- substitute(model[[WRTNAME]] <- origVals[[I]],
    list(I = i, WRTNAME = wrtName[[i]][1] ))
  }
  returnLine <- substitute(return(outVal))
  calcWrapperFunction <- function(){}
  body(calcWrapperFunction) <- nimble:::putCodeLinesInBrackets(c(argSaveLines, calcLine, callLine, argReplaceLines, returnLine))
  calcFunctionArgNames <- list('model' = NA)
  for(i in seq_along(wrtName)){
    calcFunctionArgNames[[ wrtName[[i]][1]]] <- NA
  }
  formals(calcWrapperFunction) <- calcFunctionArgNames
  return(calcWrapperFunction)
}

test_ADModelCalculate <- function(model, calcNodeNames = NULL, wrt = NULL, 
                                  testR = TRUE, testCompiled = TRUE,  verbose = TRUE){
  temporarilyAssignInGlobalEnv(model)  
  if(testR){
    for(i in seq_along(calcNodeNames)){
      for(j in seq_along(wrt)){
        wrtNames <- strsplit(wrt[[j]], '\\[')
        RCalcADTestFunction <- makeCalcWrapperFunction(calcNodeNames[[i]], wrtNames)
        argList <- list('model' = quote(model))
        for(k in seq_along(wrtNames)){
            argList[[wrtNames[[k]][1]]] <- model[[wrtNames[[k]][1]]]
        }
        argList <- c(list('RCalcADTestFunction'), argList)
        fxnCall <- as.call(argList)
        fxnCall[[1]] <- quote(RCalcADTestFunction)  
        wrapperDerivs <- eval(substitute(nimDerivs(FXNCALL, wrt = WRT),
                                         list(FXNCALL = fxnCall,
                                              WRT = wrt[[j]])))
        calcDerivs <- nimDerivs(model$calculate(calcNodeNames[[i]]), wrt = wrt[[j]])
        browser()
        expect_equal(wrapperDerivs$value, calcDerivs$value)
        expect_equal(wrapperDerivs$gradient, calcDerivs$gradient, tolerance = .0001)
        expect_equal(wrapperDerivs$hessian, calcDerivs$hessian, tolerance = .001)
      }
    }
  }
  
  if(testCompiled){
    expect_message(cModel <- compileNimble(model))
  }
}

test_ADModelCalculate(ADMod1, calcNodeNames = list(c('a', 'y'), c('a'), c(ADMod1$getDependencies('a'))),
                      wrt = list(c('a', 'y'), c('a[1]', 'y[1]'), c('a[1:2]', 'y[1:2]')), testR = TRUE,
                                      testCompiled = FALSE)
