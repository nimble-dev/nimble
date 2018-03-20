
######################################################################
######################################################################


## Testing bnp distributions

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
context('Testing NIMBLE distributions')

## dCRP nimble function

x <- c(1,1,3,1,1,3)
conc <- 1

truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
  (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))

try(test_that("dCRP nimble function calculates density correctly: ",
              expect_equal(dCRP(x, conc, size=length(x), log=FALSE),
                           truth,
                           info = paste0("incorrect dCRP nimble function calculation"))))

## test use through log scale
ltruth <- log((conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
                (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1)))

try(test_that("dCRP nimble function calculates density correctly in log scale: ",
              expect_equal(dCRP(x, conc, size=length(x), log=TRUE),
                           ltruth,
                           info = paste0("incorrect dCRP nimble function calculation in log scale"))))


## test use through compile nimble function
cdCRP <- compileNimble(dCRP)

try(test_that("Test that dCRP works correctly in compiled nimble function: ",
              expect_equal(cdCRP(x, conc, size=length(x)), (truth), 
                           info = paste0("incorrect dCRP value in compiled nimble function"))))


## test use through compile nimble function in log scale

try(test_that("Test that dCRP works correctly in compiled nimble function in log scale: ",
              expect_equal(cdCRP(x, conc, size=length(x), log=TRUE), (ltruth), 
                           info = paste0("incorrect dCRP value in compiled nimble function  in log scale"))))



## test: use in model
CRP_code <- nimbleCode({
  x[1:6] ~ dCRP(conc, size=6)
})

Consts <- list(conc = 1)
Inits <- list(x = c(1,1,3,1,1,3))
CRP_model <- nimbleModel(CRP_code, data=Inits, constants=Consts)

CRP_model$x <- x

try(test_that("Test that dCRP calculation is correct in model likelihood calculation: ",
              expect_equal(exp(CRP_model$calculate()), (truth),
                           info = paste0("incorrect likelihood value for dCRP"))))


c_CRP_model <- compileNimble(CRP_model)
c_CRP_model$x
try(test_that("Test that dCRP (compiled) calculation is correct in model likelihood calculation: ",
              expect_equal(exp(c_CRP_model$calculate()), (truth),
                           info = paste0("incorrect likelihood value for compiled dCRP"))))


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

## test of misspecification of 'conc'
try(test_that("dCRP nimble function correctly returns NaN for conc parameter < 0: ",
              expect_equal(dCRP(x, conc=-1, size=length(x), log=FALSE),
                           NaN,
                           info = paste0("incorrect parameters space allowed"))))


## wrong specification of size: size larger than length xi:
CRP_code <- nimbleCode({
  x ~ dnorm(xi[1],1)
  y ~ dnorm(xi[2],1)
  xi[1:2] ~ dCRP(1, 4)
})


data <- list(x = 0, y = 0)
inits <- list(xi=c(1,3))
CRP_model <- nimbleModel(CRP_code, data=data, inits=inits)

CRP_model$xi
CRP_model$simulate('xi')
CRP_model$xi

# length(xi) larger than size
CRP_code <- nimbleCode({
  x ~ dnorm(xi[1],1)
  y ~ dnorm(xi[2],1)
  xi[1:3] ~ dCRP(1, 2)
})


data <- list(x = 0, y = 0)
inits <- list(xi=c(1,3,3))
CRP_model <- nimbleModel(CRP_code, data=data, inits=inits)

CRP_model$xi
CRP_model$simulate('xi')
CRP_model$xi


## random sampling
set.seed(0)
n <- 6
r_samps <- t(replicate(10000, rCRP(n = n, conc)))
# K is the number of unique components in x of length 6
true_EK <- sum(conc/(conc+1:n-1))


try(test_that("Test that random samples (rCRP) have correct mean of K: ",
              expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
                           tol = 0.01,
                           info = "Difference in expected mean of K exceeds tolerance")))

## sampling via `simulate`: does not work!
#set.seed(0)
#simul_samp <- function(model) {
#  model$simulate()
#  return(model$x)
#}

#simul_samps <- t(replicate(10000, simul_samp(c_CRP_model)))

#c_CRP_model$simulate()

#try(test_that("Test that random samples (simulate) have correct mean of K: ",
#              expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
#                           tol = 0.01,
#                           info = "Difference in expected mean of K exceeds tolerance")))



## Test stick_breaking nimble function:


######################################################################
######################################################################


context("Testing of bnp conjugacy")

#RwarnLevel <- options('warn')$warn
#options(warn = -1)

#source(system.file(file.path('tests', 'dynamicIndexingTestLists.R'), package = 'nimble'))

#nimbleOptions(allowDynamicIndexing = TRUE)

## check variations on use of dynamic indexing in BUGS code including building and compiling,
## dependencies, and valid and invalid dynamic index values

#ans1 <- sapply(testsDynIndex, test_dynamic_indexing_model)
#ans2 <- sapply(testsInvalidDynIndex, test_dynamic_indexing_model)
#ans3 <- sapply(testsInvalidDynIndexValue, test_dynamic_indexing_model)

## check conjugacy detection

test_that("Testing conjugacy detection with bnp models", { 
  
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(thetatilde[xi[i]], sd = 1)
    xi[1:4] ~ dCRP(1)
    for(j in 1:5) 
      thetatilde[j] ~ dnorm(0, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_conjugate_dnorm_dnorm",
               info = "failed to detect normal-normal conjugacy")
  
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(thetatilde[xi[i]+1], sd = 1)
    xi[1:(4+1)] ~ dCRP(1)
    for(j in 1:5) 
      thetatilde[j] ~ dnorm(0, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_nonconjugate",
               info = "failed to detect non conjugacy")
  
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(thetatilde[2*xi[i]], sd = 1)
    xi[1:4] ~ dCRP(1)
    for(j in 1:5) 
      thetatilde[j] ~ dnorm(0, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_nonconjugate",
               info = "failed to detect non conjugacy")
  
  
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(thetatilde[xi[i]], sd = 1)
    xi[1:4] ~ dCRP(1)
    for(j in 1:5) 
      thetatilde[j] ~ dpois(1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rpois(5, 1)))
  conf <- configureMCMC(m)
  expect_equal(length(grep("dCRP_conjugate", conf$getSamplers()[[1]]$name)), 0,
               info = "incorrectly  conjugacy detected in poisson-normal case")
  
  # y1 and y2 are independent
  code = nimbleCode({
    for(i in 1:4){ 
      y1[i] ~ dnorm(thetatilde1[xi1[i]], sd = 1)
      y2[i] ~ dnorm(thetatilde2[xi2[i]], sd = 1)}
    xi1[1:4] ~ dCRP(1)
    xi2[1:4] ~ dCRP(1)
    for(j in 1:5){
      thetatilde1[j] ~ dnorm(0,1)
      thetatilde2[j] ~ dnorm(0,1)}
  })
  m = nimbleModel(code, data = list(y1 = rnorm(4), y2 = rnorm(4)),
                  inits = list(xi1 = rep(1,4), xi2 = rep(1,4), thetatilde1=rnorm(5), thetatilde2=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_conjugate_dnorm_dnorm",
               info = "failed to detect normal-normal conjugacy")
  
  code = nimbleCode({
    for(i in 1:4){ 
      y1[i] ~ dnorm(thetatilde[xi1[i]], sd = 1)
      y2[i] ~ dnorm(thetatilde[xi2[i]], sd = 1)}
    xi1[1:4] ~ dCRP(1)
    xi2[1:4] ~ dCRP(1)
    for(j in 1:5){
      thetatilde[j] ~ dnorm(0,1)}
  })
  m = nimbleModel(code, data = list(y1 = rnorm(4), y2 = rnorm(4)),
                  inits = list(xi1 = rep(1,4), xi2 = rep(1,4), thetatilde=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_nonconjugate",
               info = "failed to detect non conjugacy") # more than 1 dCRP distr
  
  code = nimbleCode({
    for(i in 1:4){ 
      y1[i] ~ dnorm(thetatilde[xi[i]], sd = 1)
      y2[i] ~ dnorm(thetatilde[xi[i]], sd = 1)}
    xi[1:4] ~ dCRP(1)
    for(j in 1:5){
      thetatilde[j] ~ dnorm(0,1)}
  })
  m = nimbleModel(code, data = list(y1 = rnorm(4), y2 = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(5)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_nonconjugate",
               info = "failed to detect non conjugacy")
  
  code = nimbleCode({
    for(i in 1:4)
      y[i] ~ dnorm(thetatilde[xi[i]] + mu, sd = 1)
    xi[1:4] ~ dCRP(1)
    for(j in 1:4)
      thetatilde[j] ~ dnorm(0,1)
    mu ~ dnorm(0,1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(4), mu=rnorm(1)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_nonconjugate",
               info = "failed to detect non conjugacy")
  
  code = nimbleCode({
    for(i in 1:4){
      theta[i] <- thetatilde[xi[i]]
      y[i] ~ dnorm(theta[i], sd = 1) }
    xi[1:4] ~ dCRP(1)
    for(j in 1:4)
      thetatilde[j] ~ dnorm(0,1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rnorm(4)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_conjugate_dnorm_dnorm",
               info = "failed to detect normal normal conjugacy with determ node")
  
  code = nimbleCode({
    for(i in 1:4){
      theta[i] <- thetatilde[xi[i]]
      y[i] ~ dnorm(0, var=theta[i]) }
    xi[1:4] ~ dCRP(1)
    for(j in 1:4)
      thetatilde[j] ~ dinvgamma(1,1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), thetatilde=rinvgamma(4,1,1)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[1]]$name, "dCRP_nonconjugate",
               info = "failed to detect non conjugacy for variance") # nos included yet
  
})






######################################################################
######################################################################

# testing models

testBUGSmodel(example = 'test', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

## simple test of normal normal DPM model with no deterministic nodes:

model <- function() {
  for(i in 1:3)
    thetatilde[i] ~ dnorm(mean=0, var=1) 
  xi[1:3] ~ dCRP(1)
  for(i in 1:3)
    y[i] ~ dnorm(thetatilde[xi[i]] , var=1)
}

inits <- list(thetatilde=rnorm(3), xi=rep(1,3))
data <- list(y = rnorm(3))

testBUGSmodel(example = 'testt', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

