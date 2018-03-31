
######################################################################
######################################################################


source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))


## Test stick_breaking nimble function:
context('Testing stick breaking function')

set.seed(0)
x <- rbeta(5, 1, 1)

truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
sum(truth)
ltruth <- log(truth)

##-- test: use through nimble function
try(test_that("stick_breaking nimble function calculates weights correctly: ",
              expect_equal(stick_breaking(x, log=FALSE),
                           truth,
                           info = paste0("incorrect stick_breaking nimble function calculation"))))

try(test_that("stick_breaking nimble function calculates log weights correctly: ",
              expect_equal(stick_breaking(x, log=TRUE),
                           ltruth,
                           info = paste0("incorrect stick_breaking nimble function log calculation"))))

##-- test: use through compile nimble function
cSB <- compileNimble(stick_breaking)
try(test_that("compiled stick_breaking nimble function calculates weights correctly: ",
              expect_equal(cSB(x, log=FALSE),
                           truth,
                           info = paste0("incorrect compiled stick_breaking nimble function calculation"))))

try(test_that("compiled stick_breaking nimble function calculates log weights correctly: ",
              expect_equal(cSB(x, log=TRUE),
                           ltruth,
                           info = paste0("incorrect compiled stick_breaking nimble function log calculation"))))


##-- test: use in model
SB_code <- nimbleCode({
  for(i in 1:5) z[i] ~ dbeta(1, 1)
  w[1:6] <- stick_breaking(z[1:5])
})

set.seed(1)
Inits <- list(z = rbeta(5, 1, 1))
SB_model <- nimbleModel(SB_code, data=Inits)

SB_model$z <- x
SB_model$calculate()

# Chris, this test is not passed. Don't know why.... Later the compiled version is passed (test in line 71)!
                                        # maybe diference in decimals?
# Claudia, it's awkward becasue SB_model$w is an array not a vector. You'll need to convert to a vector before comparing to 'truth' 
try(test_that("Test that SB_model calculation is correct in weights calculation: ",
              expect_equal(SB_model$w, truth,
                           info = paste0("incorrect stick breaking weigths in model"))))

# Claudia, I'm not sure we need this - if the elements are all equal in the test above, then the mean is necessarily equal. Right?
try(test_that("Test that SB_model calculation is correct in mean of weights calculation: ",
              expect_equal(mean(SB_model$w), mean(truth), tol=0.01,
                           info = paste0("incorrect stick breaking mean weigths in model"))))

##-- test: use in  compile model
c_SB_model <- compileNimble(SB_model)
c_SB_model$z <- x
c_SB_model$calculate()
c_SB_model$w

try(test_that("Test that compiled SB_model calculation is correct in weights calculation: ",
              expect_equal(c_SB_model$w, truth,
                           info = paste0("incorrect stick breaking weigths in compiled model"))))


##-- test: random sampling from a compiled model adding one more level:
#-- definition of the model
SB_code2 <- nimbleCode({
  for(i in 1:5) 
    z[i] ~ dbeta(1, 1)
  w[1:6] <- stick_breaking(z[1:5])
  for(i in 1:10){
    xi[i] ~ dcat(w[1:6])
  }
})

set.seed(1)
Inits <- list(z = rbeta(5, 1, 1))
data <- list(xi = rep(1,10))
SB_model2 <- nimbleModel(SB_code2, data=data, inits=Inits)

c_SB_model2 <- compileNimble(SB_model2)

c_SB_model2$z <- x
c_SB_model2$calculate()

#-- checking some computations
try(test_that("Test that SB_model2 calculation is correct in weights calculation: ",
              expect_equal(c_SB_model2$w, truth,
                           info = paste0("incorrect stick breaking weigths in SB_model2"))))

#-- sampling via simulate:
set.seed(0)
simul_samp <- function(model) {
  model$simulate()
  return(model$w)
}
simul_samps <- t(replicate(10000, simul_samp(c_SB_model2)))

trueE <- c(0.5^(1:5) )

#-- checking the mean of the components of a vector that has a generalized dirichelt distribution
#-- if z_i ~ beta then weights, w, defined by a SB representation have a generalized dirichlet distribution
#-- and the expectation of w_j=0.5^j (in this case a=b=1).

try(test_that("Test that expectation of w_j is correct based on simulations: ",
              expect_equal(apply(simul_samps, 2, mean)[1:5], trueE, tol=0.01,
                           info = paste0("incorrect weights (w) sampling  in SB_model2"))))


##-- test: simple models
#-- 
model <- function() {
  for(j in 1:5) 
    z[j] ~ dbeta(1, 1)
  w[1:6] <- stick_breaking(z[1:5])
  for(i in 1:10){
    xi[i] ~ dcat(w[1:6])
  }
}

Inits <- list(z = rep(0.5,5))
Data <- list(xi = rep(1,10))

testBUGSmodel(example = 'test1', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)


#--
model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    y[i] ~ dnorm( thetatilde[xi[i]], var=1)
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=rep(1, 10))
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test2', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)


#--
model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    theta[i] <- thetatilde[xi[i]]
    y[i] ~ dnorm( theta[i], var=1)
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=rep(1, 10))
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test3', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)


#--
model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    theta[i] <- thetatilde[xi[i]]
    y[i] ~ dnorm( theta[i], var=s2tilde[xi[i]])
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
    s2tilde[i] ~ dinvgamma(1, 1)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=rep(1, 10), s2tilde=rep(1,10))
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test4', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)


##-- test: conjugacy detection for stick breaking parameter ,z:
test_that("Testing conjugacy detection with bnp stick breaking models", { 
  
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dbeta(1, 1)
      }
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[6]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
  # replace beta by uniform distribution: here no conjugacy is detected.
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dunif(0,1)
      }
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[6]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
  
  # concentration parameter added in beta distribution
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dbeta(1, conc)
      }
      conc ~ dgamma(1,1)
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4), conc=1))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[7]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
})

