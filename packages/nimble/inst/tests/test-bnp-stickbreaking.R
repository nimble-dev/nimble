
######################################################################
######################################################################


source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))


## Test stick_breaking nimble function:
context('Testing stick breaking function')


##-- test: use through nimble function
test_that("stick_breaking nimble function calculation and use is correct:", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  ltruth <- log(truth)
  
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation"))
  
  expect_equal(stick_breaking(x, log=TRUE),
               ltruth,
               info = paste0("incorrect stick_breaking nimble function log calculation"))
  
  cSB <- compileNimble(stick_breaking)
  
  expect_equal(cSB(x, log=FALSE),
               truth,
               info = paste0("incorrect compiled stick_breaking nimble function calculation"))
  
  expect_equal(cSB(x, log=TRUE),
               ltruth,
               info = paste0("incorrect compiled stick_breaking nimble function log calculation"))
  
  x <- c(0.1, 0.4, -0.1, 0.3)
  aux <- stick_breaking(x, log=FALSE)
  expect_equal(aux, rep(NaN, length(x)+1),
               info = "incorrect argument use (negative component) of stick breaking function")
  
  
  x <- c(0.1, 5, 0.4, 0.3)
  aux <- stick_breaking(x, log=FALSE)
  expect_equal(aux, rep(NaN, length(x)+1),
               info = "incorrect argument use (larger than 1 component) of stick breaking function")
  
  x <- c(0.1, 0.2, 0, 0.3, 0.8)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation with one 0 component"))
  
  x <- c(0.1, 0.2, 1, 0.3, 0.8)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation with one 1 component"))
})




##-- test: use in model
test_that("Stick breaking model calculation is correct:", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  
  SB_code <- nimbleCode({
    for(i in 1:5) z[i] ~ dbeta(1, 1)
    w[1:6] <- stick_breaking(z[1:5])
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  SB_model <- nimbleModel(SB_code, data=Inits)
  
  SB_model$z <- x
  SB_model$calculate()
  
  expect_equal(c(SB_model$w), truth,
               info = paste0("incorrect stick breaking weigths in model"))
  
  c_SB_model <- compileNimble(SB_model)
  
  c_SB_model$z <- x
  c_SB_model$calculate()
  c_SB_model$w
  
  expect_equal(c(c_SB_model$w), truth,
               info = paste0("incorrect stick breaking weigths in compiled model"))

})


##-- test: random sampling from a compiled model adding one more level:
test_that("random sampling from model works fine:", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  
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
  
  expect_equal(c_SB_model2$w, truth,
               info = paste0("incorrect stick breaking weigths in SB_model2"))
  
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
  
  expect_equal(apply(simul_samps, 2, mean)[1:5], trueE, tol=0.01,
               info = paste0("incorrect weights (w) sampling  in SB_model2"))
  
  
  # wrong specification of stick variables
  SB_code3 <- nimbleCode({
    for(i in 1:5) 
      z[i] ~ dgamma(10, 10)
    w[1:6] <- stick_breaking(z[1:5])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  # wrong prior and starting values for stick variables: how to recognize this wrong specification?
  set.seed(1)
  Inits <- list(z = rgamma(5, 10, 10))
  data <- list(xi = rep(1,10))
  expect_message(nimbleModel(SB_code3, data=data, inits=Inits)) # message is sent because z >1.
  
  
  
  # good stating values for stick variables
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  data <- list(xi = rep(1,10))
  SB_model3 <- nimbleModel(SB_code3, data=data, inits=Inits)
  expect_message(SB_model3$simulate(), message="values in 'z' have to be in") # message is sent because z >1.
  #expect_equal(c(SB_model3$w), rep(NaN, length(x)+1),
  #             info = "incorrect distribution for stick variables not identified")

  
  # wrong specification of length in stick variables, should be 5
  SB_code4 <- nimbleCode({
    for(i in 1:4) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:4])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  # wrong length for stick variables, a warning is sent in the nimbleModel function
  # how to recognize this wrong specification in the test?
  set.seed(1)
  Inits <- list(z = rbeta(4, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(nimbleModel(SB_code4, data=data, inits=Inits), message = "number of items to replace")
  
  
  # wrong specification of length in stick variables, should be 5
  # no warning in nimbleModel function
  SB_code5 <- nimbleCode({
    for(i in 1:2) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:2])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(2, 10, 10))
  data <- list(xi = rep(1,10))
  SB_model5 <- nimbleModel(SB_code5, data=data, inits=Inits)
  cSB_model5 <- compileNimble(SB_model5)
  expect_output(cSB_model5$calculate('w'), "Error in mapCopy")
  
  
  
  # longer vector of stick variables
  SB_code6 <- nimbleCode({
    for(i in 1:10) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:10])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  # wrong length for stick variables, a warning is sent in the nimbleModel function
  # how to recognize this wrong specification in the test?
  set.seed(1)
  Inits <- list(z = rbeta(10, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(nimbleModel(SB_code6, data=data, inits=Inits), message = "number of items to replace")
})



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



 # test of BNP model

test_that("Testing BNP model using stick breaking representation", { 
  
Code=nimbleCode(
  {
    for(i in 1:Trunc) {
      thetatilde[i] ~ dnorm(mean=0, var=40) 
      s2tilde[i] ~ dinvgamma(shape=1, scale=0.5) 
    }
    for(i in 1:(Trunc-1)) {
      z[i] ~ dbeta(1, 1)
    }
    w[1:Trunc] <- stick_breaking(z[1:(Trunc-1)])
    
    for(i in 1:N) {
      xi[i] ~ dcat(w[1:Trunc])
      theta[i] <- thetatilde[xi[i]]
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=s2[i])
    }
  }
)

Consts <- list(N=50, Trunc=25)
set.seed(1)
Inits <- list(thetatilde = rnorm(Consts$Trunc, 0, sqrt(40)),
           s2tilde = rinvgamma(Consts$Trunc, shape=1, scale=0.5),
           z = rbeta(Consts$Trunc-1, 1, 1),
           xi = sample(1:10, size=Consts$N, replace=TRUE))
Data = list(y = c(rnorm(Consts$N/2,5,sqrt(4)), rnorm(Consts$N/2,-5,sqrt(4))))

model = nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel = compileNimble(model)


#-- MCMC configuration:
modelConf = configureMCMC(model, print=FALSE, thin=100)

expect_match(modelConf$getSamplers()[[51]]$name, "conjugate_dbeta_dcat",
             info = "failed to detect categorical-beta conjugacy in BNP model")

modelMCMC = buildMCMC(modelConf)
CmodelMCMC = compileNimble(modelMCMC, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)

CmodelMCMC$run(10000)

#-- results from the algorithm:
samples = as.matrix(CmodelMCMC$mvSamples)
s2Sam = samples[, 1:25]
thetaSam = samples[, 26:50]
zSam = samples[, 51:74]
Tr = 25
Wpost = t(apply(zSam, 1, function(x)c(x[1], x[2:(Tr-1)]*cumprod(1-x[1:(Tr-2)]), cumprod(1-x[1:(Tr-1)])[N=Tr-1])))

# grid:
ngrid = 302
grid = seq(-10, 25,len=ngrid)

# posterior samples of the density
nsave = 100
predSB = matrix(0, ncol=ngrid, nrow=nsave)
for(i in 1:nsave) {
  predSB[i, ] = sapply(1:ngrid, function(j)sum(Wpost[i, ]*dnorm(grid[j], thetaSam[i,],sqrt(s2Sam[i,]))))
}

#hist(Data$y, freq=FALSE, xlim=c(min(grid), max(grid)))
#points(grid, apply(predSB, 2, mean), col="blue", type="l", lwd=2)

f0 <- function(x) 0.5*dnorm(x,5,sqrt(4)) + 0.5*dnorm(x,-5,sqrt(4))
fhat <- apply(predSB, 2, mean)
f0grid <- sapply(grid, f0)

L1dist <- mean(abs(f0grid - fhat))

expect_equal(L1dist, 0.01, tol=0.01,
             info = "wrong estimation of density in DPM of normal distrbutions")

})





#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

## Testing CRP distribution

context('Testing NIMBLE CRP distribution')

test_that("dCRP nimble function calculates density correctly: ",{
  
  x <- c(1,1,2,1,1,2)
  conc <- 1
  
  truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
    (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))
  ltruth <- log(truth)
  
  expect_equal(dCRP(x, conc, size=length(x), log=FALSE),
               truth,
               info = paste0("incorrect dCRP nimble function calculation"))
  
  expect_equal(dCRP(x, conc, size=length(x), log=TRUE),
               ltruth,
               info = paste0("incorrect dCRP nimble function calculation in log scale"))
  
  cdCRP <- compileNimble(dCRP)
  
  expect_equal(cdCRP(x, conc, size=length(x)), (truth), 
               info = paste0("incorrect dCRP value in compiled nimble function")))
  
  expect_equal(cdCRP(x, conc, size=length(x), log=TRUE), (ltruth), 
             info = paste0("incorrect dCRP value in compiled nimble function  in log scale"))
  
  expect_equal(dCRP(x, conc=-1, size=length(x), log=FALSE),
               NaN,
               info = paste0("incorrect parameters space allowed"))
  
  expect_error(dCRP(x, conc=1, size=3, log=FALSE), "length of x has to be equal to size")
  
  expect_error(dCRP(x, conc=1, size=10, log=FALSE), "length of x has to be equal to size")
  
})


##-- test: use in model:
test_that("CRP model calculation and dimensionas are correct:", {
  
  x <- c(1,1,2,1,1,2)
  conc <- 1
  
  truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
    (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))
  ltruth <- log(truth)
  
  CRP_code <- nimbleCode({
    x[1:6] ~ dCRP(conc, size=6)
  })
  
  Consts <- list(conc = 1)
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model <- nimbleModel(CRP_code, data=Inits, constants=Consts)
  
  CRP_model$x <- x
  expect_equal(exp(CRP_model$calculate()), truth,
               info = paste0("incorrect likelihood value for dCRP"))
  
  c_CRP_model <- compileNimble(CRP_model)
  c_CRP_model$x
  expect_equal(exp(c_CRP_model$calculate()), truth,
               info = paste0("incorrect likelihood value for compiled dCRP"))
  
  
  # different length of x and size:
  CRP_code2 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=10)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model2 <- nimbleModel(CRP_code2, data=Inits)
  expect_error(CRP_model2$calculate(), "length of x has to be equal to size")
  
  # different length of x and size:
  CRP_code3 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=3)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model3 <- nimbleModel(CRP_code3, data=Inits)
  expect_error(CRP_model3$calculate(), "length of x has to be equal to size")
  
  
})

  
##-- test: random sampling from a compiled model adding one more level:
test_that("random sampling from CRP and model works fine:", {
  
  conc <- 1
  set.seed(0)
  size <- 6
  r_samps <- t(replicate(10000, rCRP(n = 1, conc, size = size)))
  # K is the number of unique components in x of length 6
  true_EK <- sum(conc/(conc+1:size-1))
  
  expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K exceeds tolerance")
  
  # sampling from the model:
  set.seed(1)
  CRP_code <- nimbleCode({
    x[1:6] ~ dCRP(conc=1, size=6)
    for(i in 1:6){
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[x[i]], 1)
    }
  })
  Inits <- list(x = c(1,1,2,1,1,2), mu = 1:6)
  Data <- list( y =  rnorm(6))
  CRP_model <- nimbleModel(CRP_code, data=Data, inits=Inits)
  c_CRP_model <- compileNimble(CRP_model)
  
  simul_samp <- function(model) {
    model$simulate()
    return(model$x)
  }
  simul_samps <- t(replicate(10000, simul_samp(c_CRP_model)))
  
  expect_equal(mean(apply(simul_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K, from compiled model, exceeds tolerance")
  
})
  



context("Testing of bnp conjugacy")

test_that("Testing conjugacy detection with models using CRP: ", { 
  
  # dnorm_dnorm
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")

  
  # dgamma_dpois
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dpois(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rpois(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  
  # dgamma_dexp
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dexp(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  
  # dgamma_dgamma
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      y[i] ~ dgamma(4, mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rgamma(4, 4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dgamma")
  
 
  # dbeta_dbern
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dbeta(1,1)
      y[i] ~ dbern(mu[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rbinom(4, size=1, prob=0.5)),
                  inits = list(xi = rep(1,4), mu=rbeta(4, 1, 1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbern")
  
  
  # dbeta_dbern
  code=nimbleCode(
    {
      for(i in 1:4){
        p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  set.seed(1)
  p0 <- matrix(0, ncol=3, nrow=4)
  y0 <- matrix(0, ncol=3, nrow=4)
  for(i in 1:4){
    p0[i,]=rdirch(1, c(1, 1, 1))
    y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
  }
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = rep(1,4), p=p0), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[5]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  # dnorm_dnorm_dinvgamma
  code = nimbleCode({
    for(i in 1:4) {
      s2[i] ~ dinvgamma(1, 1)
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), s2=rinvgamma(4, 1,1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[9]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
})

#---- OLD TESTS:
#set.seed(0)
#x <- rbeta(5, 1, 1)

#truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
#sum(truth)
#ltruth <- log(truth)

#try(test_that("stick_breaking nimble function calculates weights correctly: ",
#              expect_equal(stick_breaking(x, log=FALSE),
#                           truth,
#                           info = paste0("incorrect stick_breaking nimble function calculation"))))

#try(test_that("stick_breaking nimble function calculates log weights correctly: ",
#              expect_equal(stick_breaking(x, log=TRUE),
#                           ltruth,
#                           info = paste0("incorrect stick_breaking nimble function log calculation"))))

##-- test: use through compile nimble function

#cSB <- compileNimble(stick_breaking)
#try(test_that("compiled stick_breaking nimble function calculates weights correctly: ",
#              expect_equal(cSB(x, log=FALSE),
#                           truth,
#                           info = paste0("incorrect compiled stick_breaking nimble function calculation"))))

#try(test_that("compiled stick_breaking nimble function calculates log weights correctly: ",
#              expect_equal(cSB(x, log=TRUE),
#                          ltruth,
#                           info = paste0("incorrect compiled stick_breaking nimble function log calculation"))))


#SB_code <- nimbleCode({
#  for(i in 1:5) z[i] ~ dbeta(1, 1)
#  w[1:6] <- stick_breaking(z[1:5])
#})

#set.seed(1)
#Inits <- list(z = rbeta(5, 1, 1))
#SB_model <- nimbleModel(SB_code, data=Inits)

#SB_model$z <- x
#SB_model$calculate()

# Chris, this test is not passed. Don't know why.... Later the compiled version is passed (test in line 71)!
# maybe diference in decimals?
# Claudia, it's awkward becasue SB_model$w is an array not a vector. You'll need to convert to a vector before comparing to 'truth' 
#try(test_that("Test that SB_model calculation is correct in weights calculation: ",
#              expect_equal(SB_model$w, truth,
#                           info = paste0("incorrect stick breaking weigths in model"))))

# Claudia, I'm not sure we need this - if the elements are all equal in the test above, then the mean is necessarily equal. Right?
# If the above test works we don't need this! 
#try(test_that("Test that SB_model calculation is correct in mean of weights calculation: ",
#              expect_equal(mean(SB_model$w), mean(truth), tol=0.01,
#                           info = paste0("incorrect stick breaking mean weigths in model"))))


#try(test_that("Test that compiled SB_model calculation is correct in weights calculation: ",
#              expect_equal(c_SB_model$w, truth,
#                           info = paste0("incorrect stick breaking weigths in compiled model"))))


#SB_code2 <- nimbleCode({
#  for(i in 1:5) 
#    z[i] ~ dbeta(1, 1)
#  w[1:6] <- stick_breaking(z[1:5])
#  for(i in 1:10){
#    xi[i] ~ dcat(w[1:6])
#  }
#})

#set.seed(1)
#Inits <- list(z = rbeta(5, 1, 1))
#data <- list(xi = rep(1,10))
#SB_model2 <- nimbleModel(SB_code2, data=data, inits=Inits)

#c_SB_model2 <- compileNimble(SB_model2)

#c_SB_model2$z <- x
#c_SB_model2$calculate()

#-- checking some computations
#try(test_that("Test that SB_model2 calculation is correct in weights calculation: ",
#              expect_equal(c_SB_model2$w, truth,
#                           info = paste0("incorrect stick breaking weigths in SB_model2"))))

#-- sampling via simulate:
#set.seed(0)
#simul_samp <- function(model) {
#  model$simulate()
#  return(model$w)
#}
#simul_samps <- t(replicate(10000, simul_samp(c_SB_model2)))

#trueE <- c(0.5^(1:5) )

#-- checking the mean of the components of a vector that has a generalized dirichelt distribution
#-- if z_i ~ beta then weights, w, defined by a SB representation have a generalized dirichlet distribution
#-- and the expectation of w_j=0.5^j (in this case a=b=1).

#try(test_that("Test that expectation of w_j is correct based on simulations: ",
#              expect_equal(apply(simul_samps, 2, mean)[1:5], trueE, tol=0.01,
#                           info = paste0("incorrect weights (w) sampling  in SB_model2"))))

