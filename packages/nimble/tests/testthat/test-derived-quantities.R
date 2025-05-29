source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

code <- nimbleCode({
  b[1] ~ dnorm(0, sd = 2)
  b[2] ~ dnorm(0, sd = 2)
  sigma ~ dunif(0, 2)
  for (i in 1:n){
    mu[i] <- b[1] + b[2]*x[i]
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
})

test_that("mean derived parameter", {

  # Single derived parameter (mean)
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes="b[1]")
  expect_is(conf$derivedConfs[[1]], "derivedConf")

  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)

  expect_is(out, "list")
  expect_equal(names(out), c("samples", "derived"))
  expect_equivalent(out$derived$mean[1,1], 
                  mean(out$samples[,"b[1]"]))

  # Means for two parameters
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]", "b[2]"))
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$mean[1,],
             colMeans(out$samples[,1:2]))

  # Using range
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1:2]"))
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$mean[1,],
             colMeans(out$samples[,1:2]))

  # Combination of two separate parameters
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=list(c("b[1]", "sigma")))
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(round(out$derived$mean[1,], 4), round(apply(out$samples[,c("b[1]", "sigma")], 2, mean),4))

  # What if the parameter doesn't exist in the model?
  # result is empty matrix
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("alpha"), interval=11)
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$mean, matrix(0, 1, 0))
  
  # If parameter index is too big
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[3]"), interval=11)
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$mean, matrix(0, 1, 0))
})

test_that("interval and thinning rate for derived parameter", {

  # Changing interval
  # every other mcmc iteration will be used to calculate mean (starting with 2)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=2)
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(mean(out$samples[c(2,4,6,8,10),1]),
                  out$derived$mean[1,1])

  # Try interval of 3
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=3)
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(mean(out$samples[c(3,6,9),1]),
                  out$derived$mean[1,1])

  # What happens if interval is larger than sample?
  # output is NA
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=11)
  mcmc <- buildMCMC(conf)
  set.seed(1)
  out <- runMCMC(mcmc, niter=10)
  expect_true(is.na(out$derived$mean[1,1]))

  # How does interval interact with thinning rate?
  # First create unthinned samples for reference
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  mcmc <- buildMCMC(conf)
  out_unthinned <- runMCMC(mcmc, niter=10)

  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 2, print=FALSE)
  # Without setting interval; this should set interval to the same as thin automatically
  conf$addDerivedQuantity("mean", nodes=c("b[1]"))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(mean(out_unthinned[c(2,4,6,8,10), 1]),
                  out$derived$mean[1,1])

  # What if we manually set interval to 1?
  # then all samples (thinned and  unthinned) are used in mean calculation
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 2, print=FALSE)
  # Without setting interval; this should set interval to the same as thin automatically
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=1)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(mean(out_unthinned[, 1]),
                  out$derived$mean[1,1])

  # Setting different thin and interval rates
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 2, print=FALSE)
  # Without setting interval; this should set interval to the same as thin automatically
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=3)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(mean(out_unthinned[c(3,6,9), 1]),
                  out$derived$mean[1,1])
})

test_that("recording frequency for derived parameters", {
  # Testing recording frequency
  # The default 0 (calculated once at the end)
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 1, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=1, recordingFrequency=0)
  mcmc <- buildMCMC(conf)
  out_ref <- runMCMC(mcmc, niter=10)
  expect_equivalent(mean(out_ref$samples[,1]), out_ref$derived$mean[1,1])

  # frequency 1 (cumulative mean)
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 1, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=1, recordingFrequency=1)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(cumsum(out$samples[,1]) / seq_along(out$samples[,1]), 
                  out$derived$mean[,1])

  # Vary the thinning rate
  # This should not affect derived quantity results when we set interval and frequency to 1
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 2, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=1, recordingFrequency=1)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  # Comparing to the reference output above with no thinning
  expect_equivalent(cumsum(out_ref$samples[,1]) / seq_along(out_ref$samples[,1]), 
                  out$derived$mean[,1])

  # frequency = 2
  # cumulative mean but calculate it half as often
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 1, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=1, recordingFrequency=2)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(c(mean(out$samples[1:2]), mean(out$samples[1:4]), mean(out$samples[1:6]),
                  mean(out$samples[1:8]), mean(out$samples[1:10])),
                  out$derived$mean[,1])

  # interval = 2 and frequency = 2
  # only use every other sample and only calculate after 2 recorded samples
  # result should be length 2
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, thin = 1, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=2, recordingFrequency=2)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(out$derived$mean,
    c(mean(out$samples[c(2,4),1]), mean(out$samples[c(2,4,6,8)])))

  # What if recording frequency is greater than number of iterations?
  # Output is numeric(0)
  conf <- configureMCMC(mod, thin = 1, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=1, recordingFrequency=11)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(out$derived$mean, numeric(0))

  # Same thing in this case when combo of interval + recording frequency means 
  # there aren't any derived values
  conf <- configureMCMC(mod, thin = 1, print=FALSE)
  conf$addDerivedQuantity("mean", nodes=c("b[1]"), interval=6, recordingFrequency=2)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(out$derived$mean, numeric(0))
})

test_that("mean derived parameter with compiled mcmc", {
  # Check compiled version works as expected
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes="b[1]")
  mcmc <- buildMCMC(conf)
  Cmod <- compileNimble(mod)
  Cmcmc <- compileNimble(mcmc, project=mod)
  out <- runMCMC(Cmcmc, niter=10)
  expect_equivalent(mean(out$samples[,1]), out$derived$mean[1,1]) 
})

test_that("variance derived parameter", {
  # Variance
  # fewer tests here since it should have same behavior as mean
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("variance", nodes="b[1]")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(var(out$samples[,1]), out$derived$variance[1,1]) 

  # Check compiled version
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("variance", nodes="b[1]")
  mcmc <- buildMCMC(conf)
  Cmod <- compileNimble(mod)
  Cmcmc <- compileNimble(mcmc, project=mod)
  out <- runMCMC(Cmcmc, niter=10)
  expect_equivalent(var(out$samples[,1]), out$derived$variance[1,1]) 

  # Combination of both mean and variance
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("mean", nodes="b[1]")
  conf$addDerivedQuantity("variance", nodes="b[1]")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equivalent(out$derived$mean, mean(out$samples[,1]))
  expect_equivalent(out$derived$variance, var(out$samples[,1]))
})

test_that("logProb derived quantity", {

  set.seed(1)
  inits <- list(b=c(0, 0), sigma=1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=inits)
  conf <- configureMCMC(mod, print=FALSE)

  # Single node
  conf$addDerivedQuantity("logProb", nodes="b[1]")
  mcmc <- buildMCMC(conf)

  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$logProb), c(10,1))
  expect_equivalent(out$derived$logProb[10,"b[1]"], mod$logProb_b[1])
 
  # non-scalar parameter
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("logProb", nodes="b")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), c("b[1]", "b[2]"))
  expect_equivalent(out$derived$logProb[10,], mod$logProb_b)

  # Multiple parameters
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("logProb", nodes=c("b[1]", "sigma"))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), c("b[1]", "sigma"))
  expect_equivalent(out$derived$logProb[10,], c(mod$logProb_b[1], mod$logProb_sigma))

  # Multiple parameters including non-scalars
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("logProb", nodes=c("b", "sigma"))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), c("b[1]", "b[2]", "sigma"))
  expect_equivalent(out$derived$logProb[10,], c(mod$logProb_b, mod$logProb_sigma))

  # All parameters
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # Here .all means the sum of logprobs of all parameters
  conf$addDerivedQuantity("logProb", nodes=".all")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), c("_all_nodes_"))
  expect_equivalent(out$derived$logProb[10,1], sum(c(mod$logProb_b, mod$logProb_sigma, mod$logProb_y)))

  # Functions of parameters
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # Here taking sum logprobs of b
  conf$addDerivedQuantity("logProb", nodes=list("b", "sigma"))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), c("b", "sigma"))
  expect_equivalent(out$derived$logProb[10,], c(sum(mod$logProb_b), mod$logProb_sigma))

  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # Here taking sum of logprobs of b[1] and sigma
  conf$addDerivedQuantity("logProb", nodes=list(c("b[1]", "sigma")))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), "sum1")
  expect_equivalent(out$derived$logProb[10,1], sum(c(mod$logProb_b[1], mod$logProb_sigma)))

  # Create name for sum
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # Here taking sum of logprobs of b[1] and sigma
  conf$addDerivedQuantity("logProb", nodes=list(b_and_sigma=c("b[1]", "sigma")))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(colnames(out$derived$logProb), "b_and_sigma")
  expect_equivalent(out$derived$logProb[10,1], sum(c(mod$logProb_b[1], mod$logProb_sigma)))

  # What happens if parameter doesn't exist? An empty matrix is returned
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # alpha not in model
  conf$addDerivedQuantity("logProb", nodes="alpha")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$logProb, matrix(0, 10, 0))

  # What if an impossible index range is requested?
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # alpha not in model
  conf$addDerivedQuantity("logProb", nodes="b[1:3]")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  # The extra element of b is ignored
  expect_equal(colnames(out$derived$logProb), c("b[1]", "b[2]"))

  # Request only the non-existing index - result is empty matrix
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # alpha not in model
  conf$addDerivedQuantity("logProb", nodes="b[3]")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$logProb, matrix(0, 10, 0))

  # Change interval
  set.seed(1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE)
  # Here taking sum of logprobs of b[1] and sigma
  conf$addDerivedQuantity("logProb", nodes="b[1]", interval=5)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$logProb), c(2,1))
})

test_that("predictive derived quantity", {

  code2 <- nimbleCode({
    b[1] ~ dnorm(0, sd = 2)
    b[2] ~ dnorm(0, sd = 2)
    sigma ~ dunif(0, 2)
    for (i in 1:n){
      mu[i] <- b[1] + b[2]*x[i]
      y[i] ~ dnorm(mu[i], sd = sigma)
    }

    bsq <- b[1]^2 # derived quantity
    z ~ dnorm(bsq, sd = sigma) # posterior predictive node
  })

  # Simple example
  set.seed(1)
  mod <- nimbleModel(code2, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE, samplePredictiveNodes=FALSE)

  conf$addDerivedQuantity("predictive", nodes=c("bsq", "z"))
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$predictive), c(10,2))
  expect_equal(colnames(out$derived$predictive), c("bsq", "z"))

  # Only the predictive node
  set.seed(1)
  mod <- nimbleModel(code2, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE, samplePredictiveNodes=FALSE)

  conf$addDerivedQuantity("predictive", nodes="z")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$predictive), c(10,1))
  expect_equal(colnames(out$derived$predictive), "z")

  # Change the interval
  set.seed(1)
  mod <- nimbleModel(code2, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE, samplePredictiveNodes=FALSE)

  conf$addDerivedQuantity("predictive", nodes="z", interval=5)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$predictive), c(2,1))
  expect_equal(colnames(out$derived$predictive), "z")

  # Set a non-predictive node
  set.seed(1)
  mod <- nimbleModel(code2, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE, samplePredictiveNodes=FALSE)

  conf$addDerivedQuantity("predictive", nodes=c("b[1]", "z"))
  # warning message here
  expect_message(mcmc <- buildMCMC(conf))
  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$predictive), c(10, 2))
  expect_equal(colnames(out$derived$predictive), c("b[1]", "z"))

  # Set .all
  set.seed(1)
  mod <- nimbleModel(code2, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE, samplePredictiveNodes=FALSE)

  # Check that samples of both bsq and z are returned when setting .all
  conf$addDerivedQuantity("predictive", nodes=".all")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(dim(out$derived$predictive), c(10,2))
  expect_equal(colnames(out$derived$predictive), c("bsq", "z"))

  # Non-existing parameter - empty matrix is returned
  set.seed(1)
  mod <- nimbleModel(code2, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=list(b=c(0,0), sigma=1))
  conf <- configureMCMC(mod, print=FALSE, samplePredictiveNodes=FALSE)
  conf$addDerivedQuantity("predictive", nodes="alpha")
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)
  expect_equal(out$derived$predictive, matrix(0, 10, 0))

})

test_that("user-facing mcmc configuration methods", {

  set.seed(1)
  inits <- list(b=c(0, 0), sigma=1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=inits)
  conf <- configureMCMC(mod, print=FALSE)

  # Add node
  conf$addDerivedQuantity("logProb", nodes="b[1]")
  mcmc <- buildMCMC(conf)

  # Print
  co <- capture.output(conf$printDerivedQuantities())
  expect_equal(co, "[1] derived quantity: logProb,  execution interval: thin,  nodes: b[1]")
  
  # Get
  conf$addDerivedQuantity("logProb", nodes="sigma")
  dq <- conf$getDerivedQuantities()
  expect_equal(length(dq), 2)
  expect_is(dq[[1]], "derivedConf")

  dqf <- conf$getDerivedQuantityDefinition(1)
  expect_is(dqf, "list")

  # Remove nodes
  conf2 <- conf
  conf2$removeDerivedQuantities()
  co <- capture.output(conf2$printDerivedQuantities())
  expect_equal(co, character(0))

  # Remove nodes
  conf2 <- conf
  conf2$removeDerivedQuantity()
  co <- capture.output(conf2$printDerivedQuantities())
  expect_equal(co, character(0))
})

test_that("configureMCMC and nimbleMCMC derived quantity interface", {

  set.seed(1)
  inits <- list(b=c(0, 0), sigma=1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=inits)
  conf <- configureMCMC(mod, print=FALSE, mean="b[1]", logProb="sigma")
  co <- capture.output(conf$printDerivedQuantities())
  expect_equal(co, c("[1] derived quantity: mean,  execution interval: thin,  nodes: b[1]",    
                     "[2] derived quantity: logProb,  execution interval: thin,  nodes: sigma"))

  # logprob=TRUE
  set.seed(1)
  inits <- list(b=c(0, 0), sigma=1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=inits)
  conf <- configureMCMC(mod, print=FALSE, logProb=TRUE)
  co <- capture.output(conf$printDerivedQuantities())
  expect_equal(co, "[1] derived quantity: logProb,  execution interval: thin,  nodes: .all")    

  # nimbleMCMC
  mod <- nimbleMCMC(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                    mean="b[1]", logProb="sigma", nchain=1, niter=10, nburnin=0) 
  expect_equal(names(mod$derived), c("mean", "logProb"))

  mod <- nimbleMCMC(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                    logProb=TRUE, nchain=1, niter=10, nburnin=0) 
  expect_equal(dim(mod$derived$logProb), c(10,1))
})

test_that("custom derived quantity nimbleFunction", {
  # multiplies the specified nodes by some value multValue
  test_derived <- nimbleFunction(
    name = 'test_derived',
    contains = derived_BASE,     ## all derived functions must contain (inherit from) the derived_BASE class
    setup = function(model, mcmc, interval, control) {
      nodes <- extractControlElement(control, 'nodes', defaultValue = character())
      multValue <- extractControlElement(control, 'multValue', defaultValue = 2)
      ## node list generation
      nodes <- model$expandNodeNames(nodes)
      ## names generation
      names <- if(length(nodes) < 2) c(nodes,'','') else nodes     ## vector
      ## numeric value generation
      nResults <- length(nodes)
      vals <- numeric(max(nResults, 2))
      results <- array(0, c(1, nResults))
    },
    run = function(timesRan = double()) {
      if(nResults == 0)   return()
      vals <<- values(model, nodes)
      results[1,] <<- vals * multValue
    },
    methods = list(
      set_interval = function(newInterval = double()) {
        interval <<- newInterval
      },
      get_results = function() {
        returnType(double(2))
        return(results)
      },
      get_names = function() {
        returnType(character(1))
        return(names)
      }
    )
  )
  temporarilyAssignInGlobalEnv(test_derived)

  set.seed(1)
  inits <- list(b=c(0, 0), sigma=1)
  mod <- nimbleModel(code, constants=list(n=10, x=rnorm(10)), data = list(y=rnorm(10)),
                     inits=inits)
  conf <- configureMCMC(mod, print=FALSE)
  conf$addDerivedQuantity("test_derived", nodes=c("b[1]", "sigma"), multValue=2)
  mcmc <- buildMCMC(conf)
  out <- runMCMC(mcmc, niter=10)

  expect_equal(out$samples[10,c(1,3)]*2, out$derived$test_derived[1,])

  # Check in C++
  Cmod <- compileNimble(mod)
  Cmcmc <- compileNimble(mcmc, project=mod)
  out <- runMCMC(Cmcmc, niter=10)
  expect_equal(out$samples[10,c(1,3)]*2, out$derived$test_derived[1,])

  # check error when custom derived quantity function not defined
  expect_error(conf$addDerivedQuantity("fake", nodes=c("b[1]", "sigma")))
})
