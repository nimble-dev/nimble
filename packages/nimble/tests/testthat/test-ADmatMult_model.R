source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

nimbleOptions(useADmatMultAtomic = TRUE)

# This model does a matrix mult and uses the results.
mc <- nimbleCode({
  for(i in 1:nCov) {
    beta[i] ~ dnorm(0, sd = 100)
  }
  sigma ~ dunif(0, 100) # This serves as a parameter connected to a lot
  p ~ dnorm(0, sd = 100) # Similar role
  pred_y[1:nObs] <- p*(X[1:nObs, 1:nCov] %*% beta[1:nCov])[,1]
  for(i in 1:nObs)
    y[i] ~ dnorm(pred_y[i], sd = sigma)
})

# Set up model inputs
set.seed(1)
nCov <- 5
nObs <- 7
constants <- list(nCov = nCov, nObs = nObs)
data <- list()

X <- matrix(runif(nCov * nObs, min = 1, max = 3), nrow = nObs)
beta <- runif(nCov, min = 1, max = 3)
y <- rnorm(nObs, mean = X %*% beta, sd = 5)
inits <- list(beta = rnorm(nCov), sigma = 5, p = .8, y = y, X = X)

X2 <- matrix(runif(nCov * nObs, min = 1, max = 3), nrow = nObs)
beta2 <- runif(nCov, min = 1, max = 3)
y2 <- rnorm(nObs, mean = X2 %*% beta2, sd = 5)
inits2 <- list(beta = rnorm(nCov), sigma = 6, p = .7, y = y2, X = X2)

# build and compile model
m <- nimbleModel(mc, constants = constants, data = data, inits = inits,
                 buildDerivs = TRUE)
cm <- compileNimble(m)

# Make test function compatible with test_AD2
argTypes <- list(arg1 = "double(1)")
op <- list(
  expr = quote( {
    values(model, derivNodes) <<- arg1
    out <- model$calculate(calcNodes)
    return(out)
  }),
  args = list(arg1 = quote(double(1))),
  outputType = quote(double())
)

matMultTest_pieces <- make_AD_test2(op = op, argTypes = argTypes, includeModelArgs = TRUE)

matMultTest <- nimbleFunction(
  setup = function(model, nodesList) {
    derivNodes <- nodesList$derivNodes
    updateNodes <- nodesList$updateNodes
    constantNodes <- nodesList$constantNodes
    calcNodes <- nodesList$calcNodes
    nNodes <- length(derivNodes)
  },
  run = matMultTest_pieces$run,
  methods = matMultTest_pieces$methods,
  buildDerivs = matMultTest_pieces$buildDerivs
)

## Create test cases of interest
## Part of what atomic matrix multiplication does is
## carve out an inner square of CppAD variables from
## possible CppAD constants around it.  CppAD constants are
## values that are not being used through derivative tracking.
## The way we manipulate different situations is through changes in
## what are constantNodes and what are updateNodes

setup_update_and_constant_nodes <- function(model,
                                            derivNodes,
                                            forceConstantNodes = character(),
                                            forceUpdateNodes = character()) {
  ## "update" means "CppAD dynamic"
  ##  derivNodes <- model$expandNodeNames(derivNodes) # do not do this because do not want vector node names
  nNodes <- length(derivNodes)
  calcNodes <- model$getDependencies(derivNodes)
  ucNodes <- makeModelDerivsInfo(m, derivNodes, calcNodes, dataAsConstantNodes = TRUE)
  updateNodes <- ucNodes$updateNodes
  constantNodes <- ucNodes$constantNodes
  updateNodes <- setdiff(updateNodes, forceConstantNodes) # remove forceConstants from updates
  constantNodes <- setdiff(constantNodes, forceUpdateNodes) # remove forceUpdates from constants
  constantNodes <- union(constantNodes, forceConstantNodes) # add forceConstants to constants
  updateNodes <- union(updateNodes, forceUpdateNodes)     # add forceUpdates to updates
  list(derivNodes = derivNodes, updateNodes = updateNodes,
       constantNodes = constantNodes, calcNodes = calcNodes)
}

# X is dynamic.  beta is part dynamic, part variable.
nodesList_case1 <- setup_update_and_constant_nodes(m, c('beta[3]', 'beta[4]'))
v1_case1 <- list(arg1 = c(0.3, 0.4))
v2_case1 <- list(arg1 = c(0.5, 0.6))

# Both X and beta have part variable
nodesList_case2 <- setup_update_and_constant_nodes(m, c('X[2, 2]', 'beta[3]', 'beta[4]'))
v1_case2 <- list(arg1 = c(0.2, 0.3, 0.4))
v2_case2 <- list(arg1 = c(0.4, 0.5, 0.6))

# X is part variable. beta is dynamic
nodesList_case3 <- setup_update_and_constant_nodes(m, c('X[2, 2]', 'X[2, 3]'))
v1_case3 <- list(arg1 = c(0.2, 0.3))
v2_case3 <- list(arg1 = c(0.4, 0.5))

# p is variable.  X and beta are dynamic
nodesList_case4 <- setup_update_and_constant_nodes(m, c('p'))
v1_case4 <- list(arg1 = .5)
v2_case4 <- list(arg1 = .6)



## F unction to create nf instance with one set of
## constantNodes and updateNodes and then check inputs
## for it.
checkCase <- function(Rmodel, Cmodel,
                      deriv_nf, nodesList,
                      v1, v2,
                      order,
                      varValues = list(),
                      varValues2 = list(), ...) {
  # This sets up an instance of a checker fxn
  Rfxn <- deriv_nf( Rmodel, nodesList)
  Cfxn <- compileNimble(Rfxn, project = Rmodel)
  
  test_AD2_oneCall(Rfxn, Cfxn,
                   recordArgs = v1, testArgs = v2,
                   order = order,
                   Rmodel = Rmodel, Cmodel = Cmodel,
                   recordInits = varValues, testInits = varValues2,
                   nodesToChange = c(nodesList$updateNodes),
                   ...)
}

##
resetTols()
RCrelTol <- c(1e-15, 1e-6, 1e-2)

checkCase(m, cm, matMultTest,
          nodesList_case1, v1_case1, v2_case1,
          order = c(0, 1, 2),
          varValues = inits, varValues2 = inits2,
          RCrelTol = RCrelTol)

checkCase(m, cm, matMultTest,
          nodesList_case2, v1_case2, v2_case2,
          order = c(0, 1, 2),
          varValues = inits, varValues2 = inits2,
          RCrelTol = RCrelTol)

checkCase(m, cm, matMultTest,
          nodesList_case3, v1_case3, v2_case3,
          order = c(0, 1, 2),
          varValues = inits, varValues2 = inits2,
          RCrelTol = RCrelTol) 

checkCase(m, cm, matMultTest,
          nodesList_case4, v1_case4, v2_case4,
          order = c(0, 1, 2),
          varValues = inits, varValues2 = inits2,
          RCrelTol = RCrelTol)


