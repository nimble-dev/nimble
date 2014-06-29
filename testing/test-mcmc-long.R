# running a few BUGS examples via C++ for longer runs

# currently only works with loadAllCode() (error with compileNimbleFunction w/ R CMD INSTALL)
source('loadAllCode.R')

useInits <- TRUE

example <- "pump"; vol <- "vol1"
#example <- "blocker"; vol <- "vol1"
#example <- "bones"; vol <- "vol1"; useInits <- FALSE # for some reason bones-init overwrites 'grade'
#example <- "dugongs"; vol <- "vol2"

# assumes run from examples/demos if not using R CMD INSTALL
if(!"package:nimble" %in% search()) {
  dir <- file.path("..", "..", "packages", "nimble", "inst", "classic-bugs")
  dir <- file.path(dir, vol, example)
} else {
  dir <- file.path(system.file("classic-bugs", package = "nimble"), vol, example)
}


Rmodel <- readBUGSmodel(example, dir = dir, useInits = useInits)
if(example == "dugongs") simulate(Rmodel, 'gamma')

Cmodel <- compileBUGSmodel(Rmodel)

stochNodes <- Rmodel$getNodeNames(stochOnly = TRUE, includeData = FALSE)
stochVars <- unique(gsub("\\[.*\\]", "", stochNodes))
mcmcspec <- MCMCspec(Rmodel, control = list(adaptInterval=100))
mcmcspec$addMonitors(stochVars)

Rmcmc <- buildMCMC(mcmcspec)
Cmcmc <- compileNimbleFunction(buildMCMC)

RmvSample  <- nfVar(Rmcmc, 'mvSamples')
CmvSample <- nfVar(Cmcmc, 'mvSamples')

#calculate(Rmodel, 'lambda')
Rniter <- 20
set.seed(0);     Rmcmc(Rniter)
Cniter <- 100000
set.seed(0);     Cmcmc(Cniter)

R_samples <- modelValues2Matrix(RmvSample)
C_samples <- modelValues2Matrix(CmvSample)

if(require(testthat)) {
  context(paste0("testing ", example, " MCMC"))
  test_that(paste0("test of equality of output from R and C versions of ", example, " MCMC"), {
    expect_that(R_samples, equals(C_samples[1:Rniter,]), info = paste("R and C posterior samples are not equal"))
#    for(varName in stochVars) 
#      expect_that(unlist(RmvSample[[varName]]), equals(unlist(CmvSample[[varName]][1:Rniter])), info = paste0("R and C values of variable ", varName, " not equal"))
  })
}

burnin <- 10000

tmp <- C_samples[(burnin+1):Cniter, ]
means <- colMeans(tmp)
lower <- apply(tmp, 2, quantile, .025)
upper <- apply(tmp, 2, quantile, .975)
output <- cbind(mns, lower, upper)
dimnames(output)[[2]] <- c('mean', 'lower', 'upper')
cat("Posterior summary: mean and 95% credible interval\n")
print(output)





