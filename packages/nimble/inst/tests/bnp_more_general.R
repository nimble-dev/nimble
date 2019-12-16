## Checking find cluster nodes for more general cases
## These extend (greatly) from cases from bnp_more_general_cases.Rmd (Sept. 2019)

library(nimble, lib.loc = '/tmp/nim-bmg')

## Basic more general cases

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

### check wrapper sampler identifies clusterID correctly 
samplers <- conf$getSamplers()
wh <- which(sapply(samplers, function(x) x$name == "CRP_cluster_wrapper"))
ids <- sapply(wh, function(i) samplers[[i]]$control$clusterID)
all(ids == rep(1:5, each = 3))

## truncated: findClusterNodes gives same result
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
        }}
    for(i in 1:n2) {
        for(j in 1:J) {
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
n2 <- 4
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n2), n2, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

## intermediate nodes
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(theta[i,j], 1)
            theta[i, j] <- thetaTilde[xi[i], j]
            thetaTilde[i,j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

## truncated and intermediate nodes

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(theta[i, j], 1)
            theta[i, j] <- thetaTilde[xi[i], j]
        }}
    for(i in 1:n2)
        for(j in 1:J)
            thetaTilde[i, j] ~ dnorm(0, 1)
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
n2 <- 4
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

## column in instead of row
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[j, i] ~ dnorm(thetaTilde[j, xi[i]], 1)
            thetaTilde[j, i] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),J,n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), J, n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

## column in instead of row, different loop ordering
## FIX
code <- nimbleCode({
    for(j in 1:J) {
        for(i in 1:n) {
            y[j, i] ~ dnorm(thetaTilde[j, xi[i]], 1)
            thetaTilde[j, i] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),J,n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), J, n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()        # INCORRECT

## column instead of row, intermediate nodes, indexes swapped
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[j, i] ~ dnorm(theta[j, i], 1)
            theta[j, i] <- thetaTilde[xi[i], j]
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),J,n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

## index offset
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i]+1, j], 1)
        }}
    for(i in 2:(n+1)) {
        for(j in 1:J) {
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*(n+1)), n+1, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()

## index offset, with intermediate
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(theta[i, j], 1)
            theta[i,j] <- thetaTilde[xi[i]+1, j]
        }}
    for(i in 2:(n+1)) {
        for(j in 1:J) {
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*(n+1)), n+1, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

### check wrapper sampler identifies clusterID correctly 
samplers <- conf$getSamplers()
wh <- which(sapply(samplers, function(x) x$name == "CRP_cluster_wrapper"))
ids <- sapply(wh, function(i) samplers[[i]]$control$clusterID)
all(ids == rep(1:5, each = 3))

## detection of differing number of cluster parameters
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(thetaTilde[xi[i]], s2tilde[xi[i]])
        s2tilde[i] ~ dunif(0, 10)
    }
    for(i in 1:n2) 
        thetaTilde[i] ~ dnorm(0, 1)
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
n2 <- 4
J <- 3
constants <- list(n = n, J = J)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n2), s2tilde = runif(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # correctly errors out



## detection of differing number of cluster parameters, complicated indexing
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(thetaTilde[xi[i]], s2tilde[xi[i]+1])
        s2tilde[i] ~ dunif(0, 10)
        thetaTilde[i] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), s2tilde = runif(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # correctly errors out

## multiple obs, but only single thetaTilde in each group
## check with Claudia it is ok to do non-conjugate here
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i]], 1)
        }
        thetaTilde[i] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   # not conjugate
conf <- configureMCMC(model)
conf$printSamplers() 
mcmc <- buildMCMC(conf)


## Model A: dnorm obs, dmnorm prior
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J)
            y[i,j] ~ dnorm(thetaTilde[xi[i],j], 1)
        thetaTilde[i, 1:J] ~ dmnorm(mn[1:J], iden[1:J,1:J])
    }
    mn[1:J] <- mu*ones[1:J]
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J),
              iden = diag(3), mu = 0)
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') #  no conjugacy because dnorm vs. dmnorm
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Model B: separate prior specifications
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:2) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
        }
        thetaTilde[i, 1] ~ dnorm(0, 1)
        thetaTilde[i, 2] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 2
constants <- list(n = n, J=J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)


## Model C: dmnorm obs, dmnorm prior
code <- nimbleCode({
    for(i in 1:n) {
        y[i, 1:J] ~ dmnorm(thetaTilde[xi[i],1:J], iden[1:J,1:J])
        thetaTilde[i, 1:J] ~ dmnorm(mn[1:J], iden[1:J,1:J])
    }
    mn[1:J] <- mu*ones[1:J]
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J),
              iden = diag(3), mu = 0)
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # detects conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Model D: dmnorm obs, dnorm prior
code <- nimbleCode({
    for(i in 1:n) {
        y[i, 1:J] ~ dmnorm(thetaTilde[xi[i],1:J], iden[1:J,1:J])
        for(j in 1:J)
            thetaTilde[i, j] ~ dnorm(0,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J),
              iden = diag(J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # no conj because of dnorm vs. dmnorm
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Model E: dmnorm obs, separate prior declarations
## FIX this - the wrapped sampler doesn't set up clusterID correctly
code <- nimbleCode({
    for(i in 1:n) {
        y[i, 1:2] ~ dmnorm(thetaTilde[xi[i],1:2], iden[1:2,1:2])
        thetaTilde[i, 1] ~ dnorm(0,1)
        thetaTilde[i, 2] ~ dnorm(0,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 2
constants <- list(n = n)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J),
              iden = diag(J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # no conj because of dnorm vs. dmnorm
conf <- configureMCMC(model)
conf$printSamplers()         # INCORRECT
mcmc <- buildMCMC(conf)


## Model F: dnorm obs, separate priors
code <- nimbleCode({
    for(i in 1:n) {
        y1[i] ~ dnorm(thetaTilde1[xi[i]], 1)
        y2[i] ~ dnorm(thetaTilde2[xi[i]], 1)
        thetaTilde1[i] ~ dnorm(mu, sigma)
        thetaTilde2[i] ~ dnorm(mu, sigma)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y1 = rnorm(n), y2 = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde1 = rnorm(n), thetaTilde2 = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not detected - because two clusterVars
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Model F: non conjugate, non-identical priors
## FIX 
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
    }
    thetaTilde[1] ~ dnorm(mu, sigma)
    for(i in 2:n) 
        thetaTilde[i] ~ dgamma(1,1)
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # non-conjugate
conf <- configureMCMC(model)
conf$printSamplers()   # INCORRECT
mcmc <- buildMCMC(conf)

## Model F: non conjugate, non-identical priors, part 2
## doesn't detect presence of some conjugacies because only check first set
## weird setup causes us to say we can't handle this case
code <- nimbleCode({
    for(i in 3:n) {
        y[i] ~ dt(thetaTilde[xi[i]], 1, 1)
    }
    for(i in 1:2)
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
    for(i in 1:n) 
        thetaTilde[i] ~ dnorm(0,1)
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()     # doesn't detect some thetaTilde samplers can be conjugate
## weird model setup causes us to say we can't handle this case
mcmc <- buildMCMC(conf)  # error with number of cluster parameters not the same

## Model F: non conjugate, non-identical priors, part 3
code <- nimbleCode({
    for(i in 2:n) {
        y[i] ~ dt(thetaTilde[xi[i]], 1, 1)
    }
    y[1] ~ dnorm(thetaTilde[xi[1]], 1)
    for(i in 1:n) 
        thetaTilde[i] ~ dnorm(0,1)
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') 
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # errors based on unusual indexing

## another case of mixed distributions
code <- nimbleCode({
    for(i in 1:(n-1)) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    for(j in 1:J) {
        y[n, j] ~ dt(thetaTilde[xi[n], j], 1 ,1)
        thetaTilde[n, j] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)  # errors based on unusual indexing


## Model G: separate declarations for obs
code <- nimbleCode({
    for(i in 1:n) {
        y1[i] ~ dnorm(thetaTilde[xi[i], 1], 1)
        y2[i] ~ dnorm(thetaTilde[xi[i], 2], 1)
        for(j in 1:2)
            thetaTilde[i, j] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 2
constants <- list(n = n)
data <- list(y1 = rnorm(n), y2 = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(n*J),n,J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not detected - because two clusterVars
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Model H: separate obs and prior declarations
code <- nimbleCode({
    for(i in 1:n) {
        y1[i] ~ dnorm(thetaTilde[xi[i], 1], 1)
        y2[i] ~ dnorm(thetaTilde[xi[i], 2], 1)
        thetaTilde[i, 1] ~ dnorm(0, 1)
        thetaTilde[i, 2] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 2
constants <- list(n = n)
data <- list(y1 = rnorm(n), y2 = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(n*J),n,J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not detected - two clusterVars
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Model I: separate declartions for obs, dmnorm for prior
code <- nimbleCode({
    for(i in 1:n) {
        y1[i] ~ dnorm(thetaTilde[xi[i], 1], 1)
        y2[i] ~ dnorm(thetaTilde[xi[i], 2], 1)
        thetaTilde[i, 1:2] ~ dmnorm(mn[1:2], iden[1:2,1:2])
    }
    mn[1:2] <- mu*ones[1:2]
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 2
constants <- list(n = n)
data <- list(y1 = rnorm(n), y2 = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(n*J), n, J), ones = rep(1,2))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not conj - dmnorm and dnorm
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## model J: mixed indexing of two cluster vars
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i,j] ~ dnorm(theta[i], sd = sigma[i,j])        
            sigma[i,j] <- sigmaTilde[xi[i], j]
            sigmaTilde[i,j] ~ dinvgamma(a, b)
        }
        theta[i] <- thetaTilde[xi[i]]
        thetaTilde[i] ~ dnorm(mu, phi)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 2
constants <- list(n = n)
data <- list(y = matrix(rnorm(n*J), n, J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), sigmaTilde = matrix(rgamma(n*J, 1, 1), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)


## Model K: more complicated intermediate nodes
## We do not have normal-gamma conj set up; this is ok.
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(theta[i], sd = sigma[i])
        theta[i] <- thetaTilde[xi[i]]
        sigma[i] <- 1 / tau[i]
        tau[i] <- tauTilde[xi[i]]
        thetaTilde[i] ~ dnorm(mu, var = sigmaTilde[i])
        sigmaTilde[i] <- 1 / tauTilde[i]
        tauTilde[i] ~ dgamma(a, b)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})
n <- 5
J <- 2
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), tauTilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## another variation where we have intermediate nodes and dependency in cluster nodes
## clustering of deterministic nodes - not allowed
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(theta[i], sd = sigma[i])
        theta[i] <- thetaTilde[xi[i]]
        sigma[i] <- sigmaTilde[xi[i]]
        sigmaTilde[i] <- 1 / tauTilde[i]
        thetaTilde[i] ~ dnorm(mu, var = sigmaTilde[i])
        tauTilde[i] ~ dgamma(a, b)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})
n <- 5
J <- 2
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), tauTilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

## Chris needs to do more checking of the cases below this line.

## 3-d cases
code <- nimbleCode({
  for(i in 1:3) {
    for(j in 1:4) {
        y[i,j] ~ dnorm(thetaTilde[xi[i], 2, j] , var = 1)
    }
  }
  for(i in 1:3) {
      for(j in 1:4) {
          for(k in 1:2) {
              thetaTilde[i, k, j] ~ dnorm(0,1)
          }
      }
  }
  xi[1:3] ~ dCRP(1, size=3)
})

Inits <- list(xi = c(1, 1, 1), 
              thetaTilde = array(0, c(3,2,4)))

y <- matrix(5, nrow=3, ncol=4) 
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:3]')

## 3-d: this is ok but arbitrary ordering by k and j
code <- nimbleCode({
    for(i in 1:3) {
        for(k in 1:2) {
            for(j in 1:4) {
                y[i,k,j] ~ dnorm(thetaTilde[xi[i], k, j] , var = 1)
            }
        }
    }
    for(i in 1:3) {
        for(j in 1:4) {
            for(k in 1:2) {
                thetaTilde[i, k, j] ~ dnorm(0,1)
            }
        }
    }
    xi[1:3] ~ dCRP(1, size=3)
})

Inits <- list(xi = c(1, 1, 1), 
              thetaTilde = array(0, c(3,2,4)))

y <- array(5, c(3,2,4))
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:3]')

## more 3-d
code <- nimbleCode({
    for(i in 1:3) {
        for(k in 1:2) {
            for(j in 1:4) {
                y[k,i,j] ~ dnorm(thetaTilde[k, xi[i], j] , var = 1)
            }
        }
    }
    for(i in 1:3) {
        for(j in 1:4) {
            for(k in 1:2) {
                thetaTilde[k, i, j] ~ dnorm(0,1)
            }
        }
    }
    xi[1:3] ~ dCRP(1, size=3)
})

Inits <- list(xi = c(1, 1, 1), 
              thetaTilde = array(0, c(2,3,4)))

y <- array(5, c(2,3,4))
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:3]')


## weird indexing

code <- nimbleCode({
  for(i in 1:3) {
    for(j in 1:4) {
      y[i,j] ~ dnorm(thetaTilde[xi[3-i+1], j] , var = 1) 
      thetaTilde[i, j] ~ dnorm(0,1)
    }
  }
  xi[1:3] ~ dCRP(1, size=3)
})

Inits <- list(xi = c(1, 1, 1), 
              thetaTilde = matrix(0, nrow=3,  ncol=4))

y <- matrix(5, nrow=3, ncol=4) #rnorm(5*8, 5, 1)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:3]')

## weird indexing - findClusterNodes finds that 'i+j' is a function and errors
code <- nimbleCode({
  for(i in 1:3) {
    for(j in 1:4) {
      y[i,j] ~ dnorm(thetaTilde[xi[i+j], j] , var = 1) 
      thetaTilde[i, j] ~ dnorm(0,1)
    }
  }
  xi[1:7] ~ dCRP(1, size=7)
})

Inits <- list(xi = rep(1,7), 
              thetaTilde = matrix(0, nrow=3,  ncol=4))

y <- matrix(5, nrow=3, ncol=4) #rnorm(5*8, 5, 1)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:7]')

## Cases of multiple use of CRP node 

## use of xi[i] and xi[j] in same statement - errors out (correctly)

code <- nimbleCode({
  for(i in 1:3) {
    for(j in 1:4) {
      y[i,j] ~ dnorm(thetaTilde[xi[i], xi[j]] , var = 1) 
      thetaTilde[i, j] ~ dnorm(0,1)
    }
  }
  xi[1:4] ~ dCRP(1, size=4)
})
Inits <- list(xi = rep(1,7), 
              thetaTilde = matrix(0, nrow=3,  ncol=4))

y <- matrix(5, nrow=3, ncol=4) #rnorm(5*8, 5, 1)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:4]')

## use of xi[i] and xi[i] in same parameter - errors out (correctly)

code <- nimbleCode({
  for(i in 1:4) {
    for(j in 1:4) {
      y[i,j] ~ dnorm(thetaTilde[xi[i], xi[i]] , var = 1) 
      thetaTilde[i, j] ~ dnorm(0,1)
    }
  }
  xi[1:4] ~ dCRP(1, size=4)
})
Inits <- list(xi = rep(1,7), 
              thetaTilde = matrix(0, nrow=4,  ncol=4))

y <- matrix(5, nrow=4, ncol=4)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:4]')

## xi[i] and xi[j] in different parameters

code <- nimbleCode({
  for(i in 1:4) {
    for(j in 1:3) {
      y[i,j] ~ dnorm(thetaTilde[xi[i]], var = s2Tilde[xi[j]])
    }
  }
  for(i in 1:4)
      thetaTilde[i] ~ dnorm(0,1)
  for(i in 1:3)
      s2Tilde[i] ~ dnorm(0,1)
  xi[1:4] ~ dCRP(1, size=4)
})
Inits <- list(xi = rep(1,4))
y <- matrix(5, nrow=4, ncol=3)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:4]')
nimble:::checkCRPconjugacy(model, 'xi[1:4]')
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf) # errors as there are 3 thetatildes and 4 s2tildes

## variation on the above
code <- nimbleCode({
  for(i in 1:4) {
    for(j in 1:4) {
      y[i,j] ~ dnorm(thetaTilde[xi[i]], var = s2Tilde[xi[j]])
    }
  }
  for(i in 1:4)
      thetaTilde[i] ~ dnorm(0,1)
  for(i in 1:4)
      s2Tilde[i] ~ dnorm(0,1)
  xi[1:4] ~ dCRP(1, size=4)
})
Inits <- list(xi = rep(1,4))
y <- matrix(5, nrow=4, ncol=4)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:4]')  # errors
nimble:::checkCRPconjugacy(model, 'xi[1:4]')


## xi[i] used twice in same parameter

code <- nimbleCode({
    for(i in 1:4) {
        y[i] ~ dnorm(b0[xi[i]] + b1[xi[i]]*x[i], var = 1)
    }
  for(i in 1:4) {
      b0[i] ~ dnorm(0,1)
      b1[i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(x = rnorm(4), xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # conjugacy not detected
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)

## xi[i] used twice in same parameter, variation

code <- nimbleCode({
    for(i in 1:4) {
        y[i] ~ dnorm(beta[xi[i], 1] + beta[xi[i], 2]*x[i], var = 1)
    }
  for(i in 1:4) {
      beta[i,1] ~ dnorm(0,1)
      beta[i,2] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(x = rnorm(4), xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # conjugacy not detected
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)

## xi[i] used twice in same parameter, another variation

code <- nimbleCode({
    for(i in 1:4) {
        y[i] ~ dnorm(beta[xi[i], 1] + beta[xi[i], 2]*x[i], var = 1)
    }
    for(i in 1:4) 
        for(j in 1:2) 
      beta[i,j] ~ dnorm(0,1)
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(x = rnorm(4), xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # conjugacy not detected
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)

## use of inprod
code <- nimbleCode({
    for(i in 1:4) {
        y[i] ~ dnorm(inprod(beta[1:2, xi[i]], x[i]), var = 1)
    }
  for(i in 1:4) {
      beta[1, i] ~ dnorm(0,1)
      beta[2, i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(x = rnorm(4), xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # non-identity conjugacy detected

## Cases of trickiness in structure of observations

code <- nimbleCode({
    for(i in 1:4) {
        for(j in 2:3) {
            y[i,j] ~ dnorm(mu[xi[i]] + y[i,j-1], 1)
        }
    }
  for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*3), 4 ,3))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # non-identity conjugacy detected, but moreGeneral version not yet written
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # errors with lack of independence

## another non-indep of y's case
code <- nimbleCode({
    for(i in 1:4) {
        for(j in 2:3) {
            y[i,j] ~ dnorm(mu[xi[i]], exp(y[i,j-1]))
        }
    }
  for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*3), 4 ,3))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # non-identity conjugacy detected, but moreGeneral version not yet written
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # errors with lack of independence

## obs dependent across group

code <- nimbleCode({
    for(i in 2:4) {
        for(j in 1:3) {
            y[i,j] ~ dnorm(mu[xi[i]] + y[i-1,j], 1)
        }
    }
  for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*3), 4 ,3))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # no conj as it checks y[1,1] first (I think)
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # errors with lack of independence but error indicates within group issue

code <- nimbleCode({
    for(i in 1:3) {
        for(j in 1:3) {
            y[i,j] ~ dnorm(mu[xi[i]] + y[i+1,j], 1)
        }
    }
  for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*3), 4 ,3))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # no conj because of restrictLink
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # errors with lack of independence but error indicates within group issue

## Cases of trickiness in structure of observations, but with multiple y declarations

code <- nimbleCode({
    for(i in 1:4) {
        y[i,1] ~ dnorm(mu[xi[i]], 1)
        for(j in 2:3) {
            y[i,j] ~ dnorm(mu[xi[i]] + y[i,j-1], 1)
        }
    }
  for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*3), 4 ,3))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # no conjugacy - presumably because multiple decl so multiple clusterVars
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # errors with lack of independence

## y dependence in multivariate way
code <- nimbleCode({
    for(i in 1:4) {
        y[i,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[xi[i],1:2,1:2])
        mu[i, 1:2] ~ dmnorm(mu0[1:2], iden[1:2,1:2])
        sigma[i,1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*2), 4 ,2))
inits = list(xi = rep(1,4), iden = diag(2))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # not conjugate
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  

## should be conj once we have dmnorm-invwish-dmnorm
code <- nimbleCode({
    for(i in 1:4) {
        y[i,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[xi[i],1:2,1:2])
        mu[i, 1:2] ~ dmnorm(mu0[1:2], sigmaAux[i,1:2,1:2])
        sigmaAux[i,1:2,1:2] <- sigma[i,1:2,1:2]/kappa
        sigma[i,1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*2), 4 ,2))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # not conjugate
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  

## y dependence in multivariate way
## check back on this when bnp_dmnorm_dmnorm is merged in
code <- nimbleCode({
    for(i in 1:4) {
        y[i,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[1:2,1:2])
        mu[i, 1:2] ~ dmnorm(mu0[1:2], iden[1:2,1:2])
  }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y =matrix( rnorm(4*2), 4 ,2))
inits = list(xi = rep(1,4), iden = diag(2), sigma = diag(2))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]') # we're detecting conjugacy
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # we print it is conjugate, but we assign nonconj
mcmc$samplerFunctions[[1]]$helperFunctions

## clusters not independent

code <- nimbleCode({
    for(i in 1:4) 
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
    thetaTilde[1] ~ dnorm(a,1)
    for(i in 2:4)
        thetaTilde[i] ~ dnorm(thetaTilde[i-1], 1)
    xi[1:4] ~ dCRP(1, size=4)
    a ~ dunif(0,1)
})
data = list(y = rnorm(4))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # not conj
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)  # error - clusters not indep

## clusters not indep - alternative case
code <- nimbleCode({
    for(i in 1:4) 
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
    thetaTilde[4] ~ dnorm(0,1)
    for(i in 1:3)
        thetaTilde[i] ~ dnorm(thetaTilde[i+1], 1)
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # not conj
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf) # error - clusters not indep


## clusters not indep, with mv decl

code <- nimbleCode({
    for(i in 1:4) 
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
    thetaTilde[1:4] ~ dmnorm(z[1:4], iden[1:4,1:4]) # formally indep but not known to us
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(xi = rep(1,4), iden = diag(4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # error - not indep
conf <- configureMCMC(model) # error - not indep


## clusters indep G0 but not IID

code <- nimbleCode({
    for(i in 1:4)  {
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
        thetaTilde[i] ~ dnorm(i, 1)
    }
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # not conj
conf <- configureMCMC(model)  # non-conjugate

code <- nimbleCode({
    for(i in 1:4) 
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
    for(i in 1:3) 
    thetaTilde[i] ~ dnorm(0, 1)
    thetaTilde[4] ~ dnorm(5, 2)
  xi[1:4] ~ dCRP(1, size=4)
})
data = list(y = rnorm(4))
inits = list(xi = rep(1,4))
model <- nimbleModel(code, data = data, inits = inits)
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # not conj
conf <- configureMCMC(model)  # non-conjugate
mcmc <- buildMCMC(conf)

## clusters IID G0; IID within cluster

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)


## clusters IID G0; indep within cluster

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
            thetaTilde[i, j] ~ dnorm(j, 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # detects conj
conf <- configureMCMC(model)
conf$printSamplers()  # uses conjugate sampler
mcmc <- buildMCMC(conf)

## indep and diff dists 

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
        }
        thetaTilde[i, 1] ~ dnorm(0, 1)
        thetaTilde[i, 2] ~ dgamma(3, 1)
        thetaTilde[i, 3] ~ dnorm(5, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') 
conf <- configureMCMC(model)
conf$printSamplers()
mcmc <- buildMCMC(conf)

code <- nimbleCode({
    for(i in 1:3) {
            y[i] ~ dnorm(thetaTilde[xi[i]], 1)
        }
        thetaTilde[1] ~ dnorm(0, 1)
        thetaTilde[2] ~ dgamma(1, 1)
        thetaTilde[3] ~ dnorm(5, 1)
    xi[1:3] ~ dCRP(alpha, size = 3)
})

n <- 3
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:3]')
nimble:::checkCRPconjugacy(model, 'xi[1:3]')  # correctly says no conjugacy

code <- nimbleCode({
    for(i in 1:3) {
        y[i] ~ dnorm(theta[i], 1)
        theta[i] <- thetaTilde[xi[i]]
    }
    thetaTilde[1] ~ dnorm(0, 1)
    thetaTilde[2] ~ dgamma(1, 1)
    thetaTilde[3] ~ dnorm(5, 1)
    xi[1:3] ~ dCRP(alpha, size = 3)
})

n <- 3
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:3]')
nimble:::checkCRPconjugacy(model, 'xi[1:3]')  # correctly says no conjugacy

code <- nimbleCode({
    for(i in 1:3) {
        y[i] ~ dnorm(theta[i], 1)
        thetaTilde[i] ~ dnorm(0,1)
    }
    theta[1] <- thetaTilde[xi[1]]
    theta[2] <- exp(thetaTilde[xi[2]])
    theta[3] <- thetaTilde[xi[3]]
    xi[1:3] ~ dCRP(alpha, size = 3)
})

n <- 3
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:3]')
nimble:::checkCRPconjugacy(model, 'xi[1:3]')  # correctly says no conjugacy (because findClusterNodes gives character(0) 

code <- nimbleCode({
    for(i in 1:4) {
        y[i] ~ dnorm(theta[i], 1)
        thetaTilde[i] ~ dnorm(0,1)
    }
    for(i in 1:2)
        theta[i] <- thetaTilde[xi[i]]
    for(i in 3:4)
        theta[i] <- exp(thetaTilde[xi[i]])
    xi[1:4] ~ dCRP(alpha, size = 4)
})

n <- 4
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:4]')
nimble:::checkCRPconjugacy(model, 'xi[1:4]')  # correctly says no conjugacy (because findClusterNodes gives character(0) 


## mutilde and s2tilde; not conjugate

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i]])
            thetaTilde[i, j] ~ dnorm(j, 1)
        }
        s2tilde[i] ~ dinvgamma(1,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()   # non-conj
mcmc <- buildMCMC(conf)

## mutilde and s2tilde; not conjugate

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i]], var = s2tilde[xi[i]])
        }
        thetaTilde[i] ~ dnorm(0, 1)
        s2tilde[i] ~ dinvgamma(1,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), s2tilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()   # non-conj
mcmc <- buildMCMC(conf)

## another case
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i]], sd = s2tilde[xi[i]])
        }
        thetaTilde[i] ~ dnorm(0, 1)
        s2tilde[i] ~ dinvgamma(1,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), s2tilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()   # non-conj
mcmc <- buildMCMC(conf)

## mutilde and s2tilde, IID across cluster, indep within

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i], j])
            thetaTilde[i, j] ~ dnorm(j, var = s2tilde[i, j]/kappa)
            s2tilde[i, j] ~ dinvgamma(1,1)
        }

    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = matrix(rgamma(J*n, 1, 1),n,J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   
conf <- configureMCMC(model)  # conj
conf$printSamplers()  
mcmc <- buildMCMC(conf) 



## mutilde and s2tilde; conj
                                      
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i], j])
            thetaTilde[i, j] ~ dnorm(theta0, var = s2tilde[i, j]/kappa)
            s2tilde[i, j] ~ dinvgamma(1,1)
        }

    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = matrix(rgamma(J*n, 1, 1),n,J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   #  conjugate
conf <- configureMCMC(model)
conf$printSamplers()   # conj
mcmc <- buildMCMC(conf)  

## mutilde and s2tilde; not conj because multiple obs per cluster param set
                                      
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i]], var = s2tilde[xi[i]])
        }
        thetaTilde[i] ~ dnorm(theta0, var = s2tilde[i]/kappa)
        s2tilde[i] ~ dinvgamma(1,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), s2tilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   # not conjugate
conf <- configureMCMC(model)  
conf$printSamplers()   
mcmc <- buildMCMC(conf)

## weird case of multiple thetaTilde but not s2tilde

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i]])
            thetaTilde[i,j] ~ dnorm(theta0, var = s2tilde[i]/kappa)
        }
        s2tilde[i] ~ dinvgamma(1,1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = rgamma(n, 1, 1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')   # not conjugate because thetaTilde and s2tilde don't match 1:1
conf <- configureMCMC(model)
conf$printSamplers()  
mcmc <- buildMCMC(conf)

## clusters IID G0; dep within cluster

code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
        }
        thetaTilde[i, 1] ~ dnorm(0,1)
        thetaTilde[i, 2] ~ dnorm(thetaTilde[i,1], 1)
        thetaTilde[i, 3] ~ dnorm(thetaTilde[i,2], 1)
        }}
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()  # uses conjugate sampler
mcmc <- buildMCMC(conf)

## function in index
code <- nimbleCode({
    for(i in 1:n) {
             y[i] ~ dnorm(thetaTilde[xi[n-i+1]], 1)
        }
    for(i in 1:n) {
            thetaTilde[i] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')  # errors out

code <- nimbleCode({
    for(i in 1:n) {
             y[i] ~ dnorm(thetaTilde[xi[exp(i)]], 1)
        }
    for(i in 1:n) {
            thetaTilde[i] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')  # errors out

code <- nimbleCode({
    for(i in 1:n) {
             y[i] ~ dnorm(thetaTilde[xi[i]], s2tilde[xi[n-i+1]])
        }
    for(i in 1:n) {
            thetaTilde[i] ~ dnorm(0, 1)
    }
    for(i in 1:n)
        s2tilde[i] ~ dunif(0, 10)
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n), s2tilde = runif(n+1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')  # errors out

## use of tilde var in two separate parameters

code <- nimbleCode({
    for(i in 1:n) {
             y[i] ~ dnorm(thetaTilde[xi[i]], exp(thetaTilde[xi[i]]))
        }
    for(i in 1:n) {
            thetaTilde[i] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]') 
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()  
mcmc <- buildMCMC(conf)

code <- nimbleCode({
    for(i in 1:n) {
             y[i] ~ dnorm(theta[i], exp(thetaTilde[xi[i]]))
        }
    for(i in 1:n) {
        thetaTilde[i] ~ dnorm(0, 1)
        theta[i] <- thetaTilde[xi[i]]
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]') 
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()  
mcmc <- buildMCMC(conf)

## case where inconsistent indexing would mess up cluster creation
## and nosample-empty-clusters
code <- nimbleCode({
    for(i in 1:n) {
             y[i] ~ dnorm(thetaTilde[xi[i]+1], exp(thetaTilde[xi[i]]))
        }
    for(i in 1:(n+1)) {
        thetaTilde[i] ~ dnorm(0, 1)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n+1))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')  
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conjugate
conf <- configureMCMC(model)
conf$printSamplers()  # incorrect wrapping because of xi[i]+1
mcmc <- buildMCMC(conf) # correct error trapping 

code <- nimbleCode({
    for(i in 1:n) {
        y1[i] ~ dnorm(thetaTilde[xi[i]], 1)
        y2[i] ~ dnorm(thetaTilde[xi[i]], 1)
        thetaTilde[i] ~ dnorm(mu, sigma)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y1 = rnorm(n), y2 = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::checkCRPconjugacy(model, 'xi[1:5]') # not detected - because two clusterVars
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)

## thetaTildes indexed by multiple cluster variables
code <- nimbleCode({
    for(i in 1:n) {
        y[i] ~ dnorm(thetaTilde[xi[i]], 1)
        z[i] ~ dnorm(thetaTilde[eta[i]], 1)
        thetaTilde[i] ~ dnorm(0, 1)
    }                  
    xi[1:n] ~ dCRP(alpha, size = n)
    eta[1:n] ~ dCRP(alpha, size = n)
})

n <- 5
constants <- list(n = n)
data <- list(y = rnorm(n), z = rnorm(n))
inits <- list(alpha = 1, xi = rep(1, n), eta = rep(1,n),
              thetaTilde = rnorm(n))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')
nimble:::findClusterNodes(model, 'eta[1:5]')  # not conj
nimble:::checkCRPconjugacy(model, 'xi[1:5]')  # not conj
conf <- configureMCMC(model)
conf$printSamplers()  ## non-conjugate for xi and eta
mcmc <- buildMCMC(conf)  ## correctly errors out because tildeNodes have deps other than dataNodes

## what happens here - warning
code <- nimbleCode({
    for(i in 1:n) {
        for(j in 1:J) {
            y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
        }}
    for(i in 1:n) {
        for(j in 1:J) {
            thetaTilde[i, j] ~ dnorm(0, 1)
        }}
    xi[1:(n+1)] ~ dCRP(alpha, size = n+1)
})

n <- 5
n2 <- 4
J <- 3
constants <- list(n = n, J = J)
data <- list(y = matrix(rnorm(n*J),n,J))
inits <- list(alpha = 1, xi = rep(1, n+1),
              thetaTilde = matrix(rnorm(J*n), n, J))
model <- nimbleModel(code, data = data, constants = constants, inits = inits)
nimble:::findClusterNodes(model, 'xi[1:5]')


## Crossed clustering (model 3)

code <- nimbleCode({
  for(i in 1:3) {
    for(j in 1:4) {
      y[i,j] ~ dnorm( thetaTilde[xi[i], eta[j]] , var = 1) 
       thetaTilde[i, j] ~ dnorm(0, 1)
    }
  }
  xi[1:3] ~ dCRP(1, size=3)
  eta[1:4] ~ dCRP(1, size=4)
})
Inits <- list(xi = rep(1,3), eta = rep(1,4), 
              thetaTilde = matrix(0, nrow=3,  ncol=4))

y <- matrix(5, nrow=3, ncol=4)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:3]')
nimble:::findClusterNodes(model, 'eta[1:4]')
nimble:::checkCRPconjugacy(model, 'xi[1:3]')
conf <- configureMCMC(model)
mcmc <- buildMCMC(conf)

## crossed clustering with truncation 

code <- nimbleCode({
  for(i in 1:3) {
    for(j in 1:4) {
        y[i,j] ~ dnorm( thetaTilde[xi[i], eta[j]] , var = 1)
    }}
  for(i in 1:2) {
      for(j in 1:3) {
       thetaTilde[i, j] ~ dnorm(0, 1)
    }
  }
  xi[1:3] ~ dCRP(1, size=3)
  eta[1:4] ~ dCRP(1, size=4)
})
Inits <- list(xi = rep(1,3), eta = rep(1,4), 
              thetaTilde = matrix(0, nrow=3,  ncol=4))

y <- matrix(5, nrow=3, ncol=4)
Data <- list(y=y)
model <- nimbleModel(code, data=Data, inits=Inits)
nimble:::findClusterNodes(model, 'xi[1:3]')
nimble:::findClusterNodes(model, 'eta[1:4]')





