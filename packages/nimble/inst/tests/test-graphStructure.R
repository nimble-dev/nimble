## These tests address construction and querying of the graph structure
## test-getDependencies has some more basic tests.
source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

goldFileName <- 'graphStructureTestLog_Correct.Rout'
tempFileName <- 'graphStructureTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForGraphStructureTesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForGraphStructureTesting'), goldFileName) else tempFileName

sink_with_messages(outputFile)

cases <- list()
caseName <- 'graph structure tests case 1'
m <- nimbleModel(
    code = nimbleCode({
        a ~ dnorm(0,1)
        for(i in 1:10) x[i] ~ dnorm(a, 1)
        for(j in 1:10) y[j] ~ dnorm(x[j], 1)
    })
)
cases[[caseName]] <- list(
    m$getDependencies('a'),
    m$getDependencies('a', omit = 'a'),
    m$getDependencies('a', self = FALSE),
    m$getDependencies('a', omit = 'a', self = FALSE),
    m$getDependencies('x[2]'),
    m$getDependencies('x[2]', self = FALSE),
    m$getDependencies('x[2]', omit = 'x[2]'),
    m$getDependencies('x[2]', omit = 'y[2:3]'),
    m$getDependencies('x[2:4]'),
    m$getDependencies('x[2:4]', self = FALSE),
    m$getDependencies('x[2:4]', omit = 'x[4:6]'),
    m$getDependencies('a', downstream = TRUE),
    m$getDependencies('a', downstream = TRUE, self = FALSE),
    m$getDependencies('a', downstream = TRUE, returnScalarComponents = TRUE),
    m$getDependencies('a', downstream = TRUE, self = FALSE, returnScalarComponents = TRUE)
    )

###
caseName <- 'graph structure tests case 2 (dmnorm fully split)'
m <- nimbleModel(
    code = nimbleCode({
        x[1:5] ~ dmnorm(mu[1:5], cov[1:5, 1:5])                        
        for(i in 1:5) {
            ypred[i] <- x[i] + 1
            y[i] ~ dnorm(ypred[i], 1)
        }
    }),
    constants = list(mu = rep(0,5), cov = diag(5))
)

cases[[caseName]] <- list(
    m$getDependencies('x[1:5]'),
    m$getDependencies('x[1:5]', returnScalarComponents = TRUE),
    m$getDependencies('x[1:5]', self = FALSE),
    m$getDependencies('x[1:5]', self = FALSE, returnScalarComponents = TRUE),
    m$getDependencies('x[1:5]',  omit = c('x[2:3]','y[c(1, 5)]')),
    m$getDependencies('x[2]'),
    m$getDependencies('x[2]', returnScalarComponents = TRUE),
    m$getDependencies('x[2]', self = FALSE),
    m$getDependencies('x[2:4]'),
    m$getDependencies('x[2:4]', self = FALSE),
    m$getDependencies('mu[2:4]'),
    m$getDependencies('mu[2:4]', includeRHSonly = TRUE),
    m$getDependencies('mu[2:4]', includeRHSonly = TRUE, returnScalarComponents = TRUE),
    m$getDependencies('mu', includeRHSonly = TRUE)
)

###
caseName <- 'split in the middle of a vector node (original case of Issue #340)'
code <- nimbleCode({
    x[1:4] ~ dmnorm(mu[1:4], C[1:4,1:4])
    y[1] ~ dnorm(x[2], 1)
})
constants <- list(mu = rep(0,4), C = diag(4))
Rmodel <- nimbleModel(code, constants = constants)
cases[[caseName]] <- list(
    Rmodel$getDependencies('x'),
    Rmodel$getDependencies('x[1]'),
    Rmodel$getDependencies('x[2]'),
    Rmodel$getDependencies('x[3]')
)

##
caseName <- 'split vector node, with dependencies through LHSinferred (from Issue #734)'
code <- nimbleCode({
    for(i in 1:3) mu[i] ~ dnorm(0, 1) ## How mu[1:3] is set up is not important to the bug
    x[1:3] ~ dmnorm(mu[1:3], C[1:3,1:3])
    y[1] ~ dnorm(x[2], 1)
    z <- y[1] + 1
})
constants <- list(C = diag(3))
## constants missing in test-graphstructure
Rmodel <- nimbleModel(code, constants = constants)
cases[[caseName]] <- list(
    Rmodel$getDependencies('mu', downstream = TRUE)
)

##
caseName <- 'scalar split in the middle of a matrix node'
code <- nimbleCode({
    x[1:5, 1:5] <- 2 * inputx[1:5, 1:5]
    y[1] ~ dnorm(x[2, 3], 1)
})

constants <- list(inputx = diag(5))
Rmodel <- nimbleModel(code, constants = constants)
cases[[caseName]] <- list(
    Rmodel$getDependencies('x[1:5, 2:3]'),
    Rmodel$getDependencies('x[1:5, 2:3]', self = FALSE),
    Rmodel$getDependencies('x[2, 3]'),
    Rmodel$getDependencies('x[1:5, 1]'),
    Rmodel$getDependencies('x[1:5, 2:3]')
)

caseName <- 'vector splits of matrix node'
code <- nimbleCode({
    x[1:5, 1:5] <- 2 * inputx[1:5, 1:5]
    for(i in 1:5) y[i,1:5] ~ dmnorm(x[i, 1:5], cov[1:5, 1:5])
})

constants <- list(inputx = diag(5), cov = diag(5))
Rmodel <- nimbleModel(code, constants = constants)
cases[[caseName]] <- list(
    Rmodel$getDependencies('x[2, 1:5]'),
    Rmodel$getDependencies('x[2, 4]'),
    Rmodel$getDependencies('x[1:3, 4]'),
    Rmodel$getDependencies('x[1:3, 4]'),
    Rmodel$getDependencies('x[1:3, 3:5]')
)

caseName <- 'some double splitting'
code <- nimbleCode({
    x[1:5, 1:5] <- 2 * inputx[1:5, 1:5]
    a[1:3, 1:4] <- 2 * x[2:4, 2:5]
    for(i in 1:3) y[1:3, i + 1] ~ dmnorm(a[1:3, i], cov[1:3, 1:3])
    z[2] ~ dnorm(a[2, 3], 1)
    z[3] ~ dnorm(a[1, 3], 1)
})

constants <- list(inputx = diag(5), cov = diag(5))
Rmodel <- nimbleModel(code, constants = constants)
cases[[caseName]] <- list(
    Rmodel$getDependencies('x'),
    Rmodel$getDependencies('x[2:4, 2:5]'),
    Rmodel$getDependencies('x[2:4, 3]'),
    Rmodel$getDependencies('x[1, 1:5]'),
    Rmodel$getDependencies('a[1:3, 2]'),
    Rmodel$getDependencies('a[1:3, 3]'),
    Rmodel$getDependencies('a[2:3, 2:3]')
)

caseName <- 'some wierd double splitting'
code <- nimbleCode({
    x[1:5, 1:5] <- 2 * inputx[1:5, 1:5]
    for(i in 1:3) y[i,1:5] ~ dmnorm(x[i, 1:5], cov[1:5, 1:5])
    y[5,3] ~ dnorm(x[2, 3], 1)
    y[6, 1:5] ~ dmnorm(x[5, 1:5], cov[1:5, 1:5])
})

constants <- list(inputx = diag(5), cov = diag(5))
Rmodel <- nimbleModel(code, constants = constants)

cases[[caseName]] <- list(
    Rmodel$getDependencies('x'),
    Rmodel$getDependencies('x[2,3]'),
    Rmodel$getDependencies('x[5,1]'),
    Rmodel$getDependencies('x[2:5,1]'),
    Rmodel$getDependencies('x[2:4,1:5]'),
    Rmodel$getDependencies('x[1:2,1:5]')
)

caseName <- 'elemental tests of makeVertexNamesFromIndexArray2'
indArr <- matrix(1, nrow = 5, ncol = 5)
indArr[1:5, 1:2] <- 4
indArr[2, 1:3] <- 2
indArr[3, 2:5] <- 3

indArr2 <- matrix(1, nrow = 3, ncol = 3)
indArr2[2:3, 1:3] <- 2

cases[[caseName]] <- list(
    nimble:::makeVertexNamesFromIndexArray2(indArr, 1, 'x'),
    nimble:::makeVertexNamesFromIndexArray2(indArr2, 1, 'x')
)

sink(NULL)

if(!generatingGoldFile) {
    trialResults <- readLines(tempFileName)
    correctResults <- readLines(system.file(file.path('tests', goldFileName), package = 'nimble'))
    compareFilesByLine(trialResults, correctResults)
}

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
