source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of compareMCMCs")
## iCase <- 1

## code1 <- nimbleCode({
##     a ~ dnorm(0, 1)
##     b ~ dnorm(a, 1)
##     for(i in 1:3) y[i] ~ dnorm(b, 1)
## })

## input1 <- list(model = code1, data = list(y = 1:3), inits = list(a = 0.5))

## test1 <- try(results1 <- compareMCMCs(input1))
## expect_true(!inherits(test1, 'try-error'), paste('case', iCase))
## iCase <- iCase + 1

## make_MCMC_comparison_pages(results1, 'model1')
