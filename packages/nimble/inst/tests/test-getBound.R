source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of getBound")

code <- nimbleCode({
    y1 ~ T(dnorm(mu, sd = sig),b,d)  # non-constant trunc bounds
    y2 ~ T(dnorm(0,1), y2l, y2u)       # constant trunc bounds
    y3 ~ T(dnorm(0,1), exp(y2l), y3u)       # constant trunc bounds
    y2l ~ dnorm(0,1)
    y3u <- 15*cc
    mu ~ dnorm(0, 1)                 
    sig ~ dunif(0, 5)
    b ~ dunif(a,d)                   # non-constant unif bounds
    a ~ dunif(al,4)                  # constant unif bounds
    z[1:2] ~ dmnorm(mu0[1:2], prec0[1:2,1:2])           # mv dist
    alpha[1:3] ~ ddirch(alpha0[1:3]) # mv dist with actual bound
    for(i in 1:3)
        alpha0[i] ~ dgamma(1,1)
    y4 ~ T(dgamma(1, 1), y4l, y4u)     # invalid trunc bound
    y5 ~ T(dgamma(1, 1), mu, y5u)      # stochastically-invalid trunc bound
    y6 ~ T(dgamma(1, 1), y6l, y6u)     # reversed trunc bound
})

consts <- list(y2u = 0.8, mu0 = rep(0, 2), prec0 = diag(rep(1,2)),
               y4l = -2, y4u = 5, y5u = 0.5, y6l = 2, y6u = 1,
               al = 0.1, dl = 0.15)
inits <- list(b = 0.7, mu = 1, sig = 1, d = 3,
             z = rep(1, 2), alpha0 = 1:3, cc = 1, y2l = 0.3, y3u = 15, cc = 2)
m <- nimbleModel(code, inits = inits,
                constants = consts,
                data = list(y1 = 0, y2 = 0, y3 = 1, y4 = 1, y5 = 1, y6 = 1))
cm <- compileNimble(m)

test <- nimbleFunction(
    setup = function(model, node, bnd) {},
    run = function() {
        output <- numeric(2)
        output[1] <- model$getBound(node, bnd)
        output[2] <- getBound(model, node, bnd)
        return(output)
        returnType(double(1))
    }
)

test_getBound(m, cm, test, 'y1', 'lower', inits$b, "non-constant lower truncation bound")
test_getBound(m, cm, test, 'y1', 'upper', inits$d, "non-constant upper truncation bound")
test_getBound(m, cm, test, 'y2', 'lower', inits$y2l, "non-constant lower truncation bound")
test_getBound(m, cm, test, 'y2', 'upper', consts$y2u, "non-constant upper truncation bound")
test_getBound(m, cm, test, 'y3', 'lower', exp(m$y2l)[1], "functional lower truncation bound")
test_getBound(m, cm, test, 'y3', 'upper', m$y3u[1], "functional/deterministic upper truncation bound")
test_getBound(m, cm, test, 'b', 'lower', m$a[1], "dunif non-constant lower truncation bound")
test_getBound(m, cm, test, 'b', 'upper', m$d[1], "dunif non-constant upper truncation bound")
test_getBound(m, cm, test, 'a', 'lower', consts$al, "dunif constant lower truncation bound")
test_getBound(m, cm, test, 'a', 'upper', 4, "dunif constant upper truncation bound")
test_getBound(m, cm, test, 'z[1:2]', 'lower', -Inf, "dmnorm constant lower truncation bound")
test_getBound(m, cm, test, 'z[1:2]', 'upper', Inf, "dmnorm constant upper truncation bound")
test_getBound(m, cm, test, 'alpha[1:3]', 'lower', 0, "ddirch constant lower truncation bound")
test_getBound(m, cm, test, 'alpha[1:3]', 'upper', 1, "ddirch constant upper truncation bound")
test_getBound(m, cm, test, 'y4', 'lower', consts$y4l, "invalid constant lower truncation bound")
test_getBound(m, cm, test, 'y4', 'upper', consts$y4u, "invalid constant upper truncation bound")
test_getBound(m, cm, test, 'y5', 'lower', m$mu[1], "stochastically-invalid constant lower truncation bound")
test_getBound(m, cm, test, 'y5', 'upper', consts$y5u, "stochastically-invalid constant upper truncation bound")
test_getBound(m, cm, test, 'y6', 'lower', consts$y6l, "reversed constant lower truncation bound")
test_getBound(m, cm, test, 'y6', 'upper', consts$y6u, "reversed constant upper truncation bound")

# test after a simulate
set.seed(0)
simulate(m)
set.seed(0)
simulate(cm)
test_getBound(m, cm, test, 'y1', 'lower', m$b[1], "non-constant lower truncation bound")
test_getBound(m, cm, test, 'y1', 'upper', m$d[1], "non-constant upper truncation bound")
test_getBound(m, cm, test, 'y2', 'lower', m$y2l[1], "non-constant lower truncation bound")
test_getBound(m, cm, test, 'y2', 'upper', consts$y2u, "non-constant upper truncation bound")
test_getBound(m, cm, test, 'y3', 'lower', exp(m$y2l[1]), "functional lower truncation bound")
test_getBound(m, cm, test, 'y3', 'upper', m$y3u[1], "functional/deterministic upper truncation bound")
test_getBound(m, cm, test, 'b', 'lower', m$a[1], "dunif non-constant lower truncation bound")
test_getBound(m, cm, test, 'b', 'upper', m$d[1], "dunif non-constant upper truncation bound")
test_getBound(m, cm, test, 'a', 'lower', consts$al, "dunif constant lower truncation bound")
test_getBound(m, cm, test, 'a', 'upper', 4, "dunif constant upper truncation bound")
test_getBound(m, cm, test, 'z[1:2]', 'lower', -Inf, "dmnorm constant lower truncation bound")
test_getBound(m, cm, test, 'z[1:2]', 'upper', Inf, "dmnorm constant upper truncation bound")
test_getBound(m, cm, test, 'alpha[1:3]', 'lower', 0, "ddirch constant lower truncation bound")
test_getBound(m, cm, test, 'alpha[1:3]', 'upper', 1, "ddirch constant upper truncation bound")
test_getBound(m, cm, test, 'y4', 'lower', consts$y4l, "invalid constant lower truncation bound")
test_getBound(m, cm, test, 'y4', 'upper', consts$y4u, "invalid constant upper truncation bound")
test_getBound(m, cm, test, 'y5', 'lower', m$mu[1], "stochastically-invalid constant lower truncation bound")
test_getBound(m, cm, test, 'y5', 'upper', consts$y5u, "stochastically-invalid constant upper truncation bound")
test_getBound(m, cm, test, 'y6', 'lower', consts$y6l, "reversed constant lower truncation bound")
test_getBound(m, cm, test, 'y6', 'upper', consts$y6u, "reversed constant upper truncation bound")



