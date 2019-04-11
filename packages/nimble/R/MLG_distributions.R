

dMLG_chol <- nimbleFunction(name = 'dMLG_chol',
    run = function(x = double(1), c = double(1), cholesky = double(2),
                   shape = double(0), rate = double(0),
                   sigma = double(0, default = 1), prec_param = integer(0, default = 1),
                   identity = integer(0, default = 0),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        n <- length(x)
        logdens <- n * (shape * log(rate) - lgamma(shape))
        logdens <- logdens - n * log(sigma)
        if(identity) {
            Vqc <- (x - c) / sigma
        } else {
            det <- sum(log(diag(cholesky)))
            if(prec_param) {
                logdens <- logdens + det
                Vqc <- (cholesky %*% (x - c))[,1] / sigma
            } else {
                logdens <- logdens - det
                Vqc <- forwardsolve(cholesky, x - c) / sigma 
            } 
        }
        logdens <- logdens + shape * sum(Vqc) - rate*sum(exp(Vqc))
        if(log) return(logdens) else return(exp(logdens))
    }
)

rMLG_chol <- nimbleFunction(name = 'rMLG_chol',
    run = function(n = double(0), c = double(1), cholesky = double(2),
                   shape = double(0), rate = double(0),
                   sigma = double(0, default = 1), prec_param = integer(0, default = 1),
                   identity = integer(0, default = 0)) {
        returnType(double(1))

        if(n != 1) {
            stop("rMLG only handles n = 1 at the moment.\n")
        }

        r <- length(c) # or of cholesky?
        w <- log(rgamma(r, shape, rate = rate))
        if(identity) {
            q <- c + w
        } else {
            if(prec_param) {
                q <- c + sigma * forwardsolve(cholesky, w)
            } else {
                q <- (c + sigma * cholesky %*% w)[,1]
            }
        }
        return(q)
    }
)

dcMLG <- nimbleFunction(name = 'dcMLG',
    run = function(x = double(1), data = double(1), offset = double(1), coeff = double(2),
                   prior_cholesky = double(2), prior_sigma = double(0),
                   prior_shape = double(0), prior_rate = double(0),
                   log = integer(0, default = 0)) {
        stop("dcMLG not functional; cMLG distribution only intended for conjugate sampling via 'rcMLG'.")
    }
)


rcMLG <- nimbleFunction(name = 'rcMLG',
    run = function(n = double(0), data = double(1), offset = double(1), coeff = double(2),
                   prior_cholesky = double(2), prior_sigma = double(0),
                   prior_shape = double(0), prior_rate = double(0)) {
        returnType(double(1))        
        if(n != 1) {
            stop("rcMLG only handles n = 1 at the moment.\n")
        }

        ## note this should handle case of prior_cholesky not being precision scale

        ## add 'd' stuff
        r <- dim(cholesky)[1]
        n <- length(data)
        
        w1 <- log(rgamma(n, data, rate = exp(offset)))
        w2 <- log(rgamma(r, prior_shape, rate = prior_rate))
        # new rate should be prior_rate*exp(-(1/prior_sigma)*prior_cholesky %*% prior_mean
        
        H1t <- t(coeff)
        Htw <- H1t %*% w1 + prior_cholesky %*% w2 / prior_sigma

        H <- H1t %*% coeff + t(prior_cholesky) %*% prior_cholesky / (prior_sigma^2)
        cholesky <- chol(t(H) %*% H)
        return(backsolve(cholesky, forwardsolve(t(cholesky), Htw)))
    }
)

dvpois <- nimbleFunction(name = 'dvpois',
    run = function(x = double(1), lambda = double(1), log = integer(0, default = 0)) {
        returnType(double(0))
        logdens <- sum(dpois(x, lambda, log = TRUE))
        if(log) return(logdens) else return(exp(logdens))
    }
)


rvpois <- nimbleFunction(name = 'rvpois',
    run = function(n = double(0), lambda = double(1)) {
        returnType(double(1))
        if(n != 1) {
            stop("rvpois only handles n = 1 at the moment.\n")
        }

        m <- length(lambda)
        return(rpois(m, lambda))
    }
)

registerDistributions(list(
        dvpois = list(
            BUGSdist = 'dvpois(lambda)',
            discrete = TRUE,
            altParams = c('log_lambda = log(lambda)'),
            range = c(0, Inf),
            types    = c('value = double(1)', 'lambda = double(1)', 'log_lambda = double(1)')
    ),
    ## not clear how to set identity, constraints; should 'c' be vector?
    ## set identity = TRUE if no cholesky? would need to check V is cholesky...
    dMLG = list(BUGSdist = 'dMLG(c, cholesky, shape, rate, sigma, Q, K, prec_param, identity, constraints)',
                      Rdist    = c('dMLG_chol(c, cholesky, shape, rate, sigma = 1, prec_param = 1)',
                                   'dMLG_chol(c, cholesky = chol(Q), shape, rate, sigma, prec_param = 1)',
                                   'dMLG_chol(c, cholesky = chol(Q), shape, rate, sigma = 1, prec_param = 1)',
                                   'dMLG_chol(c, cholesky = t(chol(K)), shape, rate, sigma, prec_param = 0)',
                                   'dMLG_chol(c, cholesky = t(chol(K)), shape, rate, sigma = 1, prec_param = 0)'),
                      #altParams= c('Q = calc_dMLG_altParams(cholesky, prec_param, 1)',
                      #             'K = calc_dMLG_altParams(cholesky, prec_param, 0)'),
                types    = c('value = double(1)', 'c = double(1)', 'cholesky = double(2)'),
                conjugacy = list(prior = 'dMLG',
                                 link = 'linear_exp',
                                 dependents = list(
                                     dvpois = list(param = 'log_lambda', contribution_value = 'value',
                                                   contribution_offset = 'offset', contribution_coeff = 'coeff')),
                                 posterior = 'dcMLG(contribution_value, contribution_offset, contribution_coeff, prior_cholesky, prior_sigma, prior_shape, prior_rate)')
         ),
    dcMLG = list(BUGSdist = 'dcMLG(data, offset, coeff, cholesky, sigma, shape, rate)',
                 types    = c('value = double(1)', 'data = double(1)', 'offset = double(1)', 'coeff = double(2)', 'cholesky = double(2)'))
), verbose = FALSE)

calc_dMLG_altParams <- function() {}
