
# the format of this list should be as follows:
# the name of the element is the name of the BUGS density
# the first string in each element is the BUGS density definition with the BUGS parameters as the first parameters, in order (the parameter names do NOT need to match the names in the BUGS manual since BUGS does not use parameter names), followed by any alternative parameter names
# the second and subsequent strings in an element are any reparameterizations, with the density name as in R/NIMBLE and the canonical parameter names as used in R's math library (and NIMBLE extensions)
# if R/NIMBLE uses the same parameterization (in the same order) and same density name, the element may simply be the BUGS density definition (e.g., dichisq())
# if R/NIMBLE uses the same parameterization but in a different order, then one should define the reordering as a reparameterization (e.g., dbin())



distributionsInputList <- list(
    
    
    ############################################
    #### univariate distributions, discrete ####
    ############################################
    
    
    dbern   = list(BUGSdist = 'dbern(prob)',
                   Rdist    = 'dbinom(size = 1, prob)',
                   discrete = TRUE),
    
    dbin    = list(BUGSdist = 'dbin(prob, size)',
                   Rdist    = 'dbinom(size, prob)',
                   discrete = TRUE),
    
    dcat    = list(BUGSdist = 'dcat(prob)',
                   ##Rdist    = 'dcat(prob, K = length(prob))',
                   Rdist    = 'dcat(prob)',
                   types    = c('value = integer()', 'prob = double(1)'), 
                   discrete = TRUE),
    
    dmulti  = list(BUGSdist = 'dmulti(prob, size)',
   ##                Rdist    = 'dmulti(size, prob, K = length(prob))',
                   Rdist    = 'dmulti(size, prob)',
                   types    = c('value = double(1)', 'prob = double(1)', 'size = integer(0)'),
                   discrete = TRUE),
    
    dnegbin = list(BUGSdist = 'dnegbin(prob, size)',
                   Rdist    = 'dnbinom(size, prob)',
                   discrete = TRUE),
    
    dpois   = list(BUGSdist = 'dpois(lambda)',
                   discrete = TRUE),
    
    
    ##############################################
    #### univariate distributions, continuous ####
    ##############################################
    
    
    dbeta   = list(BUGSdist = 'dbeta(shape1, shape2, mean, sd)',
                   Rdist    = 'dbeta(shape1 = mean^2*(1-mean)/sd^2-mean, shape2 = mean*(1-mean)^2/sd^2+mean-1)',
                   altParams= c('mean = shape1/(shape1+shape2)', 'sd = sqrt(shape1*shape2/((shape1*shape2)^2*(shape1+shape2+1)))')),
    
    dchisq  = list(BUGSdist = 'dchisq(df)'),
    
    ## ddexp   = list('ddexp(location, scale, rate)'),   ## 'ddexp' function not implemented yet?  -DT
    
    dexp    = list(BUGSdist = 'dexp(rate, scale)',
                   Rdist    = 'dexp(rate = 1/scale)',
                   altParams= 'scale = 1/rate'),
    
    dgamma  = list(BUGSdist = 'dgamma(shape, rate, scale, mean, sd)',
                   Rdist    = c('dgamma(shape, scale = 1/rate)', 'dgamma(shape = mean^2/sd^2, scale = sd^2/mean)'),
                   altParams= 'rate = 1/scale'),
    
    ## gen.gamma = list(BUGSdist = 'gen.gamma(r, mu, beta)'),   ## not sure the state of this?  -DT
    
    dlnorm  = list(BUGSdist = 'dlnorm(meanlog, tau, sdlog)',
                   Rdist    = 'dlnorm(meanlog, sdlog = 1/sqrt(tau))',
                   altParams= c('tau = sdlog^-2', 'var = sdlog^2')),
    
    dlogis  = list(BUGSdist = 'dlogis(location, rate, scale)',
                   Rdist    = 'dlogis(location, scale = 1/rate)'),
    
    dnorm   = list(BUGSdist = 'dnorm(mean, tau, sd, var)',
                   Rdist    = c('dnorm(mean, sd = 1/sqrt(tau))', 'dnorm(mean, sd = sqrt(var))'),
                   altParams= c('tau = sd^-2', 'var = sd^2')),
    
    ## dpar    = list(BUGSdist = 'dpar(alpha, c)'),   ## not sure the state of this?  -DT
    
    dt      = list(BUGSdist = 'dt(mu, tau, df)',
                   Rdist    = 'dt_nonstandard(df, mu, sigma = 1/sqrt(tau))',
                   altParams = c('tau = sd^-2')),
    
    dunif   = list(BUGSdist = 'dunif(min, max)'),
    
    dweib   = list(BUGSdist = 'dweib(shape, lambda, scale, rate)',
                   Rdist    = c('dweibull(shape, scale = lambda^(-1/shape))', 'dweibull(shape, scale = 1/rate)'),
                   altParams= 'rate = 1/scale'),
    
    
    ####################################
    #### multivariate distributions ####
    ####################################
    
    
    ddirch  = list(BUGSdist = 'ddirch(alpha)',
##                   Rdist    = 'ddirch(alpha, K = length(alpha))',
                   Rdist    = 'ddirch(alpha)',
                   types    = c('value = double(1)', 'alpha = double(1)')),
    
    dmnorm  = list(BUGSdist = 'dmnorm(mean, prec, cov, chol, prec_param)',
                   Rdist    = c('dmnorm_chol(mean, chol = chol(prec), prec_param = 1)', 'dmnorm_chol(mean, chol = chol(cov), prec_param = 0)', 'dmnorm_chol(mean, chol, prec_param)'),
##                   altParams= c('prec = if(prec_param) crossprod(chol) else inverse(crossprod(chol))', 'cov = if(prec_param) inverse(crossprod(chol)) else crossprod(chol)'),
        altParams= c('prec = chol', 'cov = chol'), ## NOT CORRECT. These are placeholders to get other parts working
        types    = c('value = double(1)', 'mean = double(1)', 'chol = double(2)', 'prec_param = integer()', 'prec = double(2)', 'cov = double(2)')),
    
    ## dmt     = list(BUGSdist = 'dmt(mu, T, k)'),   ## not sure the state of this?  -DT
    
    dwish   = list(BUGSdist = 'dwish(R, df, S)',
                   ##Rdist    = c('dwish_chol(chol = chol(R), df, p = dim(R)[1], scale_param = 0)', 'dwish_chol(chol = chol(S), df, p = dim(S)[1], scale_param = 1)'),
                   Rdist    = c('dwish_chol(chol = chol(R), df, scale_param = 0)', 'dwish_chol(chol = chol(S), df, scale_param = 1)'),
                   ##types    = c('value = double(2)', 'chol = double(2)', 'p = integer()', 'scale_param = integer()'))
                   altParams = c('R = chol', 'S = chol'), ##NOT CORRECT. These are placeholders to get other parts working.
                   types    = c('value = double(2)', 'R = double(2)', 'S = double(2)', 'chol = double(2)', 'scale_param = integer()'))
)





