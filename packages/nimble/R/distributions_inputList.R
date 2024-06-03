
# the format of this list should be as follows:
# the name of the element is the name of the BUGS density
# the first string in each element is the BUGS density definition with the BUGS parameters as the first parameters, in order (the parameter names do NOT need to match the names in the BUGS manual since BUGS does not use parameter names), followed by any alternative parameter names
# the second and subsequent strings in an element are any reparameterizations, with the density name as in R/NIMBLE and the canonical parameter names as used in R's math library (and NIMBLE extensions)
# if R/NIMBLE uses the same parameterization (in the same order) and same density name, the element may simply be the BUGS density definition (e.g., dichisq())
# if R/NIMBLE uses the same parameterization but in a different order, then one should define the reordering as a reparameterization (e.g., dbin())


# See end for addition of buildDerivs field

distributionsInputList <- list(
    
    
    ############################################
    #### univariate distributions, discrete ####
    ############################################
    
    
    dbern   = list(BUGSdist = 'dbern(prob)',
                   Rdist    = 'dbinom(size = 1, prob)',
                   discrete = TRUE,
                   range    = c(0, 1),
                   pqAvail  = TRUE),
    
    dbin    = list(BUGSdist = 'dbin(prob, size)',
                   Rdist    = 'dbinom(size, prob)',
                   discrete = TRUE,
                   range    = c('lower = 0', 'upper = size'),
                   pqAvail  = TRUE,
                   alias    = 'dbinom'),
    
    dcat    = list(BUGSdist = 'dcat(prob)',
                   Rdist    = 'dcat(prob)',
                   altParams= c('k = length(prob)'),
                   types    = c('prob = double(1)'),
                   range    = c(1, Inf), 
                   discrete = TRUE),
    
    ## construct used to enforce constraints - 0/1 random variable depending on if cond is TRUE
    dconstraint = list(BUGSdist = 'dconstraint(cond)',
                       range    = c(0, 1),
                       discrete = TRUE),

    ## construct used to enforce censoring.
    ## takes values 0,1,...,len(c), depending on which interval t falls into
    dinterval     = list(BUGSdist = 'dinterval(t, c)',
                         types    = c('c = double(1)'),
                         range    = c(0, Inf),
                         discrete = TRUE),

    dnegbin = list(BUGSdist = 'dnegbin(prob, size)',
                   Rdist    = 'dnbinom(size, prob)',
                   discrete = TRUE,
                   range    = c(0, Inf),
                   pqAvail  = TRUE,
                   alias    = 'dnbinom'),
    
    dpois   = list(BUGSdist = 'dpois(lambda)',
                   discrete = TRUE,
                   range    = c(0, Inf),
                   pqAvail  = TRUE),
    
    
    ##############################################
    #### univariate distributions, continuous ####
    ##############################################
    
    
    dbeta   = list(BUGSdist = 'dbeta(shape1, shape2, mean, sd)',
                   Rdist    = 'dbeta(shape1 = mean^2*(1-mean)/sd^2-mean, shape2 = mean*(1-mean)^2/sd^2+mean-1)',
                   altParams= c('mean = shape1/(shape1+shape2)',
                                'sd = sqrt(shape1*shape2/((shape1 + shape2)^2*(shape1+shape2+1)))'),
                   range    = c(0, 1),
                   pqAvail  = TRUE),
    
    dchisq  = list(BUGSdist = 'dchisq(df)',
                   range    = c(0, Inf),
                   pqAvail  = TRUE,
                   alias    = 'dchisqr'),

    ddexp   = list(BUGSdist = 'ddexp(location, rate, scale, var)',
                   Rdist    = c('ddexp(location, scale = 1/rate)',
                                'ddexp(location, scale = sqrt(var/2))'),
                   altParams= c('rate = 1/scale',
                                'var = 2*scale^2'),
                   pqAvail  = TRUE,
                   alias    = 'dlaplace'),  
    
    dexp    = list(BUGSdist = 'dexp(rate, scale)',
                   Rdist    = 'dexp_nimble(rate = 1/scale)',
                   altParams= 'scale = 1/rate',
                   range    = c(0, Inf),
                   pqAvail  = TRUE),

    dflat   = list(BUGSdist = 'dflat()',
                   pqAvail  = FALSE),
    
    dhalfflat   = list(BUGSdist = 'dhalfflat()',
                       range    = c(0, Inf),
                       pqAvail  = FALSE),
    
    dgamma  = list(BUGSdist = 'dgamma(shape, rate, scale, mean, sd)',
                   Rdist    = c('dgamma(shape, scale = 1/rate)',
                                'dgamma(shape = mean^2/sd^2, scale = sd^2/mean)'),
                   altParams= c('rate = 1/scale',
                                'mean = scale*shape',
                                'sd = scale * sqrt(shape)'),
                   range    = c(0, Inf),
                   pqAvail  = TRUE),

    # (shape,scale) is BUGSdist as scale provides conjugacy
    # calculation of shape/scale from mean/sd not obvious
    # (solution to cubic polynomial) so not using as alternative param
    dinvgamma  = list(BUGSdist = 'dinvgamma(shape, scale, rate)',
                      Rdist    = c('dinvgamma(shape, rate = 1/scale)'),
                      altParams= c('scale = 1/rate',
                                   'mean = 1 / (rate * (max(shape,1)-1))',
                                   'sd = 1 / (rate * (max(shape,1)-1) * sqrt(max(shape,2)-2))'), # max ensures Inf moment when appropriate
                      range    = c(0, Inf),
                      pqAvail  = TRUE),
    
    # intended solely for use in dhalfflat conjugacy                           
    dsqrtinvgamma  = list(BUGSdist = 'dsqrtinvgamma(shape, scale, rate)',
                          Rdist    = c('dsqrtinvgamma(shape, rate = 1/scale)'),
                          range    = c(0, Inf),
                          pqAvail  = FALSE),
    
    ## gen.gamma = list(BUGSdist = 'gen.gamma(r, mu, beta)'),   ## not sure the state of this?  -DT
    
    dlnorm  = list(BUGSdist = 'dlnorm(meanlog, taulog, sdlog, varlog)',
                   Rdist    = c('dlnorm(meanlog, sdlog = 1/sqrt(taulog))',
                                'dlnorm(meanlog, sdlog = sqrt(varlog))'),
                   altParams= c('taulog = sdlog^-2',
                                'varlog = sdlog^2'),
                   range    = c(0, Inf),
                   pqAvail  = TRUE),
    
    dlogis  = list(BUGSdist = 'dlogis(location, rate, scale)',
                   Rdist    = 'dlogis(location, scale = 1/rate)',
                   altParams= 'rate = 1/scale',
                   pqAvail  = TRUE),
    
    dnorm   = list(BUGSdist = 'dnorm(mean, tau, sd, var)',
                   Rdist    = c('dnorm(mean, sd = 1/sqrt(tau))',
                                'dnorm(mean, sd = sqrt(var))'),
                   altParams= c('tau = sd^-2',
                                'var = sd*sd'),
                   pqAvail  = TRUE),
    
    dt      = list(BUGSdist = 'dt(mu, tau, df, sigma, sigma2)',
                   Rdist    = c('dt_nonstandard(df, mu, sigma = 1/sqrt(tau))',
                                'dt_nonstandard(df, mu, sigma = sqrt(sigma2))'),
                   altParams= c('tau = sigma^-2',
                                 'sigma2 = sigma^2'),
                   pqAvail  = TRUE),
    
    dunif   = list(BUGSdist = 'dunif(min, max, mean, sd)',
                   Rdist    = c('dunif(min = mean - sqrt(3)*sd, max = mean + sqrt(3)*sd)'),
                   altParams= c('mean = (min + max)/2',
                                 'sd = (max - min)/sqrt(12)'),
                   range    = c('lower = min', 'upper = max'),
                   pqAvail  = TRUE),
    
    dweib   = list(BUGSdist = 'dweib(shape, lambda, scale, rate)',
                   Rdist    = c('dweibull(shape, scale = lambda^(-1/shape))',
                                'dweibull(shape, scale = 1/rate)'),
                   altParams= c('rate = 1/scale',
                                'lambda = scale^(-shape)'),
                   range    = c(0, Inf),
                   pqAvail  = TRUE,
                   alias    = 'dweibull'),
    
    
    ####################################
    #### multivariate distributions ####
    ####################################
    
    
    dcar_normal = list(BUGSdist = 'dcar_normal(adj, weights, num, tau, c, zero_mean)',
                       Rdist    = c('dcar_normal(adj, weights,           num, tau, c,                                zero_mean = 0)',
                                    'dcar_normal(adj, weights,           num, tau, c = CAR_calcNumIslands(adj, num), zero_mean    )',
                                    'dcar_normal(adj, weights,           num, tau, c = CAR_calcNumIslands(adj, num), zero_mean = 0)',
                                    'dcar_normal(adj, weights = adj/adj, num, tau, c,                                zero_mean    )',
                                    'dcar_normal(adj, weights = adj/adj, num, tau, c,                                zero_mean = 0)',
                                    'dcar_normal(adj, weights = adj/adj, num, tau, c = CAR_calcNumIslands(adj, num), zero_mean    )',
                                    'dcar_normal(adj, weights = adj/adj, num, tau, c = CAR_calcNumIslands(adj, num), zero_mean = 0)'),
                       types    = c('value = double(1)', 'adj = double(1)', 'weights = double(1)', 'num = double(1)', 'tau = double(0)', 'c = double(0)', 'zero_mean = double(0)'),
                       mixedSizes = TRUE,
                       alias    = 'car.normal'),
    
    dcar_proper = list(BUGSdist = 'dcar_proper(mu, C, adj, num, M, tau, gamma, evs)',
                       Rdist    = c('dcar_proper(mu, C,                       adj, num, M,                  tau, gamma, evs = CAR_calcEVs3(C, adj, num))',
                                    'dcar_proper(mu, C = CAR_calcC(adj, num), adj, num, M = CAR_calcM(num), tau, gamma, evs = CAR_calcEVs2(   adj, num))'),
                       types    = c('value = double(1)', 'mu = double(1)', 'C = double(1)', 'adj = double(1)', 'num = double(1)', 'M = double(1)', 'tau = double(0)', 'gamma = double(0)', 'evs = double(1)'),
                       mixedSizes = TRUE,
                       alias    = 'car.proper'),

    ddirch  = list(BUGSdist = 'ddirch(alpha)',
                   Rdist    = 'ddirch(alpha)',
                   types    = c('value = double(1)', 'alpha = double(1)'),
                   range    = c(0, 1),
                   alias    = 'ddirich'),
    
    dmnorm  = list(BUGSdist = 'dmnorm(mean, prec, cov, cholesky, prec_param)',
                   Rdist    = c('dmnorm_chol(mean, cholesky = chol(prec), prec_param = 1)',
                                'dmnorm_chol(mean, cholesky = chol(cov), prec_param = 0)',
                                'dmnorm_chol(mean, cholesky, prec_param)'),
                   altParams= c('prec = calc_dmnormAltParams(cholesky, prec_param, 1)',
                                'cov = calc_dmnormAltParams(cholesky, prec_param, 0)'),
                   types    = c('value = double(1)', 'mean = double(1)', 'cholesky = double(2)', 'prec = double(2)', 'cov = double(2)')),
    
    dmulti  = list(BUGSdist = 'dmulti(prob, size)',
                   Rdist    = 'dmulti(size, prob)',
                   types    = c('value = double(1)', 'prob = double(1)'),
                   range    = c(0, Inf),
                   discrete = TRUE,
                   alias    = 'dmultinom'),
    
    dmvt  = list(BUGSdist = 'dmvt(mu, prec, df, scale, cholesky, prec_param)',
                   Rdist    = c('dmvt_chol(mu, cholesky = chol(prec), df = df, prec_param = 1)',
                                'dmvt_chol(mu, cholesky = chol(scale), df = df, prec_param = 0)',
                                'dmvt_chol(mu, cholesky, df = df, prec_param)'),
                   altParams= c('prec = calc_dmnormAltParams(cholesky, prec_param, 1)',
                                'scale = calc_dmnormAltParams(cholesky, prec_param, 0)'),
                   types    = c('value = double(1)', 'mu = double(1)', 'cholesky = double(2)', 'df = double(0)', 'prec = double(2)', 'scale = double(2)')),
    
    dlkj_corr_cholesky  = list(BUGSdist = 'dlkj_corr_cholesky(eta, p)',
                   Rdist    = 'dlkj_corr_cholesky(eta, p)',
                   types    = c('value = double(2)')),

    dwish   = list(BUGSdist = 'dwish(R, df, S, cholesky, scale_param)',
                   Rdist    = c('dwish_chol(cholesky = chol(R), df, scale_param = 0)',
                                'dwish_chol(cholesky = chol(S), df, scale_param = 1)',
                                'dwish_chol(cholesky, df, scale_param)'),
                   altParams= c('R = calc_dwishAltParams(cholesky, scale_param, 0)',
                                'S = calc_dwishAltParams(cholesky, scale_param, 1)'),
                   alias    = 'dwishart',
                   types    = c('value = double(2)', 'R = double(2)', 'S = double(2)', 'cholesky = double(2)')),

    dinvwish   = list(BUGSdist = 'dinvwish(S, df, R, cholesky, scale_param)',
                      Rdist    = c('dinvwish_chol(cholesky = chol(S), df, scale_param = 1)',
                                   'dinvwish_chol(cholesky = chol(R), df, scale_param = 0)'),
                      altParams= c('R = calc_dwishAltParams(cholesky, scale_param, 0)',
                                   'S = calc_dwishAltParams(cholesky, scale_param, 1)'),
                      alias    = 'dinvwishart',
                      types    = c('value = double(2)', 'S = double(2)', 'R = double(2)', 'cholesky = double(2)'))
)

distsNotAllowedInAD <- c(
  paste0('d', c('interval', 'constraint'))
)

# Add buildDerivs=TRUE to all distributions except distsNotAllowedInAD
distributionsInputList <-
  seq_along(distributionsInputList) |>
  lapply(\(i) {
    this_name <- names(distributionsInputList)[i]
    buildDerivs <- !(this_name %in% distsNotAllowedInAD)
    c(distributionsInputList[[i]], list(buildDerivs=buildDerivs))
  }) |> setNames(names(distributionsInputList))
