Table: (\#tab:parameterizations) Distribution parameterizations allowed in NIMBLE. The first column
indicates the supported parameterizations for
distributions given in Table \@ref(tab:distributions). The second
column indicates the
relationship to the *canonical* parameterization used in NIMBLE. 

  Parameterization          NIMBLE re-parameterization
  ------------------------- --------------------------------------------------------
  `dbern(prob)`             `dbin(size = 1, prob)`
  `dbeta(shape1, shape2)`   canonical
  `dbeta(mean, sd)`         `dbeta(shape1 = mean^2 * (1-mean) / sd^2 - mean,`
                            `shape2 = mean * (1 - mean)^2 / sd^2 + mean - 1)`
  `dbin(prob, size)`        canonical
  `dcat(prob)`              canonical
  `dchisq(df)`              canonical
  `ddexp(location, scale)`  canonical
  `ddexp(location, rate)`   `ddexp(location, scale = 1 / rate)`
  `ddexp(location, var)`    `ddexp(location, scale = sqrt(var / 2))`
  `ddirch(alpha)`           canonical
  `dexp(rate)`              canonical
  `dexp(scale)`             `dexp(rate = 1/scale)`
  `dgamma(shape, scale)`    canonical
  `dgamma(shape, rate)`     `dgamma(shape, scale = 1 / rate)`
  `dgamma(mean, sd)`        `dgamma(shape = mean^2/sd^2, scale = sd^2/mean)`
  `dinvgamma(shape, rate)`  canonical
  `dinvgamma(shape, scale)` `dgamma(shape, rate = 1 / scale)`
  `dlogis(location, scale)` canonical
  `dlogis(location, rate)`  `dlogis(location, scale = 1 / rate`
  `dlnorm(meanlog, sdlog)`  canonical
  `dlnorm(meanlog, taulog)` `dlnorm(meanlog, sdlog = 1 / sqrt(taulog)`
  `dlnorm(meanlog, varlog)` `dlnorm(meanlog, sdlog = sqrt(varlog)`
  `dmulti(prob, size)`      canonical
  `dmnorm(mean, cholesky, ` canonical (precision)
  `...prec_param=1)` 
  `dmnorm(mean, cholesky, ` canonical (covariance)
  `...prec_param=0)`
  `dmnorm(mean, prec)`      `dmnorm(mean, cholesky = chol(prec), prec_param=1)` 
  `dmnorm(mean, cov)`       `dmnorm(mean, cholesky = chol(cov), prec_param=0)`
  `dmvt(mu, cholesky, df,`  canonical (precision/inverse scale)
  `...prec_param=1)`       
  `dmvt(mu, cholesky, df,`  canonical (scale)
  `...prec_param=0)`
  `dmvt(mu, prec, df)`      `dmvt(mu, cholesky = chol(prec), df, prec_param=1)`
  `dmvt(mu, scale, df)`     `dmvt(mu, cholesky = chol(scale), df, prec_param=0)` 
  `dnegbin(prob, size)`     canonical
  `dnorm(mean, sd)`         canonical
  `dnorm(mean, tau)`        `dnorm(mean, sd = 1 / sqrt(var))`
  `dnorm(mean, var)`        `dnorm(mean, sd = sqrt(var))`
  `dpois(lambda)`           canonical
  `dt(mu, sigma, df)`       canonical
  `dt(mu, tau, df)`         `dt(mu, sigma = 1 / sqrt(tau), df)`
  `dt(mu, sigma2, df)`      `dt(mu, sigma = sqrt(sigma2), df)`
  `dunif(min, max)`         canonical
  `dweib(shape, scale)`     canonical
  `dweib(shape, rate)`      `dweib(shape, scale = 1 / rate)`
  `dweib(shape, lambda)`    `dweib(shape, scale = lambda^(- 1 / shape)`
  `dwish(cholesky, df,`     canonical (scale) 
  `...scale_param=1)`
  `dwish(cholesky, df,`     canonical (inverse scale)
  `...scale_param=0)` 
  `dwish(R, df)`            `dwish(cholesky = chol(R), df, scale_param = 0)`
  `dwish(S, df)`            `dwish(cholesky = chol(S), df, scale_param = 1)`
  `dinvwish(cholesky, df,`  canonical (scale)
  `...scale_param=1)` 
  `dinvwish(cholesky, df,`  canonical (inverse scale)
  `...scale_param=0)`
  `dinvwish(R, df)`         `dinvwish(cholesky = chol(R), df, scale_param = 0)`
  `dinvwish(S, df)`         `dinvwish(cholesky = chol(S), df, scale_param = 1)`

