ADtestCode1 <- nimbleCode({
  a[1] ~ dnorm(0, 1)
  a[2] ~ dnorm(0, 1)
  y[1] ~ dnorm(a[1], 1)
  y[2] ~ dnorm(a[2], 1)
})
ADtestMod1 <- nimbleModel(code = ADtestCode1, data = list(y = numeric(2)), dimensions = list(y = c(2)))
cADtestMod1 <- compileNimble(ADtestMod1)
