
runNim <- function(code, ## model created by nimbleCode()
                       constants,
                       data,
                       inits, ## inital values for parameters to estimate
                       monitors, ## vector of parameters to monitor
                       niter, ## number of samples
                       thin,...){ ## thinning rate
  ## build model
  R.model <- nimbleModel(code=code,
                        constants=constants,
                        data=data,
                        inits=inits,
                        check=FALSE,...)
  message('R model created')

  ## configure and build mcmc
  mcmc.spec <- configureMCMC(R.model,
                            print=FALSE,
                            monitors = monitors,
                            thin=thin)
  mcmc <- buildMCMC(mcmc.spec)
  message('MCMC built')

  ## compile model in C++
  C.model <- compileNimble(R.model)
  C.mcmc <- compileNimble(mcmc, project = R.model)
  message('NIMBLE model compiled')

  ## run model
  message('running model')
  C.mcmc$run(niter)
  return(C.mcmc)
}
