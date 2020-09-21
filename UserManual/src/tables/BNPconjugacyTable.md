Table: (#tab:BNPconjugacy) Conjugate relationships for the `dCRP` distribution  supported by NIMBLE’s MCMC engine.

  Baseline Distribution  Mixture (Dependent Node) Distribution  Parameter
  ---------------------  -------------------------------------  -----------
  Beta                   Bernoulli                              `prob`
	                     Binomial                               `prob`
                         Negative Binomial                      `prob`
  Dirichlet              Multinomial                            `prob`
  Gamma                  Poisson                                `lambda`
                         Normal                                 `tau`
                         Gamma                                  `rate`
                         Inverse Gamma                          `scale`
                         Exponential                            `rate`
                         Weibull                                `lambda`
  Inverse Gamma          Normal                                 `var`
  Normal                 Normal                                 `mean`
  Multivariate Normal    Multivariate Normal                    `mean`
  Wishart                Multivariate Normal                    `prec`
  Inverse Wishart        Multivariate Normal                    `cov`




