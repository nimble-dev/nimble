Table: (#tab:conjugacy) Conjugate relationships supported by NIMBLE’s MCMC engine.

  Prior Distribution  Sampling (Dependent Node) Distribution  Parameter
  ------------------  --------------------------------------  -----------
  Beta                Bernoulli                               `prob`
	                  Binomial                                `prob`
                      Negative Binomial                       `prob`
  Dirichlet           Multinomial                             `prob`
                      Categorical                             `prob`
  Flat                Normal                                  `mean`
                      Lognormal                               `meanlog`
  Gamma               Poisson                                 `lambda`
                      Normal                                  `tau`
                      Lognormal                               `taulog`
                      Gamma                                   `rate`
                      Inverse Gamma                           `scale`
                      Exponential                             `rate`
                      Double Exponential                      `rate`
                      Weibull                                 `lambda`
  Halfflat            Normal                                  `sd`
                      Lognormal                               `sdlog`
  Inverse Gamma       Normal                                  `var`
                      Lognormal                               `varlog`
                      Gamma                                   `scale`
                      Inverse Gamma                           `rate`
                      Exponential                             `scale`
                      Double Exponential                      `scale`
  Normal              Normal                                  `mean`
                      Lognormal                               `meanlog`
  Multivariate Normal Multivariate Normal                     `mean`
  Wishart             Multivariate Normal                     `prec`
  Inverse Wishart     Multivariate Normal                     `cov`




