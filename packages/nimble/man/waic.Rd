% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_WAIC.R
\name{waic}
\alias{waic}
\alias{getWAIC}
\alias{getWAICdetails}
\alias{buildWAIC}
\alias{WAIC}
\alias{enableWAIC}
\title{Using WAIC}
\description{
Details of the WAIC measure for comparing models. NIMBLE implements an online
WAIC algorithm, computed during the course of the MCMC iterations.
}
\details{
To obtain WAIC, set \code{WAIC = TRUE} in nimbleMCMC. If using a more
customized workflow, set \code{enableWAIC = TRUE} in \code{configureMCMC}
or (if skipping \code{configureMCMC}) in \code{buildMCMC}, followed by
setting \code{WAIC = TRUE} in \code{runMCMC}, if using runMCMC to manage
sample generation.

By default, NIMBLE calculates WAIC using an online algorithm that updates
required summary statistics at each post-burnin iteration of the MCMC.

One can also use \code{calculateWAIC} to run an offline version of the
WAIC algorithm after all MCMC sampling has been done. This allows calculation
of WAIC from a matrix (or dataframe) of posterior samples and also retains
compatibility with WAIC in versions of NIMBLE before 0.12.0. However, the
offline algorithm is less flexible than the online algorithm and only
provides conditional WAIC without the ability to group data points. See
\code{help(calculateWAIC)} for details.
}
\section{\code{controlWAIC} list}{


The \code{controlWAIC} argument is a list that controls the behavior of the
WAIC algorithm and is passed to either \code{configureMCMC} or (if not using
\code{configureMCMC}) \code{buildMCMC}. One can supply any of the following
optional components:

\code{online}: Logical value indicating whether to calculate WAIC during the
course of the MCMC. Default is \code{TRUE} and setting to \code{FALSE} is
primarily for backwards compatibility to allow use of the old
\code{calculateWAIC} method that calculates WAIC from monitored values after
the MCMC finishes.

\code{dataGroups}: Optional list specifying grouping of data nodes,
one element per group, with each list element containing the node names
for the data nodes in that group. If provided, the predictive density values
computed will be the joint density values, one joint density per group.
Defaults to one data node per 'group'. See details.

\code{marginalizeNodes}: Optional set of nodes (presumably latent nodes)
over which to marginalize to compute marginal WAIC (i.e., WAIC based on a
marginal likelihood), rather than the default conditional WAIC (i.e., WAIC
conditioning on all parent nodes of the data nodes). See details.

\code{niterMarginal}: Number of Monte Carlo iterations to use when
marginalizing (default is 1000).

\code{convergenceSet}: Optional vector of numbers between 0 and 1 that
specify a set of shorter Monte Carlo simulations for marginal WAIC
calculation as fractions of the full (\code{niterMarginal}) Monte Carlo
simulation. If not provided, NIMBLE will use 0.25, 0.50, and 0.75.
NIMBLE will report the WAIC, lppd, and pWAIC that would have been obtained
for these smaller Monte Carlo simulations, allowing assessment of the number
of Monte Carlo samples needed for stable calculation of WAIC.

\code{thin}: Logical value for specifying whether to do WAIC calculations
only on thinned samples (default is \code{FALSE}). Likely only useful for
reducing computation when using marginal WAIC.

\code{nburnin_extra}: Additional number of pre-thinning MCMC iterations
to discard before calculating online WAIC. This number is discarded in
addition to the usual MCMC burnin, \code{nburnin}. The purpose of this
option is to allow a user to retain some samples for inspection without
having those samples used for online WAIC calculation (default = 0).
}

\section{Extracting WAIC}{


The calculated WAIC and related quantities can be obtained in various ways
depending on how the MCMC is run. If using \code{nimbleMCMC} and setting
\code{WAIC = TRUE}, see the \code{WAIC} component of the output list. If using
\code{runMCMC} and setting \code{WAIC = TRUE}, either see the \code{WAIC}
component of the output list or use the \code{getWAIC} method of the MCMC
object (in the latter case \code{WAIC = TRUE} is not required). If using
the \code{run} method of the MCMC object, use the \code{getWAIC} method of
the MCMC object.

The output of running WAIC (unless one sets \code{online = FALSE}) is a list
containing the following components:

\code{WAIC}: The computed WAIC, on the deviance scale. Smaller values are
better when comparing WAIC for two models.

\code{lppd}: The log predictive density component of WAIC.

\code{pWAIC}: The pWAIC estimate of the effective number of parameters,
computed using the \emph{p}WAIC2 method of Gelman et al. (2014).

To get further information, one can use the \code{getWAICdetails} method
of the MCMC object.  The result of running \code{getWAICdetails} is a list
containing the following components:

\code{marginal}: Logical value indicating whether marginal (\code{TRUE}) or
conditional (\code{FALSE}) WAIC was calculated.

\code{niterMarginal}: Number of Monte Carlo iterations used in computing
marginal likelihoods if using marginal WAIC.

\code{thin}: Whether WAIC was calculated based only on thinned samples.

\code{online}: Whether WAIC was calculated during MCMC sampling.

\code{nburnin_extra}: Number of additional iterations discarded as burnin,
in addition to original MCMC burnin.

\code{WAIC_partialMC}, \code{lppd_partialMC}, \code{pWAIC_partialMC}: The
computed marginal WAIC, lppd, and pWAIC based on fewer Monte Carlo
simulations, for use in assessing the sensitivity of the WAIC calculation
to the number of Monte Carlo iterations.

\code{niterMarginal_partialMC}: Number of Monte Carlo iterations used for the
values in \code{WAIC_partialMC}, \code{lppd_partialMC}, \code{pWAIC_partialMC}.

\code{WAIC_elements}, \code{lppd_elements}, \code{pWAIC_elements}: Vectors of
individual WAIC, lppd, and pWAIC values, one element per data node (or group
of nodes in the case of specifying \code{dataGroups}). Of use in computing
the standard error of the difference in WAIC between two models, following
Vehtari et al. (2017).
}

\section{Online WAIC}{


As of version 0.12.0, NIMBLE provides enhanced WAIC functionality, with user
control over whether to use conditional or marginal versions of WAIC and
whether to group data nodes. In addition, users are no longer required to
carefully choose MCMC monitors. WAIC by default is now calculated in an online
manner (updating the required summary statistics at each MCMC iteration),
using all post-burnin samples. The WAIC (Watanabe, 2010) is calculated from
Equations 5, 12, and 13 in Gelman et al. (2014) (i.e., using 'pWAIC2').

Note that there is not a unique value of WAIC for a model. By default, WAIC
is calculated conditional on the parent nodes of the data nodes, and the
density values used are the individual density values of the data nodes.
However, by modifying the \code{marginalizeNodes} and \code{dataGroups}
elements of the control list, users can request a marginal WAIC (using a
marginal likelihood that integrates over user-specified latent nodes) and/or
a WAIC based on grouping observations (e.g., all observations in a cluster)
to use joint density values. See the MCMC Chapter of the NIMBLE
\href{https://r-nimble.org/html_manual/cha-mcmc.html}{User Manual}
for more details.

For more detail on the use of different predictive distributions, see Section
2.5 from Gelman et al. (2014) or Ariyo et al. (2019).

Note that based on a limited set of simulation experiments in Hug and Paciorek
(2021) our tentative recommendation is that users only use marginal WAIC if
also using grouping.
}

\examples{
code <- nimbleCode({
  for(j in 1:J) {
    for(i in 1:n) 
      y[j, i] ~ dnorm(mu[j], sd = sigma)
    mu[j] ~ dnorm(mu0, sd = tau)
  }
  sigma ~ dunif(0, 10)
  tau ~ dunif(0, 10)
})
J <- 5
n <- 10
groups <- paste0('y[', 1:J, ', 1:', n, ']') 
y <- matrix(rnorm(J*n), J, n)
Rmodel <- nimbleModel(code, constants = list(J = J, n = n), data = list(y = y),
                      inits = list(tau = 1, sigma = 1))

## Various versions of WAIC available via online calculation.
## Conditional WAIC without data grouping:
conf <- configureMCMC(Rmodel, enableWAIC = TRUE)
## Conditional WAIC with data grouping
conf <- configureMCMC(Rmodel, enableWAIC = TRUE, controlWAIC = list(dataGroups = groups))
## Marginal WAIC with data grouping:
conf <- configureMCMC(Rmodel, enableWAIC = TRUE, controlWAIC =
            list(dataGroups = groups, marginalizeNodes = 'mu'))
\dontrun{
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
output <- runMCMC(Cmcmc, niter = 1000, WAIC = TRUE)
output$WAIC              # direct access
## Alternatively call via the `getWAIC` method; this doesn't require setting
## `waic=TRUE` in `runMCMC`
Cmcmc$getWAIC()          
Cmcmc$getWAICdetails()
}
}
\references{
Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
widely applicable information criterion in singular learning theory.
\emph{Journal of Machine Learning Research} 11: 3571-3594.

Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models.
\emph{Statistics and Computing} 24(6): 997-1016.

Ariyo, O., Quintero, A., Munoz, J., Verbeke, G. and Lesaffre, E. (2019).
Bayesian model selection in linear mixed models for longitudinal data.
\emph{Journal of Applied Statistics} 47: 890-913.

Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC.
\emph{Statistics and Computing} 27: 1413-1432.

Hug, J.E.  and Paciorek, C.J. (2021). A numerically stable online
implementation and exploration of WAIC through variations of the
predictive density, using NIMBLE. \emph{arXiv e-print} <arXiv:2106.13359>.
}
\seealso{
\code{\link{calculateWAIC}} \code{\link{configureMCMC}}
\code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
}
\author{
Joshua Hug and Christopher Paciorek
}
