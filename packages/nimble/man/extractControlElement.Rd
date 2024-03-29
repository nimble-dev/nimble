% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_utils.R
\name{extractControlElement}
\alias{extractControlElement}
\title{Extract named elements from MCMC sampler control list}
\usage{
extractControlElement(controlList, elementName, defaultValue, error)
}
\arguments{
\item{controlList}{control list object, which is passed as an argument to all MCMC sampler setup functions.}

\item{elementName}{character string, giving the name of the element to be extracted from the control list.}

\item{defaultValue}{default value of the control list element, giving the value to be used when the \code{elementName} does not exactly match the name of an element in the \code{controlList}.}

\item{error}{character string, giving the error message to be printed if no \code{defaultValue} is provided and \code{elementName} does not match the name of an element in the \code{controlList}.}
}
\value{
The element of \code{controlList} whose name matches \code{elementName}. If no \code{controlList} name matches \code{elementName}, then \code{defaultValue} is returned.
}
\description{
Extract named elements from MCMC sampler control list
}
\author{
Daniel Turek
}
