% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions_implementations.R
\name{nimSvd}
\alias{nimSvd}
\alias{svd}
\title{Singular Value Decomposition of a Matrix}
\usage{
nimSvd(x, vectors = "full")
}
\arguments{
\item{x}{a symmetric numeric matrix (double or integer) whose spectral decomposition is to be computed.}

\item{vectors}{character that determines whether to calculate left and right singular vectors.  Can take values \code{'none'}, \code{'thin'} or \code{'full'}.  Defaults to \code{'full'}.  See \sQuote{Details}.}
}
\value{
The singular value decomposition of \code{x} is returned as a \code{\link{nimbleList}} with elements:
\itemize{
\item d length m vector containing the singular values of \code{x}, sorted in decreasing order.
\item v matrix with columns containing the left singular vectors of \code{x}, or an empty matrix if \code{vectors = 'none'}.
\item u matrix with columns containing the right singular vectors of \code{x}, or an empty matrix if \code{vectors = 'none'}.
}
}
\description{
Computes singular values and, optionally, left and right singular vectors of a numeric matrix.
}
\details{
Computes the singular value decomposition of a numeric matrix using the Eigen C++ template library.  

The \code{vectors} character argument determines whether to compute no left and right singular vectors (\code{'none'}), thinned left and right singular vectors (\code{'thin'}), or full left and right singular vectors (\code{'full'}).  For a
matrix \code{x} with dimensions \code{n} and \code{p}, setting \code{vectors = 'thin'} will does the following (quoted from eigen website): 
In case of a rectangular n-by-p matrix, letting m be the smaller value among n and p, there are only m singular vectors; 
the remaining columns of U and V do not correspond to actual singular vectors. 
Asking for thin U or V means asking for only their m first columns to be formed. 
So U is then a n-by-m matrix, and V is then a p-by-m matrix. 
Notice that thin U and V are all you need for (least squares) solving.

Setting \code{vectors = 'full'} will compute full matrices for U and V, so that U will be of size n-by-n, and V will be of size p-by-p.

In a \code{nimbleFunction}, \code{svd} is identical to \code{nimSvd}. 
 
 \code{returnType(svdNimbleList())} can be used within a \code{link{nimbleFunction}} to specify that the function will return a \code{\link{nimbleList}} generated by the \code{nimSvd} function.  \code{svdNimbleList()} can also be used to define a nested \code{\link{nimbleList}} element.  See the \href{https://r-nimble.org/html_manual/cha-welcome-nimble.html}{User Manual} for usage examples.
}
\examples{
 singularValuesDemoFunction <- nimbleFunction(
   setup = function(){
     demoMatrix <- diag(4) + 2
   },
   run = function(){
     singularValues <- svd(demoMatrix)$d
     returnType(double(1))
     return(singularValues)
 })
}
\seealso{
\code{\link{nimEigen}} for spectral decompositions.
}
\author{
NIMBLE development team
}
