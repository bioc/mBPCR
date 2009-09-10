\name{estGlobParam}
\alias{estGlobParam}
\title{Estimate global parameters of copy number data}


\description{
Function to estimate the global parameters of copy number data: the mean and the variance of the segment levels (called \code{nu} and \code{rhoSquare}, respectively), the variance of the noise (\code{sigmaSquare}). It is possible
to choose the estimator of \code{rhoSquare} (i.e. either \eqn{\hat{\rho}_1^2} or \eqn{\hat{\rho}^2}) and by default \eqn{\hat{\rho}_1^2} is used.
}

\usage{

  estGlobParam(y, nu=NULL, rhoSquare=NULL, sigmaSquare=NULL, typeEstRho=1)
}


\arguments{
\item{y}{array containing the log2ratio of the copy number data}
\item{nu}{mean of the segment levels. If \code{nu=NULL}, then the algorithm estimates it on the sample.}
\item{rhoSquare}{variance of the segment levels. If \code{rhoSquare=NULL}, then the algorithm estimates it on the sample.}
\item{sigmaSquare}{variance of the noise. If \code{sigmaSquare=NULL}, then the algorithm estimates it on the sample.}
\item{typeEstRho}{choice of the estimator of \code{rhoSquare}. If \code{typeEstRho=1}, then the algorithm estimates \code{rhoSquare} 
with \eqn{\hat{\rho}_1^2}, while if \code{typeEstRho=0}, it estimates \code{rhoSquare} with \eqn{\hat{\rho}^2}.}
}

\value{
A list containing:
\item{\code{nu}}{}
\item{\code{rhoSquare}}{}
\item{\code{sigmaSquare}}{}
}

\references{
  Rancoita, P. M. V., Hutter, M., Bertoni, F., Kwee, I. (2009).
  Bayesian DNA copy number analysis.  \emph{BMC Bioinformatics} 10: 10.
  \url{http://www.idsia.ch/~paola/mBPCR}  
}


\examples{
##import the 10K data of cell line REC  
data(rec10k)
##estimation of all the global parameters (the variance of the segment is estimated with \eqn{\hat{\rho}^2_1})
estGlobParam(rec10k$log2ratio)
}

\keyword{univar}
