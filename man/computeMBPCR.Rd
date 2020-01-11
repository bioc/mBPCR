\name{computeMBPCR} 
\alias{computeMBPCR} 
\title{Estimate the copy number profile}


\description{
Function to  estimate the copy number profile with a piecewise constant function using mBPCR. Eventually, it is possible to estimate the profile with a
 smoothing curve using either the Bayesian Regression Curve with \eqn{K_2} (BRC with \eqn{K_2}) or the Bayesian Regression Curve Averaging over k (BRCAk).  It is also possible
to choose the estimator of the variance of the levels \code{rhoSquare} (i.e. either \eqn{\hat{\rho}_1^2} or \eqn{\hat{\rho}^2}) and by default \eqn{\hat{\rho}_1^2} is used.
}

\usage{
 computeMBPCR(y, kMax=50, nu=NULL, rhoSquare=NULL, sigmaSquare=NULL, typeEstRho=1,
              regr=NULL)
}


\arguments{
\item{y}{array containing the log2ratio of the copy number data}
\item{kMax}{maximum number of segments}
\item{nu}{mean of the segment levels. If \code{nu=NULL}, then the algorithm estimates it on the sample.}
\item{rhoSquare}{variance of the segment levels. If \code{rhoSquare=NULL}, then the algorithm estimates it on the sample.}
\item{sigmaSquare}{variance of the noise. If \code{sigmaSquare=NULL}, then the algorithm estimates it on the sample.}
\item{typeEstRho}{choice of the estimator of \code{rhoSquare}. If \code{typeEstRho=1}, then the algorithm estimates \code{rhoSquare} 
with \eqn{\hat{\rho}_1^2}, while if \code{typeEstRho=0}, it estimates \code{rhoSquare}  with \eqn{\hat{\rho}^2}.}
\item{regr}{choice of the computation of the regression curve. If \code{regr=NULL}, then the regression curve is not computed, 
if \code{regr="BRC"} the Bayesian Regression Curve with \eqn{K_2} is computed (BRC with \eqn{K_2}), if \code{regr="BRCAk"} the Bayesian 
Regression Curve Averaging over k is computed (BRCAk).}
}

\details{
By default, the function estimates the copy number profile with mBPCR and estimating rhoSquare on the sample, using \eqn{\hat{\rho}_1^2}. It is
also possible to use \eqn{\hat{\rho}^2} as estimator of \code{rhoSquare}, by setting \code{typeEstRho=0}, or to directly set the value of the parameter.
 
The function gives also the possibility to estimate the profile with a Bayesian regression curve: if \code{regr="BRC"} the Bayesian Regression Curve with \eqn{K_2} is computed (BRC with \eqn{K_2}), if \code{regr="BRCAk"} the Bayesian 
Regression Curve Averaging over k is computed (BRCAk).    

}

\value{
A list containing: 
\item{\code{estK}}{the estimated number of segments}
\item{\code{estBoundaries}}{the estimated boundaries}
\item{\code{estPC}}{the estimated profile with mBPCR}
\item{\code{regrCurve}}{the estimated bayesian regression curve. It is returned only if \code{regr!=NULL}.}
\item{\code{nu}}{}
\item{\code{rhoSquare}}{}
\item{\code{sigmaSquare}}{}
\item{\code{postProbT}}{for each probe, the posterior probablity to be a breakpoint}
}

\references{
  Rancoita, P. M. V., Hutter, M., Bertoni, F., Kwee, I. (2009).
  Bayesian DNA copy number analysis.  \emph{BMC Bioinformatics} 10: 10.
  \url{http://www.idsia.ch/~paola/mBPCR}  
}

\seealso{\code{\link{estProfileWithMBPCR}}, \code{\link{plotEstProfile}}, \code{\link{writeEstProfile}}, \code{\link{estGlobParam}}}

\examples{
##import the 250K NSP data of chromosome 11 of cell line JEKO-1 
data(jekoChr11Array250Knsp)


##first example 
## we select a part of chromosome 11
y <- jekoChr11Array250Knsp$log2ratio[6400:6900]
p <- jekoChr11Array250Knsp$PhysicalPosition[6400:6900]
##we estimate the profile using the global parameters estimated on the whole genome
##the profile is estimated with mBPCR and with the Bayesian Regression Curve
results <- computeMBPCR(y, nu=-3.012772e-10, rhoSquare=0.0479, sigmaSquare=0.0699, regr="BRC")
plot(p, y)
points(p, results$estPC, type='l', col='red')
points(p, results$regrCurve,type='l', col='green')

###second example 
### we select a part of chromosome 11
#y <- jekoChr11Array250Knsp$log2ratio[10600:11600]
#p <- jekoChr11Array250Knsp$PhysicalPosition[10600:11600]
###we estimate the profile using the global parameters estimated on the whole genome
###the profile is estimated with mBPCR and with the Bayesian Regression Curve Ak
#results <- computeMBPCR(y, nu=-3.012772e-10, rhoSquare=0.0479, sigmaSquare=0.0699, regr="BRCAk")
#plot(p,y)
#points(p, results$estPC, type='l', col='red')
#points(p, results$regrCurve, type='l', col='green')

}


\keyword{regression}
\keyword{smooth}



