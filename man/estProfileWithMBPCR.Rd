\name{estProfileWithMBPCR}
\alias{estProfileWithMBPCR}
\title{Estimate and print the copy number profile of some chromosomes of a sample}

\description{
Function to estimate the copy number profile with a piecewise constant function using mBPCR. Eventually, it is possible to estimate the profile with a
 smoothing curve, using either the Bayesian Regression Curve with \eqn{K_2} (BRC with \eqn{K_2}) or the Bayesian Regression Curve Averaging over k (BRCAk). It is also possible
to choose the estimator of the variance of the levels \code{rhoSquare} (i.e. either \eqn{\hat{\rho}_1^2} or \eqn{\hat{\rho}^2}) and by default \eqn{\hat{\rho}_1^2} is used.

}

\usage{
  estProfileWithMBPCR(snpName, chr, position, logratio, chrToBeAnalyzed, maxProbeNumber, 
                      rhoSquare=NULL, kMax=50, nu=NULL, sigmaSquare=NULL, typeEstRho=1, 
                      regr=NULL, hg='hg18')
}


\arguments{
\item{snpName}{array containing the name of each probe}
\item{chr}{array containing the name of the chromosome to which each of the probes belongs. The possible values of the elements of \code{chr} are: the integers from 1 to 22, 'X' and 'Y'.}
\item{position}{array containing the physical position of each probe}
\item{logratio}{array containing the log2ratio of the raw copy number data}
\item{chrToBeAnalyzed}{array containing the name of the chromosomes that the user wants to analyze. The possible values of the chromosomes are: the integers from 1 to 22, 'X' and 'Y'.}
\item{maxProbeNumber}{maximum number of probes that a chromosome (or arm of a chromosome) can have to be analyzed. The procedure of profile estimation
needs the computation of an array of length \eqn{(length(chromosome)+1)*(length(chromosome)+2)/2}. To be sure to have set this parameter
 correctly, try to create the array \code{A <- array(1, dim=(maxProbeNumber+1)*(maxProbeNumber+2)/2)}, before starting with the estimation procedure.}
\item{rhoSquare}{variance of the segment levels. If \code{rhoSquare=NULL}, then the algorithm estimates it on the sample.}
\item{kMax}{maximum number of segments}
\item{nu}{mean of the segment levels. If \code{nu=NULL}, then the algorithm estimates it on the sample.}
\item{sigmaSquare}{variance of the noise. If \code{sigmaSquare=NULL}, then the algorithm estimates it on the sample.}
\item{typeEstRho}{choice of the estimator of \code{rhoSquare}. If \code{typeEstRho=1}, then the algorithm estimates \code{rhoSquare} 
with \eqn{\hat{\rho}_1^2}, while if \code{typeEstRho=0}, it estimates \code{rhoSquare} with \eqn{\hat{\rho}^2}.}
\item{regr}{choice of the computation of the regression curve. If \code{regr=NULL}, then the regression curve is not computed, 
if \code{regr="BRC"} the Bayesian Regression Curve is computed (BRC with \eqn{K_2}), if \code{regr="BRCAk"} the Bayesian 
Regression Curve Averaging over k is computed (BRCAk).}
\item{hg}{genome build used for retrieving the base positions of the centromeres in case the chromosomes need to be divided into two parts for the estimation (see explanation of \code{maxProbeNumber}). Current available options are: \code{'hg18'}, \code{'hg19'} and \code{'hg38'}.}
}

\details{
By default, the function estimates the copy number profile with mBPCR and estimating rhoSquare on the sample, using \eqn{\hat{\rho}_1^2}. It is
also possible to use \eqn{\hat{\rho}^2} as estimator of \code{rhoSquare}, by setting \code{typeEstRho=0}, or to directly set the value of the parameter.
 
The function gives also the possibility to estimate the profile with a Bayesian regression curve: if \code{regr="BRC"} the Bayesian Regression Curve with \eqn{K_2} is computed (BRC with \eqn{K_2}), if \code{regr="BRCAk"} the Bayesian 
Regression Curve Averaging over k is computed (BRCAk).    

See function \code{writeEstProfile}, to have the results in nicer tables or to write them on files. 
}


\value{
A list containing:
\item{\code{estPC}}{an array containing the estimated profile with mBPCR}
\item{\code{estBoundaries}}{the list of estimated breakpoints for each of the analyzed chomosomes}
 \item{\code{postProbT}}{the list of the posterior probablity to be a breakpoint for each estimated breakpoint of the analyzed chomosomes}
\item{\code{regrCurve}}{an array containing the estimated bayesian regression curve}

\code{estPC} and \code{regrCurve} have the same length of \code{logratio}, hence their components, 
corresponding to the not analyzed chromosomes, are equal to \code{NA}.
}

\references{
  Rancoita, P. M. V., Hutter, M., Bertoni, F., Kwee, I. (2009).
  Bayesian DNA copy number analysis.  \emph{BMC Bioinformatics} 10: 10.
  \url{http://www.idsia.ch/~paola/mBPCR}  
}

\seealso{\code{\link{plotEstProfile}}, \code{\link{writeEstProfile}}, \code{\link{computeMBPCR}}}


\examples{
##import the 10K data of cell line REC  
data(rec10k)
##estimation of the profile of chromosome 5
results <- estProfileWithMBPCR(rec10k$SNPname, rec10k$Chromosome, rec10k$PhysicalPosition, 
rec10k$log2ratio, chrToBeAnalyzed=5, maxProbeNumber=2000)
##plot the estimated profile of chromosome 5
y <- rec10k$log2ratio[rec10k$Chromosome == 5]
p <- rec10k$PhysicalPosition[rec10k$Chromosome == 5]
plot(p, y)
points(p, results$estPC[rec10k$Chromosome == 5], type='l', col='red')

###for the estimation of the profile of all chromosomes
#results <- estProfileWithMBPCR(rec10k$SNPname, rec10k$Chromosome, rec10k$PhysicalPosition, 
#rec10k$log2ratio, chrToBeAnalyzed=c(1:22,'X'), maxProbeNumber=2000)

}

\keyword{regression}
\keyword{smooth}

