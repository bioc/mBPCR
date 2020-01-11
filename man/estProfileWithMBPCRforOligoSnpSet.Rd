\name{estProfileWithMBPCRforOligoSnpSet}
\alias{estProfileWithMBPCRforOligoSnpSet}
\title{Estimate and print the copy number profile of some chromosomes of samples in an oligoSnpSet object}

\description{
Function to estimate the copy number profile with a piecewise constant function using mBPCR. Eventually, it is possible to estimate the profile with a
 smoothing curve, using either the Bayesian Regression Curve with \eqn{K_2} (BRC with \eqn{K_2}) or the Bayesian Regression Curve Averaging over k (BRCAk). It is also possible
to choose the estimator of the variance of the levels \code{rhoSquare} (i.e. either \eqn{\hat{\rho}_1^2} or \eqn{\hat{\rho}^2}) and by default \eqn{\hat{\rho}_1^2} is used.

}

\usage{
  estProfileWithMBPCRforOligoSnpSet(sampleData, sampleToBeAnalyzed, chrToBeAnalyzed, 
                      maxProbeNumber, ifLogRatio=1, rhoSquare=NULL, kMax=50, nu=NULL, 
                      sigmaSquare=NULL, typeEstRho=1, regr=NULL, hg='hg18')

}


\arguments{
\item{sampleData}{object of type oligoSnpSet. The following fields must not be empty: \code{assayData(sampleData)$copyNumber} (it contains the raw copy number values with scale log2 multiplied by 100 and transformed as integers), \code{featureNames(featureData(sampleData))} (it contains the names of the SNPs), \code{featureData(sampleData)$chromosome}  (it contains the names of the chromosomes to which each of the SNPs belongs), \code{featureData(sampleData)$position}  (it contains the physical positions of the SNPs).}
\item{sampleToBeAnalyzed}{vector containing the number of the columns corresponding to the samples the user wants to analyze.}
\item{chrToBeAnalyzed}{array containing the name of the chromosomes that the user wants to analyze. The possible values of the chromosomes are: the integers from 1 to 22, 'X' and 'Y'.}
\item{maxProbeNumber}{maximum number of probes that a chromosome (or arm of a chromosome) can have to be analyzed. The procedure of profile estimation
needs the computation of an array of length \eqn{(length(chromosome)+1)*(length(chromosome)+2)/2}. To be sure to have set this parameter
 correctly, try to create the array \code{A <- array(1, dim=(maxProbeNumber+1)*(maxProbeNumber+2)/2)}, before starting with the estimation procedure.}
\item{ifLogRatio}{denotes whether the original log2 data were centered at zero (i.e. they were in log2ratio scale) or not. By default, they are considered as derived by log2ratio data (\code{ifLogRatio=1}), otherwise (\code{ifLogRatio=0}) they are transformed in order to be derived from log2ratio data.}
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


}


\value{
A list containing:
\item{estPC}{an oligoSnpSet equal to sampleData apart from the field \code{assayData(estPC)$copyNumber}, which contains the estimated profile with mBPCR in scale log2ratio multiplied by 100}
\item{regrCurve}{an oligoSnpSet equal to sampleData apart from the field \code{assayData(regrCurve)$copyNumber}, which contains the estimated bayesian regression curve in scale log2ratio multiplied by 100. This object is returned only if \code{regr!=NULL}.}

The matrices \code{assayData(estPC)$copyNumber} and \code{assayData(regrCurve)$copyNumber} have the same dimension of \code{assayData(sampleData)$copyNumber}, hence their elements, 
corresponding to the not analyzed chromosomes and samples, are equal to \code{NA}.
}

\references{
  Rancoita, P. M. V., Hutter, M., Bertoni, F., Kwee, I. (2009).
  Bayesian DNA copy number analysis.  \emph{BMC Bioinformatics} 10: 10.
  \url{http://www.idsia.ch/~paola/mBPCR}  
}

\seealso{\code{\link{estProfileWithMBPCR}}, \code{\link{computeMBPCR}}}


\examples{
###import an example of oligoSnpSet data
#data(oligoSetExample, package="oligoClasses")
##estimation of chromosome 2 in sample 1 
#r <-estProfileWithMBPCRforOligoSnpSet(oligoSet, sampleToBeAnalyzed=1, chrToBeAnalyzed=2, 
#maxProbeNumber=1000, ifLogRatio=0, rhoSquare=0.0889637)
##plot of the estimated chromosome
#library(SNPchip)
#cc <- r$estPC
#cc1 <- cc[chromosome(cc) == "2",1]
#par(las=1)
#plot(position(cc1), copyNumber(cc1)/100, ylim=c(-0.23, 0.1), ylab="copy number",
#xlab="base position")
}

\keyword{regression}
\keyword{smooth}
