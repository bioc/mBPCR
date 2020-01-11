\name{plotEstProfile}
\alias{plotEstProfile}
\title{Plot the estimated profile of copy number data}


\description{
Function to plot the estimated profiles of copy number data.

}

\usage{
  plotEstProfile(sampleName='', chr, position, logratio, chrToBePlotted, estPC, 
                 maxProbeNumber, legendPosition='bottomleft', regrCurve=NULL, 
                 regr=NULL, hg='hg18')

}



\arguments{
\item{sampleName}{name of the sample, if the user wants to put it in the title of the graph}
\item{chr}{array containing the name of the chromosome to which each probe belongs. The possible values of the elements of \code{chr} are: the integers from 1 to 22, 'X' and 'Y'.}
\item{position}{array containing the physical position of each probe}
\item{logratio}{array containing the log2ratio of the raw copy number data}
\item{chrToBePlotted}{array containing the name of the estimated chromosomes, that the user wants to plot. The possible values of the chromosomes are: the integers from 1 to 22, 'X' and 'Y'.}
\item{estPC}{array containing the estimated copy number profile as a piecewise constant function. If \code{estPC=NULL}, only the estimated Bayesian regression curve is plotted.}
\item{maxProbeNumber}{maximum number of probes that a chromosome (or arm of a chromosome) can have to be analyzed. The procedure of profile estimation
needs the computation of an array of length \eqn{(length(chromosome)+1)*(length(chromosome)+2)/2}. To be sure to have set this parameter
 correctly, try to create the array \code{A <- array(1, dim=(maxProbeNumber+1)*(maxProbeNumber+2)/2)}, before starting with the estimation procedure.}
\item{legendPosition}{string containing the position of the legend in the plot. The possible values are the same used in the function \code{plot}.}
\item{regrCurve}{array containing the estimated regression curve. If \code{regrCurve=NULL}, then the estimated Bayesian regression curve is not plotted. If \code{regrCurve!=NULL} and also \code{estPC!=NULL} both estimated profiles are plotted on the same graph.}
\item{regr}{choice of the computation of the regression curve. If \code{regr=NULL}, then the regression curve was not computed (then the estimated Bayesian regression curve is not plotted),
 if \code{regr="BRC"} the Bayesian Regression Curve was computed (mBRC with \eqn{K_2}), 
if \code{regr="BRCAk"} the Bayesian Regression Curve Averaging over k was computed (BRCAk).}
\item{hg}{genome build used for retrieving the base positions of the centromeres in case the chromosomes need to be divided into two parts for the estimation (see explanation of \code{maxProbeNumber}). Current available options are: \code{'hg18'}, \code{'hg19'} and \code{'hg38'}.}
}

\details{
The function plots the estimated profiles of the chromosomes of \code{chrToBePlotted}, separately. 
}

\examples{
##import the 10K data of cell line REC
data(rec10k)
##estimation of chromosomes 3 and 5
results <- estProfileWithMBPCR(rec10k$SNPname, rec10k$Chromosome, rec10k$PhysicalPosition, 
rec10k$log2ratio, chrToBeAnalyzed=c(3,5), maxProbeNumber=2000)
##plot the corresponding estimated profiles
plotEstProfile(sampleName='rec10k', rec10k$Chromosome, rec10k$PhysicalPosition, rec10k$log2ratio, 
chrToBePlotted=c(3,5), results$estPC, maxProbeNumber=2000)
}

\keyword{hplot}
\keyword{print}

