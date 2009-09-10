\name{writeEstProfile}
\alias{writeEstProfile}
\title{Write the estimated profile of copy number data}


\description{

Function to write nicely the results of the copy number profile estimation. The function either writes the tables
 directly on a tab delimited file or returns the corresponding tables. 
}

\usage{

  writeEstProfile(path='', sampleName='', snpName, chr, position, logratio,
                  chrToBeWritten, estPC, estBoundaries=NULL, postProbT=NULL,
                  regrCurve=NULL, regr=NULL)
}


\arguments{
\item{path}{path of the folder where the user wants to write the results of the estimation (it must end with '\\' in windows, or '//' in linux). If \code{path=''}, they will be written in the working directory. If \code{path=NULL}, the tables will not be written on a file, but only returned by the function.}
\item{sampleName}{name of the sample. If the name of the sample if provided, it is used to named the files.}
\item{snpName}{array containing the name of each probe}
\item{chr}{array containing the name of the chromosome to which each probe belongs. The possible values of the elements of \code{chr} are: the integers from 1 to 22, 'X' and 'Y'.}
\item{position}{array containing the physical position of each probe}
\item{logratio}{array containing the log2ratio of the raw copy number data}
\item{chrToBeWritten}{array containing the name of the estimated chromosomes, of which the user wants to write the results. The possible values of the chromosomes are: the integers from 1 to 22, 'X' and 'Y'.}
\item{estPC}{array containing the estimated copy number profile as a piecewise constant function}
\item{estBoundaries}{list containing the vectors of the estimated breakpoints, for each of the chromosomes mentioned in \code{chrToBeWritten}. If \code{estBoundaries=NULL}, then this information is not written.}
\item{postProbT}{list containing the vectors of the posterior probabilities to be a breakpoint of the estimated breakpoints, for each of the chromosomes mentioned in \code{chrToBeWritten}. If \code{postProbT=NULL}, 
then this information is not written in the file containing the estimated breakpoints.}
\item{regrCurve}{array containing the estimated regression curve. If \code{regrCurve=NULL}, then the file containing this information is not written.}
\item{regr}{choice of the computation of the regression curve. If \code{regr=NULL}, then the regression curve was not computed (then the file containing this information is not written),
 if \code{regr="BRC"} the Bayesian Regression Curve with \eqn{K_2} was computed (BRC with \eqn{K_2}), 
if \code{regr="BRCAk"} the Bayesian Regression Curve Averaging over k was computed (BRCAk).}
}

\details{
The function writes or returns at maximum three tables: 

-one containing the estimated profile with mBPCR (the columns are: 'SNPname', 'chromosome', 'position', 'rawLog2ratio', 'mBPCRestimate')

-one containing a summary about the estimated profile with mBPCR  (the columns are: 'SNPname(start)', 'SNPname(end)', 'chromosome', 'position(start)', 'position(end)', 'nProbes', 'mBPCRestimate' and, eventually, 'breakpointPostProb'). This table is not created if \code{estBoundaries=NULL}.

-one containing the estimated profile with a regression curve (the columns are: 'SNPname', 'chromosome', 'position', 'rawLog2ratio' and the name of the regression curve used). This table is not created if \code{regrCurve=NULL}.
}

\examples{

##import the 10K data of cell line REC
data(rec10k)
##estimation of chromosome 5
results <- estProfileWithMBPCR(rec10k$SNPname, rec10k$Chromosome, rec10k$PhysicalPosition, rec10k$log2ratio, chrToBeAnalyzed=5, maxProbeNumber=2000)
##write the estimated profile of chromosome 5 in a file in  the working directory
writeEstProfile(path='', sampleName='rec10k', rec10k$SNPname, rec10k$Chromosome, rec10k$PhysicalPosition, rec10k$log2ratio, chrToBeWritten=5, results$estPC, results$estBoundaries, results$postProbT)


#### the same result can be obtained in the following way, by using the function computeMBPCR for the estimation
#
##estimation of the global parameters
#param <- estGlobParam(rec10k$log2ratio)
##estimation of chromosome 5
#results <- computeMBPCR(rec10k$log2ratio[rec10k$Chromosome == 5], nu=param$nu, rhoSquare=param$rhoSquare, sigmaSquare=param$sigmaSquare)
##write the estimated profile of chromosome 5 in a file in the working directory
#estPC <- array(dim=length(rec10k$SNPname))
#estBoundaries <- list(dim=1)
#postProbT <- list(dim=1)
#estPC[rec10k$Chromosome == 5] <- results$estPC
#estBoundaries[[1]] <- results$estBoundaries
#postProbT[[1]] <- c(results$postProbT[results$estBoundaries[-results$estK]],1)
#writeEstProfile(path='', sampleName='rec10k', rec10k$SNPname, rec10k$Chromosome, rec10k$PhysicalPosition, rec10k$log2ratio, chrToBeWritten=5, estPC, estBoundaries, postProbT)
}

\keyword{file}
