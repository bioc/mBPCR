\name{importCNData} 
\alias{importCNData} 
\title{Import the copy number data}


\description{
Function to  import the raw copy number data from a tab delimited file.
}

\usage{
  importCNData(path, NRowSkip, ifLogRatio=1)
}

\arguments{
\item{path}{path of the tab delimited file containing the copy number data. The file must contain a table, where in the first column there are the names of the probes (snpName),
in the second one, the chromosome to which each probe belongs (the possible values of the chromosomes are: the integers from 1 to 22, 'X' and 'Y'), in the third one, the phisical positions of the probes and in the forth one, the copy number data.}
\item{NRowSkip}{number of row to skip in the file, before the table. The names of the columns are to be skipped.}
\item{ifLogRatio}{denotes if the data are either the log2ratio of raw copy number data or raw copy number data. By default, they are considered as log2ratio data, otherwise (\code{ifLogRatio=0}) they are transformed in log2ratio data.}
}

\value{
A list containing:
\item{\code{snpName}}{an array containing the names of the probes}
\item{\code{chr}}{an array containing the name of the chromosome to which each probe belongs}
\item{\code{position}}{an array containing the physical position of each probe}
\item{\code{logratio}}{an array containing the log2ratio of the raw copy number data}
}

\examples{

###import the 10K data of cell line REC
path <- system.file("extdata", "rec10k.txt", package = "mBPCR")
rec10k <- importCNData(path, NRowSkip=1)
plot(rec10k$position[rec10k$chr == 3], rec10k$logratio[rec10k$chr == 3])	    
}

\keyword{file}
