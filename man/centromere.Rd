\name{centromere}
\alias{centromere}
\title{Retrieve base positions of centromeres}


\description{
Function to retrieve base positions of the centromere of a specific chromosome.
}

\usage{
  centromere(chr, hg='hg18')
}


\arguments{
\item{chr}{chromosome of which you want to retrieve the centromere base positions.}
\item{hg}{genome build used for retrieving the centromere base positions of the selected chromosome. Current available options are: \code{'hg18'}, \code{'hg19'} and \code{'hg38'}.}
}

\value{
The function returns the start and end base positions of the centromere of the selected chromosome, by using the specified genome build. The function is based on the annotation provided in the package \code{GWASTools}.
}

\examples{
#centromere base positions of chromosome 1 in genome build hg18
centromere(1, hg='hg18')
}

\keyword{manip}
