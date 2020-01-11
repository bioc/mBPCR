\name{mBPCR-internal} 
\alias{computeA10}
\alias{computeLA0Vect}
\alias{computePCReg}
\alias{computeRecursions}
\alias{computeRegrCurve}
\alias{indexLA0}
\title{Internal mBPCR functions}


\description{
Internal functions of package mBPCR.
}


\usage{
computeA10(i, j, y, nu, rhoSquare, sigmaSquare)
computeLA0Vect(y, nu, rhoSquare, sigmaSquare)
computePCReg(y, lA0, lL, lR, nu, rhoSquare, sigmaSquare, kMax=50, regr=NULL)
computeRecursions(lA0, n, kMax=50)
computeRegrCurve(y, typeRegr="BRC", n, kMax=50, lL, lR, lA0, nu, rhoSquare, 
                 sigmaSquare, option)
indexLA0(r, c, n)
}


\details{
These functions are not to be called directly by the user
}

\keyword{internal}
