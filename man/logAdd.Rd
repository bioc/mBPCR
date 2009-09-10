\name{logAdd}
\alias{logAdd}
\title{Overflow-safe computation of the logarithm of a sum}


\description{
Function to compute the logarithm of a sum of small numbers, avoiding overflow.
}

\usage{
  logAdd(x)
}


\arguments{
\item{x}{array or matrix containing the logarithm of the terms of the sum. If \code{x} is a matrix, the function return the results by column.}
}

\value{
If \code{x} is an array, the function returns \eqn{log(sum_i(e^x[i]))}, otherwise it returns an array containing the results by column.
}

\examples{
x <- log(c(0.0001, 0.0003, 0.000006))
y <- logAdd(x) 
##verification that the computation is correct
z <- sum(c(0.0001, 0.0003, 0.000006))
z
exp(y)
}

\keyword{math}
