% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_xvar.R
\name{check.xvar}
\alias{check.xvar}
\title{Check covariance of exposure values are well defined.}
\usage{
check.xvar(Xvar)
}
\arguments{
\item{Xvar}{Input of covariance matrix of exposure values}
}
\value{
Xvar
}
\description{
Function that verifies \code{Xvar} are well defined for
  \code{paf.confidence} and \code{pif.confidence}
}
\examples{

#Example 1
Xvar <- 0.2
check.xvar(Xvar)

}
\seealso{
\code{\link{check.confidence}}, \code{\link{check.thetas}}, 
  \code{\link{check.cft}}, \code{\link{check.rr}}, 
  \code{\link{check.exposure}}, \code{\link{check.integrals}}
}
\author{
Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}

Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
}
\keyword{internal}
