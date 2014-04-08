\name{providerClasses}

\alias{providerClasses}
\alias{providerClasses,RefNet-method}
\alias{\S4method{providerClasses}{RefNet}}

\title{providerClasses}

\description{

  RefNet providerClasses are the names of the groups of
  available data sources, currently "PSICQUIC" and "native"
}

\usage{
\S4method{providerClasses}{RefNet}(object)
}

\arguments{
  \item{object}{a \code{RefNet} object.}
}

\value{
This method returns a list of the current groups.
}
\author{Paul Shannon}

\examples{
  refnet <- RefNet()
  providerClasses(refnet)

}
\seealso{
  RefNet,
  interactions,
  addStandardNames,
  IDMapper
}


\keyword{utilities}
