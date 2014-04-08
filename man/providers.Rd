\name{providers}

\alias{providers}
\alias{providers,RefNet-method}
\alias{\S4method{providers}{RefNet}}

\title{providers}

\description{

  RefNet providers are the names of data sources, many from PSICQUIC,
  some from RefNet.db.  
}

\usage{
\S4method{providers}{RefNet}(object)
}

\arguments{
  \item{object}{a \code{RefNet} object.}
}

\value{
This method returns a named list of these
  data sources, arrayed in two named groups, "PSICQUIC" and "native".
}
\author{Paul Shannon}

\examples{
  refnet <- RefNet()
  providers(refnet)

}
\seealso{
  RefNet,
  providerClasses,
  interactions,
  addStandardNames,
  IDMapper
}


\keyword{utilities}
