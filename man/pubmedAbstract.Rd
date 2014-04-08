\name{pubmedAbstract}

\alias{pubmedAbstract}

\title{pubmedAbstract}

\description{

  Get the text of a journal article's abstract from PubMed.
}

\usage{
  pubmedAbstract(pmid, split=TRUE)
}

\arguments{
  \item{pmid}{A character string, the pubmedID.}
  \item{split}{Whether or not to split the text on embedded newlines.}
}

\value{
  A single character string, or a vector of strings split on
  the embedded newline characters.
}


\author{Paul Shannon}

\examples{
  pubmedAbstract("22959076")

}
\seealso{
  RefNet,
  providerClasses,
  interactions,
  addStandardNames
}


\keyword{utilities}
