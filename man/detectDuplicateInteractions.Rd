\name{detectDuplicateInteractions}

\alias{detectDuplicateInteractions}

\title{detectDuplicateInteractions}

\description{
  Assign a shared "dupGroup" number to each duplicated interaction,
  here defined as A.canonical, type, B.canonical.   All other
  information about the interaction is ignored.  An A-B interaction
  is treated the same as a B-A interaction.   This assessment
  prepares a possibly large interaction set for scrutiny by eye,
  or programmatically, for filtering, in which a single preferred
  interaction can be selected, out of each dupGroup, for further use.
}

\usage{
  detectDuplicateInteractions(tbl)
}

\arguments{
  \item{tbl}{a data.frame, returned by interactions method}
  }


\value{
  A data.frame identical to the one passed as an argument, but with one
  additional column "dupGroup".  dupGroup 0 (zero) contains all of the
  interactions which have NO duplicates.
}


\author{Paul Shannon}

\examples{
  filename <- system.file(package="RefNet", "extdata", "tbl.28g2.RData")
  load(filename)
  tbl.withDupsColumn <- detectDuplicateInteractions(tbl.28g2)
  table(tbl.withDupsColumn$dupGroup)

}
\seealso{
  RefNet,
  providerClasses,
  interactions,
  addStandardNames,
  pickBestFromDupGroup
}


\keyword{utilities}
