\name{pickBestFromDupGroup}

\alias{pickBestFromDupGroup}

\title{pickBestFromDupGroup}

\description{
  Assign a shared "dupGroup" number to each duplicated interaction,
  here defined as A.canonical, type, B.canonical.   All other
  information about the interaction is ignored.  An A-B interaction
  is treated the same as a B-A interaction.   This assessment
  prepares a possibly large interaction set for scrutiny by eye,
  or programmatically, for filtering, in which a single preferred
  interaction can be selected, out of each dupGroup, for further use.

  Note that preferred.interaction.types elements may be substrings
  of the full (and often unwieldy) interaction type names used
  in the underlying data. 
}

\usage{
  pickBestFromDupGroup(dupGrp, tbl.dups,
                       preferred.interaction.types)
}

\arguments{
  \item{dupGrp}{an integer}
  \item{tbl.dups}{data.frame, returned by detectDuplicateInteractions.}
  \item{preferred.interaction.types}{list of character strings}
}


\value{
  The row name (or names, for dupGroup 0) of the tbl.dups row
  which is the best match to the ordered list of preferred.interaction.types.
}


\author{Paul Shannon}

\examples{
    filename <- system.file(package="RefNet", "extdata", "tbl.dups.RData")
    load(filename)
    preferred.types <- c("direct", "physical", "aggregation")

      # get the best from dupGroup 1
    best.1 <- pickBestFromDupGroup(1, tbl.dups, preferred.types)
    tbl.dups[best.1, c("A.common", "B.common", "type", "provider", "publicationID")]

      # get all of the best.  not every dupGroup will pass muster
    dupGroups <- sort(unique(tbl.dups$dupGroup))
    bestOfDups <- unlist(lapply(dupGroups, function(dupGroup) 
                       pickBestFromDupGroup(dupGroup, tbl.dups, preferred.types)))
    deleters <- which(is.na(bestOfDups))
    if(length(deleters) > 0)
        bestOfDups <- bestOfDups[-deleters]
    length(bestOfDups)
    tbl.dups[bestOfDups, c("A.common", "B.common", "type", "provider", "publicationID")]
}
\seealso{
  RefNet,
  providerClasses,
  interactions,
  addStandardNames,
  detectDuplicateInteractions
}


\keyword{utilities}
