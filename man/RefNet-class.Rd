\name{RefNet-class}
\docType{class}

\alias{class:RefNet}
\alias{RefNet}
\alias{RefNet-class}
\alias{show,RefNet-method}

\title{RefNet}

\description{
   A query interface to a large collecion of molecular interactions.

}


\section{Constructor}{
  \describe{
    \item{}{
      \code{RefNet}: establishes connectin to the central PSICQUIC web server,
      and loads all "native" data sources.
    }
  }
}


\section{Methods}{
  \describe{
    \item{}{
      \code{providers(x)}: lists short names of the data providers
     }

   \item{}{
      \code{providerClasses(x)}: lists the names of logical groupings of
      the providers.
       }

    \item{}{
      \code{interactions(x,id,
         species, speciesExclusive,
          type, provider,detectionMethod,
         publicationID, quiet)}: retrieves all interactions matching the specified pattern.
         }

   \item{}{
      \code{providerClasses(x)}: lists the names of logical groupings of
      the providers.
       }

   \item{}{
      \code{show(x)}: displays current providers and related data
       }
    }
  }


\section{Functions}{
  \describe{
    \item{}{
      \code{detectDuplicateInteractions(tbl)}: adds a column grouping
      all interactions by their two participants.   The 0th dupGroup
      is all singletons.
      }
      \item{}{
	\code{pickBestFromDupGroup(dupGrr, tbl.dups, preferred.interaction.types)}:
           matches, in order, the preferred types against the
           interactions types, returning tbl.dups rowname for best match.
        }
      \item{}{
	\code{pubmedAbstract(pmid, split=FALSE)}:
           returns the text of the pubmed abstract
        }
      }
    }


\author{Paul Shannon}

\seealso{
  PSICQUIC,
  interactions,
  providerClasses,
  providers,
  pubmedAbstract,
  detectDuplicateInteractions,
  pickBestFromDupGroup
}

  
\examples{
    # List the sources
  refNet <- RefNet()
  show(refNet)

}
\keyword{methods}
\keyword{classes}
