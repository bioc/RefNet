\name{interactions}
\alias{interactions}
\alias{interactions,RefNet-method}
\alias{\S4method{interactions}{RefNet}}

\title{interactions}

\description{
  This is the primary interface to the RefNet services, in
  which approximately thirty interaction providers return annotated
  molecular (often protein-protein) interactions in response to
  a standard query language ("MIQL" -- the molecular interaction query
  language).  This method is thus an R interface to MIQL, with several
  conveniences added.  
  
  In the most typical and basic use, one specifies a species and a gene symbol,
  and RefNet returns a data.frame listing all of the interactions
  for that gene (or the protein it encodes) from all RefNet
  providers, of all interaction types and all detection methods.

  More focused queries are easily accomplished: all arguments (execpt
  for the RefNet object the method dispatches on) are optional, and each
  specified argument limits the search space: to one or more genes, a publication
  identifier (typically a pubmed id), interaction types, a detection
  methods, and/or one or more species.

  In principle you could obtain all interactions, of all types, from all
  providers, by using default values for all of the search-constraining
  arugments.  But this would be very unwise!  Highly specific or
  exploratory queries of modest scope are   encouraged, as are
  bulk downloads using \emph{ftp} directed at the providers' web sites
  for acquiring large interaction datasets.
  
}
\usage{
\S4method{interactions}{RefNet}(object,
                                  id=NA,
                                  species=NA,
                                  speciesExclusive=TRUE,
                                  type=NA,
                                  provider=NA,
                                  detectionMethod=NA,
                                  publicationID=NA,
                                  quiet=TRUE)
}

\arguments{
  \item{object}{
    a \code{RefNet} object.
  }

  \item{id}{

    one or more \code{character} strings.  These are most often gene
    symbols (e.g., "RAD17"), to which most of the RefNet providers
    respond, translating these into the native identifiers of their data
    store, which is typically a protein identifier of one sort or
    another.  You may also provide such a protein identifier directly
    (e.g., uniprot "O75943", or refseq "NP_002684").  Protein
    identifiers typically return fewer interactions, presumably because
    providers do not translate them into their native identifiers, and
    only a subset of the providers find matches.

    If multiple ids are supplied, then only interactions which occur
    between any two of the ids is returned.   If only one id is
    supplied, all interactions are reported in which that id participates.
  }
  
  
   \item{species}{
     One or more \code{character} strings.  Use an NCBI taxon code, for
     which the \code{speciesIds} method is useful.  Note the next argument.
   }

   \item{speciesExclusive}{
     a \code{logical} variable, default TRUE. If
     set to FALSE than interactions will be included with proteins from species other than
     those explicitly named. Infection and transgenic experiments are
     examples of cross-species protein-protein interactions.  We
     anticipate that the most common query seeks interactions among
     proteins in the same cell of the same organism.  For this reason, 
     this argument defaults to TRUE.
     }
     
   \item{type}{
     a \code{list} of \code{character} strings.  See the method
     \code{interactionTypes}.
   }
   
   \item{provider}{
     a \code{list} of \code{character} strings.  See the
     method \code{providers} for the currently available set.
     This argument defaults to NA, and all currently live
     providers are queried.
   }

   \item{detectionMethod}{a \code{list} of \code{character} strings.  See
     the method \code{detectionMethods}
   }

   \item{publicationID}{a \code{list} of \code{character} strings.
     Usually a pubmedID, but sometimes an OMIM reference number.
   }

   \item{quiet}{
     a \code{logical} of \code{character} strings, default
     TRUE.  Set to FALSE if you wish to trace the exection of your query
     against all providers, and explore the details of the REST API.
   }
}

\value{
  A data.frame with 16 columns, and one row for every interaction,
  described as an annotated A/B relationship, with self-describing
  column names.
}

\author{Paul Shannon}

\examples{
  refnet <- RefNet()
  providers(refnet)
      # query all providers for all known Myc interactions
  tbl.1 <- interactions(refnet, id="Myc", species="9606",
                        provider="gerstein-2012")

      # all Myc/EP400 interactions known to BioGrid
      # make sure that BioGrid is currently available

  if("BioGrid" \%in\% providers(refnet))
     tbl.2 <- interactions(refnet, c("Myc", "EP400"), species="9606",
                           provider="BioGrid")

      # or those between Myc and any other molecule, detected by
      # the "pull down" proteomics technique and judged to be
      # a "direct interaction"

  # recon2 temporarily unavailable.  will be added to AnnotationHub soon
  # tbl.3 <- interactions(refnet, "L-alanine transaminase",  provider="recon2")
  # tbl.3[, c("A.common", "type", "B.common", "A.canonical", "B.canonical")]


}
\seealso{
  providers,
  addStandardNames,
  IDMapper
}



\keyword{utilities}
