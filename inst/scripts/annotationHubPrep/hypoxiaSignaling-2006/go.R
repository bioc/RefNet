# go.R
#------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library (org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
run = function (levels)
{
  if (0 %in% levels) {
      tbl <<- read.table("hypoxiaSignaling-2006-orig.tsv", sep="\t", header=TRUE, as.is=TRUE)
    } # 0

  if (1 %in% levels) {
     a <- tbl$a.canonical
       # allow for, e.g., geneID::post-translational-modification::cellularLocation
       # get just the id
     base.names <- unlist(lapply(strsplit(a, "::"), "[", 1))
     genes <<- mget(base.names, org.Hs.egSYMBOL, ifnotfound=NA)
     NAs <- which(is.na(genes))
     if(length(NAs) > 0)
         genes[NAs] <<- names(genes)[NAs]
   
     commonNames <- as.character(genes)
     ids.with.suffixes <- grep ("::", a)
     suffixes <- unlist(lapply(a[4], function(n) {x <- strsplit(n, "::")[[1]]; tokenCount <- length(x); x[2:tokenCount]}))
     commonNames[ids.with.suffixes] <- paste(commonNames[ids.with.suffixes], suffixes, sep="::")
     tbl$a.common <<- commonNames
     
     entrez.canonicals <- grep("^[0-9]+$", base.names)
     a.canonicalIdType <- rep(NA, nrow(tbl))
     a.canonicalIdType [entrez.canonicals] <- "entrezGeneID"
     tbl$a.canonicalIdType <<- a.canonicalIdType
     tbl$a.organism <<- rep("9606", nrow(tbl))
     } # 1

  if (2 %in% levels) {
     b <- tbl$b.canonical
       # allow for, e.g., geneID::post-translational-modification::cellularLocation
       # get just the id
     base.names <- unlist(lapply(strsplit(b, "::"), "[", 1))
     genes <<- mget(base.names, org.Hs.egSYMBOL, ifnotfound=NA)
     NAs <- which(is.na(genes))
     if(length(NAs) > 0)
         genes[NAs] <<- names(genes)[NAs]
   
     commonNames <- as.character(genes)
     ids.with.suffixes <- grep ("::", b)
     browser()
     suffixes <- unlist(lapply(a[4], function(n) {x <- strsplit(n, "::")[[1]]; tokenCount <- length(x); x[2:tokenCount]}))
     commonNames[ids.with.suffixes] <- paste(commonNames[ids.with.suffixes], suffixes, sep="::")
     tbl$b.common <<- commonNames
     
     entrez.canonicals <- grep("^[0-9]+$", base.names)
     b.canonicalIdType <- rep(NA, nrow(tbl))
     b.canonicalIdType [entrez.canonicals] <- "entrezGeneID"
     tbl$b.canonicalIdType <<- b.canonicalIdType
     tbl$b.organism <<- rep("9606", nrow(tbl))
     } # 2

  if (3 %in% levels) {
     tbl$cellType <- rep(NA, nrow(tbl))
     tbl$a.modification <- rep(NA, nrow(tbl))
     tbl$a.cellularComponent <- rep(NA, nrow(tbl))
     tbl$b.modification <- rep(NA, nrow(tbl))
     tbl$b.cellularComponent <- rep(NA, nrow(tbl))
     preferred.column.order <- c("a.canonical",
                                 "b.canonical",
                                 "relation",
                                 "pmid",
                                 "a.organism",
                                 "b.organism",
                                 "a.common",
                                 "a.canonicalIdType",
                                 "b.common",
                                 "b.canonicalIdType",
                                 "cellType",
                                 "a.modification",
                                 "a.cellularComponent",
                                 "b.modification",
                                 "b.cellularComponent",
                                 "note")
     tbl <<- tbl[, preferred.column.order]
     colnames(tbl)[16] <<- "comment"
     } # 3

  if (4 %in% levels) {
      write.table(tbl, file="hypoxiaSignaling-2006.tsv", sep="\t",
                  row.names=FALSE, quote=FALSE)
                 
    } # 4

  if (5 %in% levels) {
    } # 5

  if (6 %in% levels) {
    } # 6

  if (7 %in% levels) {
    } # 7

  if (8 %in% levels) {
    } # 8

  if (9 %in% levels) {
    } # 9

  if (10 %in% levels) {
    } # 10

  if (11 %in% levels) {
    } # 11

  if (12 %in% levels) {
    } # 12

  if (13 %in% levels) {
    } # 13

  if (14 %in% levels) {
    } # 14

  if (15 %in% levels) {
    } # 15

  if (16 %in% levels) {
    } # 16

  if (17 %in% levels) {
    } # 17

  if (18 %in% levels) {
    } # 18

  if (19 %in% levels) {
    } # 19

  if (20 %in% levels) {
    } # 20


} # run
#------------------------------------------------------------------------------------------------------------------------
