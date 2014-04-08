# go.R
#------------------------------------------------------------------------------------------------------------------------
source("../standardColumns.R")
#------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library (org.Hs.eg.db)

if(!exists("assignGeneIDs"))
    source("~/s/data/public/human/symToGeneID.R")

DATA_ROOT <- "~/s/data/public/human/networks/gersteinEncode"

DESTINATION_DIR <- "../../extdata"
stopifnot(file.exists(DESTINATION_DIR))

.printf <- function(...) print(noquote(sprintf(...)))

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_readFile()
    test_addStandardColumns()
    test_createIdMapping()
    test_addStandardIdentifiers()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
runAll = function (levels)
{
    runTests()
    tbl <- addStandardColumns(readFile())
    lookup <- createIdMapping(tbl)
    tbl.gerstein <- addStandardIdentifiers(tbl, lookup)
    checkEquals(length(setdiff(standardColumns, colnames(tbl.gerstein))), 0)
    destinationFile <- "gerstein-2012.tsv"
    write.table(tbl.gerstein, file=destinationFile, row.names=FALSE, sep="\t", quote=FALSE)
    .printf("wrote tbl.gerstein (%d interactions) as %s", nrow(tbl.gerstein), destinationFile)

} # runAll
#------------------------------------------------------------------------------------------------------------------------
readFile <- function()
{
    filename <- file.path(DATA_ROOT, "Hs_Tr.txt")
    tbl <- read.table(filename, sep="\t", header=FALSE, as.is=TRUE)
    colnames(tbl) <- c("a.common", "b.common")

    invisible(as.data.frame(tbl))

} # readFile
#-------------------------------------------------------------------------------
test_readFile <- function()
{
    print("--- test_readFile")
    tbl <- readFile()
    checkTrue(is(readFile(), "data.frame"))
    checkEquals(dim(tbl), c (6895, 2))
    checkEquals(colnames(tbl), c("a.common", "b.common"))
    
} # test_readFile
#-------------------------------------------------------------------------------
#  "a.canonical",
#  "a.canonicalIdType",
#  "a.cellularComponent",
#  "a.common",
#  "a.modification",
#  "a.organism",
#  "b.canonical",
#  "b.canonicalIdType",
#  "b.cellularComponent",
#  "b.common",
#  "b.modification",
#  "b.organism",
#  "cellType",
#  "comment",
#  "detectionMethod",
#  "pmid",
#  "relation"
#  "reciprocal"
addStandardColumns <- function(tbl)
{
    stopifnot(colnames(tbl) == c("a.common", "b.common"))
    count <- nrow(tbl)

    tbl2 <- cbind(a.canonical=rep(NA, count),
                  b.canonical=rep(NA, count),
                  relation=rep("psi-mi:MI:0407(direct interaction)", count),
                  bidirectional=rep(FALSE, count),
                  detectionMethod=rep(
                     "psi-mi:MI:0402(chromatin immunoprecipitation assay)", count),
                  pmid=rep("22955619", count),
                  a.organism=rep("9606", count),
                  b.organism=rep("9606", count),
                  a.common=tbl$a.common,
                  a.canonicalIdType=rep(NA, count),
                  b.common=tbl$b.common,
                  b.canonicalIdType=rep(NA, count),
                  cellType=rep(NA, count),
                  a.modification=rep(NA, count),
                  a.cellularComponent=rep(NA, count),
                  b.modification=rep(NA, count),
                  b.cellularComponent=rep(NA, count),
                  provider=rep("gerstein.2012", count),
                  comment=rep(NA, count))

    invisible(as.data.frame(tbl2))

} # addStandardColumns
#-------------------------------------------------------------------------------
test_addStandardColumns <- function()
{
    print("--- test_addStandardColumns")
    tbl <- readFile()
    tbl <- addStandardColumns(tbl)
    checkTrue(is(tbl, "data.frame"))
    checkEquals(dim(tbl), c(6895, 19))
    
} # test_addStandardColumns
#-------------------------------------------------------------------------------
createIdMapping <- function(tbl)
{
    id.map <- assignGeneIDs(unique(c(tbl$a.common, tbl$b.common)))
    lookup <- id.map$mapped
    multiples <- id.map$mapped[id.map$multiples]
    multiples.trimmed <- sapply(multiples, "[", 1)
    lookup[names(multiples.trimmed)] <- as.character(multiples.trimmed)
    deleters <- unique(id.map$failures)
    lookup <- lookup[!names(lookup) %in% deleters]

    lookup

} # createIdMapping
#-------------------------------------------------------------------------------
test_createIdMapping <- function()
{
    print("--- test_createIdMapping")
    tbl <- addStandardColumns(readFile())
    lookup <- createIdMapping(tbl)

    expected.size <- 2864  # empirically obtained
    
    checkEquals(length(lookup), expected.size)
    geneIDs <- as.character(lookup)
    geneIDs.duplicated <- which(duplicated(geneIDs))
       # sometimes multiple symbols map to the same geneID
    checkEquals(length(geneIDs.duplicated), 14)

    checkEquals(length(unique(names(lookup))), expected.size)
    checkEquals(length(unique(as.character(lookup))), expected.size-length(geneIDs.duplicated))
    

} # test_createIdMapping
#-------------------------------------------------------------------------------
addStandardIdentifiers <- function(tbl, lookup)
{
    tbl$a.canonical <- as.character(lookup[tbl$a.common])
    tbl$b.canonical <- as.character(lookup[tbl$b.common])
    idType <- rep("entrezGeneID", nrow(tbl))
    tbl$a.canonicalIdType <- idType
    tbl$b.canonicalIdType <- idType

    invisible(tbl[, standardColumns])  # ensure latest ordering of columns is used

} # addStandardIdentifiers
#-------------------------------------------------------------------------------
test_addStandardIdentifiers <- function()
{
    print("--- test_addStandardIdentifiers")
    tbl <- addStandardColumns(readFile())
    lookup <- createIdMapping(tbl)
    tbl <- addStandardIdentifiers(tbl, lookup)

      # check for proper columns, then proper order
    checkEquals(length(setdiff(standardColumns, colnames(tbl))), 0)
    checkEquals(standardColumns, colnames(tbl))

    x <- as.list(tbl[1,])
    checkEquals(x$a.common, "AR")
    checkEquals(x$b.common, "ACPP")
    checkEquals(x$a.organism, "9606")
    checkEquals(x$b.organism, "9606")
    checkEquals(x$relation, "psi-mi:MI:0407(direct interaction)")
    checkEquals(x$detectionMethod, "psi-mi:MI:0402(chromatin immunoprecipitation assay)")
    checkEquals(x$provider, "gerstein.2012")
    checkEquals(x$pmid, "22955619")
    checkEquals(x$a.canonical, "367")
    checkEquals(x$b.canonical, "55")
    checkEquals(x$a.canonicalIdType, "entrezGeneID")
    checkEquals(x$a.common, "AR")
    checkEquals(x$b.common, "ACPP")
    
} # test_addStandardIdentifiers
#-------------------------------------------------------------------------------
