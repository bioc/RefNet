# go.R
#------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library (org.Hs.eg.db)
source("../standardColumns.R")

if(!exists("assignGeneIDs"))
    source("~/s/data/public/human/symToGeneID.R")

DATA_DIR <- "."

DESTINATION_DIR <- "."
stopifnot(file.exists(DESTINATION_DIR))

.printf <- function(...) print(noquote(sprintf(...)))

#total.count <- 574122    # number of interactions
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_readData()
    test_convertComplexInteractionsToStandardColumns()
    test_convertModifiesInteractionsToStandardColumns()
    test_convertProductOfInteractionsToStandardColumns()
    test_convertSubstrateOfInteractionsToStandardColumns()
    test_convertAll()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
run = function ()
{
   tbl <- readData()
   tbl.out <- convertAll(tbl)
   write.table(tbl.out, file="recon202.tsv", row.names=FALSE, sep="\t", quote=FALSE)

    # check the line count, upload to amazon s3, which requires credentials, typically
    #   a pair of access keys in ~/.aws/config
    # wc -l recon202.tsv  #  48938
    # aws s3 ls s3://refnet-networks
    # aws s3 cp recon202.tsv s3://refnet-networks/ --region us-west-1


   browser()
   x <- 99

} # run
#------------------------------------------------------------------------------------------------------------------------
readData <- function(quiet=TRUE)
{
    filename <- file.path(DATA_DIR, "tbl.recon202.RData")
    stopifnot(file.exists(filename))
    table.name <- load(filename)

    if(!quiet)
        printf("loaded %s, %d rows, %d cols", table.name, nrow(tbl.recon202), ncol(tbl.recon202))

    invisible(tbl.recon202)

} # readData
#-------------------------------------------------------------------------------
test_readData <- function()
{
    print("--- test_readData")
    tbl <- readData(quiet=TRUE)

    checkEquals(dim(tbl), c (48937, 18))
    checkEquals(colnames(tbl),
                c("a.orig",             "interaction",        "b.orig",            
                  "a.name",             "b.name",             "a.type",            
                  "a.geneID",           "a.uniprot",          "a.sboTerm",         
                  "a.chebi",            "a.kegg.compound",    "a.kegg.genes",      
                  "a.kegg.drug",        "a.hmdb",             "a.pubchem.substance",
                  "reversible",         "pubmed",             "compartment"))        
    
} # test_readData
#-------------------------------------------------------------------------------
convertComplexInteractionsToStandardColumns <- function(tbl)
{
   tokens <- strsplit(tbl$a.orig, "_")
   a.canonical <- unlist(lapply(tokens, "[", 2))
   a.canonicalIdType <- rep("entrezGeneID", nrow(tbl))
   
   a.common <- unlist(mget(a.canonical, org.Hs.egSYMBOL), use.names=FALSE)
   relation <- rep("complexMember", nrow(tbl))
   a.cellularComponent <- unlist(lapply(tokens, "[", 4))

   tokens <- strsplit(tbl$b.orig, "_")
   as.geneIDs <- lapply(tokens, function(toks) toks[seq(from=2, to=length(toks), by=3)])
   as.colonSeparatedGeneIDs <- lapply(as.geneIDs, function(geneIDs) paste(geneIDs, collapse=":"))
   b.canonical <- unlist(as.colonSeparatedGeneIDs)
   b.canonicalIdType <- rep("entrezGeneID", nrow(tbl))

   as.symbols <- lapply(as.geneIDs,
                       function(cplxMembers) as.character(mget(cplxMembers, org.Hs.egSYMBOL)))
   as.colonSeparatedSymobls <- lapply(as.symbols, function(symbols) paste(symbols, collapse=":"))

   b.common <- unlist(as.colonSeparatedSymobls)
   
   compartment <- unlist(lapply(tokens,
                                function(toks) unique(toks[seq(from=4, to=length(toks), by=3)])))
   b.cellularComponent <- compartment

   a.organism <- "9606"
   b.organism <- "9606"

   bidirectional <- rep(FALSE, nrow(tbl))
   detectionMethod <- rep(NA, nrow(tbl))
   provider <- rep("recon2", nrow(tbl))
   comment <- rep(NA, nrow(tbl))
   
   cellType <- rep(NA, nrow(tbl))
   detectionMethod <- rep(NA, nrow(tbl))
   a.modification <- rep(NA, nrow(tbl))
   b.modification <- rep(NA, nrow(tbl))
   pmid <- tbl$pubmed

   return(data.frame(a.canonical=a.canonical,
                     b.canonical=b.canonical,
                     relation=relation,
                     bidirectional=bidirectional,
                     detectionMethod=detectionMethod,
                     pmid=pmid,
                     a.organism=a.organism,
                     b.organism=b.organism,
                     a.common=a.common,
                     a.canonicalIdType=a.canonicalIdType,
                     b.common=b.common,
                     b.canonicalIdType=b.canonicalIdType,
                     cellType=cellType,
                     a.modification=a.modification,
                     a.cellularComponent=a.cellularComponent,
                     b.modification=b.modification,
                     b.cellularComponent=b.cellularComponent,
                     provider=provider,
                     comment=comment,
                     stringsAsFactors=FALSE))
 
} # convertComplexInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
test_convertComplexInteractionsToStandardColumns <- function()
{
    print("--- test_convertComplexInteractionsToStandardColumns")
    tbl <- readData()

       # test just 3 complexMember rows to start
    case.1 <- subset(tbl, interaction=="complexMember")[1:3,]
    tbl.1 <- convertComplexInteractionsToStandardColumns(case.1)
    checkEquals(dim(tbl.1), c(3, 19))
    checkEquals(colnames(tbl.1), standardColumns)
    checkEquals(unique(tbl.1$relation), "complexMember")
    checkEquals(unique(tbl.1$provider), "recon2")
    checkEquals(unique(tbl.1$bidirectional), FALSE)
    checkEquals(tbl.1$a.common, c("PIGK", "PIGT", "PIGU"))
    checkEquals(tbl.1$a.canonical, c("10026", "51604", "128869"))
       # all 3 are the same complex
    checkEquals(unique(tbl.1$b.canonical), "10026:128869:51604:8733:94005")
    checkEquals(unique(tbl.1$b.common), "PIGK:PIGU:PIGT:GPAA1:PIGS")
    
      # now run the whole table
    case.2 <- subset(tbl, interaction=="complexMember")
    tbl.2 <- convertComplexInteractionsToStandardColumns(case.2)
    checkEquals(dim(tbl.2), c(8116, 19))
    checkEquals(unique(tbl.2$relation), "complexMember")
    checkEquals(unique(tbl.2$provider), "recon2")
    checkEquals(unique(tbl.2$bidirectional), FALSE)
    checkEquals(colnames(tbl.2), standardColumns)


} # test_convertComplexInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
convertModifiesInteractionsToStandardColumns <- function(tbl)
{
   tokens <- strsplit(tbl$a.orig, "_")
   a.canonical <- unlist(lapply(tokens, "[", 2))
   
   a.common <- unlist(mget(a.canonical, org.Hs.egSYMBOL, ifnotfound=NA), use.names=FALSE)
   failures <- which(is.na(a.common))
       # may as well use the unmapped canonical name for common names as well
   a.common[failures] <- a.canonical[failures]
   a.canonicalIdType <- rep("entrezGeneID", nrow(tbl))
   a.canonicalIdType[failures] <- NA
   relation <- rep("modifies", nrow(tbl))
   a.cellularComponent <- unlist(lapply(tokens, "[", 4))

   b.canonical <- tbl$b.orig
   b.canonicalIdType <- rep("recon2ReactionID", nrow(tbl))
   b.common <- tbl$b.name
   
   b.cellularComponent <- rep(NA, nrow(tbl))

   a.organism <- "9606"
   b.organism <- "9606"

   bidirectional <- rep(FALSE, nrow(tbl))
   detectionMethod <- rep(NA, nrow(tbl))
   provider <- rep("recon2", nrow(tbl))
   comment <- rep(NA, nrow(tbl))
   
   cellType <- rep(NA, nrow(tbl))
   detectionMethod <- rep(NA, nrow(tbl))
   a.modification <- rep(NA, nrow(tbl))
   b.modification <- rep(NA, nrow(tbl))
   pmid <- tbl$pubmed

   return(data.frame(a.canonical=a.canonical,
                     b.canonical=b.canonical,
                     relation=relation,
                     bidirectional=bidirectional,
                     detectionMethod=detectionMethod,
                     pmid=pmid,
                     a.organism=a.organism,
                     b.organism=b.organism,
                     a.common=a.common,
                     a.canonicalIdType=a.canonicalIdType,
                     b.common=b.common,
                     b.canonicalIdType=b.canonicalIdType,
                     cellType=cellType,
                     a.modification=a.modification,
                     a.cellularComponent=a.cellularComponent,
                     b.modification=b.modification,
                     b.cellularComponent=b.cellularComponent,
                     provider=provider,
                     comment=comment,
                     stringsAsFactors=FALSE))


} # convertModifiesInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
test_convertModifiesInteractionsToStandardColumns <- function()
{
    print("--- test_convertModifiesInteractionsToStandardColumns")
    tbl <- readData()

       # test just 3 complexMember rows to start
    case.1 <- subset(tbl, interaction=="modifies")[1:3,]
    tbl.1  <- convertModifiesInteractionsToStandardColumns(case.1)
    checkEquals(dim(tbl.1), c(3,19))
    checkEquals(colnames(tbl.1), standardColumns)
    checkEquals(unique(tbl.1$a.canonicalIdType), "entrezGeneID")
    checkEquals(unique(tbl.1$b.canonicalIdType), "recon2ReactionID")
    checkEquals(tbl.1$a.common, c("AOC2", "AOC1", "AOC3"))
    checkEquals(unique(tbl.1$b.common), "1,3-Diaminopropane:oxygen oxidoreductase (deaminating)")
    checkEquals(unique(tbl.1$b.canonical), "R_13DAMPPOX")    

      # test the whole table
    case.2 <- subset(tbl, interaction=="modifies")
    tbl.2  <- convertModifiesInteractionsToStandardColumns(case.2)
    checkEquals(dim(tbl.2), c(nrow(case.2), length(standardColumns)))
    checkEquals(as.list(table(tbl.2$a.canonicalIdType)),
                list(entrezGeneID=9728))
       # nrow (case.2) - 9728 which have been successfully mapped
    checkEquals(length(which(is.na(tbl.2$a.canonicalIdType))), 254)

} # test_convertModifiesInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
convertProductOfInteractionsToStandardColumns <- function(tbl)
{
   tokens <- strsplit(tbl$a.orig, "_")
   a.canonical <- unlist(lapply(tokens, "[", 2))
   
   a.common <- tbl$a.name
   a.canonicalIdType <- rep("reconMetabolite", nrow(tbl))
   relation <- rep("productOf", nrow(tbl))
   a.cellularComponent <- unlist(lapply(tokens, "[", 3))

   b.canonical <- tbl$b.orig
   b.canonicalIdType <- rep("recon2ReactionID", nrow(tbl))
   b.common <- tbl$b.name
   
   b.cellularComponent <- rep(NA, nrow(tbl))

   a.organism <- "9606"
   b.organism <- "9606"

   bidirectional <- tbl$reversible
   detectionMethod <- rep(NA, nrow(tbl))
   provider <- rep("recon2", nrow(tbl))
   comment <- rep(NA, nrow(tbl))
   
   cellType <- rep(NA, nrow(tbl))
   detectionMethod <- rep(NA, nrow(tbl))
   a.modification <- rep(NA, nrow(tbl))
   b.modification <- rep(NA, nrow(tbl))
   pmid <- tbl$pubmed

   return(data.frame(a.canonical=a.canonical,
                     b.canonical=b.canonical,
                     relation=relation,
                     bidirectional=bidirectional,
                     detectionMethod=detectionMethod,
                     pmid=pmid,
                     a.organism=a.organism,
                     b.organism=b.organism,
                     a.common=a.common,
                     a.canonicalIdType=a.canonicalIdType,
                     b.common=b.common,
                     b.canonicalIdType=b.canonicalIdType,
                     cellType=cellType,
                     a.modification=a.modification,
                     a.cellularComponent=a.cellularComponent,
                     b.modification=b.modification,
                     b.cellularComponent=b.cellularComponent,
                     provider=provider,
                     comment=comment,
                     stringsAsFactors=FALSE))

} # convertProductOfInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
test_convertProductOfInteractionsToStandardColumns <- function()
{
    print("--- test_convertProductOfInteractionsToStandardColumns")
    tbl <- readData()

       # test just 3 productOf rows to start
    case.1 <- subset(tbl, interaction=="productOf")[1:3,]
    tbl.1  <- convertProductOfInteractionsToStandardColumns(case.1)

    checkEquals(dim(tbl.1), c(3,19))
    checkEquals(colnames(tbl.1), standardColumns)
    checkEquals(unique(tbl.1$a.canonicalIdType), "reconMetabolite")
    checkEquals(unique(tbl.1$b.canonicalIdType), "recon2ReactionID")
    checkEquals(tbl.1$a.common,
                c("10-formyltetrahydrofolate-[Glu](5)",
                  "10-formyltetrahydrofolate-[Glu](5)",
                  "10-formyltetrahydrofolate-[Glu](6)"))
    checkEquals(tbl.1$b.common, c("5-glutamyl-10FTHF transport, lysosomal",
                                  "5-glutamyl-10FTHF transport, mitochondrial",
                                  "6-glutamyl-10FTHF transport, lysosomal"))
    checkEquals(tbl.1$b.canonical, c("R_10FTHF5GLUtl", "R_10FTHF5GLUtm", "R_10FTHF6GLUtl"))

      # test the whole table
    case.2 <- subset(tbl, interaction=="productOf")
    tbl.2  <- convertProductOfInteractionsToStandardColumns(case.2)
    checkEquals(dim(tbl.2), c(nrow(case.2), length(standardColumns)))
    checkEquals(as.list(table(tbl.2$a.canonicalIdType)), list(reconMetabolite=15863))
    checkEquals(as.list(table(tbl.2$b.canonicalIdType)), list(recon2ReactionID=15863))


} # test_convertProductOfInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
convertSubstrateOfInteractionsToStandardColumns <- function(tbl)
{
   tokens <- strsplit(tbl$a.orig, "_")
   a.canonical <- unlist(lapply(tokens, "[", 2))
   
   a.common <- tbl$a.name
   a.canonicalIdType <- rep("reconMetabolite", nrow(tbl))
   relation <- rep("substrateOf", nrow(tbl))
   a.cellularComponent <- unlist(lapply(tokens, "[", 3))

   b.canonical <- tbl$b.orig
   b.canonicalIdType <- rep("recon2ReactionID", nrow(tbl))
   b.common <- tbl$b.name
   
   b.cellularComponent <- rep(NA, nrow(tbl))

   a.organism <- "9606"
   b.organism <- "9606"

   bidirectional <- tbl$reversible
   detectionMethod <- rep(NA, nrow(tbl))
   provider <- rep("recon2", nrow(tbl))
   comment <- rep(NA, nrow(tbl))
   
   cellType <- rep(NA, nrow(tbl))
   a.modification <- rep(NA, nrow(tbl))
   b.modification <- rep(NA, nrow(tbl))
   pmid <- tbl$pubmed

   return(data.frame(a.canonical=a.canonical,
                     b.canonical=b.canonical,
                     relation=relation,
                     bidirectional=bidirectional,
                     detectionMethod=detectionMethod,
                     pmid=pmid,
                     a.organism=a.organism,
                     b.organism=b.organism,
                     a.common=a.common,
                     a.canonicalIdType=a.canonicalIdType,
                     b.common=b.common,
                     b.canonicalIdType=b.canonicalIdType,
                     cellType=cellType,
                     a.modification=a.modification,
                     a.cellularComponent=a.cellularComponent,
                     b.modification=b.modification,
                     b.cellularComponent=b.cellularComponent,
                     provider=provider,
                     comment=comment,
                     stringsAsFactors=FALSE))

} # convertSubstrateOfInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
test_convertSubstrateOfInteractionsToStandardColumns <- function()
{
    print("--- test_convertSubstrateOfInteractionsToStandardColumns")
    tbl <- readData()

       # test just 3 productOf rows to start
    case.1 <- subset(tbl, interaction=="substrateOf")[1:3,]
    tbl.1  <- convertSubstrateOfInteractionsToStandardColumns(case.1)

    checkEquals(dim(tbl.1), c(3,19))
    checkEquals(colnames(tbl.1), standardColumns)
    checkEquals(unique(tbl.1$a.canonicalIdType), "reconMetabolite")
    checkEquals(unique(tbl.1$b.canonicalIdType), "recon2ReactionID")
    checkEquals(tbl.1$a.common,
                c("10-formyltetrahydrofolate-[Glu](5)",
                  "10-formyltetrahydrofolate-[Glu](5)",
                  "10-formyltetrahydrofolate-[Glu](6)"))

    checkEquals(tbl.1$b.common, c("5-glutamyl-10FTHF transport, lysosomal",
                                  "5-glutamyl-10FTHF transport, mitochondrial",
                                  "6-glutamyl-10FTHF transport, lysosomal"))

    checkEquals(tbl.1$b.canonical, c("R_10FTHF5GLUtl", "R_10FTHF5GLUtm","R_10FTHF6GLUtl"))

      # test the whole table
    case.2 <- subset(tbl, interaction=="substrateOf")
    tbl.2  <- convertSubstrateOfInteractionsToStandardColumns(case.2)
    checkEquals(dim(tbl.2), c(nrow(case.2), length(standardColumns)))
    checkEquals(as.list(table(tbl.2$a.canonicalIdType)), list(reconMetabolite=14976))
    checkEquals(as.list(table(tbl.2$b.canonicalIdType)), list(recon2ReactionID=14976))

} # test_convertSubstrateOfInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
convertAll <- function(tbl)
{
    tbl.1 <- subset(tbl, interaction=="substrateOf")
    tbl.2 <- subset(tbl, interaction=="productOf")
    tbl.3 <- subset(tbl, interaction=="modifies")
    tbl.4 <- subset(tbl, interaction=="complexMember")

    tbl.1c <- convertSubstrateOfInteractionsToStandardColumns(tbl.1)
    tbl.2c <- convertProductOfInteractionsToStandardColumns(tbl.2)
    tbl.3c <- convertModifiesInteractionsToStandardColumns(tbl.3)
    tbl.4c <- convertComplexInteractionsToStandardColumns(tbl.4)

    invisible(rbind(tbl.1c, tbl.2c, tbl.3c, tbl.4c))
    
} # convertAll
#----------------------------------------------------------------------------------------------------
test_convertAll <- function()
{
    print("--- test_convertAll")
    tbl <- readData()

       # test just 3 productOf rows to start
    tbl.1 <- subset(tbl, interaction=="substrateOf")[1:3,]
    tbl.2 <- subset(tbl, interaction=="productOfOf")[1:3,]
    tbl.3 <- subset(tbl, interaction=="modifies")[1:3,]
    tbl.4 <- subset(tbl, interaction=="complexMember")[1:3,]

    tbl.small <- rbind(tbl.1, tbl.2, tbl.3, tbl.4)

    tbl.out <- convertAll(tbl.small)
    checkEquals(dim(tbl.out), c(nrow(tbl.small), length(standardColumns)))
    checkEquals(as.list(table(tbl.out$relation)),
                list(complexMember=3,
                     modifies=3,
                     productOf=3,
                     substrateOf=3))

} # test_convertAll
#----------------------------------------------------------------------------------------------------
