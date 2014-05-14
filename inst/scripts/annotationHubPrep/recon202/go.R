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
    #test_convertToStandardColumns()        # canonical geneIDs added also
    #test_fixNonStandardGeneNames()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
run = function ()
{
   runTests()
   tbl <- readData()
   tbl <- convertToStandardColumns(tbl)
   tbl <- fixNonStandardGeneNames(tbl)
   write.table(tbl, file="stamlabTFs-2012.tsv", row.names=FALSE, sep="\t", quote=FALSE)
    # take a quick look:
    #   t(tbl[1:2,])
    #                      1                                      2                                     
    #  a.canonical         "196"                                  "196"                                 
    #  b.canonical         "79365"                                "4849"                                
    #  relation            "transcription factor binding"         "transcription factor binding"        
    #  bidirectional       "FALSE"                                "FALSE"                               
    #  detectionMethod     "psi-mi:MI:0606(DNase I footprinting)" "psi-mi:MI:0606(DNase I footprinting)"
    #  pmid                "22959076"                             "22959076"                            
    #  a.organism          "9606"                                 "9606"                                
    #  b.organism          "9606"                                 "9606"                                
    #  a.common            "AHR"                                  "AHR"                                 
    #  a.canonicalIdType   "entrezGeneID"                         "entrezGeneID"                        
    #  b.common            "BHLHE41"                              "CNOT3"                               
    #  b.canonicalIdType   "entrezGeneID"                         "entrezGeneID"                        
    #  cellType            "AG10803"                              "AG10803"                             
    #  a.modification      NA                                     NA                                    
    #  a.cellularComponent NA                                     NA                                    
    #  b.modification      NA                                     NA                                    
    #  b.cellularComponent NA                                     NA                                    
    #  provider            "stamlabTF"                            "stamlabTF"                           
    #  comment             NA                                     NA                                    

    # check the line count, upload to amazon s3, which requires credentials, typically
    #   a pair of access keys in ~/.aws/config
    # wc -l stamlabTFs-2012.tsv  #  574123
    # aws s3 ls s3://refnet-networks
    # aws s3 cp stamlabTFs-2012.tsv s3://refnet-networks/ --region us-west-1

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
convertToStandardColumns <- function(tbl)
{
    browser()
    tbl.0 <- convertComplexInteractionsToStandardColumns(subset(tbl, interaction=="complexMember"))
    tbl.1 <- convertModifiesInteractionsToStandardColumns(subset(tbl, interaction=="modifies"))
    tbl.2 <- convertProductOfInteractionsToStandardColumns(subset(tbl, interaction=="productOf"))
    tbl.3 <- convertSubstrateOfInteractionsToStandardColumns(subset(tbl, interaction=="substrateOf"))
    
    a.organism <- rep("9606", nrow(tbl))
    b.organism <- rep("9606", nrow(tbl))
    
    relation <- rep("transcription factor binding", nrow(tbl))
    detectionMethod <- rep("psi-mi:MI:0606(DNase I footprinting)", nrow(tbl))
    pmid <- rep("22959076", nrow(tbl))
    bidirectional <- rep("FALSE", nrow(tbl))
    a.modification <- rep(NA, nrow(tbl))
    b.modification <- rep(NA, nrow(tbl))
    a.cellularComponent <- rep(NA, nrow(tbl))
    b.cellularComponent <- rep(NA, nrow(tbl))
    provider <- rep("stamlabTF", nrow(tbl))
    comment <- rep(NA, nrow(tbl))
    a.canonical <- as.character(mget(tbl$a.common, org.Hs.egSYMBOL2EG, ifnotfound=NA))
    b.canonical <- as.character(mget(tbl$b.common, org.Hs.egSYMBOL2EG, ifnotfound=NA))

    a.canonicalIdType <- rep("entrezGeneID", nrow(tbl))
    b.canonicalIdType <- rep("entrezGeneID", nrow(tbl))

    tbl <- cbind(tbl,
                 a.canonical=a.canonical,
                 b.canonical=b.canonical,
                 relation=relation,
                 bidirectional=bidirectional,
                 detectionMethod=detectionMethod,
                 pmid=pmid,
                 a.organism=a.organism,
                 b.organism=b.organism,
                 #a.common=a.common,
                 a.canonicalIdType=a.canonicalIdType,
                 #b.common=b.common,
                 b.canonicalIdType=b.canonicalIdType,
                 #cellType=cellType,
                 a.modification=a.modification,
                 a.cellularComponent=a.cellularComponent,
                 b.modification=b.modification,
                 b.cellularComponent=b.cellularComponent,
                 provider=provider,
                 comment=comment)

     tbl[, standardColumns]

} # convertToStandardColumns
#-------------------------------------------------------------------------------
test_convertToStandardColumns <- function()
{
    print("--- test_convertToStandardColumns")
    tbl <- readData()
    tbl <- convertToStandardColumns(tbl)

    checkEquals(colnames(tbl), standardColumns);
#    checkTrue(nrow(tbl) == 27277)
#
#    checkTrue(all(tbl$a.canonicalIdType %in% recognized.canonical.idTypes))
#    checkTrue(all(tbl$species=="9606"))
#    checkTrue(all(tbl$provider=="stamlabTF"))
#    checkTrue(all(tbl$publicationID=="22959076"))
#    checkTrue(all(tbl$detectionMethod=="psi-mi:MI:0606(DNase I footprinting)"))
#
#    checkEquals(length(which(tbl$a.canonical=="NA")), 222)
#    checkEquals(length(which(tbl$b.canonical=="NA")), 131)
    
} # test_convertToStandardColumns
#----------------------------------------------------------------------------------------------------
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

}
#----------------------------------------------------------------------------------------------------
convertProductOfInteractionsToStandardColumns <- function(tbl)
{

}
#----------------------------------------------------------------------------------------------------
convertSubstrateOfInteractionsToStandardColumns <- function(tbl)
{

}
#----------------------------------------------------------------------------------------------------
test_convertModifiesInteractionsToStandardColumns <- function()
{
    print("--- test_convertModifiesInteractionsToStandardColumns")

} # test_convertModifiesInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
test_convertProductOfInteractionsToStandardColumns <- function()
{
    print("--- test_convertProductOfInteractionsToStandardColumns")

} # test_convertProductOfInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
test_convertSubstrateOfInteractionsToStandardColumns <- function()
{
    print("--- test_convertSubstrateOfInteractionsToStandardColumns")

} # test_convertSubstrateOfInteractionsToStandardColumns
#----------------------------------------------------------------------------------------------------
# four stamlab gene symbols are not (or are no longer) HUGO standard.
#  HTLF="3344",    # FOXN2
#  ZFP161="7541",  # ZBTB14
#  ZNF238="10472", # ZBTB18
#  INSAF="3637"   # record withdrawn by ncbi, 2008
fixNonStandardGeneNames <- function(tbl)
{
    
   htlf.a.hits <- grep("HTLF", tbl$a.common)
   htlf.b.hits <- grep("HTLF", tbl$b.common)

   if(length(htlf.a.hits) > 0){
       tbl$a.common[htlf.a.hits] <- "FOXN2"
       tbl$a.canonical[htlf.a.hits] <- "3344"
       }

   if(length(htlf.b.hits) > 0){
       tbl$b.common[htlf.b.hits] <- "FOXN2"
       tbl$b.canonical[htlf.b.hits] <- "3344"
       }

   zfp161.a.hits <- grep("ZFP161", tbl$a.common)
   zfp161.b.hits <- grep("ZFP161", tbl$b.common)

   if(length(zfp161.a.hits) > 0){
       tbl$a.common[zfp161.a.hits] <- "ZBTB14"
       tbl$a.canonical[zfp161.a.hits] <- "7541"
       }

   if(length(zfp161.b.hits) > 0){
       tbl$b.common[zfp161.b.hits] <- "ZBTB14"
       tbl$b.canonical[zfp161.b.hits] <- "7541"
       }

     #  ZNF238="10472", # ZBTB18

   znf238.a.hits <- grep("ZNF238", tbl$a.common)
   znf238.b.hits <- grep("ZNF238", tbl$b.common)

   if(length(znf238.a.hits) > 0){
       tbl$a.common[znf238.a.hits] <- "ZBTB18"
       tbl$a.canonical[znf238.a.hits] <- "10472"
       }

   if(length(znf238.b.hits) > 0){
       tbl$b.common[znf238.b.hits] <- "ZBTB18"
       tbl$b.canonical[znf238.b.hits] <- "10472"
       }

     #  INSAF="3637" record withdrawn by ncbi, 2008
   insaf.a.hits <- grep("INSAF", tbl$a.common)
   insaf.b.hits <- grep("INSAF", tbl$b.common)

   if(length(insaf.a.hits) > 0){
       tbl$a.canonical[insaf.a.hits] <- "3637"
       }

   if(length(insaf.b.hits) > 0){
       tbl$b.canonical[insaf.b.hits] <- "3637"
       }

   invisible(tbl)

} # fixNonStandardGeneNames
#--------------------------------------------------------------------------------
test_fixNonStandardGeneNames <- function(tbl)
{
   print("--- test_fixNonStandardGeneNames")

   tbl <- readData(max.directories.for.testing=2)
   tbl <- addStandardColumns(tbl)
   tbl <- fixNonStandardGeneNames(tbl)

   checkEquals(length(which(tbl$a.canonical=="NA")), 0)
   checkEquals(length(which(tbl$b.canonical=="NA")), 0)

      # make sure that mapping gene symbols in those 2 directories
      # catches all nonstandard gene symbols in all directories
   
   tbl <- readData()
   tbl <- addStandardColumns(tbl)
   tbl <- fixNonStandardGeneNames(tbl)
   checkEquals(length(which(tbl$a.canonical=="NA")), 0)
   checkEquals(length(which(tbl$b.canonical=="NA")), 0)
   checkEquals(length(which(tbl$a.canonical=="NA")), 0)
   checkEquals(length(which(tbl$b.canonical=="NA")), 0)

} # test_fixNonStandardGeneNames
#--------------------------------------------------------------------------------
