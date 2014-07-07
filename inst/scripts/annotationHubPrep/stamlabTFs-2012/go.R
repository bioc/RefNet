# go.R
#------------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library (org.Hs.eg.db)
source("../standardColumns.R")

if(!exists("assignGeneIDs"))
    source("~/s/data/public/human/symToGeneID.R")

DATA_ROOT <- "~/s/data/public/human/networks/stamlabTF"

DESTINATION_DIR <- "."
stopifnot(file.exists(DESTINATION_DIR))

.printf <- function(...) print(noquote(sprintf(...)))

total.count <- 574122    # number of interactions
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_readData()
    test_addStandardColumns()        # canonical geneIDs added also
    test_fixNonStandardGeneNames()
    
} # runTests
#------------------------------------------------------------------------------------------------------------------------
run = function ()
{
   runTests()
   tbl <- readData()
   tbl <- addStandardColumns(tbl)
   tbl <- fixNonStandardGeneNames(tbl)
   write.table(tbl, file="stamlabTFs-2012.tsv", row.names=FALSE, sep="\t", quote=FALSE)
    # take a quick look:
    #   t(tbl[1:2,])
    #                      1                                      2                                     
    #  a.id         "196"                                  "196"                                 
    #  b.id         "79365"                                "4849"                                
    #  relation            "transcription factor binding"         "transcription factor binding"        
    #  bidirectional       "FALSE"                                "FALSE"                               
    #  detectionMethod     "psi-mi:MI:0606(DNase I footprinting)" "psi-mi:MI:0606(DNase I footprinting)"
    #  pmid                "22959076"                             "22959076"                            
    #  a.organism          "9606"                                 "9606"                                
    #  b.organism          "9606"                                 "9606"                                
    #  a.name            "AHR"                                  "AHR"                                 
    #  a.idType   "entrezGeneID"                         "entrezGeneID"                        
    #  b.name            "BHLHE41"                              "CNOT3"                               
    #  b.idType   "entrezGeneID"                         "entrezGeneID"                        
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
readData <- function(max.directories.for.testing=NA, quiet=TRUE)
{
    cellLineDirectories <- list.files(DATA_ROOT)
    if(!is.na(max.directories.for.testing))
        max <- max.directories.for.testing
    else
        max <- length(cellLineDirectories)

    mtx <- matrix(data=NA_character_, nrow=total.count, ncol=3)
    next.free.row <- 1
    
    for(dir in cellLineDirectories[1:max]){
        file.name <- file.path(DATA_ROOT, dir, "genes-regulate-genes.txt")
        if(!file.exists(file.name))
            next
        x <- scan(file.name, what=character(), quiet=quiet)
        m.tmp <- matrix(x, ncol=2, byrow=TRUE)
        last.row <- next.free.row + (length(x)/2) -1
        mtx [next.free.row:last.row, 1:2] <- m.tmp
        mtx [next.free.row:last.row, 3] <- rep(dir, nrow(m.tmp))
        next.free.row <- last.row + 1
        if(!quiet)
           .printf("new cellLine: %s, rowCount %d/%d", dir, nrow(m.tmp), nrow(mtx))
        }# for
    
     tbl <- data.frame(mtx[1:last.row,])
     colnames(tbl) <- c("a.name", "b.name", "cellType")

     invisible(tbl)

} # readData
#-------------------------------------------------------------------------------
test_readData <- function()
{
    print("--- test_readData")
    tbl <- readData(max.directories.for.testing=1)
          # 45 directories, one per celltype, check just one to start
    checkEquals(dim(tbl), c (12482, 3))
    checkEquals(as.list(table(tbl$cellType))$AG10803, 12482)
    checkEquals(colnames(tbl), c("a.name","b.name", "cellType"))
                
    tbl <- readData(max.directories.for.testing=2)
    checkEquals(dim(tbl), c(27277, 3))
    checkEquals(as.list(table(tbl$cellType))$AG10803, 12482)
    checkEquals(colnames(tbl), c("a.name","b.name", "cellType"))
                
    tbl <- readData()
    checkEquals(dim(tbl), c(574122, 3))
    checkEquals(length(which(is.na(tbl$a.name))), 0)
    checkEquals(length(which(is.na(tbl$b.name))), 0)
    checkEquals(length(which(is.na(tbl$cellType))), 0)
    
} # test_readData
#-------------------------------------------------------------------------------
addStandardColumns <- function(tbl)
{
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
    a.name <- tbl$a.name
    b.name <- tbl$b.name
    a.id <- as.character(mget(a.name, org.Hs.egSYMBOL2EG, ifnotfound=NA))
    b.id <- as.character(mget(b.name, org.Hs.egSYMBOL2EG, ifnotfound=NA))

    cellType <- tbl$cellType

    a.idType <- rep("entrezGeneID", nrow(tbl))
    b.idType <- rep("entrezGeneID", nrow(tbl))

    tbl <- cbind(tbl,
                 a.id=a.id,
                 b.id=b.id,
                 relation=relation,
                 bidirectional=bidirectional,
                 detectionMethod=detectionMethod,
                 pmid=pmid,
                 a.organism=a.organism,
                 b.organism=b.organism,
                 a.name=a.name,
                 a.idType=a.idType,
                 b.name=b.name,
                 b.idType=b.idType,
                 cellType=cellType,
                 a.modification=a.modification,
                 a.cellularComponent=a.cellularComponent,
                 b.modification=b.modification,
                 b.cellularComponent=b.cellularComponent,
                 provider=provider,
                 comment=comment)

     tbl[, standardColumns]

} # addStandardColumns
#-------------------------------------------------------------------------------
test_addStandardColumns <- function()
{
    print("--- test_addStandardColumns")
    tbl <- readData(max.directories.for.testing=2)
    tbl <- addStandardColumns(tbl)

    checkEquals(colnames(tbl), standardColumns);
    checkTrue(nrow(tbl) == 27277)

    checkTrue(all(tbl$a.idType %in% recognized.canonical.idTypes))
    checkTrue(all(tbl$species=="9606"))
    checkTrue(all(tbl$provider=="stamlabTF"))
    checkTrue(all(tbl$publicationID=="22959076"))
    checkTrue(all(tbl$detectionMethod=="psi-mi:MI:0606(DNase I footprinting)"))

    checkEquals(length(which(tbl$a.id=="NA")), 222)
    checkEquals(length(which(tbl$b.id=="NA")), 131)
    
} # test_addStandardColumns
#-------------------------------------------------------------------------------
# four stamlab gene symbols are not (or are no longer) HUGO standard.
#  HTLF="3344",    # FOXN2
#  ZFP161="7541",  # ZBTB14
#  ZNF238="10472", # ZBTB18
#  INSAF="3637"   # record withdrawn by ncbi, 2008
fixNonStandardGeneNames <- function(tbl)
{
    
   htlf.a.hits <- grep("HTLF", tbl$a.name)
   htlf.b.hits <- grep("HTLF", tbl$b.name)

   if(length(htlf.a.hits) > 0){
       tbl$a.name[htlf.a.hits] <- "FOXN2"
       tbl$a.id[htlf.a.hits] <- "3344"
       }

   if(length(htlf.b.hits) > 0){
       tbl$b.name[htlf.b.hits] <- "FOXN2"
       tbl$b.id[htlf.b.hits] <- "3344"
       }

   zfp161.a.hits <- grep("ZFP161", tbl$a.name)
   zfp161.b.hits <- grep("ZFP161", tbl$b.name)

   if(length(zfp161.a.hits) > 0){
       tbl$a.name[zfp161.a.hits] <- "ZBTB14"
       tbl$a.id[zfp161.a.hits] <- "7541"
       }

   if(length(zfp161.b.hits) > 0){
       tbl$b.name[zfp161.b.hits] <- "ZBTB14"
       tbl$b.id[zfp161.b.hits] <- "7541"
       }

     #  ZNF238="10472", # ZBTB18

   znf238.a.hits <- grep("ZNF238", tbl$a.name)
   znf238.b.hits <- grep("ZNF238", tbl$b.name)

   if(length(znf238.a.hits) > 0){
       tbl$a.name[znf238.a.hits] <- "ZBTB18"
       tbl$a.id[znf238.a.hits] <- "10472"
       }

   if(length(znf238.b.hits) > 0){
       tbl$b.name[znf238.b.hits] <- "ZBTB18"
       tbl$b.id[znf238.b.hits] <- "10472"
       }

     #  INSAF="3637" record withdrawn by ncbi, 2008
   insaf.a.hits <- grep("INSAF", tbl$a.name)
   insaf.b.hits <- grep("INSAF", tbl$b.name)

   if(length(insaf.a.hits) > 0){
       tbl$a.id[insaf.a.hits] <- "3637"
       }

   if(length(insaf.b.hits) > 0){
       tbl$b.id[insaf.b.hits] <- "3637"
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

   checkEquals(length(which(tbl$a.id=="NA")), 0)
   checkEquals(length(which(tbl$b.id=="NA")), 0)

      # make sure that mapping gene symbols in those 2 directories
      # catches all nonstandard gene symbols in all directories
   
   tbl <- readData()
   tbl <- addStandardColumns(tbl)
   tbl <- fixNonStandardGeneNames(tbl)
   checkEquals(length(which(tbl$a.id=="NA")), 0)
   checkEquals(length(which(tbl$b.id=="NA")), 0)
   checkEquals(length(which(tbl$a.id=="NA")), 0)
   checkEquals(length(which(tbl$b.id=="NA")), 0)

} # test_fixNonStandardGeneNames
#--------------------------------------------------------------------------------
