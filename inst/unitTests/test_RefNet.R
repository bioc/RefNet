# test_RefNet.R
#-------------------------------------------------------------------------------
# move these to paulsTests when development slows
library(RefNet)
library(RUnit)
library(PSICQUIC)
library(org.Hs.eg.db)

if(!exists("refnet")){
    refnet <- RefNet()
    object <- refnet   # for easy debugging
    }

# PSICQUIC is included in RefNet and -- being the worthy result of an ongoing
# community effort -- sets some useful standards for us to emulate.
# make sure it is loaded, and that a sample result is available to compare
# our results against.

if(!exists("psicquic")){
    psicquic <- PSICQUIC()
    idMapper <- IDMapper("9606")
    }

if(exists("psicquic"))
   if (!exists("tbl.pqDemo"))
      if("IntAct" %in% providers(psicquic))
         tbl.pqDemo <- interactions(psicquic, id="CCNG1", species="9606",
                                    provider="IntAct")

#-------------------------------------------------------------------------------
paulsTests <- function()
{
    test_ctor()

    test_.combinations()
    test_.findHits()
    test_.filterOnColumnValue()

    test_.smartRbind()
    test_interactions()
    test_recon2Interactions()

    test_providerMix_PSICQUIC_and_native()
    test_pubmedAbstract()

    test_detectDuplicateInteractions()
    test_pickBestFromDupGroup()

} # paulsTests
#-------------------------------------------------------------------------------
test_ctor <- function()
{
    print("--- test_ctor")
    refnet <- RefNet()
    providerClasses <-  c("native", "PSICQUIC")
    checkEquals(providerClasses(refnet), providerClasses)
    providers <- providers(refnet)
    checkEquals(names(providers), providerClasses)
       # gerstein-2012  and  hypoxiaSignaling-2006  are the first two entries
       # we made into the  AnnotationHub, so let's check just those.
    checkTrue(all(c("gerstein-2012", "hypoxiaSignaling-2006") %in% providers(refnet)$native))
    if(class(refnet@psicquic) == "PSICQUIC")
       checkTrue(length(providers$PSICQUIC) > 10)

} # test_ctor
#-------------------------------------------------------------------------------
test_.combinations <- function()
{
    print("--- test_.combinations")
    checkEquals(RefNet:::.combinations(99), 99)
    checkEquals(RefNet:::.combinations(list()), list())
    checkEquals(RefNet:::.combinations(c(1,8)), list(c(1,8)))
    c148 <- RefNet:::.combinations(c(1,4,8))
    checkTrue(list(c(1,4)) %in% c148)
    checkTrue(list(c(1,8)) %in% c148)
    checkTrue(list(c(4,8)) %in% c148)

    cabc <- RefNet:::.combinations(c("a", "b", "c"), as.character)
    checkTrue(list(c("a", "b")) %in% cabc)
    checkTrue(list(c("a", "c")) %in% cabc)
    checkTrue(list(c("b", "c")) %in% cabc)


} # test_.combinations
#-------------------------------------------------------------------------------
test_.findHits <- function()
{
    print("--- test_.findHits")

    tbl <- refnet@sources[["gerstein-2012"]]
    gerstein.colnames <- list("A", "B", "A.common", "B.common", "altA", "altB",
                              "A.canonical", "B.canonical")

      # make sure that case is ignored

    upper.case.count <- length(RefNet:::.findHits(tbl, gerstein.colnames, "MYC"))
    checkTrue(upper.case.count > 700)
    checkEquals(length(RefNet:::.findHits(tbl, gerstein.colnames, "myc")),
                upper.case.count)
    checkEquals(length(RefNet:::.findHits(tbl, gerstein.colnames, "mYc")),
                upper.case.count)

      # gerstein reports just one interaction between esr1 and avp
    id <- c("ESR1", "AVP")
    hits <- RefNet:::.findHits(tbl, gerstein.colnames, id)
    checkEquals(length(hits), 1)
    checkEquals(paste(tbl[hits, 1:2]), id)

       # now just AVP, which is mentioned only as a TF target in gerstein,
       # not as a TF, thus not found in tbl$a.sym, only in tbl$b.sym

    hits.avp <- RefNet:::.findHits(tbl, gerstein.colnames, "AVP")
    checkEquals(length(hits.avp), 9)
    checkEquals(unique(tbl$B[hits.avp]), "AVP")
    checkEquals(length(grep("AVP", tbl$B.common)), 9)
    checkEquals(length(grep("AVP", tbl$A.common)), 0)

       # now just ESR1, which is mentioned as both a TF and a target in gerstein

    hits.esr1 <- RefNet:::.findHits(tbl, gerstein.colnames, "ESR1")
    checkEquals(length(hits.esr1), 127)
       # we know by inspection that there are no self-loops reported for ESR1
       # so the a.sym and b.sym hits do not overlap, and sum to the
       # total number of hits
    esr1.in.a <- grep("ESR1", tbl$A[hits.esr1], value=TRUE)
    esr1.in.b <- grep("ESR1", tbl$B[hits.esr1], value=TRUE)
    checkEquals(length(esr1.in.a), 111)
    checkEquals(length(esr1.in.b), 16)
    checkEquals(unique(esr1.in.a), "ESR1")
    checkEquals(unique(esr1.in.b), "ESR1")

       # now submit 3 gene symbols, to obtain all interactions between
       # any two of them.  manual extraction finds only two:
       # subset(tbl, a.sym %in% ids & b.sym %in% ids)[, 1:2]
       #       a.sym b.sym
       # 2895   ESR1   AVP
       # 4648 POU3F2   AVP

    ids <- c("POU3F2", "ESR1", "AVP")
    hits.3 <- RefNet:::.findHits(tbl, gerstein.colnames, ids)
    checkEquals(sort(rownames(tbl[hits.3,])),
                sort(rownames(subset(tbl, A %in% ids & B %in% ids))))

    tbl <- refnet@sources[["gerstein-2012"]]
    checkEquals(length(RefNet:::.findHits(tbl, "A.common", ids=NA)), nrow(tbl))


} # test_.findHits
#-------------------------------------------------------------------------------
test_.smartRbind <- function()
{
    print("--- test_.smartRbind")
    tbl.a <- data.frame(a=rep("aaa", 3), b=rep("bbb", 3), stringsAsFactors=FALSE)
    tbl.b <- data.frame(b=rep("BBB", 3), c=rep("CCC", 3), d=rep("DDD", 3),
                        stringsAsFactors=FALSE)
    tbl.ab <- RefNet:::.smartRbind(tbl.a, tbl.b)
    checkEquals(dim(tbl.ab), c(6,4))


} # test_.smartRbind
#-------------------------------------------------------------------------------
test_interactions <- function()
{
    print("--- test_interactions")
    tbl <- interactions(refnet, id="JUN",
                        provider=c("gerstein-2012"),
                        species="9606")
    gerstein.rows <- nrow(tbl)
    checkEquals(gerstein.rows, 259)

       # get the entire gerstein table
    tbl <- interactions(refnet, species="9606", provider=c("gerstein-2012"))
    all.gerstein.rows <- nrow(tbl)
    checkTrue(all.gerstein.rows > 6000)

       # ask for a non-existent provider
    suppressMessages(checkException(
            interactions(refnet, id="JUN", species="9606",
                         provider="bogusProvider"), silent=TRUE))

    tbl <- interactions(refnet, id="JUN",
                        provider=c("hypoxiaSignaling-2006"),
                        species="9606")
    hypoxia.rows <- nrow(tbl)
    checkTrue(hypoxia.rows == 1)

    tbl <- interactions(refnet, id="JUN",
                        provider=c("gerstein-2012", "hypoxiaSignaling-2006"),
                        species="9606")
    checkEquals(nrow(tbl), hypoxia.rows + gerstein.rows)

} # test_interactions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
test_recon2Interactions <- function()
{
    print("--- test_recon2Interactions")
    if(!"recon2" %in% providers(refnet))
        return()

    tbl <- interactions(refnet, id="SHMT1", species="9606", type="modifies",
                        provider="recon2")
    checkEquals(tbl$A, c("_6470_1_c","_6470_2_c","_6470_2_c","_6470_1_c"))
    checkEquals(tbl$B, c("R_GHMT2r","R_GHMT2r","R_GHMT3","R_GHMT3"))
    checkEquals(tbl$aliasA, c("SHMT1","SHMT1","SHMT1","SHMT1"))
    checkEquals(tbl$aliasB, c("glycine hydroxymethyltransferase, reversible",
                              "glycine hydroxymethyltransferase, reversible",
                              "glycine hydroxymethyltransferase",
                              "glycine hydroxymethyltransferase"))
    checkEquals(tbl$altA, c("6470", "6470", "6470", "6470"))
    checkEquals(tbl$altB, c("-", "-", "-", "-"))
    checkEquals(tbl$type, c("modifies", "modifies","modifies","modifies"))

        # now test the type argument.  the above reaction, catalzyed by SHMT1
        # has 2 substrates, 2 products, and 2 modifiers
    tbl <- interactions(refnet, id="R_GHMT3", type=NA, provider="recon2")
    checkEquals(nrow(tbl), 6)  # all interactions for this reaction

    tbl <- interactions(refnet, id="R_GHMT3", type="modifies", provider="recon2")
    checkEquals(nrow(tbl), 2)
    checkEquals(tbl$type, c("modifies", "modifies"))

    tbl <- interactions(refnet, id="R_GHMT3", type="substrateOf",
                        provider="recon2")
    checkEquals(nrow(tbl), 2)
    checkEquals(tbl$type, c("substrateOf", "substrateOf"))

    tbl <- interactions(refnet, id="R_GHMT3", type="productOf",
                        provider="recon2")
    checkEquals(nrow(tbl), 2)
    checkEquals(tbl$type, c("productOf", "productOf"))

    tbl <- interactions(refnet, id="R_GHMT3", type=c("productOf", "substrateOf"),
                        provider="recon2")
    checkEquals(nrow(tbl), 4)
    checkEquals(sort(tbl$type), c("productOf", "productOf",
                                  "substrateOf", "substrateOf"))

} # test_recon2Interactions
#-------------------------------------------------------------------------------
# demonstrate mixing PSICQUIC sources with RefNet "native"
test_providerMix_PSICQUIC_and_native <- function()
{
    print("--- test_providerMix_PSICQUIC_and_native")

    if(class(refnet@psicquic) == "logical")  # detects NA
        return()

    if(!"BioGrid" %in% providers(refnet)$PSICQUIC){
        print("  BioGrid not available, skipping test")
        return()
        }

    tbl.1 <- interactions(refnet, "TERT", provider="gerstein-2012")
    tbl.2 <- interactions(refnet, "TERT", provider="BioGrid")
    tbl.3 <- interactions(refnet, "TERT", provider=c("BioGrid", "gerstein-2012"))
    checkEquals(nrow(tbl.3), nrow(tbl.1) + nrow(tbl.2))
    checkEquals(sort(unique(tbl.3$provider)), c("BioGrid", "gerstein.2012"))

    tbl.4 <- interactions(refnet, "TERT", provider=c("BioGrid", "gerstein-2012"))
    checkEquals(nrow(tbl.4), nrow(tbl.1) + nrow(tbl.2))
    checkEquals(sort(unique(tbl.4$provider)), c("BioGrid", "gerstein.2012"))

} # test_providerMix_PSICQUIC_and_native
#-------------------------------------------------------------------------------
test_pubmedAbstract <- function()
{
    print("--- test_pubmedAbstract")
    text.lines <- pubmedAbstract("22959076")   # split is TRUE by default
    checkTrue(nchar(paste(text.lines, collapse="")) > 1000)
    checkTrue(length(text.lines) > 20)
    title <- paste("Circuitry and dynamics of human transcription factor",
                   "regulatory networks")
    checkEquals(grep(title, text.lines), 5)


} # test_pubmedAbstract
#-------------------------------------------------------------------------------
test_detectDuplicateInteractions <- function()
{
    print("--- test_detectDuplicateInteractions")
    filename <- system.file(package="RefNet", "extdata", "tbl.28g2.RData")
    load(filename, envir=.GlobalEnv)
    tbl.ap <- detectDuplicateInteractions(head(tbl.28g2))
    checkEquals(tbl.ap$dupGroup, c(0,0,1,1,2,2))
    checkEquals(sort(unique(tbl.ap$pairSig)),
                c("23304:4609","4609:7157","55643:7157","56254:7157"))

    tbl.ap <- detectDuplicateInteractions(tbl.28g2)
    checkEquals(dim(tbl.ap), c(62, 22))
       # sample subset, dupGroup1, with 11 members
    tbl.dg1 <- subset(tbl.ap, dupGroup==1)[, c("A.common", "B.common",
                          "detectionMethod","type", "provider", "publicationID")]
    checkEquals(nrow(tbl.dg1), 11)

       # now try a more complex case, one which requires that
       # rows lacking canonical names be eliminated first

    # sample data file, with some rows' canonical names unassigned, created
    # this way:
    #    tbl <- interactions(refnet, id="E2F3", species="9606")
    #    tbl <- addStandardNames(idMapper, tbl)
    #     save(tbl, file="../inst/extdata/tbl.e2f3.RData")

    filename <- system.file(package="RefNet", "extdata", "tbl.e2f3.RData")
    load(filename)

    removers <- with(tbl, unique(c(grep("^-$", A.canonical),
                                   grep("^-$", B.canonical))))
    if(length(removers) > 0)
        tbl <- tbl[-removers,]

    tbl.dups <- detectDuplicateInteractions(tbl)
    checkTrue(ncol(tbl.dups) > ncol(tbl))
    checkTrue("dupGroup" %in% colnames(tbl.dups))
    distinct.groups <- unique(tbl.dups$dupGroup)
    checkEquals(length(distinct.groups), 136)
    checkEquals(sum(as.integer(table(tbl.dups$dupGroup))), nrow(tbl))
       # 309 interactions are singletons
    checkEquals(nrow(subset(tbl.dups, dupGroup==0)), 309)

} # test_detectDuplicateInteractions
#-------------------------------------------------------------------------------
test_pickBestFromDupGroup <- function()
{
    print("--- test_pickBestFromDupGroup")
    filename <- system.file(package="RefNet", "extdata", "tbl.dups.RData")
    load(filename, envir=.GlobalEnv)

       # dupGroup zero is a special case.  it contains all
       # of the a-b/b-a interactions which have NO duplicates
       # if the are of the correct type, we want to keep them all
       # this tbl.dups has two singletons, both of approved types

    dupGroups <- sort(unique(tbl.dups$dupGroup))   # 0, 1, 2, 3, 4
    preferred.types <- c("direct", "physical", "aggregation")
    singleton.rows <- pickBestFromDupGroup(0, tbl.dups, preferred.types)
    checkEquals(length(singleton.rows), 2)
    checkEquals(sort(singleton.rows), c("15","5"))
    checkEquals(tbl.dups["5", "dupGroup"], 0)
    checkEquals(tbl.dups["15", "dupGroup"], 0)

       # now check dupGroup 3, which has 23 (!) entries
    checkEquals(nrow(subset(tbl.dups, dupGroup==3)), 23)

    rowName.best <- pickBestFromDupGroup(3, tbl.dups, preferred.types)
       # "physical" matches "physical interaction" and is the best of the 23
       # there are 3 with that interaction type: "12", "13", "42".
       # we pick the first
    checkEquals(rowName.best, "12")

} # test_pickBestFromDupGroup
#-------------------------------------------------------------------------------
# 'type' is a standard column in all of our interaction sources.  here we
# test our ability to exclude rows in an interaction table which do not
# match specified types
test_.filterOnColumnValue <- function()
{
    print("--- test_.filterOnColumnValue")

    provider <- "hypoxiaSignaling-2006"
    checkTrue(provider %in% unlist(providers(refnet), use.names=FALSE))
    tbl <- refnet@sources[[provider]]

    checkEquals(nrow(RefNet:::.filterOnColumnValue(tbl, "b.modification",
                                                   "hydroxylated")), 1)
    checkEquals(nrow(RefNet:::.filterOnColumnValue(tbl, "b.modification",
                                                   "proteolyzed")), 1)
    checkEquals(nrow(RefNet:::.filterOnColumnValue(tbl, "b.modification",
                                                   "bogus")), 0)

    checkEquals(nrow(RefNet:::.filterOnColumnValue(tbl, "B.common", "VEGFA")), 5)
    checkEquals(nrow(RefNet:::.filterOnColumnValue(tbl, "B.common", "HIF1A")), 5)
    checkEquals(nrow(RefNet:::.filterOnColumnValue(tbl, "B.common", "KDR")), 1)

} # test_.filterOnColumnValue
#-------------------------------------------------------------------------------
