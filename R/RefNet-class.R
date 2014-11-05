#-------------------------------------------------------------------------------
.RefNet <- setClass("RefNet",
                    slots=c(providers="list",
                            psicquic="ANY",
                            sources="list"))

#-------------------------------------------------------------------------------
.printf <- function(...) print(noquote(sprintf(...)))
#-------------------------------------------------------------------------------
RefNet <- function()
{
    print(noquote("initializing PSICQUIC..."))
    psicquic <- PSICQUIC()

    print(noquote("initializing RefNet from AnnotationHub..."))
    ah <- AnnotationHub()
    ah <- subset(ah, ah$dataprovider=="RefNet")
    tbl.refnet <- mcols(ah)
    pathnames <- sub("refnet/", "refnet.", tbl.refnet$rdatapath)
    pathnames <- sub("-", ".", pathnames)
    titles <- sub("interactions from ", "", tbl.refnet$title)
   
    refnet.tables <- vector('list', length=length(pathnames))
    for(i in seq_len(length(pathnames))){
        pathname <- pathnames[[i]]
        tbl <- ah[[pathname]]
        refnet.tables[[i]] <- tbl
        } # for i
    names(refnet.tables) <- titles
   
    object <- .RefNet()
    object@psicquic <- psicquic

    object@sources <- refnet.tables;

  
    suppressWarnings(psicquic.na <- is.na(psicquic))
    
    if(psicquic.na)
        psicquic.providers <- list()
    else
        psicquic.providers <- providers(psicquic)

    object@providers <- list(native=names(object@sources),
                             PSICQUIC=psicquic.providers)
   
    print(noquote("RefNet ready."))
    
    object
 
} # RefNet ctor
#-------------------------------------------------------------------------------
setGeneric("providerClasses", signature="object", function(object)
           standardGeneric ("providerClasses"))

setMethod ("providerClasses", "RefNet",

   function (object){
       names(object@providers)
       })

#-------------------------------------------------------------------------------
setMethod("providers", "RefNet",

   function (object){
       object@providers
       })

#-------------------------------------------------------------------------------
setMethod('show', 'RefNet',

    function(object) {
        providers <- providers(object)
        all.providers <- unlist(providers, use.names=FALSE)
        msg <- sprintf("RefNet object with %d providers in %d classes",
                       length(all.providers),
                       length(names(providers)))
        cat (msg, '\n', sep='')
        if(length(all.providers) == 0)
            return()
        classes <- providerClasses(object)
        for(class in classes){
            msg <- sprintf("| provider class '%s':", class)
            cat(msg, "\n", sep="")
            for(provider in providers[[class]]){
                msg <- sprintf("|     %s", provider)
                cat(msg, "\n", sep="")
                } # for provider
            } # for class

        }) # show

#-------------------------------------------------------------------------------
.interactions <- function(object, id, species, speciesExclusive, type,
                          provider, detectionMethod, publicationID, quiet)
{
    result <- data.frame()
    available.providers <- unlist(providers(object), use.names=FALSE)
    
    if(all(is.na(provider)))
        provider <- available.providers

    unrecognized.providers <- setdiff(provider, available.providers)
                                      
    if(length(unrecognized.providers) > 0){
       msg <- sprintf("RefNet::interactions does not recognize provider '%s'",
                      paste(unrecognized.providers, collapse=","))
       stop(msg)
       }

    providers <- provider  # for clarity, make this explicitly plural
   

    native.providers <- intersect(providers, providers(object)$native)

    for(i in seq_len(length(native.providers))) {
        native.provider <- native.providers[i]
        tbl <- object@sources[[native.provider]]
        standard.colnames <- list("A", "B", "A.common", "B.common", 
                                  "A.canonical", "B.canonical", "publicationID")
        hits <- .findHits(tbl, standard.colnames, id)
        if(length(hits) > 0) {
            tbl.sub  <- tbl[hits,]
            if(all(!is.na(type)))
                tbl.sub <- .filterOnColumnValue(tbl.sub, "type", type)
            if(all(!is.na(detectionMethod)))
                tbl.sub <- .filterOnColumnValue(tbl.sub, "detectionMethod",
                                                detectionMethod)
            if(nrow(tbl.sub) > 0)
                result <- .smartRbind(result, tbl.sub)
           } # hits > 0
        } # for i

    psicquic.providers <- intersect(providers, providers(object)$PSICQUIC)
    if(length(psicquic.providers) > 0){
           # dispatch to the interactions method in the PSICQUIC class
        new.rows <- interactions(object@psicquic, id, species, speciesExclusive,
                                 type, psicquic.providers, detectionMethod,
                                 publicationID, quiet)
        result <- .smartRbind(result, new.rows)
        }

    result

} # .interactions
#-------------------------------------------------------------------------------
setMethod ("interactions", "RefNet",

    function(object, id, species, speciesExclusive, type,
             provider, detectionMethod, publicationID, quiet) {

        .interactions(object, id, species, speciesExclusive, type,
                      provider, detectionMethod, publicationID, quiet)
        })

#-------------------------------------------------------------------------------
# given a list (of any scalar type), return all combinations, not distinguishing
# order (i.e., c(1,2) is the same as c(2,1), cast to the specified type.
# i am surprised (and may be wrong): this function is not provided by base R!
.combinations <- function(vec, casting.function=as.numeric)
{
    if(length(vec) <= 1)
        return(vec)
    
    mtx <- outer(vec, vec, paste, sep=",")
    pair.strings <- mtx[upper.tri(mtx)]
    tmp <- strsplit(pair.strings, ",")
    lapply(strsplit(pair.strings, ","), casting.function)

} # .combinations
#-------------------------------------------------------------------------------
# in our convention, hits of an identifier against a refnet data source is
# defined as 
#  1) if one id only, then a hit is any interaction which mentions the id, in
#     any of the candidate columns (which are typically any identifier column in
#     the refnet data source
#  2) if two or more ids, then a hit is any interaction between any two of the
#     specified ids.
# definition 2 requires that all combinations of pairs be tested, for which the
# .combinations function is useful.

.findHits <- function(tbl, candidate.colnames, ids=NA)
{
    stopifnot(all(candidate.colnames %in% colnames(tbl)))

       # prevent partial matches, e.g., CAD and "ACADM", "CAD","SMARCAD1"

    if(all(is.na(ids)))
        return(seq_len(nrow(tbl)))
    
    ids.delimited <- sprintf("^%s$", ids)

       # for each delimited id, search in each of the candidate columns
    
    raw.hits <- with(tbl,
                     lapply(ids.delimited,
                        function(id)
                           lapply(candidate.colnames,
                                  function(col) grep(id, tbl[,col], ignore.case=TRUE))))
    names(raw.hits) <- ids

        # collapse the multiple column info for each id into a list of
        # rownumbers for each id.  rows.by.id becomes a named list, one element
        # per id.  we will use this simple list to identify interactions between
        # possibly multiple ids.  if just one id, then all of its interactions
        # will be returned
    
    
    rows.by.id <- lapply(raw.hits,
                         function(raw.hit) unlist(raw.hit, use.names=FALSE))

    if(length(ids) == 1){
        real.hits <- unique(rows.by.id[[1]])
    } else {
        possible.pairs <- .combinations(1:length(ids), as.numeric)
        real.hits <- c()
        for(pair in possible.pairs){
            real.hits <- c(real.hits, intersect(rows.by.id[[pair[1]]],
                                                rows.by.id[[pair[2]]]))
            } # for pair
        } # else: length(ids) > 1
                            
    real.hits

}# .findHits
#-------------------------------------------------------------------------------
.smartRbind <- function(tbl.a, tbl.b, fill="-")
{
    orig.a.cols <- colnames(tbl.a)
    orig.b.cols <- colnames(tbl.b)
   
    extra.a.cols <- setdiff(orig.a.cols, orig.b.cols)
    extra.b.cols <- setdiff(orig.b.cols, orig.a.cols)

    empty.a.col <- rep(fill, nrow(tbl.a))
    empty.b.col <- rep(fill, nrow(tbl.b))

    for(extra.a.col in extra.a.cols)
        tbl.b <- cbind(tbl.b, empty.b.col, stringsAsFactors=FALSE)

    for(extra.b.col in extra.b.cols)
        tbl.a <- cbind(tbl.a, empty.a.col, stringsAsFactors=FALSE)

    colnames(tbl.a) <- c(orig.a.cols, extra.b.cols)
    colnames(tbl.b) <- c(orig.b.cols, extra.a.cols)

    tbl.b <- tbl.b[, colnames(tbl.a)]
    result <- rbind(tbl.a, tbl.b)
   
} # .smartRbind
#-------------------------------------------------------------------------------
pubmedAbstract <- function (pmid, split=TRUE)
{
    p2 = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="
    p3 = pmid
    p4 = "&retmode=text&rettype=abstract"
    p5 = "&email=pshannon@fhcrc.org"          # TODO -- generalize this
    url = paste (p2, p3, p4,p5, sep='')
    text <- getURL(url)
    if(split)
        text <- strsplit(text, "\n")[[1]]

    text

} # pubmedAbstract
#-------------------------------------------------------------------------------
extractPubmedIDs <- function(rawIDs)
{
    hits <- gregexpr("pubmed:([0-9]+)", rawIDs)
    pmids.list <- regmatches(rawIDs, hits)
    pubmedIDs <- lapply(pmids.list, function(pmids) gsub("pubmed:", "", pmids))
    names(pubmedIDs) <- rawIDs

    pubmedIDs

} # extractPubmedIDs
#-------------------------------------------------------------------------------
# duplicate means:  a->b is equivalent to b->a
# duplicate means:  a->b is equivalent to b->a
detectDuplicateInteractions <- function(tbl)
{
   stopifnot(all(c("A.canonical", "B.canonical") %in% colnames(tbl)))
   a <- tbl$A.canonical
   b <- tbl$B.canonical

     # we detect duplicates by comparing sorted canonical names
     # make sure no unassigned rows are in the tbl
  
   unassigned.id.a <-length(grep("^-$", a))
   unassigned.id.b <-length(grep("^-$", b))
  
   if(unassigned.id.a | unassigned.id.b){
       msg.0 <- sprintf("found unassigned identifiers")
       msg.1 <- sprintf("     unassigned A.id count: %d/%d",
                              unassigned.id.a, length(a))
       msg.2 <- sprintf("     unassigned B.id count: %d/%d",
                              unassigned.id.b, length(b))
       warning(paste("", msg.0, msg.1, msg.2, sep="\n"))
       }
      
   max <- nrow(tbl)
   sigs <- vector("character", max)

   for (i in 1:max){ 
       ab.sorted <- sort(c(a[i], b[i]))
       sig <- paste(ab.sorted, collapse=":")
       sigs[i] <- sig
      } # for i

   dups <- which(duplicated(sigs))
   dup.sigs <- unique(sigs[dups])
   dup.group <- rep(0, nrow(tbl))
   if(length(dup.sigs > 0)){
       for(i in 1:length(dup.sigs)){
           member.indices <- which(sigs == dup.sigs[i])
           dup.group[member.indices] <- i
           } # for i
       } # if dup.sigs


   tbl <- cbind(tbl, pairSig=sigs, dupGroup=dup.group, stringsAsFactors=FALSE)
   tbl

}# detectDuplicateInteractions
#-------------------------------------------------------------------------------
old.detectDuplicateInteractions <- function(tbl)
{
   stopifnot(all(c("A.canonical", "B.canonical") %in% colnames(tbl)))
   a <- tbl$A.canonical
   b <- tbl$B.canonical

     # we detect duplicates by comparing sorted canonical names
     # make sure no unassigned rows are in the tbl
  
   unassigned.canonical.a <-length(grep("^-$", a))
   unassigned.canonical.b <-length(grep("^-$", b))
  
   if(unassigned.canonical.a | unassigned.canonical.b){
       msg.0 <- sprintf("found unassigned identifiers")
       msg.1 <- sprintf("     unassigned A.canonical count: %d",
                              unassigned.canonical.a)
       msg.2 <- sprintf("     unassigned B.canonical count: %d",
                              unassigned.canonical.b)
       stop(paste("", msg.0, msg.1, msg.2, sep="\n"))
       }
      
   max <- nrow(tbl)
   sigs <- vector("character", max)

   for (i in 1:max){ 
       ab.sorted <- sort(c(a[i], b[i]))
       sig <- paste(ab.sorted, collapse=":")
       sigs[i] <- sig
      } # for i

   dups <- which(duplicated(sigs))
   dup.sigs <- unique(sigs[dups])
   dup.group <- rep(0, nrow(tbl))
   if(length(dup.sigs > 0)){
       for(i in 1:length(dup.sigs)){
           member.indices <- which(sigs == dup.sigs[i])
           dup.group[member.indices] <- i
           } # for i
       } # if dup.sigs


   tbl <- cbind(tbl, pairSig=sigs, dupGroup=dup.group, stringsAsFactors=FALSE)
   tbl

}# old.detectDuplicateInteractions
#-------------------------------------------------------------------------------
# a hasty solution to an important problem
# TODO:
#   1) vectorize dupGrp argument
#   2) allow preferred detectionMethods, providers
#   3) support an "order" like precedence scheme
#
pickBestFromDupGroup <- function(dupGrp, tbl.dups, preferred.interaction.types)
{
    rows.of.interest <- which(tbl.dups$dupGroup==dupGrp)

    tbl.sub <- tbl.dups[rows.of.interest, 
                       c("A.canonical", "B.canonical", "detectionMethod",
                         "type", "provider", "publicationID")]

    
    all.types.found <- sort(unique(tbl.sub$type))

    approved.types.found <- sapply(preferred.interaction.types,
                              function(type) length(grep(type, tbl.sub$type)))
    if(!any(approved.types.found) > 0)
        return(NA)
    
    if(dupGrp==0){ # special case, no duplicates, consider keeping all of these
        best.types <- names(which(approved.types.found > 0))
        rows.by.type <- sapply(best.types,
                                  function(best.type)
                                      grep(best.type, tbl.sub$type,
                                           ignore.case=TRUE))
        best.rows <- sort(unlist(rows.by.type, use.names=FALSE))
        return(rownames(tbl.sub)[best.rows])
        }

      # dupGrp != 0, and therefore multiple rows for each interaction. choose
      # the first of the best, as specified in the order of
      # preferred.interaction.types above
         
    best.type <- names(which(approved.types.found > 0))[1]
    result <- tbl.sub[grep(best.type, tbl.sub$type),]
    return(rownames(result)[1])   # choose the first among equals

} # pickBestFromDupGroup
#-------------------------------------------------------------------------------
# use this for detectionMethod and type (that is, "interaction type")
.filterOnColumnValue <- function(tbl, colname, value)
{
    if(nrow(tbl) == 0)
        return(tbl)
   
    if(!colname %in% colnames(tbl)){
        msg <- sprintf("'%s' not in colnames(tbl): %s",
                       colname, paste(colnames(tbl), collapse=","))
        stop(msg)
        }

    if(any(is.na(value)))
        return(tbl)

    keeper.rows <- unique(unlist(lapply(value,
                                        function(t) grep(t, tbl[, colname]))))

    if(length(keeper.rows) == 0)
        return(data.frame())

    return(tbl[keeper.rows,])

} # .filterOnColumnValue
#-------------------------------------------------------------------------------
