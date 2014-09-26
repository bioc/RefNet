library(shiny)
#library(shinyIncubator)
library(RefNet)

refnet <- RefNet()
idMapper <- IDMapper(species="9606")

empty.data.frame <- data.frame(A.name=c(""), B.name=(""), type="", publicationID="", A="", B="",
                               detectionMethod="", provider="", stringsAsFactors=FALSE) # [-1,]

printf <- function(...)print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
simplify.psicquic.strings <- function(tbl)
{
   column.names <- colnames(tbl)
   for(colname in column.names){
      column <- tbl[, colname]
      fixed <- sub("\\).*$", "", sub("^.*\\(", "", column))
      tbl[, colname] <- fixed
      }

  tbl

} # simplify.psicquic.stings
#----------------------------------------------------------------------------------------------------
simplify.pubmed.ids <- function(tbl)
{
    column.number <- grep("publicationID", colnames(tbl))

    if(column.number > 0){
       column.values <- tbl[, column.number]
       new.column.values <- gsub(".*pubmed:([0-9]+).*", "\\1", column.values)
       tbl[, column.number] <- new.column.values
       colnames(tbl)[column.number] <- "pubmedID"
       }

   tbl
                          
} # simplify.psicquic.stings
#----------------------------------------------------------------------------------------------------
cleanGenes <- function(s)
{
   genes.raw <- sub("[[:space:]]+$","", s)
   toupper(strsplit(genes.raw, " +")[[1]])

} # cleanGenes
#----------------------------------------------------------------------------------------------------
shinyServer(function(input, output, session) {

    pubmedURL <- reactive({
        printf("entering pubmedURL reactive function");
        result <- paste0("http://www.ncbi.nlm.nih.gov/pubmed/?term=", input$hiddenPmidDiv)
        print(result)
        result
        })

    geneAURL <-  reactive({
        printf("entering geneAURL reactive function");
        result <- paste0("http://www.ncbi.nlm.nih.gov/gene/?term=", input$hiddenGeneA_id_Div)
        print(result)
        result
        })

    geneBURL <-  reactive({
        printf("entering geneBURL reactive function");
        #result <- paste0("http://www.ncbi.nlm.nih.gov/gene/?term=", input$hiddenGeneBDiv)
        result <- paste0("http://www.ncbi.nlm.nih.gov/gene/?term=", input$hiddenGeneB_id_Div)
        print(result)
        result
        })

    geneInput <- reactive({
        printf("entering inputGenes reactive function");
        result <- input$genes
        print(result)
        result
        })


    #observe({
    #   if(!is.na(initial.genes))
    #      updateTextInput(session, inputId="genes", value=initial.genes)
    #   })
        
    findInteractions <- reactive({
       input$goButton
       printf("reactive findInteractions running")
       #genes <- cleanGenes(input$genes)
       genes <- cleanGenes(isolate(input$genes))
       if(length(genes) == 0){
          printf("length(genes) == 0, returning");
          return(empty.data.frame)
          }
       if(all(nchar(genes) == 0)){
          printf("all(nchar(genes) == 0)), returning");
          return(empty.data.frame)
          }
       current.providers <- isolate(input$providers)
       if(length(current.providers) == 0){
          printf("length(current.providers) == 0, returning");
          return(empty.data.frame)
           }
       if(all(nchar(current.providers) == 0)){
          printf("all(nchar(current.providers) == 0)), returning");
          return(empty.data.frame);
          }

       printf("current.providers: %s", paste(current.providers, collapse=","))
       printf("   providers: %s", paste(current.providers, collapse=","))
       printf("querying for %s from %s", paste(genes, collapse=","), paste(current.providers, collapse=","))
       if("ALL providers" %in% current.providers) {
          current.providers = NA;
          printf("querying ALL providers")
          }
       tbl <- interactions(refnet, id=genes, provider=current.providers, quiet=FALSE, species="9606")
       printf("dim(tbl): %d x %d", nrow(tbl), ncol(tbl))
       if(nrow(tbl) == 0){
          printf("nrow 0, returning empty data frame")
          return(empty.data.frame)
          }
       printf("addGeneInfo...")
          # the RefNet native data tables have been updated, using
          #    A|B.name rather than A|B.common
          #    A|B.id   rather than A|B.canonical
          # but these new versions and their new column names are not yet
          # avaialable via the AnnotationHub
          # catch them here so that they can be recognized as "already name-mapped"
          # and thus so that the idMapper will let them stand
       colnames(tbl) <- sub(".common$", ".name", colnames(tbl))
       colnames(tbl) <- sub(".canonical$", ".id", colnames(tbl))

       tbl2 <- addGeneInfo(idMapper, tbl)
       xxxx <<- tbl2
       #tbl2 <- tbl
       desired.columns <- c("A.name", "B.name", "type", "publicationID", "detectionMethod", "provider", "A.id", "B.id")
       actual.columns <- intersect(desired.columns, colnames(tbl2))
       if(length(actual.columns) == length(desired.columns))
          actual.columns <- desired.columns
       tbl3 <- unique(tbl2[, actual.columns])
       printf("simplifying...")
       tbl4 <- simplify.psicquic.strings(tbl3)
       printf("colnames: %s", paste(colnames(tbl4), collapse=","))
       simplify.pubmed.ids(tbl4)
       }) # reactive


     output$pubmedAbstract <- renderUI({
        printf("pmid? %s", input$pmid);
        printf("url:  %s", pubmedURL());
        tags$iframe(
          height="500px",
          width="800px",
          seamless="seamless",
          src=pubmedURL())})

     output$geneA <- renderUI({
        tags$iframe(
          height="500px",
          width="800px",
          seamless="seamless",
          src=geneAURL())})

     output$geneB <- renderUI({
        input$hiddenGeneB_id_Div
        tags$iframe(
          height="500px",
          width="800px",
          seamless="seamless",
          src=geneBURL())})

      datatable.options <- list() # list(bFilter=TRUE, bSortClasses = TRUE)
      output$table <- renderDataTable(findInteractions(),options=datatable.options)

      output$selectedInteractionsTable <- renderDataTable(empty.data.frame, options=datatable.options)


}) # shinyServer
#----------------------------------------------------------------------------------------------------
