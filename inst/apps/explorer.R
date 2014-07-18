library(shiny)
library(RefNet)
if(!exists("refnet"))
    refnet <- RefNet()

if(!exists("idMapper"))
    idMapper <- IDMapper("9606")

providers <- unlist(providers(refnet), use.names=FALSE)
initial.genes <- "ABCB10"
initial.genes <- ""
empty.data.frame <- data.frame(A=c("a"), B=("b"), type="type", publicationID="", A="", B="",
                               detectionMethod="", provider="", stringsAsFactors=FALSE)[-1,]

#----------------------------------------------------------------------------------------------------
cleanGenes <- function(s)
{
   genes.raw <- sub("[[:space:]]+$","", s)
   toupper(strsplit(genes.raw, " +")[[1]])
}
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
uiWidgets <- fluidPage(
   tags$head(
      tags$style(HTML("th, td { white-space: nowrap; }"))),
   headerPanel("RefNet (Homo sapiens)"),
       sidebarPanel(width=2,
                    
       selectizeInput(inputId='providers', label='providers',
                      choices = providers, multiple = TRUE), # , selected=providers[1]),
                      textInput("genes", "Genes:", initial.genes),
                      submitButton("Find Interactions")
                      ),
       mainPanel(
          verbatimTextOutput("summary"),
          dataTableOutput(outputId="table")
          )
       ) # uiWidgets

#----------------------------------------------------------------------------------------------------
serverFunction <- function(input, output)
{
    findInteractions <- reactive({
       desired.columns <- c("A", "B", "type", "publicationID", "A", "B", "detectionMethod", "provider")
       genes <- cleanGenes(input$genes)
       if(length(genes) == 0)
           return(empty.data.frame)
       current.providers <- input$providers
       printf("querying for %s from %s", paste(genes, collapse=","), paste(current.providers, collapse=","))
       tbl <- interactions(refnet, id=genes, provider=current.providers, quiet=FALSE, species="9606")
       printf("dim(tbl): %d x %d", nrow(tbl), ncol(tbl))
       if(nrow(tbl) == 0){
          printf("nrow 0, returning empty data frame")
          return(empty.data.frame)
          }
       tbl2 <- addGeneInfo(idMapper, tbl)
       #tbl2 <- tbl
       desired.columns <- c("A.name", "B.name", "type",  "detectionMethod", "provider", "publicationID", "A", "B")
       actual.columns <- intersect(desired.columns, colnames(tbl2))
       tbl3 <- unique(tbl2[, actual.columns])
       tbl4 <- simplify.psicquic.strings(tbl3)
       simplify.pubmed.ids(tbl4)
       }) # reactive

    output$summary <- renderPrint({
          #dataset <- datasetInput()
          #summary(dataset)
          output$value <- renderPrint({ input$action })
          cleanGenes(input$genes)
          }) # renderPrint
      datatable.options <- list(bFilter=TRUE)
      output$table <- renderDataTable(findInteractions(),options=datatable.options)

} # serverFunction
#----------------------------------------------------------------------------------------------------
.explore <- function(
             id=NA,
             species=NA,
             speciesExclusive=TRUE,
             type=NA,
             provider=NA,
             detectionMethod=NA,
             publicationID=NA,
             quiet=TRUE,
             port=9999)
{  
    app <- list(ui=uiWidgets, server=serverFunction)
    runApp(app)

} # .explore
#----------------------------------------------------------------------------------------------------
