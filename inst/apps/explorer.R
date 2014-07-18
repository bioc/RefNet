library(shiny)
library(RefNet)
if(!exists("refnet"))
    refnet <- RefNet()

if(!exists("idMapper"))
    idMapper <- IDMapper()

providers <- unlist(providers(refnet), use.names=FALSE)

#----------------------------------------------------------------------------------------------------
cleanGenes <- function(s)
{
   genes.raw <- sub("[[:space:]]+$","", s)
   toupper(strsplit(genes.raw, " +")[[1]])
}
#----------------------------------------------------------------------------------------------------
uiWidgets <- pageWithSidebar(
       headerPanel("RefNet"),
       sidebarPanel(width=2,
       selectizeInput(inputId='providers', label='providers',
                      choices = providers, multiple = TRUE, selected=providers[1]),
                      textInput("genes", "Genes:", "ABCB10"),
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
       genes <- cleanGenes(input$genes)
       current.providers <- input$providers
       printf("querying for %s from %s", paste(genes, collapse=","), paste(current.providers, collapse=","))
       tbl <- interactions(refnet, id=genes, provider=current.providers, quiet=FALSE)
       printf("dim(tbl): %d x %d", nrow(tbl), ncol(tbl))
       #tbl2 <- addGeneInfo(idMapper, tbl)
       tbl2 <- tbl
       desired.columns <- c("A.name", "B.name", "type", "publicationID", "A", "B", "detectionMethod", "provider")
       actual.columns <- intersect(desired.columns, colnames(tbl2))
       unique(tbl2[, actual.columns])
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
