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
scriptsAndStyles <- function()
{
   return(c( 
         '<script src="shared/datatables/js/jquery.dataTables.min.js"></script>',
         '<script class="shiny-html-output" 
                  src= "/js/DTbinding.js"></script>',
         '<link rel = "stylesheet", 
                type = "text/css", 
                href = "shared/datatables/css/DT_bootstrap.css"></link>',
         '<style type="text/css">
                .rowsSelected td{
                background-color: rgba(112,164,255,0.2) 
                !important})  </style>',
         '<style type="text/css"> .selectable div table tbody tr{
                cursor: hand; cursor: pointer;}</style>',
         '<style type="text/css"> .selectable div table tbody tr td{
                -webkit-touch-callout: none;
                -webkit-user-select: none;
                -khtml-user-select: none;
                -moz-user-select: none;
                -ms-user-select: none;
                user-select: none;} </style>',
         '<style type="text/css">
                #myTable tfoot {display:table-header-group;}</style>',
         '<style type="text/css">
               th, td { white-space: nowrap; }</style>'))
          

} # scriptsAndStyles
#----------------------------------------------------------------------------------------------------
rowSelectableDataTableOutput <- function(outputId, ...)
{
     origStyle<- c( 
         '<script src="shared/datatables/js/jquery.dataTables.min.js"></script>',
         '<script class="shiny-html-output" 
                  src= "/js/DTbinding.js"></script>',
         '<link rel = "stylesheet", 
                type = "text/css", 
                href = "shared/datatables/css/DT_bootstrap.css"></link>',
         '<style type="text/css">
                .rowsSelected td{
                background-color: rgba(250,80,90,0.3) 
                !important})  </style>',
         '<style type="text/css"> .selectable div table tbody tr{
                cursor: hand; cursor: pointer;}</style>',
         '<style type="text/css"> .selectable div table tbody tr td{
                -webkit-touch-callout: none;
                -webkit-user-select: none;
                -khtml-user-select: none;
                -moz-user-select: none;
                -ms-user-select: none;
                user-select: none;} </style>',
         '<style type="text/css">
                #myTable tfoot {display:table-header-group;}
               th, td { white-space: nowrap; }
         </style>')
     
     tagList(
         singleton(
             tags$head(HTML(origStyle))
             ),
         div(id = outputId, class = "shiny-datatable-output selectable")
         )

} # rowSelectableDataTableOutput
#----------------------------------------------------------------------------------------------------
uiWidgets <- fluidPage(
   #tags$head(
      #tags$style(HTML("th, td { white-space: nowrap; }"))
      #HTML(scriptsAndStyles()),
      #div(id="rowSelTable", class = "shiny-datatable-output selectable")
      #),
   headerPanel("RefNet (Homo sapiens)"),
       sidebarPanel(width=2,
                    
       selectizeInput(inputId='providers', label='providers',
                      choices = providers, multiple = TRUE), # , selected=providers[1]),
                      textInput("genes", "Genes:", initial.genes),
                      submitButton("Find Interactions")
                      ),
       mainPanel(
          tabsetPanel(
               tabPanel('interactions', rowSelectableDataTableOutput(outputId="table")),
               #tabPanel('rowSelTbl', selectableDataTableOutput(outputId="rowSelTable",...)),
               tabPanel("Pubmed Abstract", htmlOutput("pubmedAbstract")),
               tabPanel("A", htmlOutput("geneA")),
               tabPanel("B", htmlOutput("geneB"))
             )
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

     output$pubmedAbstract <- renderUI({
        tags$iframe(
          height="500px",
          width="800px",
          seamless="seamless",
          src="http://www.ncbi.nlm.nih.gov/pubmed/?term=15115758")})

     output$geneA <- renderUI({
        tags$iframe(
          height="500px",
          width="800px",
          seamless="seamless",
          src="http://www.ncbi.nlm.nih.gov/gene/675")})

     output$geneB <- renderUI({
        tags$iframe(
          height="500px",
          width="800px",
          seamless="seamless",
          src="http://www.ncbi.nlm.nih.gov/gene/2177")})

      datatable.options <- list(bFilter=TRUE, bSortClasses = TRUE)
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
