# source("explorer.R"); 
# 
library(httr)
#----------------------------------------------------------------------------------------------------
library(shiny)
library(RefNet)

if(!exists("refnet"))
    refnet <- RefNet()

if(!exists("idMapper")){
   idMapper <- IDMapper("9606")
   }

providers <- unlist(providers(refnet), use.names=FALSE)
initial.genes <- "ABCB10"
initial.genes <- ""
initial.providers <- ""

empty.data.frame <- data.frame(A.name=c(""), B.name=(""), type="", publicationID="", A="", B="",
                               detectionMethod="", provider="", stringsAsFactors=FALSE)

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
         '<script src="shared/datatables/js/jquery.dataTables.js"></script>',
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
javascript <- function()
{
    text <-  '<script> $(document).ready(function() {$("#table").on("click", "tr", function () {
                     console.log("table clicked!");
                     var Aname = $("td", this).eq(0).text();
                     var Bname = $("td", this).eq(1).text();
                     var type  = $("td", this).eq(2).text();
                     var pmid =  $("td", this).eq(3).text();
                     var A_id = $("td", this).eq(6).text();
                     var B_id = $("td", this).eq(7).text();
                     console.log("A_id: " + A_id);
                     window.pmid = pmid;
                     console.log("pmid ? " + pmid)

                     $("#hiddenPmidDiv").html(pmid)
                     Shiny.onInputChange("hiddenPmidDiv", pmid);

                     $("#hiddenGeneADiv").html(Aname)
                     Shiny.onInputChange("hiddenGeneADiv", Aname);

                     $("#hiddenGeneBDiv").html(Bname)
                     Shiny.onInputChange("hiddenGeneBDiv", Bname);

                     $("#hiddenGeneA_id_Div").html(A_id)
                     Shiny.onInputChange("hiddenGeneA_id_Div", A_id);

                     $("#hiddenGeneB_id_Div").html(B_id)
                     Shiny.onInputChange("hiddenGeneB_id_Div", B_id);
                     })})
               </script>'

    return(text)

} # javascript
#----------------------------------------------------------------------------------------------------
uiWidgets <- fluidPage(
   #tags$head(
      #tags$style(HTML("th, td { white-space: nowrap; }"))
      #HTML(scriptsAndStyles()),
      #div(id="rowSeolTable", class = "shiny-datatable-output selectable")
      #),
   headerPanel("RefNet (Homo sapiens)"),
       sidebarPanel(width=2,
          selectizeInput(inputId='providers', label='providers',
                         choices = providers, multiple = TRUE,
                         #selected = initial.providers),
                         selected="mentha"),
          textInput("genes", "Genes:", initial.genes),
          #submitButton("Find Interactions"),
          actionButton("goButton", "Find interactions"),
          htmlOutput("hiddenPmidDiv"),
          htmlOutput("hiddenGeneADiv"),
          htmlOutput("hiddenGeneBDiv"),
          htmlOutput("hiddenGeneA_id_Div"),
          htmlOutput("hiddenGeneB_id_Div")
          ),
       mainPanel(
          tabsetPanel(
               tabPanel('interactions', rowSelectableDataTableOutput(outputId="table")),
               #tabPanel('rowSelTbl', selectableDataTableOutput(outputId="rowSelTable",...)),
               tabPanel("Pubmed Abstract", htmlOutput("pubmedAbstract")),
               tabPanel("A", htmlOutput("geneA")),
               tabPanel("B", htmlOutput("geneB"))
             ),
          HTML(javascript())
          )
       ) # uiWidgets

#----------------------------------------------------------------------------------------------------
serverFunction <- function(input, output, session)
{
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


    observe({
       if(!is.na(initial.genes))
          updateTextInput(session, inputId="genes", value=initial.genes)
       })

    findInteractions <- reactive({
       input$goButton
       printf("reactive findInteractions running")
       genes <- cleanGenes(isolate(input$genes))
       printf("       genes: %s", paste(genes, collapse=","))
       if(length(genes) == 0)
           return(empty.data.frame)
       if(nchar(genes) == 0)
           return(empty.data.frame)
       current.providers <- isolate(input$providers)
       printf("current.providers: %s", paste(current.providers, collapse=","))
       printf("   providers: %s", paste(current.providers, collapse=","))
       printf("querying for %s from %s", paste(genes, collapse=","), paste(current.providers, collapse=","))
       tbl <- interactions(refnet, id=genes, provider=current.providers, quiet=FALSE, species="9606")
       printf("dim(tbl): %d x %d", nrow(tbl), ncol(tbl))
       if(nrow(tbl) == 0){
          printf("nrow 0, returning empty data frame")
          return(empty.data.frame)
          }
       printf("addGeneInfo...")
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
    if(!all(is.na(id)))
        initial.genes <<- id
    
    if(!all(is.na(provider)))
        initial.providers <<- provider
    
    app <- list(ui=uiWidgets, server=serverFunction)
    printf("about to runApp, initial.genes: %s", paste(initial.genes, collapse=","))


    runApp(app)

} # .explore
#----------------------------------------------------------------------------------------------------
