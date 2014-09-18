printf("entering server.R");
library(shiny)
library(datasets)
# library(RefNet)

refnet <- RefNet()
idMapper <- IDMapper(species="9606")
#load("tbl.alk.RData")

cleanGenes <- function(s){
   genes.raw <- sub("[[:space:]]+$","", s)
   toupper(strsplit(genes.raw, " +")[[1]])
   }


shinyServer(function(input, output) {

  findInteractions <- reactive({
      genes <- cleanGenes(input$genes)
      #tbl <- interactions
      if(length(genes) == 0)
            # preserve the column structure, even when no rows
          df <- subset(tbl, A.name == "completely bogus entry")
      else if(length(genes) == 1)
         df <- subset(tbl, A.name %in% genes | B.name %in% genes)
      else if(length(genes) > 1)
         df <- subset(tbl, A.name %in% genes & B.name %in% genes & A.name != B.name)
      if(nrow(df) > 0)
         if(! "ANY" %in% input$providers)
            df <- subset(df, provider %in% input$providers)
      df
      })
  
  output$summary <- renderPrint({
     #dataset <- datasetInput()
     #summary(dataset)
     output$value <- renderPrint({ input$action })
     cleanGenes(input$genes)
     })
  
  datatable.options <- list(bFilter=TRUE)
  output$table <- renderDataTable(findInteractions(),options=datatable.options)
})
