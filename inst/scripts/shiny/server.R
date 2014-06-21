library(shiny)

load("tbl.alk.RData")

# Define a server for the Shiny app
shinyServer(function(input, output) {
  
  # Filter data based on selections
  output$table <- renderDataTable({
    if(input$provider != "All")         {tbl <- tbl[tbl$provider == input$provider,]}
    if(input$A.name != "All")           {tbl <- tbl[tbl$A.name == input$A.name,]}
    if(input$B.name != "All")           {tbl <- tbl[tbl$B.name == input$B.name,]}
    if(input$type != "All")             {tbl <- tbl[tbl$type == input$type,]}
    if(input$detectionMethod != "All")  {tbl <- tbl[tbl$detectionMethod == input$detectionMethod,]}
    if(input$publicationID != "All")    {tbl <- tbl[tbl$publicationID == input$publicationID,]}
    tbl
  })
  
})
