library(shiny)

load("tbl.alk.RData")

# Define the overall UI
shinyUI(
  fluidPage(
    theme="bootstrap.css",
    titlePanel("RefNet"),
          
    # Create a new Row in the UI for selectInputs
    fluidRow(
       column(2, selectInput("provider", "Provider:", 
                             c("All", sort(unique(as.character(tbl$provider)))))),
       column(2, selectInput("A.name", "A.name:",
                             c("All", sort(unique(as.character(tbl$A.name)))))),
       column(2, selectInput("B.name", "B.name:",
                             c("All", sort(unique(as.character(tbl$B.name)))))),
       column(4, selectInput("type", "Type:",
                             c("All", sort(unique(as.character(tbl$type)))))),
       column(4, selectInput("detectionMethod", "Detection Method",
                             c("All", sort(unique(as.character(tbl$detectionMethod)))))),
       column(4, selectInput("publicationID", "PublicationID",
                             c("All", sort(unique(as.character(tbl$publicationID))))))
      
    ),
    # Create a new row for the table.
    fluidRow(dataTableOutput(outputId="table"))
   ) # fluidPage
) # shinyUI
