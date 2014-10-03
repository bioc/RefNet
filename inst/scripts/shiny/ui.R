library(shiny)

load("tbl.alk.RData")

# Define the overall UI
data.table.options=list(iDisplayLength=-1,                    # initial number of records
                        aLengthMenu=c(5,10),                  # records/page options
                        bLengthChange=0,                       # show/hide records per page dropdown
                        bFilter=0,                                    # global search box on/off
                        bInfo=0,                                      # information on/off (how many records filtered, etc)
                        bAutoWidth=0                            # automatic column width calculation, disable if passing column width via aoColumnDefs
                        )

data.table <- dataTableOutput(outputId="table)") # , options=data.table.options)
shinyUI(
  fluidPage(
      titlePanel("RefNet"),
      #theme="bootstrap.css",
      mainPanel(tabsetPanel(
          tabPanel("Interactions", data.table),
          tabPanel("Genes", verbatimTextOutput("genes")),
          tabPanel("Providers", verbatimTextOutput("providers")),
          tabPanel("Interaction Types", verbatimTextOutput("types")),
          tabPanel("Detection Methods", verbatimTextOutput("methods"))
          ))
      ))
          
#    # Create a new Row in the UI for selectInputs
#    fluidRow(
#       column(2, selectInput("provider", "Provider:", 
#                             c("All", sort(unique(as.character(tbl$provider)))))),
#       column(2, selectInput("A.name", "A.name:",
#                             c("All", sort(unique(as.character(tbl$A.name)))))),
#       column(2, selectInput("B.name", "B.name:",
#                             c("All", sort(unique(as.character(tbl$B.name)))))),
#       column(4, selectInput("type", "Type:",
#                             c("All", sort(unique(as.character(tbl$type)))))),
#       column(4, selectInput("detectionMethod", "Detection Method",
#                             c("All", sort(unique(as.character(tbl$detectionMethod)))))),
#       column(4, selectInput("publicationID", "PublicationID",
#                             c("All", sort(unique(as.character(tbl$publicationID))))))
#      
#    ),
#    # Create a new row for the table.
#    fluidRow(dataTableOutput(outputId="table"))
#   ) # fluidPage
 # shinyUI
