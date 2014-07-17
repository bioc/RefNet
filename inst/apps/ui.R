library(shiny)
load("tbl.alk.RData")
providers <- c("ANY", sort(unique(tbl$provider)))

shinyUI(fluidPage(theme="bootstrap.css",
  titlePanel("RefNet"),
  sidebarLayout(
     sidebarPanel(width=2,
        #actionButton("initRefNetButton", "Initialize RefNet"),
        selectizeInput(inputId='providers', label='providers',
                       choices = providers, multiple = TRUE, selected=providers[1]),
        textInput("genes", "Genes:", "ALK"),
        submitButton("Find Interactions")
        ),
    mainPanel(
        verbatimTextOutput("summary"),
        dataTableOutput(outputId="table")
    )
  )
))
