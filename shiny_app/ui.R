library(shiny)
library(plotly)
library(shinyWidgets)

df <- readRDS( "journal_ORs.rds")
topics_list <- readRDS("topics_list.rds")

shinyUI(
  fluidPage(
    
    titlePanel("eFigure 8: Estimated journal-specific odds ratio vs journal Cite Score"),
    
    fluidRow(
    
    #sidebarLayout(
      #sidebarPanel(
      column(4,
        wellPanel(   
        selectInput("filter",
                    "Filter by:",
                    choices = c("No filter","Journal","Topic"),
                    multiple = FALSE,
                    selected= "No filter"
        ),
        conditionalPanel(
          condition = "input.filter == 'Journal'",
          pickerInput("which_journals",
                      "Select journals to plot:",
                      choices = sort(as.character(unique(df$sourcetitle))),
                      multiple = TRUE,
                      selected= unique(as.character(df$sourcetitle)),
                      options = pickerOptions(actionsBox = TRUE,
                                     liveSearch = TRUE)
                    )
        ),
        conditionalPanel(
        condition = "input.filter == 'Topic'",
          pickerInput("which_topics",
                      "Select topics to plot:",
                      choices = topics_list,
                      multiple = TRUE,
                      selected= topics_list,
                      options = pickerOptions(actionsBox = TRUE,
                                            liveSearch = TRUE)
          )
        ),
        numericInput("n_min", "Minimum sample size (number of matched sets):", min(df$n_cases), 
                     min = min(df$n_cases), max = max(df$n_cases)),
        selectInput("adjusted",
                    "Adjust for years active, h-index, and number of publications?",
                    choices = c("Yes","No"),
                    multiple = FALSE,
                    selected = "No")
      ),
      downloadButton(
        "metaData",
        "Download all data (eTable 6)"
      )
      ),
      
      column(8,
        plotlyOutput("scatterplot", height=600)
      )
    )
  )
)