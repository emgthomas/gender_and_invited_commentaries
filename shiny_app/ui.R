library(shiny)
library(plotly)
library(shinyWidgets)

# df <- readRDS( "/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/shiny_app/journal_ORs.rds")
df <- readRDS( "journal_ORs.rds")
topics_list <- readRDS("topics_list.rds")

shinyUI(
  fluidPage(
    
    titlePanel("Odds ratio by journal"),
    
    fluidRow(
    
    #sidebarLayout(
      #sidebarPanel(
      column(4,
        wellPanel(   
        # selectInput("which_journals",
        #             "Select journals to plot:",
        #             choices = sort(unique(df$journal)),
        #             multiple = TRUE,
        #             selected = unique(df$journal)),
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
        # sliderInput("n_range",
        #             "Select minimum sample size (number of matched sets)",
        #             min = 10,
        #             max = max(df$n_cases),
        #             value = c(10, max(df$n_cases)),
        #             sep = ""),
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
                    "Adjust for years active, h-index and number of publications?",
                    choices = c("Yes","No"),
                    multiple = FALSE,
                    selected = "No")
        # selectInput("citescore",
        #             "Fit meta-regression on journal cite score?",
        #             choices = c("Yes","No"),
        #             multiple = FALSE,
        #             selected = "No"),
        # # Only show this panel if including meta-regression
        # conditionalPanel(
        #   condition = "input.citescore == 'Yes'",
        #   numericInput(
        #     "numberOfKnots", "Number of internal knots for B-spline basis",
        #     3,min=0,max=10))
      )
      ),
      
      column(8,
      #mainPanel(
        plotlyOutput("scatterplot", height=600)
      )
    )
  )
)