library(shiny)
library(plotly)
library(shinyWidgets)

df <- readRDS( "journal_ORs.rds")
topics_list <- readRDS("topics_list.rds")

shinyUI(
  fluidPage(
    
    titlePanel("eAppendix: Estimated journal-specific odds ratio vs journal Cite Score"),
    
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
      h5("eAppendix Data (Excel format)"),
      downloadButton(
        "metaData",
        "Download"
      ),
      h5("Figure caption"),
      p("Each circle represents the odds ratio estimated for a single journal. Circle diameter is inversely proportional to the standard error of the log odds ratio estimate.",
        "Circle color represents the p-value for the null hypothesis of no association between gender and invited commentary authorship for that journal.",
        "Hover over the circle for more information about that journal.",
        "All odds ratios are adjusted for field of expertise through matching.",
        "The number of matched sets (sample size) for each journal is given by the number of unique intra-citing commentary (ICC) authors.",
        "The number of ICCs included is the number of ICCs with a corresponding author whose gender could be inferred.",
        "There may be more ICCs than unique ICC authors due to multiple ICCs by the same author.",
        "Data for all journals can be downloaded via the link above.",
        style = "font-size: 10pt")
      ),
      
      column(8,
        plotlyOutput("scatterplot", height=600)
      )
    )
  )
)