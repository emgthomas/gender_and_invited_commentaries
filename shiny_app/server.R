library(shiny)
library(plotly)
library(dplyr)
library(metafor)
library(splines)

df <- readRDS( "journal_ORs.rds")

shinyServer(
  function(input, output) {
    
    scatterplot_data <- reactive({
      
      out <- df
      
      # Adjust for author seniority?
      if(input$adjusted=="No"){
        out$OR_plot <- out$OR
        out$n_cases_plot <- out$n_cases
        out$pval_plot <- out$pval
        out$ci_lower_plot <- out$ci_lower
        out$ci_upper_plot <- out$ci_upper
        out$node_size_plot <- out$node_size
        out$effect_plot <- out$effect
        out$sd_plot <- out$sd
      } else {
        out$OR_plot <- out$OR_adj
        out$n_cases_plot <- out$n_cases_adj
        out$pval_plot <- out$pval_adj
        out$ci_lower_plot <- out$ci_lower_adj
        out$ci_upper_plot <- out$ci_upper_adj
        out$node_size_plot <- out$node_size_adj
        out$effect_plot <- out$effect_adj
        out$sd_plot <- out$sd_adj
      }
      
      # Filter on min number of matched sets
      out <- filter(out,n_cases_plot >= input$n_min)
      
      # Filter on journal *or* journal topic
      if(input$filter == "Journal"){
        out <- filter(out,sourcetitle %in% input$which_journals)
      }
      
      if(input$filter == "Topic"){
        out <- out[rowSums(out[,input$which_topics,drop=F]) > 0,]
      }
      
      if(nrow(out) == 0){
        return(out)
      }
      
      out$hovertext <- paste("<i>",out$sourcetitle,"</i>",
                             '<br>OR: ',sprintf(out$OR_plot, fmt="%.2f"),
                             ' (',sprintf(out$ci_lower_plot, fmt="%.2f"),
                             ',',sprintf(out$ci_upper_plot, fmt="%.2f"),')',
                             '<br>Cite Score: ', sprintf(out$citescore, fmt="%.2f"),
                             "<br>Unique ICC authors: ", out$n_cases_plot,
                             "<br>Number of ICCs ",
                             "<br>2013: ", out$npubs.2013,
                             "<br>2014: ", out$npubs.2014,
                             "<br>2015: ", out$npubs.2015,
                             "<br>2016: ", out$npubs.2016,
                             "<br>2017: ", out$npubs.2017,
                             sep="")
      
      return(out)
      
    })
    
    output$scatterplot <- renderPlotly({
      
      out <- scatterplot_data()
      
      # basic plot (just axes)
      plt <- plot_ly(out) %>%
        # add horizontal line for null value
        add_lines(x = c(0,18), y= c(1,1),
                  color=I("black"), 
                  hoverinfo="none") %>%
        # annotations
        add_trace(x = 16, y = 9, mode="text",text="Favors women",
                  hoverinfo="none",textfont=list(size=20,color=1),
                  marker=list(opacity=0)) %>%
        add_trace(x = 16, y = 1/9, mode="text",text="Favors men",
                  hoverinfo="none",textfont=list(size=20,color=1),
                  marker=list(opacity=0)) %>%
        # layout
        layout(yaxis = list(title="Odds Ratio (log scale)",range=c(-log(3),log(3)),type="log",
                            tickvals=c(1/8,1/4,1/2,1,2,4,8),
                            ticktext=c("1/8","1/4","1/2",
                                       as.character(c(1,2,4,8)))),
               xaxis = list(title="Journal Cite Score",range=c(0,18),tickmode="array"),
               showlegend=F,
               font=list(size=16),
               hovermode="closest"
        )
      
      
      if(nrow(out)>1){
        
        # add circles for each journal
        plt <- add_trace(plt,
                         x=~citescore,y=~OR_plot,
                         type='scatter',
                         mode='markers',
                         size = ~node_size_plot*2, 
                         color= ~pval_plot,
                         marker=list(sizeref=0.2,
                                     opacity=0.6),
                         hoverinfo="text",
                         hovertext=~hovertext
        ) 
      } else if(nrow(out)==1) {
        # add circle for single journal
        plt <- add_trace(plt,
                         x=c(out$citescore,18),y=c(out$OR_plot,8),
                         type='scatter',
                         mode='markers',
                         size=c(out$node_size,1), 
                         color=c(out$pval_plot,0.5),
                         marker=list(opacity=c(0.6,0),
                                     sizeref=0.2),
                         hoverinfo=c("text","skip"),
                         hovertext=c(out$hovertext,"")
        )
      } else {
        # add fake data to get colorbar
        plt <- add_trace(plt,
                         x=c(0,18),y=c(0.5,2),
                         type='scatter',
                         mode='markers',
                         color=c(0,1),
                         marker=list(opacity=0),
                         hoverinfo="none"
        )
      }
      
      plt <- colorbar(plt,title="P-value",
                      limits=c(0,1),
                      cmin=0,cmax=1,
                      colorscale="Viridis")
      
      plt
      
    })
    
  }
)

