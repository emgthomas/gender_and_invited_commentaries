library(shiny)
library(plotly)
library(dplyr)
library(metafor)
library(splines)

# df <- readRDS( "/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/shiny_app/journal_ORs.rds")
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
        out <- out[rowSums(out[,input$which_topics]) > 0,]
      }
      
      # # meta-regression of log ORs on citescore
      # if(input$citescore == "Yes"){
      #   # choose internal knots
      #   n_knots <- input$numberOfKnots
      #   if(n_knots > 0){
      #     citescore_by_pub <- unlist(sapply(1:nrow(out),function(idx,n_cases,citescore) rep(citescore[idx],n_cases[idx]),
      #                                       n_cases=out$n_cases_plot,citescore=out$citescore))
      #     knot_placement <- quantile(citescore_by_pub,probs = seq(0,1,length.out=n_knots+2)[2:(n_knots+1)])
      #   } else {
      #     knot_placement <- NULL
      #   }
      #   # run meta regression
      #   meta_analysis_cs <- rma.mv(effect_plot,sd_plot^2,random=~1| journal,
      #                              mods= ~ ns(citescore,knots=knot_placement),
      #                              data=out)
      #   # get predicted values for plotting
      #   citescore_newmods <- seq(0,max(out$citescore),0.1)
      #   citescore_bs <- ns(citescore_newmods,knots=knot_placement)
      #   plot_citescore <- predict(meta_analysis_cs,newmods=citescore_bs[1:nrow(citescore_bs),1:ncol(citescore_bs)],
      #                             transf=exp)
      #   plot_citescore <- plot_citescore[,c("pred","cr.lb","cr.ub")]
      #   plot_citescore$citescore <- citescore_newmods
      # } else {
      #   plot_citescore <- NA
      # }
      
      # return(list(plot_df=out,plot_citescore=plot_citescore))
      
      return(list(plot_df=out))
      
    })
    
    output$scatterplot <- renderPlotly({
      
      out <- scatterplot_data()
      
      OR_plot <- plot_ly(out$plot_df, x = ~citescore, y= ~OR_plot) %>%
      # add horizontal line for null value
      add_trace(x = c(0,18), y= c(1,1), mode = "lines") %>%
      # add point for OR of each journal
      add_trace(y=~OR_plot, type='scatter',mode='markers',
                size = ~node_size_plot*2, 
                color= ~pval_plot,
                marker=list(sizeref=0.2),
                text = ~paste("Journal: ", sourcetitle,
                            '<br>Cite Score: ', sprintf(citescore, fmt="%.2f"),
                            "<br>Number of cases: ", n_cases_plot,
                            '<br>OR: ',sprintf(OR_plot, fmt="%.2f"),
                            ' (',sprintf(ci_lower_plot, fmt="%.2f"),
                            ',',sprintf(ci_upper_plot, fmt="%.2f"),')',sep=""),
                hoverinfo='text'
                ) %>%
        colorbar(title="P-value",limits=c(0,1),cmin=0,cmax=1,
                 colorscale="Viridis") %>%
        # layout
        layout(yaxis = list(title="Odds Ratio (log scale)",range=c(-log(3),log(3)),type="log",
                            tickvals=c(1/8,1/4,1/2,1,2,4,8),
                            ticktext=c("1/8","1/4","1/2",
                                       as.character(c(1,2,4,8)))),
               xaxis = list(title="Journal Cite Score",range=c(0,18),tickmode="array"),
               showlegend=F,
               font=list(size=16)
        ) %>%
        # annotations
        add_trace(x = 16, y = 9, mode="text",text="Favors women",
                  hoverinfo="none",textfont=list(size=20,color=1)) %>%
        add_trace(x = 16, y = 1/9, mode="text",text="Favors men",
                  hoverinfo="none",textfont=list(size=20,color=1))
      
      # # add trace for citescore, if requested
      # if(!is.na(out$plot_citescore)){
      #   OR_plot <-  add_ribbons(OR_plot,x =out$plot_citescore$`citescore,y=out$plot_citescore$pred,
      #                 ymin =out$plot_citescore$cr.lb, ymax =out$plot_citescore$cr.ub,
      #                 color = I("gray80"), 
      #                 hoverinfo="text", text="95% prediction interval") %>%
      #     add_lines(y=out$plot_citescore$pred, x=out$plot_citescore$citescore,
      #               hoverinfo="text", text="Predicted OR as function of citescore",
      #               line = list(color="orange"))
      # }
      
      # plot
      OR_plot
      
    })
    
  }
)

