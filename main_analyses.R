# ------------------------------------------------ #
# --------------- Main analyses ------------------ #
# ------------------------------------------------ #

#######################################################################
sink(file="./results/main_analyses.txt")
#######################################################################

# packages
require(survival)
require(splines)
require(forestplot)
require(reshape2)
require(ggplot2)

# global functions
source("./code/functions.R")

## load data
icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

cat("\n\n------------ One-stage meta-analysis, all journals ----------------\n\n")

# icc_df_case <- icc_df[icc_df$case==1,]
# pubs.incl <- unique(icc_df_case$pub_id[!duplicated(icc_df_case$auth_id)])
# icc_df_dedup <- icc_df[icc_df$pub_id %in% pubs.incl]
# icc_df_dedup$pub_id <- factor(icc_df_dedup$pub_id,levels=unique(icc_df_dedup$pub_id))
all_1stage <- clogit(case ~ Gender + strata(pub_id), data = icc_df)
cat("---Unadjusted analysis---\n")
summary(all_1stage)

cat("\n\n---Adjusted for measures of seniority using natural cubic splines---\n")
knots <- c(2.5,5,7.5)
all_1stage_adj <- clogit(case ~ Gender + ns(years_in_scopus_ptile,knots=knots) + 
                           ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                           strata(pub_id), data = icc_df)
summary(all_1stage_adj)

cat("\n\n---Adjusted for measures of seniority as linear terms---\n")

all_1stage_adj <- clogit(case ~ Gender + years_in_scopus_ptile + 
                           h_index_ptile + n_pubs_ptile + 
                           strata(pub_id), data = icc_df)
summary(all_1stage_adj)

## Save key results for use in plotting ##
save(all_1stage,all_1stage_adj,file="./results/main_analyses.Rdata")

###### Plots of continuous variable effects ######

# ------------------- Years active ----------------------#
# set up data for prediction
n.points <- 1000
ptile <- seq(0,10,length.out = n.points)
ptile_ns <- ns(ptile,knots=knots)
dat_pred <- ptile_ns[1:nrow(ptile_ns),1:ncol(ptile_ns)]
# Get predicted values
varnames <- c("ns(years_in_scopus_ptile, knots = knots)")
pred_years_in_scopus <- spline_predictions(all_1stage_adj,dat_pred,varnames)

# ------------------- H-index ----------------------#
# Get predicted values
varnames <- c("ns(h_index_ptile, knots = knots)")
pred_h_index <- spline_predictions(all_1stage_adj,dat_pred,varnames)

# ------------------- Number of publications ----------------------#
# Get predicted values
varnames <- "ns(n_pubs_ptile, knots = knots)"
pred_n_pubs <- spline_predictions(all_1stage_adj,dat_pred,varnames)

# Put into data frame for plotting
plot_effects <- cbind(pred_years_in_scopus,pred_h_index,pred_n_pubs)
names(plot_effects) <- c("pred1","ci.lb1","ci.ub1",
                         "pred2","ci.lb2","ci.ub2",
                         "pred3","ci.lb3","ci.ub3")
plot_effects$ptile <- ptile*10

# Plot
cols3 <- brewer.pal(3,name="Accent")
# Plot
OR_plot <- plot_ly(plot_effects, x = ~ptile, y= ~pred1, height=600,width=1000) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred1, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[1]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred1,
              ymin=plot_effects$ci.lb1, ymax=plot_effects$ci.ub1,
              hoverinfo="none",
              color=I(cols3[1]),
              name="Years active") %>%
  add_lines(y=plot_effects$pred2, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[2]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred2,
              ymin=plot_effects$ci.lb2, ymax=plot_effects$ci.ub2,
              hoverinfo="none",
              color=I(cols3[2]),
              name="H-index") %>%
  add_lines(y=plot_effects$pred3, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[3]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred3,
              ymin=plot_effects$ci.lb3, ymax=plot_effects$ci.ub3,
              hoverinfo="none",
              color=I(cols3[3]),
              name="Number of publications") %>%
  add_trace(x = 50, y = 3, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 50, y = 1/3, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.1),log(2.1)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile",range=c(0,100),tickmode="array"),
         showlegend=T,
         font=list(size=16)
  )

export(OR_plot, "./results/OR_by_vars.pdf")

cat("\n\n------------ Effect modification by author-level variables ----------------\n\n")

# -------------------- Years active --------------------- #
cat("\n\n***Effect modification by years active***\n\n")
knots <- c(2.5,5,7.5)
icc_df$gender <- as.numeric(icc_df$Gender == "female")
icc_df$gender_years_in_scopus_ptile <- icc_df$gender*icc_df$years_in_scopus_ptile
all_1stage_EM_YiS <- clogit(case ~ gender + gender_years_in_scopus_ptile + 
                              ns(years_in_scopus_ptile,knots=knots) + 
                              ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                              strata(pub_id), data = icc_df)
summary(all_1stage_EM_YiS)

# set up data for prediction
n.points <- 1000
ptile <- seq(0,10,length.out = n.points)
ptile_ns <- ns(ptile,knots=knots)
dat_pred_f <- cbind(ptile,ptile_ns[1:nrow(ptile_ns),1:ncol(ptile_ns)])
# Get effects of years active for women
varnames_f <- c("gender_years_in_scopus_ptile",
                "ns(years_in_scopus_ptile, knots = knots)")
pred_f <- spline_predictions(all_1stage_EM_YiS,dat_pred_f,varnames_f)
# Get effect of years active for men
dat_pred_m <- cbind(ptile_ns[1:nrow(ptile_ns),1:ncol(ptile_ns)])
varnames_m <- c("ns(years_in_scopus_ptile, knots = knots)")
pred_m <- spline_predictions(all_1stage_EM_YiS,dat_pred_m,varnames_m)
# Get effect of h-index
varnames <- "ns(h_index_ptile, knots = knots)"
pred_h_index <- spline_predictions(all_1stage_EM_YiS,dat_pred,varnames)
# Get effect of number of pubs
varnames <- "ns(n_pubs_ptile, knots = knots)"
pred_n_pubs <- spline_predictions(all_1stage_EM_YiS,dat_pred,varnames)
# put in data frame
plot_effects <- cbind(pred_f,pred_m,pred_h_index,pred_n_pubs)
names(plot_effects) <- c("pred_f","ci.lb_f","ci.ub_f",
                         "pred_m","ci.lb_m","ci.ub_m",
                         "pred_h_index","ci.lb_h_index","ci.ub_h_index",
                         "pred_n_pubs","ci.lb_n_pubs","ci.ub_n_pubs")
plot_effects$ptile <- ptile*10

cols4 <- brewer.pal(4,name="Set1")
m <- list(
  l = 100,
  r = 50,
  b = 100,
  t = 100,
  pad = 4
)
OR_plot_int <- plot_ly(plot_effects, x = ~ptile, y= ~pred_f, height=700,width=600) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred_f, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols4[1]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_f,
              ymin=plot_effects$ci.lb_f, ymax=plot_effects$ci.ub_f,
              hoverinfo="none",
              color=I(cols4[1]),
              name="Female") %>%
  add_lines(y=plot_effects$pred_m, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols4[2]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_m,
              ymin=plot_effects$ci.lb_m, ymax=plot_effects$ci.ub_m,
              hoverinfo="none",
              color=I(cols4[2]),
              name="Male") %>%
  add_trace(x = 35, y = 3.5, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 35, y = 1/3.5, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile of Years Active",range=c(0,100),tickmode="array"),
         showlegend=T,
         legend = list(x = 0.3, y = -0.3),
         margin=m,
         font=list(size=16)
  )

export(OR_plot_int, "./results/OR_years_active_int1.pdf")

# Plot
OR_plot <- plot_ly(plot_effects, x = ~ptile, y= ~pred1, height=700,width=600) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred_h_index, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[2]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_h_index,
              ymin=plot_effects$ci.lb_h_index, ymax=plot_effects$ci.ub_h_index,
              hoverinfo="none",
              color=I(cols3[2]),
              name="H-index") %>%
  add_lines(y=plot_effects$pred_n_pubs, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[3]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_n_pubs,
              ymin=plot_effects$ci.lb_n_pubs, ymax=plot_effects$ci.ub_n_pubs,
              hoverinfo="none",
              color=I(cols3[3]),
              name="Number of publications") %>%
  add_trace(x = 35, y = 3.5, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 35, y = 1/3.5, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile",range=c(0,100),tickmode="array"),
         showlegend=T,
         legend = list(x = 0.3, y = -0.3),
         margin=m,
         font=list(size=16)
  )

export(OR_plot, "./results/OR_years_active_int2.pdf")


# -------------------- H-index --------------------- #
icc_df$gender_h_index_ptile <- icc_df$gender*icc_df$h_index_ptile
all_1stage_EM_HI <- clogit(case ~ gender + gender_h_index_ptile + 
                             ns(years_in_scopus_ptile,knots=knots) + 
                             ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) +
                             strata(pub_id), data = icc_df)

cat("\n\n***Effect modification by-H Index***\n\n")
summary(all_1stage_EM_HI)

# Get effects of h-index for women
varnames_f <- c("gender_h_index_ptile",
                "ns(h_index_ptile, knots = knots)")
pred_f <- spline_predictions(all_1stage_EM_HI,dat_pred_f,varnames_f)
# Get effect of h-index for men
varnames_m <- c("ns(h_index_ptile, knots = knots)")
pred_m <- spline_predictions(all_1stage_EM_HI,dat_pred_m,varnames_m)
# Get effect of years in scopus
varnames <- "ns(years_in_scopus_ptile, knots = knots)"
pred_years_in_scopus <- spline_predictions(all_1stage_EM_HI,dat_pred,varnames)
# Get effect of number of pubs
varnames <- "ns(n_pubs_ptile, knots = knots)"
pred_n_pubs <- spline_predictions(all_1stage_EM_HI,dat_pred,varnames)
# put in data frame
plot_effects <- cbind(pred_f,pred_m,pred_years_in_scopus,pred_n_pubs)
names(plot_effects) <- c("pred_f","ci.lb_f","ci.ub_f",
                         "pred_m","ci.lb_m","ci.ub_m",
                         "pred_years_in_scopus","ci.lb_years_in_scopus","ci.ub_years_in_scopus",
                         "pred_n_pubs","ci.lb_n_pubs","ci.ub_n_pubs")
plot_effects$ptile <- ptile*10

OR_plot_int <- plot_ly(plot_effects, x = ~ptile, y= ~pred_f, height=700,width=600) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred_f, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols4[1]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_f,
              ymin=plot_effects$ci.lb_f, ymax=plot_effects$ci.ub_f,
              hoverinfo="none",
              color=I(cols4[1]),
              name="Female") %>%
  add_lines(y=plot_effects$pred_m, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols4[2]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_m,
              ymin=plot_effects$ci.lb_m, ymax=plot_effects$ci.ub_m,
              hoverinfo="none",
              color=I(cols4[2]),
              name="Male") %>%
  add_trace(x = 35, y = 3.5, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 35, y = 1/3.5, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile of H-Index",range=c(0,100),tickmode="array"),
         showlegend=T,
         legend = list(x = 0.3, y = -0.3),
         margin=m,
         font=list(size=16)
  )

export(OR_plot_int, "./results/OR_h_index_int1.pdf")

# Plot
OR_plot <- plot_ly(plot_effects, x = ~ptile, y= ~pred_years_in_scopus, height=700,width=600) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred_h_index, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[1]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_years_in_scopus,
              ymin=plot_effects$ci.lb_years_in_scopus, ymax=plot_effects$ci.ub_years_in_scopus,
              hoverinfo="none",
              color=I(cols3[1]),
              name="Years active") %>%
  add_lines(y=plot_effects$pred_n_pubs, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[3]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_n_pubs,
              ymin=plot_effects$ci.lb_n_pubs, ymax=plot_effects$ci.ub_n_pubs,
              hoverinfo="none",
              color=I(cols3[3]),
              name="Number of publications") %>%
  add_trace(x = 35, y = 3.5, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 35, y = 1/3.5, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile",range=c(0,100),tickmode="array"),
         showlegend=T,
         legend = list(x = 0.3, y = -0.3),
         margin=m,
         font=list(size=16)
  )

export(OR_plot, "./results/OR_h_index_int2.pdf")

# -------------------- Number of pubs --------------------- #
icc_df$gender_n_pubs_ptile <- icc_df$gender*icc_df$n_pubs_ptile
all_1stage_EM_npubs <- clogit(case ~ gender + gender_n_pubs_ptile + 
                                ns(years_in_scopus_ptile,knots=knots) + 
                                ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) +
                                strata(pub_id), data = icc_df)


cat("\n\n***Effect modification by number of publications***\n\n")
summary(all_1stage_EM_npubs)

# Get effects of number of pubs for women
varnames_f <- c("gender_n_pubs_ptile",
                "ns(n_pubs_ptile, knots = knots)")
pred_f <- spline_predictions(all_1stage_EM_npubs,dat_pred_f,varnames_f)
# Get effect of number of pubs for men
varnames_m <- c("ns(n_pubs_ptile, knots = knots)")
pred_m <- spline_predictions(all_1stage_EM_npubs,dat_pred_m,varnames_m)
# Get effect of years in scopus
varnames <- "ns(years_in_scopus_ptile, knots = knots)"
pred_years_in_scopus <- spline_predictions(all_1stage_EM_npubs,dat_pred,varnames)
# Get effect of h-index
varnames <- "ns(h_index_ptile, knots = knots)"
pred_h_index <- spline_predictions(all_1stage_EM_npubs,dat_pred,varnames)
# put in data frame
plot_effects <- cbind(pred_f,pred_m,pred_years_in_scopus,pred_h_index)
names(plot_effects) <- c("pred_f","ci.lb_f","ci.ub_f",
                         "pred_m","ci.lb_m","ci.ub_m",
                         "pred_years_in_scopus","ci.lb_years_in_scopus","ci.ub_years_in_scopus",
                         "pred_h_index","ci.lb_h_index","ci.ub_h_index")
plot_effects$ptile <- ptile*10

OR_plot_int <- plot_ly(plot_effects, x = ~ptile, y= ~pred_f, height=700,width=600) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred_f, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols4[1]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_f,
              ymin=plot_effects$ci.lb_f, ymax=plot_effects$ci.ub_f,
              hoverinfo="none",
              color=I(cols4[1]),
              name="Female") %>%
  add_lines(y=plot_effects$pred_m, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols4[2]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_m,
              ymin=plot_effects$ci.lb_m, ymax=plot_effects$ci.ub_m,
              hoverinfo="none",
              color=I(cols4[2]),
              name="Male") %>%
  add_trace(x = 35, y = 3.5, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 35, y = 1/3.5, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile of number of publications",range=c(0,100),tickmode="array"),
         showlegend=T,
         legend = list(x = 0.3, y = -0.3),
         margin=m,
         font=list(size=16)
  )

export(OR_plot_int, "./results/OR_n_pubs_int1.pdf")

# Plot
OR_plot <- plot_ly(plot_effects, x = ~ptile, y= ~pred_years_in_scopus, height=700,width=600) %>%
  # add horizontal line for null value
  add_trace(x = c(0,100), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  add_lines(y=plot_effects$pred_n_pubs, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[1]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_years_in_scopus,
              ymin=plot_effects$ci.lb_years_in_scopus, ymax=plot_effects$ci.ub_years_in_scopus,
              hoverinfo="none",
              color=I(cols3[1]),
              name="Years active") %>%
  add_lines(y=plot_effects$pred_h_index, x=plot_effects$ptile,
            hoverinfo="none",
            color=I(cols3[2]),
            showlegend=F) %>%
  add_ribbons(x=plot_effects$ptile,y=plot_effects$pred_h_index,
              ymin=plot_effects$ci.lb_h_index, ymax=plot_effects$ci.ub_h_index,
              hoverinfo="none",
              color=I(cols3[2]),
              name="H-index") %>%
  add_trace(x = 35, y = 3, mode="text",text="More invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 35, y = 1/3, mode="text",text="Fewer invited commentaries",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  # layout
  layout(yaxis = list(title="Odds ratio relative to 0th percentile (log scale)",range=c(-log(2.1),log(2.1)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Percentile",range=c(0,100),tickmode="array"),
         showlegend=T,
         legend = list(x = 0.3, y = -0.3),
         margin=m,
         font=list(size=16)
  )

export(OR_plot, "./results/OR_n_pubs_int2.pdf")

# -------------------- Plot gender effect by years active --------------------- #
deciles <- c(1,3,5,7,9)

OR_df <- data.frame(ptile_YiS=c("10th","30th","50th","70th","90th"),
                    OR=numeric(length=5),ci.lb=numeric(length=5),ci.ub=numeric(length=5))

OR_df$OR <- exp(all_1stage_EM_YiS$coefficients[1] + deciles*all_1stage_EM_YiS$coefficients[14])

logOR_se <- sqrt(all_1stage_EM_YiS$var[1,1] + deciles^2*all_1stage_EM_YiS$var[14,14] + 2*deciles*all_1stage_EM_YiS$var[1,14])

OR_df$ci.lb <- OR_df$OR*exp(-1.96*logOR_se)
OR_df$ci.ub <- OR_df$OR*exp(1.96*logOR_se)

tablehead <- rbind(c("Percentile","Years Active","OR (95%CI)"),
                   rep(NA,3))

unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>0.9 & icc_df$years_in_scopus_ptile<1.1])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>2.8 & icc_df$years_in_scopus_ptile<3.2])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>3.8 & icc_df$years_in_scopus_ptile<4.2])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>6.9 & icc_df$years_in_scopus_ptile<7.1])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>8.94 & icc_df$years_in_scopus_ptile<9.06])
years <- c(8,14,16,27,38)
tablenum <- cbind(as.character(OR_df$ptile_YiS),
                  years,
                  sapply(1:nrow(OR_df),formatting_fun,or=OR_df$OR,ci.lb=OR_df$ci.lb,ci.ub=OR_df$ci.ub)
)

tabletext <- rbind(tablehead,tablenum)

means <- c(NA,NA,OR_df$OR)
lowers <- c(NA,NA,OR_df$ci.lb)
uppers <- c(NA,NA,OR_df$ci.ub)

# make plot/table
my_ticks <- c(2/3,1,3/2)
attr(my_ticks,"labels") <- c("2/3","1","3/2")
pdf(file="./results/ORs_interaction.pdf",width=7,height=3)
forestplot(tabletext,mean=means,lower=lowers,upper=uppers,
           #align=c("l",rep("r",ncol(tabletext)-1)),
           align=rep("c",3),
           zero=1,
           is.summary=c(TRUE,TRUE,rep(FALSE,nrow(tabletext)-2)),
           col=fpColors(box=c(cols2[5])),
           xlab="      Favors Men    Favors Women",
           graphwidth=unit(100,units="points"),
           lineheight=unit(22,units="points"),
           colgap=unit(6,"mm"),
           line.margin=0.2,
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex=0.9),
                            xlab  = gpar(cex = 1.2)),
           xlog=TRUE,xticks=my_ticks,xticks.digits=5,
           grid=T,
           boxsize=0.3,
           fn.ci_norm = fpDrawDiamondCI,
           new_page=F
)
dev.off()

cat("\n\n------------ Case control analysis by journal ----------------\n\n")

journals <- icc_df[!duplicated(icc_df$pub_sourceid),c("pub_sourceid","citescore","include.journal")]
outputs_df <- data.frame(journal=journals$pub_sourceid,
                         n_cases=numeric(length=nrow(journals))+NA,
                         effect=numeric(length=nrow(journals))+NA,
                         sd=numeric(length=nrow(journals))+NA,
                         pval=numeric(length=nrow(journals))+NA,
                         n_cases_adj=numeric(length=nrow(journals))+NA,
                         effect_adj=numeric(length=nrow(journals))+NA,
                         sd_adj=numeric(length=nrow(journals))+NA,
                         pval_adj=numeric(length=nrow(journals))+NA,
                         effect_int=numeric(length=nrow(journals))+NA,
                         sd_int=numeric(length=nrow(journals))+NA,
                         pval_int=numeric(length=nrow(journals))+NA,
                         effect_int2=numeric(length=nrow(journals))+NA,
                         sd_int2=numeric(length=nrow(journals))+NA,
                         pval_int2=numeric(length=nrow(journals))+NA,
                         citescore=journals$citescore,
                         included=journals$include.journal)

for(i in 1:nrow(journals)){
  # Unadjusted model
  mod <- tryCatch(clogit(case ~ Gender + strata(pub_id), 
                         icc_df, 
                         subset= pub_sourceid == outputs_df$journal[i]),
                  error=function(err) NA) # If the model won't run, return NA
  if(!is.na(mod)){
    outputs_df$effect[i] <- mod$coefficients[1]
    outputs_df$sd[i] <- sqrt(mod$var[1])
    outputs_df$pval[i] <- summary(mod)$coefficients[1,5]
    outputs_df$n_cases[i] <- summary(mod)$nevent
  } else {
    outputs_df$n_cases[i] <- sum(icc_df$pub_sourceid == outputs_df$journal[i])
  }
  # Adjusted model
  mod_adj <- tryCatch(clogit(case ~ Gender + 
                               ns(years_in_scopus_ptile,knots=knots) + 
                               ns(h_index_ptile,knots=knots) + 
                               ns(n_pubs_ptile,knots=knots) + 
                               strata(pub_id), 
                             icc_df, 
                             subset= pub_sourceid == outputs_df$journal[i]),
                      error=function(err) NA) # If the model won't run, return NA
  if(!is.na(mod_adj)){
    outputs_df$effect_adj[i] <- mod_adj$coefficients[1]
    outputs_df$sd_adj[i] <- sqrt(mod_adj$var[1])
    outputs_df$pval_adj[i] <- summary(mod_adj)$coefficients[1,5]
    outputs_df$n_cases_adj[i] <- summary(mod_adj)$nevent
  }
  # Adjusted model with interaction
  mod_int <- tryCatch(clogit(case ~ gender + gender_years_in_scopus_ptile + 
                               ns(years_in_scopus_ptile,knots=knots) +
                               ns(h_index_ptile,knots=knots) + 
                               ns(n_pubs_ptile,knots=knots) + 
                               strata(pub_id), 
                             icc_df, 
                             subset= pub_sourceid == outputs_df$journal[i]),
                      error=function(err) NA) # If the model won't run, return NA
  if(!is.na(mod_int)){
    outputs_df$effect_int[i] <- mod_int$coefficients[1]
    outputs_df$sd_int[i] <- sqrt(mod_int$var[1,1])
    outputs_df$pval_int[i] <- summary(mod_int)$coefficients[1,5]
    outputs_df$effect_int2[i] <- mod_int$coefficients[2]
    outputs_df$sd_int2[i] <- sqrt(mod_int$var[2,2])
    outputs_df$pval_int2[i] <- summary(mod_int)$coefficients[2,5]
  }
  if((i %% 100)==0) print(i)
}

cat(sum(!outputs_df$included),"journals did not meet criteria for obtaining individual estimate.\n")
cat("The criteria are:\n (1) must have at least two matched sets with not all the same gender\n")
cat("(2) among these matched sets, must have no zero cells when cross-tabulate gender against case status\n")

cat("---Summary of abs(log OR) for journals that did not meet criteria:---\n")
summary(abs(outputs_df$effect[!outputs_df$included]))
cat("\n\n---Summary of sd for journals that did not meet criteria:---\n")
summary(abs(outputs_df$sd[!outputs_df$included]))

cat("\n\n---Summary of abs(log OR) for journals that did meet criteria:---\n")
summary(abs(outputs_df$effect[outputs_df$included]))
cat("\n\n---Summary of sd for journals that did meet criteria:---\n")
summary(abs(outputs_df$sd[outputs_df$included]))

# Data frame with only valid ORs
outputs_select <- outputs_df[outputs_df$included,]

# Get ORs and CIs
outputs_select$node_size <- 1/outputs_select$sd
outputs_select$OR <- exp(outputs_select$effect)
outputs_select$ci_lower <- exp(outputs_select$effect-1.96*outputs_select$sd)
outputs_select$ci_upper <- exp(outputs_select$effect+1.96*outputs_select$sd)
outputs_select$node_size_adj <- 1/outputs_select$sd_adj
outputs_select$OR_adj <- exp(outputs_select$effect_adj)
outputs_select$ci_lower_adj <- exp(outputs_select$effect_adj-1.96*outputs_select$sd_adj)
outputs_select$ci_upper_adj <- exp(outputs_select$effect_adj+1.96*outputs_select$sd_adj)

# set to NA invalid estimates for adjusted OR
outputs_select$effect_adj[outputs_select$sd_adj > 2] <- NA

# Merge in journal topics
journal_topics <- readRDS(file = "./data/journal_topics.rds")
journal_topics <- journal_topics[journal_topics$pub_sourceid %in% unique(outputs_select$journal),]
outputs_select <- merge(outputs_select,journal_topics,
                        by.x="journal",by.y="pub_sourceid",all.x=T,all.y=F)

# remove empty topics
topic_counts <- colSums(outputs_select[,names(outputs_select) %in% names(journal_topics)[3:ncol(journal_topics)]],na.rm = T)
outputs_select[,names(outputs_select) %in% names(topic_counts)[topic_counts==0]] <- NULL
outputs_select[is.na(outputs_select) & col(outputs_select) > 20] <- 0
topics_list <- intersect(names(journal_topics)[-c(1,2)],names(outputs_select))

# save results for use in plotting and random effects meta-analysis
saveRDS(outputs_select, "./shiny_app/journal_ORs.rds")
saveRDS(topics_list,"./shiny_app/topics_list.rds")

#######################################################################
sink()
#######################################################################



