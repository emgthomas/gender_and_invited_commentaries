# ---------------------------------------------------------- #
# ------------------ Secondary analyses -------------------- #
# ---------------------------------------------------------- #

#######################################################################
sink(file="./results/secondary_analyses.txt")
#######################################################################

## packages
require(survival)
require(metafor)
require(splines)
require(forestplot)
require(plotly)
require(RColorBrewer)
require(lmtest)

## global functions
source("./code/functions.R")

## load data
outputs_select <- readRDS(file="./shiny_app/journal_ORs.rds")
icc_df <- readRDS(file="./data/processed_data_no_missing.rds")
topics_list <- readRDS(file="./shiny_app/topics_list.rds")
journal_topics <- readRDS(file = "./data/journal_topics.rds")
journal_topics <- journal_topics[journal_topics$pub_sourceid %in% unique(icc_df$pub_sourceid),]

cat("\n\n------------ Effect modification by journal citescore ----------------\n\n")

cat("********** Undjusted model ************\n\n")

# exclude one journal with high outlier citescore
icc_df_cs <- subset(icc_df,citescore<20)
outputs_select2 <- subset(outputs_select,citescore<20)
cat("Journal with outlier citescore of",max(icc_df$citescore),":",
    as.character(unique(icc_df$pub_source_title[icc_df$citescore>20])))

# define knots
n_knots <- 3 # number of *internal* knots
citescore_by_pub <- icc_df_cs$citescore[icc_df_cs$case==1]
knot_placement <- quantile(citescore_by_pub,probs = seq(0,1,length.out=n_knots+2)[2:(n_knots+1)])
icc_df_cs$gender <- as.numeric(icc_df_cs$Gender == "female")
citescore_ns <- ns(icc_df_cs$citescore,knots=knot_placement)
icc_df_cs$citescore_gender_1 <- citescore_ns[,1]*icc_df_cs$gender
icc_df_cs$citescore_gender_2 <- citescore_ns[,2]*icc_df_cs$gender
icc_df_cs$citescore_gender_3 <- citescore_ns[,3]*icc_df_cs$gender
icc_df_cs$citescore_gender_4 <- citescore_ns[,4]*icc_df_cs$gender
icc_df_cs$citescore_gender <- icc_df_cs$citescore*icc_df_cs$gender

# run regression with interaction by journal cite score
knots <- c(2.5,5,7.5)
all_1stage_cs <- clogit(case ~ gender + 
                          citescore_gender_1 + citescore_gender_2 + 
                          citescore_gender_3 + citescore_gender_4 +
                          strata(pub_id), data = icc_df_cs)
summary(all_1stage_cs)

# model without effect modifcation by cite score
all_1stage <- clogit(case ~ gender + strata(pub_id), data = icc_df_cs)

cat("\n\nTest for null hypothesis of no effect of Cite Score on odds ratio\n\n")
lrtest(all_1stage_cs,all_1stage)

###### Plot of journal-specific ORs with meta-regression on citescore ######
# set up data for prediction
n.points <- 1000
citescore_new <- seq(min(icc_df_cs$citescore),max(icc_df_cs$citescore),length.out = n.points)
citescore_ns <- ns(citescore_new,knots=knot_placement)
dat_pred <- cbind(gender=1,citescore_ns)
# Get predicted values
varnames <- c("gender","citescore_gender_1","citescore_gender_2",
              "citescore_gender_3","citescore_gender_4")
vars <- all_1stage_cs$assign
vars_num <- as.integer(vars[varnames])
beta <- all_1stage_cs$coefficients
vars_beta <- beta[vars_num]
lp <- dat_pred %*% vars_beta
# Get standard errors
v <- all_1stage_cs$var
v_vars <- v[vars_num,vars_num]
var_lp <- dat_pred %*% v_vars %*% t(dat_pred)
se_lp <- sqrt(diag(var_lp))
# Put into data frame for plotting
plot_citescore <- data.frame(pred=exp(lp),citescore=citescore_new)
plot_citescore$ci.lb <- as.numeric(exp(lp - 1.96*se_lp))
plot_citescore$ci.ub <- as.numeric(exp(lp + 1.96*se_lp))

# Plot
cols1 <- brewer.pal(9,name="BuGn")
cols2 <- brewer.pal(9,name="Oranges")
OR_plot <- plot_ly(subset(outputs_select2,n_cases>50), x = ~citescore, y= ~OR,height=600,width=1000, type="scatter", mode="markers") %>%
  # add horizontal line for null value
  add_trace(x = c(0,18), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  # add point for OR of each journal
  add_trace(y=~OR, type='scatter',mode='markers',
            size = ~node_size*2, 
            marker=list(sizeref=0.2),
            hoverinfo='none',
            alpha=0.8,
            color=I(cols1[4]),
            name="Journal-specific odds ratio"
  )  %>%
  # annotations
  add_trace(x = 16, y = 5, mode="text",text="Favors women",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 16, y = 1/5, mode="text",text="Favors men",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_ribbons(x=plot_citescore$citescore,y=plot_citescore$pred,
              ymin=plot_citescore$ci.lb, ymax=plot_citescore$ci.ub,
              hoverinfo="none",
              color=I(cols2[4]),
              name="95% confidence interval") %>%
  add_lines(y=plot_citescore$pred, x=plot_citescore$citescore,
            hoverinfo="none",
            color=I(cols2[5]),
            # line=list(color=I(cols2[1])),
            name="Predicted odds ratio") %>%
  # layout
  layout(yaxis = list(title="Odds Ratio (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4,8),
                      ticktext=c("1/8","1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Journal Cite Score",range=c(0,18),tickmode="array"),
         showlegend=T,
         font=list(size=16)
  )

# Save as pdf
export(OR_plot, "./results/OR_citescore.pdf")

cat("\n\n********** Adjusted model--- ************\n\n")

# run regression
knots <- c(2.5,5,7.5)
all_1stage_cs_adj <- clogit(case ~ gender +
                          citescore_gender_1 + citescore_gender_2 + citescore_gender_3 + citescore_gender_4 +
                          ns(years_in_scopus_ptile,knots=knots) +
                          ns(h_index_ptile,knots=knots) + 
                          ns(n_pubs_ptile,knots=knots) +
                          strata(pub_id), data = icc_df_cs)

summary(all_1stage_cs_adj)

# model without effect modifcation by cite score
all_1stage_adj <- clogit(case ~ gender + 
                           ns(years_in_scopus_ptile,knots=knots) +
                           ns(h_index_ptile,knots=knots) + 
                           ns(n_pubs_ptile,knots=knots) +
                           strata(pub_id), data = icc_df_cs)

cat("\n\nTest for null hypothesis of no effect of Cite Score on odds ratio\n\n")
lrtest(all_1stage_cs_adj,all_1stage_adj)

###### Plot of journal-specific ORs with meta-regression on citescore ######
# Get predicted values
vars <- all_1stage_cs_adj$assign
vars_num <- as.integer(vars[varnames])
beta <- all_1stage_cs_adj$coefficients
vars_beta <- beta[vars_num]
lp <- dat_pred %*% vars_beta
# Get standard errors
v <- all_1stage_cs_adj$var
v_vars <- v[vars_num,vars_num]
var_lp <- dat_pred %*% v_vars %*% t(dat_pred)
se_lp <- sqrt(diag(var_lp))
# Put into data frame for plotting
plot_citescore_adj <- data.frame(pred=exp(lp),citescore=citescore_new)
plot_citescore_adj$ci.lb <- as.numeric(exp(lp - 1.96*se_lp))
plot_citescore_adj$ci.ub <- as.numeric(exp(lp + 1.96*se_lp))

# Make plot
OR_adj_plot <- plot_ly(subset(outputs_select2,n_cases>50), x = ~citescore, y= ~OR_adj,height=600,width=1000, type="scatter", mode="markers") %>%
  # add horizontal line for null value
  add_trace(x = c(0,18), y= c(1,1), mode = "lines", color=I("black"),
            showlegend=F) %>%
  # add point for OR of each journal
  add_trace(y=~OR_adj, type='scatter',mode='markers',
            size = ~node_size*2, 
            marker=list(sizeref=0.2),
            hoverinfo='none',
            color=I(cols1[4]),
            alpha=0.8,
            name="Journal-specific odds ratio"
  )  %>%
  # annotations
  add_trace(x = 16, y = 5, mode="text",text="Favors women",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_trace(x = 16, y = 1/5, mode="text",text="Favors men",
            hoverinfo="none",textfont=list(size=20,color=1),
            showlegend=F) %>%
  add_ribbons(x=plot_citescore_adj$citescore,y=plot_citescore_adj$pred,
              ymin=plot_citescore_adj$ci.lb, ymax=plot_citescore_adj$ci.ub, 
              hoverinfo="none",
              color=I(cols2[4]),
              name="95% confidence interval") %>%
  add_lines(y=plot_citescore_adj$pred, x=plot_citescore_adj$citescore,
            hoverinfo="none",
            color=I(cols2[5]),
            name="Predicted odds ratio")%>%
  # layout
  layout(yaxis = list(title="Odds Ratio (log scale)",range=c(-log(2.2),log(2.2)),type="log",
                      tickvals=c(1/4,1/2,1,2,4),
                      ticktext=c("1/4","1/2",
                                 as.character(c(1,2,4)))),
         xaxis = list(title="Journal Cite Score",range=c(0,18),tickmode="array"),
         showlegend=T,
         font=list(size=16)
  )

# Save as pdf
export(OR_adj_plot, "./results/OR_adj_citescore.pdf")

cat("\n\n------------ Sub-group analyses by journal topic ----------------\n\n")
cat("Note: here we perform sub-group analyses by topic using both one- and two-stage
meta-analysis. The two-stage approach is treated as a sensitivity analysis, but
we include it here for ease of coding.\n")

# Total number of cases by topic (among journals included in 2-stage meta-analysis)
topic_case_counts <- sapply(topics_list,function(topic,df) sum(df$n_cases[df[,names(df) == topic]==1],na.rm=T),
                            df=outputs_select)
topic_journal_counts <- colSums(outputs_select[,names(outputs_select) %in% topics_list],na.rm=T)
topic_counts <- data.frame(topic=names(topic_case_counts), case_counts=topic_case_counts)
topic_counts <- merge(topic_counts,
                      data.frame(topic=names(topic_journal_counts), journal_counts=topic_journal_counts))
topic_names <- read.csv(file="./data/ASJC Codes with levels.csv",
                        sep=";") # cloned from github.com/plreyes/Scopus.git
topic_counts <- merge(topic_counts,topic_names[,c("Low","Code")],by.x="topic",by.y="Low")
topic_counts <- topic_counts[order(topic_counts$case_counts,decreasing = T),]

# Number of topics
n_topics <- nrow(topic_counts)

# Set up data frame for results
outputs_topics <- data.frame(ASJC=c(as.character(topic_counts$Code[1:n_topics])),
                             topic=c(as.character(topic_counts$topic[1:n_topics])),
                             n_cases_1stage=integer(n_topics),
                             n_journals_1stage=integer(n_topics),
                             OR_1stage=numeric(n_topics),ci.lb_1stage=numeric(n_topics),ci.ub_1stage=numeric(n_topics),
                             OR_adj_1stage=numeric(n_topics),ci.lb_adj_1stage=numeric(n_topics),ci.ub_adj_1stage=numeric(n_topics),
                             n_cases_2stage=integer(n_topics),
                             n_journals_2stage=integer(n_topics),
                             OR_2stage=numeric(n_topics),ci.lb_2stage=numeric(n_topics),ci.ub_2stage=numeric(n_topics),
                             pi.lb_2stage=numeric(n_topics),pi.ub_2stage=numeric(n_topics),
                             OR_adj_2stage=numeric(n_topics),ci.lb_adj_2stage=numeric(n_topics),ci.ub_adj_2stage=numeric(n_topics),
                             pi.lb_adj_2stage=numeric(n_topics),pi.ub_adj_2stage=numeric(n_topics)
  ,stringsAsFactors = F)
outputs_topics$n_cases_2stage[1:n_topics] <- topic_counts$case_counts[1:n_topics]
outputs_topics$n_journals_2stage[1:n_topics] <- topic_counts$journal_counts[1:n_topics]
journals_analysed <- c()

# Some journals are omitted from two-stage analyses due to large standard errors
outputs_select$OR_adj[outputs_select$sd_adj>2500] <- NA
outputs_select$sd_adj[outputs_select$sd_adj>2500] <- NA

# Run analyses
for(i in 1:n_topics){
  cat("\n\n******** Analysing topic:",as.character(topic_counts$topic[i]),"********** \n\n")
  cat("\n\n Topic number",i,"\n\n")
  
  cat("Total cases:",topic_counts$case_counts[i],"\n")
  cat("Number of journals:",topic_counts$journal_counts[i],"\n")
  
  # Which journals?
  journals <- journal_topics$pub_sourceid[journal_topics[,names(journal_topics) == topic_counts$topic[i]]==1]
  journals_analysed <- union(journals,journals_analysed)
  icc_topic_i <- icc_df[icc_df$pub_sourceid %in% journals,]
  which_journal <- outputs_select[,names(outputs_select)==topic_counts$topic[i]]==1
  outputs_topic_i <- outputs_select[which_journal,]
  
  # One-stage meta-analysis
  meta1 <- clogit(case ~ Gender + strata(pub_id), data = icc_topic_i)
  meta1_summ <- summary(meta1)
  cat("\n\n--- Unadjusted meta analysis, one-stage:\n")
  print(meta1_summ)
  outputs_topics[i,c("OR_1stage","ci.lb_1stage","ci.ub_1stage")] <- meta1_summ$conf.int[1,c(1,3,4)]
  
  cat("\n\n--- Adjusted for author seniority, one-stage:\n")
  meta1_adj <- clogit(case ~ Gender +  
                        ns(years_in_scopus_ptile,knots=knots) + 
                        ns(h_index_ptile,knots=knots) + 
                        ns(n_pubs_ptile,knots=knots) + 
                        strata(pub_id), 
                      data = icc_topic_i)
  meta1_adj_summ <- summary(meta1_adj)
  print(meta1_adj_summ)
  outputs_topics[i,c("OR_adj_1stage","ci.lb_adj_1stage","ci.ub_adj_1stage")] <- meta1_adj_summ$conf.int[1,c(1,3,4)]
  
  outputs_topics$n_cases_1stage[i] <- sum(icc_topic_i$case)
  outputs_topics$n_journals_1stage[i] <- length(unique(icc_topic_i$pub_source_title))
  
  # Two-stage meta-analysis
  meta2 <- rma.mv(effect,sd^2,random=~1|journal,data=outputs_topic_i)
  cat("\n\n--- Unadjusted meta analysis, two-stage:\n\n")
  print(summary(meta2))
  cat("\n\nMean and confidence (ci)/prediction (cr) intervals on odds ratio scale:\n")
  meta2_pred <- predict(meta2, transf=exp, digits=2)
  print(meta2_pred)
  outputs_topics[i,c("OR_2stage","ci.lb_2stage","ci.ub_2stage",
                    "pi.lb_2stage","pi.ub_2stage")] <- c(meta2_pred$pred,meta2_pred$ci.lb,meta2_pred$ci.ub,
                                                                 meta2_pred$cr.lb,meta2_pred$cr.ub)
  
  meta2_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,data=outputs_topic_i)
  cat("\n\n--- Adjusted meta analysis, two-stage:\n\n")
  print(summary(meta2_adj))
  cat("\n\nMean and confidence (ci)/prediction (cr) intervals on odds ratio scale:\n")
  meta2_adj_pred <- predict(meta2_adj, transf=exp, digits=2)
  print(meta2_adj_pred)
  outputs_topics[i,c("OR_adj_2stage","ci.lb_adj_2stage","ci.ub_adj_2stage",
                    "pi.lb_adj_2stage","pi.ub_adj_2stage")] <- c(meta2_adj_pred$pred,meta2_adj_pred$ci.lb,meta2_adj_pred$ci.ub,
                                                                 meta2_adj_pred$cr.lb,meta2_adj_pred$cr.ub)
}

# Sort by ASJC code
outputs_topics <- outputs_topics[order(outputs_topics$ASJC),]

# Get rid of some bad results
outputs_topics[outputs_topics$ci.lb_adj_2stage < 0.1,c("OR_adj_2stage","ci.lb_adj_2stage","ci.ub_adj_2stage")] <- NA

# How many significant effects?
cat(sum(outputs_topics$ci.ub_1stage < 1),"of",nrow(outputs_topics),"topic-specific ORs showed significant bias against women.")
cat(sum(outputs_topics$ci.ub_adj_1stage < 1),"of",nrow(outputs_topics),"topic-specific fully adjusted ORs showed significant bias against women.")

save(outputs_topics,file="./results/secondary_analyses.Rdata")

#### Produce forest plot3 #####

load("./results/main_analyses.Rdata")
load("./results/two_stage_analyses.Rdata")
load("./results/secondary_analyses.Rdata")

## one-stage meta-analysis
tablehead <- rbind(c("Topic","ASJC","Journals","Cases","Model 1 OR (95%CI)","Model 2 OR (95%CI)"),
                   rep(NA,6)
                   )
tablenum <- cbind(as.matrix(outputs_topics[,c("topic","ASJC","n_journals_1stage","n_cases_1stage")]),
                  sapply(1:nrow(outputs_topics),formatting_fun,
                         or=outputs_topics$OR_1stage,ci.lb=outputs_topics$ci.lb_1stage,
                         ci.ub=outputs_topics$ci.ub_1stage),
                  sapply(1:nrow(outputs_topics),formatting_fun,
                         or=outputs_topics$OR_adj_1stage,ci.lb=outputs_topics$ci.lb_adj_1stage,
                         ci.ub=outputs_topics$ci.ub_adj_1stage)
                   )
tablefoot <- rbind(rep(NA,6),
                   c("All","",as.character(length(unique(icc_df$pub_source_title))),as.character(sum(icc_df$case)),
                   formatting_fun(1,summary(all_1stage)$conf.int[1,1],
                                  summary(all_1stage)$conf.int[1,3],
                                  summary(all_1stage)$conf.int[1,4]),
                   formatting_fun(1,summary(all_1stage_adj)$conf.int[1,1],
                                  summary(all_1stage_adj)$conf.int[1,3],
                                  summary(all_1stage_adj)$conf.int[1,4]))
                   )
tabletext <- rbind(tablehead,tablenum,tablefoot)

means <- rbind(rep(NA,2),rep(NA,2),
               outputs_topics[,c("OR_1stage","OR_adj_1stage")],
               rep(NA,2),
          c(summary(all_1stage)$conf.int[1,1],summary(all_1stage_adj)$conf.int[1,1])
          )
lowers <- rbind(rep(NA,2),rep(NA,2),
                outputs_topics[,c("ci.lb_1stage","ci.lb_adj_1stage")],
               rep(NA,2),
               c(summary(all_1stage)$conf.int[1,3],summary(all_1stage_adj)$conf.int[1,3])
          )
uppers <- rbind(rep(NA,2),rep(NA,2),
                outputs_topics[,c("ci.ub_1stage","ci.ub_adj_1stage")],
                rep(NA,2),
                c(summary(all_1stage)$conf.int[1,4],summary(all_1stage_adj)$conf.int[1,4])
)

# make plot/table
my_ticks <- c(1/4,1/2,1,2,4)
attr(my_ticks,"labels") <- c("1/4","1/2","1","2","4")
pdf(file="./results/forestplot_1stage.pdf",width=14,height=15)
forestplot(tabletext,mean=means,lower=lowers,upper=uppers,
           is.summary=c(TRUE,TRUE,rep(FALSE,nrow(outputs_topics)),TRUE),
           align=c("l",rep("r",ncol(tabletext)-1)),
           col=fpColors(box=c(cols1[5], cols2[5]),summary=c(cols1[5], cols2[5])),
           zero=1,
           xlab="      Favors Men    Favors Women",
           graphwidth=unit(120,units="points"),
           lineheight=unit(22,units="points"),
           colgap=unit(6,"mm"),
           line.margin=0.2,
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex=0.9),
                            xlab  = gpar(cex = 1.2)),
           xlog=TRUE,xticks=my_ticks,xticks.digits=5,
           grid=T,
           boxsize=0.3,
           fn.ci_norm = c(fpDrawDiamondCI, fpDrawCircleCI),
           legend = c("Model 1 OR", "Model 2 OR"),
           legend_args = fpLegend(pos = list("top")),
           new_page=F
)
dev.off()

## 2-stage meta-analysis
tablenum <- cbind(as.matrix(outputs_topics[,c("topic","ASJC","n_journals_2stage","n_cases_2stage")]),
                  sapply(1:nrow(outputs_topics),formatting_fun,
                         or=outputs_topics$OR_2stage,ci.lb=outputs_topics$ci.lb_2stage,
                         ci.ub=outputs_topics$ci.ub_2stage),
                  sapply(1:nrow(outputs_topics),formatting_fun,
                         or=outputs_topics$OR_adj_2stage,ci.lb=outputs_topics$ci.lb_adj_2stage,
                         ci.ub=outputs_topics$ci.ub_adj_2stage)
)
tablefoot <- rbind(rep(NA,6),
                   c("All","",as.character(length(unique(icc_df$pub_source_title))),as.character(sum(icc_df$case)),
                     formatting_fun(1,predict(all_2stage, transf=exp)[[1]],
                                    predict(all_2stage, transf=exp)[[3]],
                                    predict(all_2stage, transf=exp)[[4]]),
                     formatting_fun(1,predict(all_2stage_adj, transf=exp)[[1]],
                                    predict(all_2stage_adj, transf=exp)[[3]],
                                    predict(all_2stage_adj, transf=exp)[[4]]))
)
tabletext <- rbind(tablehead,tablenum,tablefoot)

means <- rbind(rep(NA,2),rep(NA,2),
               outputs_topics[,c("OR_2stage","OR_adj_2stage")],
               rep(NA,2),
               c(predict(all_2stage, transf=exp)[[1]],predict(all_2stage_adj, transf=exp)[[1]])
)
lowers <- rbind(rep(NA,2),rep(NA,2),
                outputs_topics[,c("ci.lb_2stage","ci.lb_adj_2stage")],
                rep(NA,2),
                c(predict(all_2stage, transf=exp)[[3]],predict(all_2stage_adj, transf=exp)[[3]])
)
uppers <- rbind(rep(NA,2),rep(NA,2),
                outputs_topics[,c("ci.ub_2stage","ci.ub_adj_2stage")],
                rep(NA,2),
                c(predict(all_2stage, transf=exp)[[4]],predict(all_2stage_adj, transf=exp)[[4]])
)

# make plot/table
pdf(file="./results/forestplot_2stage.pdf",width=14,height=15)
forestplot(tabletext,mean=means,lower=lowers,upper=uppers,
           is.summary=c(TRUE,TRUE,rep(FALSE,nrow(outputs_topics)),TRUE),
           align=c("l",rep("r",ncol(tabletext)-1)),
           col=fpColors(box=c(cols1[5], cols2[5]),summary=c(cols1[5], cols2[5])),
           zero=1,
           xlab="      Favors Men    Favors Women",
           graphwidth=unit(120,units="points"),
           lineheight=unit(22,units="points"),
           colgap=unit(6,"mm"),
           line.margin=0.2,
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex=0.9),
                            xlab  = gpar(cex = 1.2)),
           xlog=TRUE,xticks=my_ticks,xticks.digits=5,
           grid=T,
           boxsize=0.3,
           fn.ci_norm = c(fpDrawDiamondCI, fpDrawCircleCI),
           legend = c("Model 1 OR", "Model 2 OR"),
           legend_args = fpLegend(pos = list("top")),
           new_page=F
)
dev.off()

#######################################################################
sink()
#######################################################################
