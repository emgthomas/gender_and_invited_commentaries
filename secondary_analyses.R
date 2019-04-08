# ---------------------------------------------------------- #
# ------------------ Secondary analyses -------------------- #
# ---------------------------------------------------------- #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

# packages
require(survival)
require(metafor)
require(splines)
require(forestplot)
require(plotly)
require(RColorBrewer)

# functions
formatting_fun <- function(idx,or,ci.lb,ci.ub){
  if(!is.na(or[idx])){
    return(paste0(sprintf(or[idx], fmt="%.2f")," (",sprintf(ci.lb[idx], fmt="%.2f"),",",sprintf(ci.ub[idx], fmt="%.2f"),")"))
  } else {
    return("N/A")
  }
}

## load data
outputs_select <- readRDS(file="./shiny_app/journal_ORs.rds")
icc_df <- readRDS(file="./data/processed_data_no_missing.rds")
topics_list <- readRDS(file="./shiny_app/topics_list.rds")
journal_topics <- readRDS(file = "./data/journal_topics.rds")
journal_topics <- journal_topics[journal_topics$pub_sourceid %in% unique(icc_df$pub_sourceid),]

# Save output
sink(file="./results/secondary_analyses2.txt")

cat("\n\n------------ Effect modification by author seniority ----------------\n\n")

cat("\n\n***Effect modification by years active***\n\n")
knots <- c(2.5,5,7.5)
all_1stage_EM_YiS <- clogit(case ~ Gender + years_in_scopus_ptile:factor(Gender,levels=c("female","male")) + 
                              ns(years_in_scopus_ptile,knots=knots) + 
                              ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                           strata(pub_id), data = icc_df)
summary(all_1stage_EM_YiS)

all_1stage_EM_HI <- clogit(case ~ Gender + h_index_ptile:factor(Gender,levels=c("female","male")) + 
                             ns(years_in_scopus_ptile,knots=knots) + 
                             ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) +
                              strata(pub_id), data = icc_df)

cat("\n\n***Effect modification by H Index***\n\n")
summary(all_1stage_EM_HI)

all_1stage_EM_npubs <- clogit(case ~ Gender + n_pubs_ptile:factor(Gender,levels=c("female","male")) + 
                                ns(years_in_scopus_ptile,knots=knots) + 
                                ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) +
                             strata(pub_id), data = icc_df)

cat("\n\n***Effect modification by number of publications***\n\n")
summary(all_1stage_EM_npubs)

# Plot main results
load("./results/main_analyses.Rdata")
OR_df <- data.frame(mod=c("2-stage unadjusted","2-stage adjusted",
                          "1-stage unadjusted","1-stage adjusted",
                          "1-stage adjusted","with effect modfication","by years active","",""),
                    ptile_YiS=c(rep("",4),c("10th percentile","30th percentile","50th percentile","70th percentile","90th percentile")),
                    OR=numeric(length=9),ci.lb=numeric(length=9),ci.ub=numeric(length=9))

OR_df$OR <- exp(c(all_2stage$beta[1,1],all_2stage_adj$beta[1,1],
              all_1stage$coefficients[1],all_1stage_adj$coefficients[1],
              all_1stage_EM_YiS$coefficients[1] + all_1stage_EM_YiS$coefficients[14],
              all_1stage_EM_YiS$coefficients[1] + 3*all_1stage_EM_YiS$coefficients[14],
              all_1stage_EM_YiS$coefficients[1] + 5*all_1stage_EM_YiS$coefficients[14],
              all_1stage_EM_YiS$coefficients[1] + 7*all_1stage_EM_YiS$coefficients[14],
              all_1stage_EM_YiS$coefficients[1] + 9*all_1stage_EM_YiS$coefficients[14]))

logOR_se <- c(all_2stage$se[1],all_2stage_adj$se[1],
              sqrt(all_1stage$var[1]),sqrt(all_1stage_adj$var[1]),
              sqrt(all_1stage_EM_YiS$var[1,1] + all_1stage_EM_YiS$var[14,14] + 2*all_1stage_EM_YiS$var[1,14]),
              sqrt(all_1stage_EM_YiS$var[1,1] + 9*all_1stage_EM_YiS$var[14,14] + 2*3*all_1stage_EM_YiS$var[1,14]),
              sqrt(all_1stage_EM_YiS$var[1,1] + 25*all_1stage_EM_YiS$var[14,14] + 2*5*all_1stage_EM_YiS$var[1,14]),
              sqrt(all_1stage_EM_YiS$var[1,1] + 49*all_1stage_EM_YiS$var[14,14] + 2*7*all_1stage_EM_YiS$var[1,14]),
              sqrt(all_1stage_EM_YiS$var[1,1] + 81*all_1stage_EM_YiS$var[14,14] + 2*9*all_1stage_EM_YiS$var[1,14]))

OR_df$ci.lb <- OR_df$OR*exp(-1.96*logOR_se)
OR_df$ci.ub <- OR_df$OR*exp(1.96*logOR_se)

tablehead <- rbind(c("Model","Years active","OR (95%CI)"),
                   rep(NA,3))
tablenum <- cbind(as.character(OR_df$mod),as.character(OR_df$ptile_YiS),
                  sapply(1:nrow(OR_df),formatting_fun,
                         or=OR_df$OR,ci.lb=OR_df$ci.lb,ci.ub=OR_df$ci.ub)
)

tabletext <- rbind(tablehead,tablenum)
tabletext <- rbind(tabletext[1:4,],rep("",ncol(tabletext)),
                   tabletext[5:6,],rep("",ncol(tabletext)),
                   tabletext[7:nrow(tabletext),])

means <- c(NA,NA,OR_df$OR[1:2],NA,OR_df$OR[3:4],NA,OR_df$OR[5:9])
lowers <- c(NA,NA,OR_df$ci.lb[1:2],NA,OR_df$ci.lb[3:4],NA,OR_df$ci.lb[5:9])
uppers <- c(NA,NA,OR_df$ci.ub[1:2],NA,OR_df$ci.ub[3:4],NA,OR_df$ci.ub[5:9])

# make plot/table
my_ticks <- c(2/3,1,3/2)
attr(my_ticks,"labels") <- c("2/3","1","3/2")
pdf(file="./results/ORs.pdf",width=8,height=4.5)
forestplot(tabletext,mean=means,lower=lowers,upper=uppers,
           align=c("l",rep("r",ncol(tabletext)-1)),
           zero=1,
           is.summary=c(TRUE,TRUE,rep(FALSE,nrow(tabletext)-2)),
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


cat("\n\n------------ Effect modification by journal citescore ----------------\n\n")

cat("********** Undjusted model ************\n\n")

# exclude one journal with high outlier citescore
outputs_select2 <- subset(outputs_select,citescore<20)
cat("Journal with outlier citescore of",max(outputs_select$citescore),":",
    outputs_select$sourcetitle[outputs_select$citescore>20])

# define knots
n_knots <- 3 # number of *internal* knots
citescore_by_pub <- unlist(sapply(1:nrow(outputs_select2),function(idx,n_cases,citescore) rep(citescore[idx],n_cases[idx]),
                                  n_cases=outputs_select2$n_cases,citescore=outputs_select2$citescore))
knot_placement <- quantile(citescore_by_pub,probs = seq(0,1,length.out=n_knots+2)[2:(n_knots+1)])
# run meta-regression
meta_analysis_cs <- rma.mv(effect,sd^2,random=~1|journal,
                           mods= ~ ns(citescore,knots=knot_placement),
                           data=outputs_select2)
cat("\n\n---Meta-regression on citescore, excluding one journal as outlier---\n\n")
summary(meta_analysis_cs)

###### Plot of journal-specific ORs with meta-regression on citescore ######
citescore_newmods <- seq(min(outputs_select2$citescore),max(outputs_select2$citescore),0.1)
citescore_bs <- ns(citescore_newmods,knots=knot_placement)
plot_citescore <- predict(meta_analysis_cs,newmods=citescore_bs[1:nrow(citescore_bs),1:ncol(citescore_bs)],
                          transf=exp)
plot_citescore <- plot_citescore[,c("pred","ci.lb","ci.ub","cr.lb","cr.ub")]
plot_citescore$citescore <- citescore_newmods

# Plot
cols1 <- brewer.pal(9,name="BuGn")
cols2 <- brewer.pal(9,name="Oranges")
OR_plot <- plot_ly(subset(outputs_select2,n_cases>50), x = ~citescore, y= ~OR,height=600,width=1000) %>%
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
              ymin=plot_citescore$cr.lb, ymax=plot_citescore$cr.ub, 
              hoverinfo="none",
              color=I(cols2[3]),
              name="95% prediction interval")%>%
  add_ribbons(x=plot_citescore$citescore,y=plot_citescore$pred,
              ymin=plot_citescore$ci.lb, ymax=plot_citescore$ci.ub, 
              hoverinfo="none",
              color=I(cols2[4]),
              name="95% confidence interval") %>%
  add_lines(y=plot_citescore$pred, x=plot_citescore$citescore,
            hoverinfo="none",
            color=I(cols2[5]),
            # line=list(color=I(cols2[1])),
            name="Predicted odds ratio")%>%
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
export(OR_plot, "./results/OR_metareg.pdf")

cat("\n\n********** Adjusted model--- ************\n\n")

# Adjusted model
meta_analysis_cs_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,
                           mods= ~ ns(citescore,knots=knot_placement),
                           data=subset(outputs_select2,!is.na(outputs_select2$effect_adj)))
cat("\n\nMeta-regression on citescore\n\n")
summary(meta_analysis_cs_adj)

###### Plot of journal-specific ORs with meta-regression on citescore ######
plot_citescore_adj <- predict(meta_analysis_cs_adj,newmods=citescore_bs[1:nrow(citescore_bs),1:ncol(citescore_bs)],
                          transf=exp)
plot_citescore_adj <- plot_citescore_adj[,c("pred","ci.lb","ci.ub","cr.lb","cr.ub")]
plot_citescore_adj$citescore <- citescore_newmods

OR_adj_plot <- plot_ly(subset(outputs_select2,n_cases>50), x = ~citescore, y= ~OR_adj,height=600,width=1000) %>%
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
              ymin=plot_citescore_adj$cr.lb, ymax=plot_citescore_adj$cr.ub, 
              hoverinfo="none",
              color=I(cols2[3]),
              name="95% prediction interval")%>%
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
export(OR_adj_plot, "./results/OR_adj_metareg.pdf")

cat("\n\n------------ Sub-group analyses by journal topic ----------------\n\n")

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
  meta1_adj <- clogit(case ~ Gender +  ns(years_in_scopus_ptile,knots=knots) + 
                        ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + strata(pub_id), data = icc_topic_i)
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

# # Repeat analyses for remaining topics (i.e., group all journals not in top 20 topics)
# cat("\n\n******** Analysing topic all journals not in top",n_topics,"topics ********** \n\n")
# all_journals <- unique(icc_df$pub_source_title)
# remaining_journals <- all_journals[!(all_journals %in% journals_analysed)]
# 
# # Which journals?
# icc_topic_i <- icc_df[icc_df$pub_source_title %in% remaining_journals,]
# outputs_topic_i <- outputs_select[outputs_select$journal %in% remaining_journals,]
# 
# # One-stage meta-analysis
# meta1 <- clogit(case ~ Gender + strata(pub_id), data = icc_topic_i)
# meta1_summ <- summary(meta1)
# cat("\n\n--- Unadjusted meta analysis, one-stage:\n")
# print(meta1_summ)
# outputs_topics[n_topics+1,c("OR_1stage","ci.lb_1stage","ci.ub_1stage")] <- meta1_summ$conf.int[1,c(1,3,4)]
# 
# cat("\n\n--- Adjusted for author seniority, one-stage:\n")
# meta1_adj <- clogit(case ~ Gender + years_in_scopus_quintile + 
#                       h_index_quintile + n_pubs_quintile + strata(pub_id), data = icc_topic_i)
# meta1_adj_summ <- summary(meta1_adj)
# print(meta1_adj_summ)
# outputs_topics[n_topics+1,c("OR_adj_1stage","ci.lb_adj_1stage","ci.ub_adj_1stage")] <- meta1_adj_summ$conf.int[1,c(1,3,4)]
# 
# outputs_topics$n_cases_1stage[n_topics+1] <- sum(icc_topic_i$case)
# outputs_topics$n_journals_1stage[n_topics+1] <- length(unique(icc_topic_i$pub_source_title))
# 
# # Two-stage meta-analysis
# meta2 <- rma.mv(effect,sd^2,random=~1|journal,data=outputs_topic_i)
# cat("\n\n--- Unadjusted meta analysis, two-stage:\n\n")
# print(summary(meta2))
# cat("\n\nMean and confidence (ci)/prediction (cr) intervals on odds ratio scale:\n")
# meta2_pred <- predict(meta2, transf=exp, digits=2)
# print(meta2_pred)
# outputs_topics[n_topics+1,c("OR_2stage","ci.lb_2stage","ci.ub_2stage",
#                    "pi.lb_2stage","pi.ub_2stage")] <- c(meta2_pred$pred,meta2_pred$ci.lb,meta2_pred$ci.ub,
#                                                         meta2_pred$cr.lb,meta2_pred$cr.ub)
# 
# meta2_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,data=outputs_topic_i)
# cat("\n\n--- Adjusted meta analysis, two-stage:\n\n")
# print(summary(meta2_adj))
# cat("\n\nMean and confidence (ci)/prediction (cr) intervals on odds ratio scale:\n")
# meta2_adj_pred <- predict(meta2_adj, transf=exp, digits=2)
# print(meta2_adj_pred)
# outputs_topics[n_topics+1,c("OR_adj_2stage","ci.lb_adj_2stage","ci.ub_adj_2stage",
#                    "pi.lb_adj_2stage","pi.ub_adj_2stage")] <- c(meta2_adj_pred$pred,meta2_adj_pred$ci.lb,meta2_adj_pred$ci.ub,
#                                                                 meta2_adj_pred$cr.lb,meta2_adj_pred$cr.ub)
# 
# outputs_topics$n_cases_2stage[n_topics+1] <- sum(outputs_topic_i$n_cases)
# outputs_topics$n_journals_2stage[n_topics+1] <- nrow(outputs_topic_i)

# Sort by ASJC code
outputs_topics <- outputs_topics[c(order(outputs_topics$ASJC[1:(nrow(outputs_topics)-1)]),nrow(outputs_topics)),]

# Get rid of some bad results
outputs_topics[outputs_topics$ci.lb_adj_2stage < 0.1,c("OR_adj_2stage","ci.lb_adj_2stage","ci.ub_adj_2stage")] <- NA

#### Produce forest plot3 #####

load("./results/main_analyses.Rdata")

## one-stage meta-analysis
tablehead <- rbind(c("Topic","ASJC","Journals","Cases","OR (95%CI)","AOR (95%CI)"),
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
pdf(file="./results/forestplot_1stage.pdf",width=13,height=17)
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
           #fn.ci_sum=c(fpDrawDiamondCI, fpDrawCircleCI),
           legend = c("Unadjusted OR", "Adjusted OR"),
           legend_args = fpLegend(pos = list("topright")),
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
pdf(file="./results/forestplot_2stage.pdf",width=13,height=17)
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
           legend = c("Unadjusted OR", "Adjusted OR"),
           legend_args = fpLegend(pos = list("topright")),
           new_page=F
)
dev.off()


sink()

