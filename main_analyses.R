# ------------------------------------------------ #
# --------------- Main analyses ------------------ #
# ------------------------------------------------ #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

# packages
require(survival)
require(metafor)
require(splines)

#######################################################################
sink(file="./results/main_analyses2.txt")
#######################################################################

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
deciles <- c(1,3,5,7,9)

OR_df <- data.frame(ptile_YiS=c("10th","30th","50th","70th","90th"),
                    OR=numeric(length=5),ci.lb=numeric(length=5),ci.ub=numeric(length=5))

OR_df$OR <- exp(all_1stage_EM_YiS$coefficients[1] + deciles*all_1stage_EM_YiS$coefficients[14])

logOR_se <- sqrt(all_1stage_EM_YiS$var[1,1] + deciles^2*all_1stage_EM_YiS$var[14,14] + 2*deciles*all_1stage_EM_YiS$var[1,14])

OR_df$ci.lb <- OR_df$OR*exp(-1.96*logOR_se)
OR_df$ci.ub <- OR_df$OR*exp(1.96*logOR_se)

tablehead <- rbind(c("Years\nActive","Percentile","Odds ratio"),
                   rep(NA,3))
# tablenum <- cbind(as.character(OR_df$ptile_YiS),
#                   sapply(1:nrow(OR_df),formatting_fun,
#                          or=OR_df$OR,ci.lb=OR_df$ci.lb,ci.ub=OR_df$ci.ub)
# )
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>0.9 & icc_df$years_in_scopus_ptile<1.1])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>2.8 & icc_df$years_in_scopus_ptile<3.2])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>3.8 & icc_df$years_in_scopus_ptile<4.2])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>6.9 & icc_df$years_in_scopus_ptile<7.1])
unique(icc_df$years_in_scopus[icc_df$years_in_scopus_ptile>8.94 & icc_df$years_in_scopus_ptile<9.06])
years <- c(8,14,16,27,38)
tablenum <- cbind(years,as.character(OR_df$ptile_YiS),
                  sprintf(OR_df$OR, fmt="%.2f")
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
  mod_adj <- tryCatch(clogit(case ~ Gender + ns(years_in_scopus_ptile,knots=knots) + 
                               ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + strata(pub_id), 
                             icc_df, 
                             subset= pub_sourceid == outputs_df$journal[i]),
                      error=function(err) NA) # If the model won't run, return NA
  if(!is.na(mod_adj)){
    outputs_df$effect_adj[i] <- mod_adj$coefficients[1]
    outputs_df$sd_adj[i] <- sqrt(mod_adj$var[1])
    outputs_df$pval_adj[i] <- summary(mod_adj)$coefficients[1,5]
    outputs_df$n_cases_adj[i] <- summary(mod_adj)$nevent
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



