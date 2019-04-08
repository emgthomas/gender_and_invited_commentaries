# ------------------------------------------------ #
# --------------- Main analyses ------------------ #
# ------------------------------------------------ #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

# packages
require(survival)
require(metafor)
require(splines)

# record output
sink(file="./results/main_analyses2.txt")

## load data
icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

cat("\n\n------------ One-stage meta-analysis, all journals ----------------\n\n")
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

## Progressively excluding more controls based on match quality quantiles
n_quantiles <- max(icc_df$match_quantile)
resControls <- as.data.frame(matrix(nrow=n_quantiles-1,ncol=9))
names(resControls) <- c("lowest_decile","n_case","OR","OR_ci.lb","OR_ci.ub",
                                 "n_case.2","OR.2","OR_ci.lb.2","OR_ci.ub.2")
resControls$lowest_decile <- 1:(n_quantiles-1)
for(i in 1:(n_quantiles-1)){
  # Prepare dataset
  icc_df_i <- subset(icc_df,match_quantile>=i)
  which_pubs <- tapply(icc_df_i$case,icc_df_i$pub_id,length)
  which_pubs <- names(which_pubs[which_pubs>1])
  icc_df_i <- icc_df_i[icc_df_i$pub_id %in% which_pubs,]
  # Run models
  mod <- summary(clogit(case ~ Gender + strata(pub_id), data = icc_df_i))
  mod_adj <- summary(clogit(case ~ Gender + ns(years_in_scopus_ptile,knots=knots) + 
                              ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                              strata(pub_id), data = icc_df_i))
  # Store in data frame
  resControls[i,2:9] <- c(mod$nevent,mod$coefficients[2],mod$conf.int[3:4],
                          mod_adj$nevent,mod_adj$coefficients[1,2],mod_adj$conf.int[1,3:4])
}

require(ggplot2)
require(reshape2)
resControls_df <- reshape(resControls,direction="long",varying=list(c("n_case","n_case.2"),
                                                                    c("OR","OR.2"),
                                                                    c("OR_ci.lb","OR_ci.lb.2"),
                                                                    c("OR_ci.ub","OR_ci.ub.2")),
                          times=c("unadjusted","adjusted"))
names(resControls_df)[2] <- "model"
resControls_df$lowest_decile <- as.factor(resControls_df$lowest_decile)
plot1 <- ggplot(resControls_df,aes(x=lowest_decile,y=OR,color=model)) + 
  geom_pointrange(aes(ymin=OR_ci.lb,ymax=OR_ci.ub)) +
  ylab("OR (95%CI)") + xlab("Lowest decile of included controls based on match score") +
  theme(panel.grid.major.x = element_blank())

plot2 <- ggplot(subset(icc_df,case==0),aes(y=Match_Score,x=Gender)) + 
  geom_boxplot(outlier.shape=NA) + scale_y_continuous(limits=c(0,110))

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

# save results for use in plotting
saveRDS(outputs_select, "./shiny_app/journal_ORs.rds")
saveRDS(topics_list,"./shiny_app/topics_list.rds")

cat("\n\n------------ Meta-analysis, all journals ----------------\n\n")

### Using the rma function should give the same result as rma.mv below.
# all_2stage <- rma(effect,sd^2,data=outputs_select)
# summary(all_2stage)
# predict(all_2stage, transf=exp, digits=2)

cat("********** Undjusted model ************\n\n")

# Meta-analysis for unadjusted model
all_2stage <- rma.mv(effect,sd^2,random=~1|journal,data=outputs_select)
cat("---Summary of random effects meta-analysis for unadjusted model---\n\n")
summary(all_2stage)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage, transf=exp, digits=2)

cat("\n\n********** Adjusted model ************\n\n")

# Meta-analysis for adjusted model
all_2stage_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,data=outputs_select)
cat("\n\n---Summary of random effects meta-analysis for adjusted model---\n\n")
summary(all_2stage_adj)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage_adj, transf=exp, digits=2)

sink()

## Save key results for use in plotting ##
save(all_1stage,all_1stage_adj,all_2stage,all_2stage_adj,file="./results/main_analyses.Rdata")


