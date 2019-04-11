# ---------------------------------------------------------- #
# ----------------- Sensitivity analyses ------------------- #
# ---------------------------------------------------------- #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

#######################################################################
sink(file="./results/sensitivity_analyses.txt")
#######################################################################

cat("\n\n--------------------------------------------------------------------\n\n")
cat("\n\n--------------------- Two-stage meta-analysis ----------------------\n\n")
cat("\n\n--------------------------------------------------------------------\n\n")

outputs_select <- readRDS("./shiny_app/journal_ORs.rds")

cat("********** Undjusted model ************\n\n")

# Meta-analysis for unadjusted model
all_2stage <- rma.mv(effect,sd^2,random=~1|journal,data=outputs_select)
cat("---Summary of random effects meta-analysis for unadjusted model---\n\n")
summary(all_2stage)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage, transf=exp, digits=2)

### Check: using the rma function should give the same result as rma.mv above.
# all_2stage <- rma(effect,sd^2,data=outputs_select)
# summary(all_2stage)
# predict(all_2stage, transf=exp, digits=2)

cat("\n\n********** Adjusted model ************\n\n")

# Meta-analysis for adjusted model
all_2stage_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,data=outputs_select)
cat("\n\n---Summary of random effects meta-analysis for adjusted model---\n\n")
summary(all_2stage_adj)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage_adj, transf=exp, digits=2)

## Save key results for use in plotting ##
save(all_2stage,all_2stage_adj,file="./results/two_stage_analyses.Rdata")

cat("\n\n---------------------------------------------------------\n\n")
cat("\n\n--------------------- Missing data ----------------------\n\n")
cat("\n\n---------------------------------------------------------\n\n")

icc_df_all <- readRDS(file="./data/processed_data_all.rds")
require(dplyr)

# Assume all missing cases are female, all missing control are men
# What impact will this have on results?s
icc_df_all$Gender[icc_df_all$case==1 & icc_df_all$Gender == "unknown"] <- "female"
icc_df_all$Gender[icc_df_all$case==0 & icc_df_all$Gender == "unknown"] <- "male"
icc_df_all$Gender <- factor(icc_df_all$Gender,levels=c("male","female"))

# Run analyses
cat("----------------- Results if we assume all missing cases are female, all missing controls are male ---------------\n\n")
all_1stage_miss <- clogit(case ~ Gender + strata(pub_id), data = icc_df_all)
cat("---Unadjusted analysis---\n")
summary(all_1stage_miss)

cat("\n\n---Adjusted analysis---\n")
knots <- c(2.5,5,7.5)
all_1stage_miss_adj <- clogit(case ~ Gender + ns(years_in_scopus_ptile,knots=knots) + 
                           ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                           strata(pub_id), data = icc_df_all)
summary(all_1stage_miss_adj)

cat("\n\n---------------------------------------------------------------------------------\n\n")
cat("\n\n--------------------- Stricter match criteria for controls ----------------------\n\n")
cat("\n\n---------------------------------------------------------------------------------\n\n")

icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

n_quantiles <- max(icc_df$match_quantile)
resControls <- as.data.frame(matrix(nrow=n_quantiles-1,ncol=16))
names(resControls) <- c("lowest_decile","gender_cases","h_index_female_cases","h_index_male_cases",
                        "n_pubs_female_cases","n_pubs_male_cases",
                        "years_active_female_cases","years_active_male_cases",
                        "n_case","OR","OR_ci.lb","OR_ci.ub",
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
                              #ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                              strata(pub_id), data = icc_df_i))
  # Store in data frame
  resControls$gender_cases[i] <- mean(icc_df_i$Gender[icc_df_i$case==1] == "female")
  resControls$h_index_female_cases[i] <- median(icc_df_i$H_Index[icc_df_i$case==1 & icc_df_i$Gender=="female"])
  resControls$h_index_male_cases[i] <- median(icc_df_i$H_Index[icc_df_i$case==1 & icc_df_i$Gender=="male"])
  resControls$n_pubs_female_cases[i] <- median(icc_df_i$Total_Publications_In_Scopus[icc_df_i$case==1 & icc_df_i$Gender=="female"])
  resControls$n_pubs_male_cases[i] <- median(icc_df_i$Total_Publications_In_Scopus[icc_df_i$case==1 & icc_df_i$Gender=="male"])
  resControls$years_active_female_cases[i] <- median(icc_df_i$years_in_scopus[icc_df_i$case==1 & icc_df_i$Gender=="female"])
  resControls$years_active_male_cases[i] <- median(icc_df_i$years_in_scopus[icc_df_i$case==1 & icc_df_i$Gender=="male"])
  resControls[i,9:16] <- c(mod$nevent,mod$coefficients[2],mod$conf.int[3:4],
                           mod_adj$nevent,mod_adj$coefficients[1,2],mod_adj$conf.int[1,3:4])
}

require(ggplot2)
require(reshape2)
resControls_df <- reshape(resControls,direction="long",varying=list(c("n_case","n_case.2"),
                                                                    c("OR","OR.2"),
                                                                    c("OR_ci.lb","OR_ci.lb.2"),
                                                                    c("OR_ci.ub","OR_ci.ub.2")),
                          times=c("unadjusted","adjusted"))
names(resControls_df)[names(resControls_df)=="time"] <- "model"
resControls_df$lowest_decile <- as.factor(resControls_df$lowest_decile)
plot1 <- ggplot(resControls_df,aes(x=lowest_decile,y=OR,color=model)) + 
  geom_pointrange(aes(ymin=OR_ci.lb,ymax=OR_ci.ub)) +
  ylab("OR (95%CI)") + xlab("Lowest decile of included controls based on match score") +
  theme(panel.grid.major.x = element_blank())

plot2 <- ggplot(subset(icc_df,case==0),aes(y=Match_Score,x=Gender)) + 
  geom_boxplot(outlier.shape=NA) + scale_y_continuous(limits=c(0,110))

#######################################################################
sink()
#######################################################################
