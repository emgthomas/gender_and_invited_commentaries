# ---------------------------------------------------------- #
# ----------------- Sensitivity analyses ------------------- #
# ---------------------------------------------------------- #

#######################################################################
sink(file="./results/sensitivity_analyses.txt")
#######################################################################

# Packages
require(survival)
require(splines)
require(metafor)
require(dplyr)
require(ggplot2)
require(reshape2)
require(lme4)
require(mice)

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

cat("\n\n---Meta-regression on journal citescore---\n\n")
# Meta-regression on journal citescore
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
              name="95% prediction interval") %>%
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
export(OR_plot, "./results/OR_citescore_metareg.pdf")

cat("\n\n********** Adjusted model ************\n\n")

# Meta-analysis for adjusted model
all_2stage_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,data=outputs_select)
cat("\n\n---Summary of random effects meta-analysis for adjusted model---\n\n")
summary(all_2stage_adj)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage_adj, transf=exp, digits=2)

## Save key results for use in plotting ##
save(all_2stage,all_2stage_adj,file="./results/two_stage_analyses.Rdata")

# Meta-regression by citescore
cat("\n\n---Meta-regression on journal citescore---\n\n")
meta_analysis_cs_adj <- rma.mv(effect_adj,sd_adj^2,random=~1|journal,
                           mods= ~ ns(citescore,knots=knot_placement),
                           data=subset(outputs_select2,!is.na(outputs_select2$effect_adj)))
cat("\n\nMeta-regression on citescore\n\n")
summary(meta_analysis_cs_adj)

cat("\n\n---Sub-group analyses by journal topic---\n\n")
cat("See secondary_analyses.R (for ease of coding, sub-group analyses using one-stage 
and two-stage meta-analysis are performed together)\n\n")

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
export(OR_adj_plot, "./results/OR_adj_citescore_metareg.pdf")

cat("\n\n********** Model with interaction ************\n\n")

# Main effect
all_2stage_int <- rma.mv(effect_int,sd_int^2,random=~1|journal,data=outputs_select)
cat("\n\n---Summary of random effects meta-analysis for interaction with years active---\n\n")
cat("\n\n-- Main effect--\n\n")
summary(all_2stage_int)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage_int, transf=exp, digits=2)

# Interaction term
all_2stage_int2 <- rma.mv(effect_int2,sd_int2^2,random=~1|journal,data=outputs_select)
cat("\n\n-- Interaction term--\n\n")
summary(all_2stage_int2)
cat("\n\n---Mean and confidence (ci)/prediction (cr) intervals on odds ratio scale---\n\n")
predict(all_2stage_int2, transf=exp, digits=2)

cat("\n\n---------------------------------------------------------\n\n")
cat("\n\n--------------------- Missing data ----------------------\n\n")
cat("\n\n---------------------------------------------------------\n\n")

icc_df_all <- readRDS(file="./data/processed_data_all.rds")

# Assume all missing cases are female, all missing controls are men
# What impact will this have on results?s
icc_df_all$Gender[icc_df_all$case==1 & icc_df_all$Gender == "unknown"] <- "female"
icc_df_all$Gender[icc_df_all$case==0 & icc_df_all$Gender == "unknown"] <- "male"
icc_df_all$Gender <- factor(icc_df_all$Gender,levels=c("male","female"))

# Run analyses
cat("\n\n----------------- Results if we assume all missing cases are female, all missing controls are male ---------------\n\n")
all_1stage_miss <- clogit(case ~ Gender + strata(pub_id), data = icc_df_all)
cat("\n\n---Unadjusted analysis---\n")
summary(all_1stage_miss)

cat("\n\n---Adjusted analysis---\n")
knots <- c(2.5,5,7.5)
all_1stage_miss_adj <- clogit(case ~ Gender + ns(years_in_scopus_ptile,knots=knots) + 
                           ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                           strata(pub_id), data = icc_df_all)
summary(all_1stage_miss_adj)

cat("\n\n----------------- Including unknown gender in the regression model ---------------\n\n")

icc_df_all <- readRDS(file="./data/processed_data_all.rds")

all_1stage_miss2 <- clogit(case ~ Gender + strata(pub_id), data = icc_df_all)
cat("\n\n---Unadjusted analysis---\n")
summary(all_1stage_miss2)

cat("\n\n---Adjusted analysis---\n")
knots <- c(2.5,5,7.5)
all_1stage_miss_adj2 <- clogit(case ~ Gender + ns(years_in_scopus_ptile,knots=knots) + 
                                ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                                strata(pub_id), data = icc_df_all)
summary(all_1stage_miss_adj2)

cat("\n\n----------------- Multiple imputation ---------------\n\n")

icc_df_all <- readRDS(file="./data/processed_data_all.rds")

# How many imputed datasets?
n.impute <- 10

# Creata dataframe for running imputation model
icc_df_to_impute <- icc_df_all
icc_df_to_impute$gender <- 0
icc_df_to_impute$gender[icc_df_all$Gender=="female"] <- 1
icc_df_to_impute$gender[icc_df_to_impute$Gender=="unknown"] <- NA
icc_df_to_impute <- icc_df_to_impute[!(is.na(icc_df_to_impute$H_Index) |
                                         is.na(icc_df_to_impute$Total_Publications_In_Scopus) |
                                         is.na(icc_df_to_impute$years_in_scopus_ptile)),]
icc_df_to_impute <- icc_df_to_impute[,c("pub_sourceid","pub_id","case","gender","years_in_scopus_ptile","h_index_ptile","n_pubs_ptile")]
icc_df_to_impute$case_years_in_scopus <- icc_df_to_impute$case*icc_df_to_impute$years_in_scopus_ptile
icc_df_to_impute$case_h_index <- icc_df_to_impute$case*icc_df_to_impute$h_index_ptile
icc_df_to_impute$case_n_pubs <- icc_df_to_impute$case*icc_df_to_impute$n_pubs_ptile

# run imputation model
knots <- c(2.5,5,7.5)
mod_impute <- glmer(gender ~ case + 
               case_years_in_scopus +
               case_h_index +
               case_n_pubs +
               ns(years_in_scopus_ptile,knots=knots) + 
               ns(h_index_ptile,knots=knots) + 
               ns(n_pubs_ptile,knots=knots) +
               (1 | pub_id), 
               data = icc_df_to_impute, family = binomial,
               verbose = 2, nAGQ=0)
summary(mod_impute)

# Get predicted probabilities of being female for those with missing gender
p_female <- predict(mod_impute,icc_df_to_impute,
                    type="response",
                    allow.new.levels=TRUE)
# # plot predicted probabilities by gender
# ggplot(data.frame(x=factor(icc_df_to_impute$gender),y=p_female),aes(y=y,x=x)) + geom_boxplot()
# # plot predicted probabilities by missingness
# ggplot(data.frame(x=factor(is.na(icc_df_to_impute$gender)),y=p_female),aes(y=y,x=x)) + geom_boxplot()

# generate imputed datasets
imputed_dfs <- list()
set.seed(3298639)
missing_gender <- is.na(icc_df_to_impute$gender)
for(i in 1:n.impute){
  imputed_df <- icc_df_to_impute
  imputed_df$gender[missing_gender] <- sapply(p_female[missing_gender],rbinom,n=1,size=2)
  imputed_dfs[[i]] <- imputed_df
}
saveRDS(imputed_dfs,file="./data/imputed_data.RDS")

# combine results
imputed_mods <- lapply(imputed_dfs,
                       FUN=clogit,
                       formula=case ~ gender + 
                       strata(pub_id)
                      ) %>% as.mira

cat("\n\n-- Unadjusted results after multiple imputation --\n\n")
res_imp <- summary(pool(imputed_mods))
c(OR=exp(res_imp$estimate),
      ci.lb=exp(res_imp$estimate - 1.96*res_imp$std.error),
      ci.ub=exp(res_imp$estimate + 1.96*res_imp$std.error),
  p.value=res_imp$p.value)

imputed_mods_adj <- lapply(imputed_dfs,
            FUN=clogit,
            formula=case ~ gender +
                     ns(years_in_scopus_ptile,knots=knots) + 
                     ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                     strata(pub_id)
              ) %>% as.mira

cat("\n\n-- Adjusted results after multiple imputation --\n\n")
res_imp_adj <- summary(pool(imputed_mods_adj))
res_imp_adj_df <- data.frame(
  OR=exp(res_imp_adj$estimate),
  ci.lb=exp(res_imp_adj$estimate - 1.96*res_imp_adj$std.error),
  ci.ub=exp(res_imp_adj$estimate + 1.96*res_imp_adj$std.error),
  p.value=res_imp_adj$p.value)
row.names(res_imp_adj_df) <- row.names(res_imp_adj)
res_imp_adj_df

imputed_mods_int <- lapply(imputed_dfs,
                           FUN=clogit,
                           formula=case ~ gender + years_in_scopus_ptile:gender + 
                             ns(years_in_scopus_ptile,knots=knots) + 
                             ns(h_index_ptile,knots=knots) + ns(n_pubs_ptile,knots=knots) + 
                             strata(pub_id)
                    ) %>% as.mira

cat("\n\n-- Adjusted results with interaction after multiple imputation --\n\n")
res_imp_int <- summary(pool(imputed_mods_int))
res_imp_int_df <- data.frame(
  OR=exp(res_imp_int$estimate),
  ci.lb=exp(res_imp_int$estimate - 1.96*res_imp_int$std.error),
  ci.ub=exp(res_imp_int$estimate + 1.96*res_imp_int$std.error),
  p.value=res_imp_int$p.value)
row.names(res_imp_int_df) <- row.names(res_imp_int)
res_imp_int_df

cat("\n\n---------------------------------------------------------------------------------\n\n")
cat("\n\n--------------------- Stricter match criteria for controls ----------------------\n\n")
cat("\n\n---------------------------------------------------------------------------------\n\n")

icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

######### Removing controls with low match scores ##########

all_1stage_int <- clogit(case ~ Gender + 
                           years_in_scopus_ptile:factor(Gender,levels=c("female","male")) + 
                           ns(years_in_scopus_ptile,knots=knots) + 
                           ns(h_index_ptile,knots=knots) + 
                           ns(n_pubs_ptile,knots=knots) + 
                           strata(pub_id), data = icc_df)
summary(all_1stage_int)

######### Removing controls with low match scores ##########

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

cat("\n\n--------------------------------------------------------------------------------------\n\n")
cat("\n\n------------- De-duplicating case authors present in multiple journals ---------------\n\n")
cat("\n\n--------------------------------------------------------------------------------------\n\n")

icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

# de-duplicate case authors
icc_df_case <- icc_df[icc_df$case==1,]
pubs.incl <- unique(icc_df_case$pub_id[!duplicated(icc_df_case$auth_id)])
icc_df_dedup <- icc_df[icc_df$pub_id %in% pubs.incl,]
icc_df_dedup$pub_id <- factor(icc_df_dedup$pub_id,levels=unique(icc_df_dedup$pub_id))

# keep only top two controls
icc_df_dedup <- icc_df_dedup[icc_df_dedup$match_score_rank <= 3,]

cat("---Unadjusted analysis---\n")
all_1stage_dedup <- clogit(case ~ Gender + strata(pub_id), data = icc_df_dedup)
summary(all_1stage_dedup)

cat("\n\n---Adjusted for measures of seniority using natural cubic splines---\n")
knots <- c(2.5,5,7.5)
all_1stage_dedup_adj <- clogit(case ~ Gender + 
                           ns(years_in_scopus_ptile,knots=knots) + 
                           ns(h_index_ptile,knots=knots) + 
                           ns(n_pubs_ptile,knots=knots) + 
                           strata(pub_id), data = icc_df_dedup)
summary(all_1stage_dedup_adj)

cat("\n\n---Including interaction term between gender and years active---\n")
all_1stage_dedup_int <- clogit(case ~ Gender + 
                                 years_in_scopus_ptile:factor(Gender,levels=c("female","male")) +
                                 ns(years_in_scopus_ptile,knots=knots) + 
                                 ns(h_index_ptile,knots=knots) + 
                                 ns(n_pubs_ptile,knots=knots) + 
                                 strata(pub_id), data = icc_df_dedup)
summary(all_1stage_dedup_int)

#######################################################################
sink()
#######################################################################
