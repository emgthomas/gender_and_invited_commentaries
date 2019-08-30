# ------------------------------------------------------- #
# --------------- Descriptive analyses ------------------ #
# ------------------------------------------------------- #

#######################################################################
sink(file="./results/descriptive_analyses.txt")
#######################################################################

## packages
require(gmodels)

## load data
authors <- readRDS(file="./data/authors_raw.rds")
publications <- readRDS(file="./data/publications_raw.rds")
icc_df_all <- readRDS(file="./data/processed_data_all.rds")
icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

cat("--------------------------------------------------\n\n")
cat("--------------------- Table 1 --------------------\n\n")
cat("--------------------------------------------------\n\n")

cat("-----------------Author-level variables for *unique* authors included in analysis----------------\n\n")

unique_authors <- icc_df[,c("auth_id","Gender","case",
                                "years_in_scopus","Total_Publications_In_Scopus",
                                "H_Index")]
unique_authors <- unique_authors[!duplicated(unique_authors),]
controls <- unique(unique_authors$auth_id[unique_authors$case==0])
cases <- unique(unique_authors$auth_id[unique_authors$case==1])
both <- intersect(controls,cases)

# actual unique authors, de-duplicating those that act as both case and control
unique_authors2 <- unique_authors[,c("auth_id","Gender",
                                     "years_in_scopus","Total_Publications_In_Scopus",
                                     "H_Index")]
unique_authors2 <- unique_authors2[!duplicated(unique_authors2),]
unique_authors2$status <- "Control"
unique_authors2$status[unique_authors2$auth_id %in% cases] <- "Case"
unique_authors2$status[unique_authors2$auth_id %in% both] <- "Both"
unique_authors2$status <- factor(unique_authors2$status,levels=c("Case","Control","Both"))

cat("\n\n----Gender----\n\n")
CrossTable(unique_authors2$Gender,unique_authors2$status,
           prop.r=F,prop.c=T,prop.t=F,prop.chisq=F)

cat("\n\n----Quartiles of years since first publication----\n\n")
tapply(unique_authors2$years_in_scopus,unique_authors2$status,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
cat("\nAll:\n\n")
quantile(unique_authors2$years_in_scopus,probs=c(0.25,0.5,0.75),na.rm=T)

cat("\n\n----Number of publications----\n\n")
tapply(unique_authors2$Total_Publications_In_Scopus,unique_authors2$status,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
cat("\nAll:\n\n")
quantile(unique_authors2$Total_Publications_In_Scopus,probs=c(0.25,0.5,0.75),na.rm=T)

cat("\n\n----H-Index----\n\n")
tapply(unique_authors2$H_Index,unique_authors2$status,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
cat("\nAll:\n\n")
quantile(unique_authors2$H_Index,probs=c(0.25,0.5,0.75),na.rm=T)

cat("\n\n-----------------Missing covariate data----------------\n\n")
missing_cov <- is.na(unique_authors2$years_in_scopus) | is.na(unique_authors2$Total_Publications_In_Scopus) | is.na(unique_authors2$H_Index)
CrossTable(unique_authors2$status,missing_cov,
           prop.r=T,prop.c=F,prop.t=F,prop.chisq=F)

icc_df_all <- readRDS(file="./data/processed_data_all.rds")

cat("--------------------------------------------------\n\n")
cat("--------------------- eTable 1 -------------------\n\n")
cat("--------------------------------------------------\n\n")

cat("-----------------Author-level variables for *unique* authors by gender, including those with missing gender----------------\n\n")

unique_authors <- icc_df_all[,c("auth_id","Gender","case",
                                "years_in_scopus","Total_Publications_In_Scopus",
                                "H_Index","asia")]
unique_authors <- unique_authors[!duplicated(unique_authors),]
controls <- unique(unique_authors$auth_id[unique_authors$case==0])
cases <- unique(unique_authors$auth_id[unique_authors$case==1])
both <- intersect(controls,cases)

# actual unique authors, de-duplicating those that act as both case and control
unique_authors2 <- unique_authors[,c("auth_id","Gender",
                                     "years_in_scopus","Total_Publications_In_Scopus",
                                     "H_Index","asia")]
unique_authors2 <- unique_authors2[!duplicated(unique_authors2),]
unique_authors2$status <- "Control"
unique_authors2$status[unique_authors2$auth_id %in% cases] <- "Case"
unique_authors2$status[unique_authors2$auth_id %in% both] <- "Both"
unique_authors2$status <- factor(unique_authors2$status,levels=c("Case","Control","Both"))

cat("\n\n----Gender----\n\n")
CrossTable(unique_authors2$status,unique_authors2$Gender,
           prop.r=F,prop.c=T,prop.t=F,prop.chisq=F)

cat("\n\n----Asian country of origin----\n\n")
unique_authors2$country_of_origin <- "Asian"
unique_authors2$country_of_origin[unique_authors2$asia==0] <- "Not Asian"
unique_authors2$country_of_origin[is.na(unique_authors2$asia)] <- "Unknown"
CrossTable(unique_authors2$country_of_origin,unique_authors2$Gender,
           prop.r=F,prop.c=T,prop.t=F,chisq=F,prop.chisq=F)

cat("\n\n----Quartiles of years since first publication----\n\n")
tapply(unique_authors2$years_in_scopus,unique_authors2$Gender,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
cat("\nAll:\n\n")
quantile(unique_authors2$years_in_scopus,probs=c(0.25,0.5,0.75),na.rm=T)

cat("\n\n----Number of publications----\n\n")
tapply(unique_authors2$Total_Publications_In_Scopus,unique_authors2$Gender,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
cat("\nAll:\n\n")
quantile(unique_authors2$Total_Publications_In_Scopus,probs=c(0.25,0.5,0.75),na.rm=T)

cat("\n\n----H-Index----\n\n")
tapply(unique_authors2$H_Index,unique_authors2$Gender,quantile,probs=c(0.25,0.5,0.75),na.rm=T)
cat("\nAll:\n\n")
quantile(unique_authors2$H_Index,probs=c(0.25,0.5,0.75),na.rm=T)

cat("\n\n------------------------------------------------------\n\n")
cat("--------------------- Other stats --------------------\n\n")
cat("------------------------------------------------------\n\n")

cat("\n\n-----------------Gender distribution of ICC authors----------------\n\n")
authors <- readRDS(file="./data/authors_raw.rds")
publications <- readRDS(file="./data/publications_raw.rds")
cases <- authors[authors$case==1,]
publications <- merge(publications,cases,by="pub_id",all.x=T,all.y=F)
CrossTable(publications$Gender)

cat("\n\n---Number of female authors for all ICCs, excluding unknown gender---\n\n")
publications <- publications[publications$Gender != "unknown",]
CrossTable(publications$Gender)

cat("\n\n-----------------Number of female authors among unique case authors, excluding unknown gender----------------\n\n")
icc_df <- readRDS("./data/processed_data_no_missing.rds")
cases <- icc_df[icc_df$case==1,]
cases <- cases[!duplicated(cases$auth_id),]
CrossTable(cases$Gender)

cat("\n\n---Proportion of cases who authored ICCs in multiple journals, by gender---\n\n")
cases <- icc_df[icc_df$case==1,c("Gender","auth_id")]
cases$auth_id <- factor(cases$auth_id)
cases$count <- 1
cases_n_journals <- tapply(X=cases$count,INDEX=cases$auth_id,FUN=sum)
cases_n_journals <- data.frame(auth_id=names(cases_n_journals),n_journals=cases_n_journals)
cases_n_journals <- merge(cases_n_journals,cases[!duplicated(cases$auth_id),c("Gender","auth_id")],by="auth_id",all.X=T,all.y=F)
tapply(cases_n_journals$n_journals>1,cases_n_journals$Gender,mean)
cat("\n\nAll\n")
mean(cases_n_journals$n_journals>1)
sum(cases_n_journals$n_journals>1)

cat("\n\n---Fraction and number of articles that had replies, responds, etc in the title---\n")
publications <- readRDS(file="./data/publications_raw.rds")
mean(publications$replies)
sum(publications$replies)

cat("\n\n-----------------Gender of controls by gender of case, excluding unknown gender----------------\n\n")
female_case_pubs <- icc_df[icc_df$case==1 & icc_df$Gender=="female",]$pub_id
female_cases_controls <- icc_df[icc_df$case==0 & icc_df$pub_id %in% female_case_pubs,]
female_cases_av_controls <- tapply(female_cases_controls$Gender=="female",
                                 factor(female_cases_controls$pub_id),
                                 mean)
cat("\n\n---Number of matched sets with female cases---\n\n")
length(female_cases_av_controls)
cat("\n\n---Average number of female controls for female cases---\n\n")
mean(female_cases_av_controls)

male_case_pubs <- icc_df[icc_df$case==1 & icc_df$Gender=="male",]$pub_id
male_cases_controls <- icc_df[icc_df$case==0 & icc_df$pub_id %in% male_case_pubs,]
male_cases_av_controls <- tapply(male_cases_controls$Gender=="female",
                                 factor(male_cases_controls$pub_id),
                                        mean)
cat("\n\n---Number of matched sets with male cases---\n\n")
length(male_cases_av_controls)
cat("\n\n---Average number of female controls for male cases---\n\n")
mean(male_cases_av_controls)

cat("\n\n-----------------Number of controls per case included in analysis----------------\n\n")
n_controls <- tapply(1-icc_df$case,icc_df$pub_id,sum)
summary(n_controls)

cat("\n\n---Summary of number of journals, by gender---\n")
tapply(cases_n_journals$n_journals,cases_n_journals$Gender,summary)
cat("\n\nAll\n")
summary(cases_n_journals$n_journals)

cat("\n\n-----------------Correlation of years active, h-index and number of pubs----------------\n\n")
cat("\n\n---Years in scopus percentile vs. h-index percentile---\n\n")
cor(icc_df$years_in_scopus_ptile,icc_df$h_index_ptile,use="complete.obs",method="spearman")
cat("\n\n---Years in scopus percentile vs. number of pubs percentile---\n\n")
cor(icc_df$years_in_scopus_ptile,icc_df$n_pubs_ptile,use="complete.obs",method="spearman")
cat("\n\n---H-index percentile vs. number of pubs percentile---\n\n")
cor(icc_df$h_index_ptile,icc_df$n_pubs_ptile,use="complete.obs",method="spearman")


#######################################################################
sink()
#######################################################################


