# ------------------------------------------------------- #
# --------------- Descriptive analyses ------------------ #
# ------------------------------------------------------- #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

sink(file="./results/descriptive_analyses.txt")

## load data
authors <- readRDS(file="./data/authors_raw.rds")
publications <- readRDS(file="./data/publications_raw.rds")
icc_df_all <- readRDS(file="./data/processed_data_all.rds")
icc_df <- readRDS(file="./data/processed_data_no_missing.rds")

cat("--------------------------------------------------\n\n")
cat("--------------------- Table 1 --------------------\n\n")
cat("--------------------------------------------------\n\n")


cat("-----------------Number of unique case and control authors of each gender----------------\n\n")

unique_authors <- icc_df_all[,c("auth_id","Gender","case",
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

require(gmodels)
cat("Missing gender variable by case/control status:\n")
CrossTable(unique_authors2$Gender,unique_authors2$status)

cat("-----------------Seniority metrics for unique authors by case status----------------\n\n")

cat("\n\n----Years since first publication----\n\n")
tapply(unique_authors2$years_in_scopus,unique_authors2$status,summary)
cat("All:\n\n")
summary(unique_authors2$years_in_scopus)

cat("\n\n----Number of publications----\n\n")
tapply(unique_authors2$Total_Publications_In_Scopus,unique_authors2$status,summary)
cat("All:\n\n")
summary(unique_authors2$Total_Publications_In_Scopus)

cat("\n\n----H-Index----\n\n")
tapply(unique_authors2$H_Index,unique_authors2$status,summary)
cat("All:\n\n")
summary(unique_authors2$H_Index)

cat("------------------------------------------------------\n\n")
cat("--------------------- Other stats --------------------\n\n")
cat("------------------------------------------------------\n\n")

cat("\n\n-----------------Gender distribution of cases and controls----------------\n\n")
cases <- icc_df[icc_df$case==1,c("Gender","auth_id")]
cases <- cases[!duplicated(cases$auth_id),]
cat("\n\n---Cases---\n\n")
CrossTable(cases$Gender)

controls <- icc_df[icc_df$case==0,c("Gender","auth_id")]
controls <- controls[!duplicated(controls$auth_id),]
cat("\n\n---Controls---\n\n")
CrossTable(controls$Gender)

cat("\n\n-----------------Gender of controls by gender of case----------------\n\n")


cat("\n\n-----------------Number of controls per case included in analysis----------------\n\n")
n_controls <- tapply(1-icc_df$case,icc_df$pub_id,sum)
summary(n_controls)

sink()



