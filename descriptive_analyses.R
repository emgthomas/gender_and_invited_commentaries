# ------------------------------------------------------- #
# --------------- Descriptive analyses ------------------ #
# ------------------------------------------------------- #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

sink(file="./results/descriptive_analyses.txt")

## load data
authors <- readRDS(file="./data/authors_raw.rds")
publications <- readRDS(file="./data/publications_raw.rds")

cat("-----------------Number of unique case and control authors of each gender----------------\n\n")
unique_authors <- authors[,c("AuthorID","Gender","Designation")]
unique_authors <- unique_authors[!duplicated(unique_authors),]
controls <- unique(unique_authors$AuthorID[unique_authors$Designation=="Control"])
cases <- unique(unique_authors$AuthorID[unique_authors$Designation=="Case"])
both <- intersect(controls,cases)
cat(length(both),"authors acted as both cases and controls\n\n")

# actual unique authors, de-duplicating those that act as both case and control
unique_authors2 <- unique_authors[,c("AuthorID","Gender")]
unique_authors2 <- unique_authors2[!duplicated(unique_authors2),]
unique_authors2$status <- "Control"
unique_authors2$status[unique_authors2$AuthorID %in% cases] <- "Case"
unique_authors2$status[unique_authors2$AuthorID %in% both] <- "Both"
unique_authors2$status <- factor(unique_authors2$status,levels=c("Case","Control","Both"))

require(gmodels)
cat("Missing gender variable by case/control status:\n")
CrossTable(unique_authors2$status,unique_authors2$Gender)

cat("----------------------Number of cases and controls--------------------\n\n")
cat(sum(authors$case),"ICCs were matched to at least one control\n\n")
cat("Found",sum(authors$case==0),"controls\n\n")
n.controls <- aggregate(authors[,"case"],by=list(authors$Publication),FUN = function(x) sum(1-x))
cat(sum(n.controls$x == 10)/nrow(n.controls),"% of cases had 10 matched controls\n\n\n\n")
print("Number of case and control authors:")
table(unique_authors$Designation)

cat(length(unique(publications$citing)),"potential case publications\n\n")
cat("Found at least one matched control for",length(unique(authors$Publication)),"of these, or",
    length(unique(authors$Publication))/length(unique(publications$citing)),"%")

cat("\n",length(unique(publications$pub_source_title)),"journals included before processing")

cat("\n\n------------Descriptive analyses with missing data removed (processed data)----------------")

icc_df <- readRDS("./data/processed_data_no_missing.rds")

cat("\n\nFinal numbers of cases and controls:\n")
table(icc_df$Designation)
cat("\n")
n.controls <- aggregate(icc_df[,"case"],by=list(icc_df$Publication),FUN = function(x) sum(1-x))
cat("Summary of number of controls per case:\n")
summary(n.controls$x)
cat("\n")
cat("Fraction of cases with 10 controls:",sum(n.controls$x == 10)/nrow(n.controls),"\n\n")
cat("Fraction of cases with >8 controls:",sum(n.controls$x >= 9)/nrow(n.controls),"\n\n")
cat("Fraction of cases with >7 controls:",sum(n.controls$x >= 8)/nrow(n.controls),"\n\n")

# How many unique authors?
cat(length(unique(icc_df$AuthorID)),"unique authors in study (cases and controls)\n\n")
cases <- unique(icc_df$AuthorID[icc_df$case==1])
controls <- unique(icc_df$AuthorID[icc_df$case==0])
both <- intersect(cases,controls)
cases.only <- cases[!(cases %in% controls)]
cat(length(cases.only),"authors were cases only\n\n")
controls.only <- controls[!(controls %in% cases)]
cat(length(controls.only),"authors were controls only\n\n")
cat(length(both),"authors were both cases and controls\n\n")

cat("\n",length(unique(publications$pub_source_title[publications$citing %in% icc_df$Publication])),"journals included after processing")

sink()



