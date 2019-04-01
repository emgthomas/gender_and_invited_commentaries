# ---------------------------------------------------------- #
# --------------------- Data cleaning ---------------------- #
# ---------------------------------------------------------- #

setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/jama_paper/")

# Read in author data
authors <- read.csv(file="./data/Author_Gender_MatchScores_SENIORITY.csv")
authors$case <- as.numeric(authors$Designation == "Case")

# Categorize match score into deciles
require(gtools)
authors$match_quantile <- numeric(length=nrow(authors)) + NA
# find quantiles of match score for controls
authors$match_quantile[authors$case==0] <- quantcut(authors$Match_Score[authors$case==0],q=10)
# set "quantile" to 10 for cases, so they will always be included
authors$match_quantile[authors$case==1] <- 10

# # Which controls are also focal article authors
# controls_focal <- read.csv(file="./data/Control_Equals_Cited.csv")
# controls_focal$focal_author <- 1
# controls_focal <- controls_focal[!duplicated(controls_focal[,c("Cited_Author","Publication")]),]
# controls_focal <- subset(controls_focal,Cited_Author %in% unique(authors$AuthorID))
# authors <- merge(authors,controls_focal[,c("Cited_Author","Publication","focal_author")],
#                  by.x=c("AuthorID","Publication"),by.y=c("Cited_Author","Publication"),
#                  all.x=T)
# authors$focal_author[is.na(authors$focal_author) & authors$case ==0] <- 0

# Read in publication data
publications <- read.csv(file="./data/Citing_Cited_Info_total.csv")
# publication field in authors can be linked to citing field in publications

# Read in journal topics
journal_topics <- read.csv(file="./data/Journal_Subjects.csv")
topic_names <- read.csv(file="./data/ASJC Codes with levels.csv",
                        sep=";") # cloned from github.com/plreyes/Scopus.git
journal_topics <- merge(x=journal_topics,y=topic_names,
      by.x="Subject",by.y="Code",
      all.x=T)

# Create dataframe of topics by journal
require(reshape2)
journal_topics_low <- dcast(journal_topics,pub_source_title ~ Low,fun.aggregate = length, value.var="Low")
journal_topics_middle <- dcast(journal_topics,pub_source_title ~ Middle,fun.aggregate = length, value.var="Low")
journal_topics_top <- dcast(journal_topics,pub_source_title ~ Top,fun.aggregate = length, value.var="Low")

# Create data frame of unique intra-citing pubs and their journals
# This removes duplicate "citing" articles that cited >1 other article in the same journal/issue
publications2 <- publications[!duplicated(publications$citing),]

# Merge journal titles into authors dataset
authors2 <- merge(x=authors,y=publications2,by.x="Publication",by.y="citing",all.x=T,all.y=F,sort=F)

# Create factor variables
authors2$Gender <- factor(authors2$Gender, levels=c("male","female","unknown"))
authors2$Publication <- factor(authors2$Publication, levels=unique(authors2$Publication))

# Remove controls that are also cases (within the same journal)
# and remove duplicate case authors within the same journal
# So, within journals, cases cannot also be controls, and each author can be a case only once
controls <- authors2[authors2$case==0,] # separate out controls
cases <- authors2[authors2$case==1,] # separate out cases
controls.keep <- logical(length=nrow(controls))
cases.keep <- logical(length=nrow(cases))
for(pub in unique(controls$pub_source_title)){
  which_controls <- controls$pub_source_title == pub
  which_cases <- cases$pub_source_title == pub
  # within journals, remove any controls that are also cases
  controls.keep[which_controls] <- !(controls$AuthorID[which_controls] %in% unique(cases$AuthorID[which_cases]))
  # within journals, remove any duplicate cases
  cases.keep[which_cases] <- !duplicated(cases$AuthorID[which_cases])
}
# put cases and controls back together
authors3 <- rbind(cases[cases.keep,],controls[controls.keep,])

# Include only matched sets that have a case and at least one control with gender designation
# i.e., remove invalid matched sets
checkFun <- function(status){
  if(length(status)==1) return(FALSE)
  if(sum(status) == 1 & length(status) > 1){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
included.pubs <- tapply(authors3$case,authors3$Publication,FUN=checkFun,simplify = T)
included.pubs <- data.frame(Publication=names(included.pubs),include=included.pubs)
authors4 <- merge(authors3,included.pubs,by="Publication",all.x=T)

### put this together in final dataset
# "icc" stands for "intra-citing commentary"
icc_df_all <- authors4[authors4$include,]
icc_df_all$include <- NULL
icc_df_all$Publication <- factor(icc_df_all$Publication,levels=unique(icc_df_all$Publication))
icc_df_all$AuthorID <- factor(icc_df_all$AuthorID,levels=unique(icc_df_all$AuthorID))

# Remove authors with no gender designation
icc_df <- icc_df_all[icc_df_all$Gender != "unknown",]
icc_df$Gender <- factor(icc_df$Gender,levels=c("male","female"))
included.pubs <- tapply(icc_df$case,icc_df$Publication,FUN=checkFun,simplify = T)
included.pubs <- data.frame(Publication=names(included.pubs),include=included.pubs)
icc_df <- merge(icc_df,included.pubs,by="Publication",all.x=T)
icc_df <- icc_df[icc_df$include,]
icc_df$include <- NULL

###### Identify journals for which we can get a valid journal-level OR estimate ########
min_num_matched_sets <- 2
check_genders <- function(gender) length(unique(gender))==2
check_matched_sets <- function(idx,publications,gender,case,min_num_matched_sets){
  publications <- publications[idx]
  gender <- gender[idx]
  case <- case[idx]
  # which matched sets are not all the same gender?
  out1 <- tapply(gender,as.character(publications),check_genders)
  # after excluding matched sets that are all the same gender, there any "zero cells"?
  out2 <- sum(table(gender[publications %in% names(out1)[out1]],case[publications %in% names(out1)[out1]])==0)
  return(sum(out1)>=min_num_matched_sets & out2==0)
}

included.journals <- tapply(1:nrow(icc_df),icc_df$pub_source_title,check_matched_sets,
                            publications=icc_df$Publication,
                            gender=icc_df$Gender,
                            case=icc_df$case,min_num_matched_sets=min_num_matched_sets)
included.journals <- data.frame(pub_source_title=names(included.journals),include.journal=included.journals)
icc_df <- merge(icc_df,included.journals,by="pub_source_title")
icc_df_all <- merge(icc_df_all,included.journals,by="pub_source_title")

# "First year in Scopus" covariate- a proxy for author seniority
icc_df$seniority <- 2019 - icc_df$First_Year_in_Scopus
seniority_quintiles <- quantile(icc_df$seniority[!duplicated(icc_df$AuthorID)],c(0.2,0.4,0.6,0.8),na.rm=T)
icc_df$seniority_quintile <- NA
icc_df$seniority_quintile[!is.na(icc_df$seniority)] <- 0
icc_df$seniority_quintile[icc_df$seniority > seniority_quintiles[1] & icc_df$seniority <= seniority_quintiles[2]] <- 1
icc_df$seniority_quintile[icc_df$seniority > seniority_quintiles[2] & icc_df$seniority <= seniority_quintiles[3]] <- 2
icc_df$seniority_quintile[icc_df$seniority > seniority_quintiles[3] & icc_df$seniority <= seniority_quintiles[4]] <- 3
icc_df$seniority_quintile[icc_df$seniority > seniority_quintiles[4]] <- 4

############## Save data #############
saveRDS(authors,file="./data/authors_raw.rds")
saveRDS(publications,file="./data/publications_raw.rds")
saveRDS(icc_df_all,file="./data/processed_data_all.rds")
saveRDS(icc_df,file="./data/processed_data_no_missing.rds")
saveRDS(journal_topics_low,file="./data/journal_topics_low.rds")
saveRDS(journal_topics_middle,file="./data/journal_topics_middle.rds")
saveRDS(journal_topics_top,file="./data/journal_topics_top.rds")
