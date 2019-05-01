# ---------------------------------------------------------- #
# --------------------- Data cleaning ---------------------- #
# ---------------------------------------------------------- #

#######################################################################
sink(file="./results/data_cleaning.txt")
#######################################################################

# packages
require(data.table)
require(dplyr)
require(gtools)
require(reshape2)

# Read in publications data
publications <- readRDS(file="./data/publications_raw.rds")

cat("------- Flow chart Step 1 (Included)-------\n\n")
cat("Number of included journals was determined in Scopus database (not shown here)")

cat("\n\n------- Flow chart Step 1 (Excluded)-------\n\n")
cat("Number of journals that did not contain ICC articles = ",4412 + 56 - length(unique(publications$pub_sourceid)))

cat("\n\n------- Flow chart Step 2 (Included)------- \n\n")
cat("Number of ICC articles = ",nrow(publications))

# Read in author data
authors <- readRDS(file="./data/authors_raw.rds")

cat("\n\n------- Flow chart Step 2 (Excluded)------- \n\n")
cat("Number of ICC articles where corresponding/single author could not be identified = ",nrow(publications)-length(unique(authors$pub_id[authors$case==1])))

cat("\n\n------- Flow chart Step 3 (Included)------- \n\n")
cat("Number of ICC articles with corresponding/single author = ",length(unique(authors$pub_id[authors$case==1])))

# Merge journal source IDs into authors dataset
authors2 <- merge(x=authors,y=publications,by.x="pub_id",all.x=T,all.y=F,sort=F)

cat("\n\n------- Flow chart Step 4 (Included)------- \n\n")
cat("Number of unique corresponding author-journal pairs (these will form cases) = ",sum(!duplicated(authors2[authors2$case==1,c("auth_id","pub_sourceid")])))

cat("\n\n------- Flow chart Step 5 (Included)------- \n\n")
cat("Number of controls = ",sum(1-authors2$case))

# Discard controls in bottom quartile of Match Score
ms_1st_quartile <- quantile(authors2$Match_Score[authors2$case==0],0.25)
authors3 <- authors2[authors2$Match_Score > ms_1st_quartile,]

# Remove controls that are also cases (within the same journal)
# and remove duplicate case authors within the same journal
# So, within journals, cases cannot also be controls, and each author can be a case only once
controls <- authors3[authors3$case==0,] # separate out controls
cases <- authors3[authors3$case==1,] # separate out cases
controls.keep <- logical(length=nrow(controls))
cases.keep <- logical(length=nrow(cases))
i <- 0
for(pub in unique(controls$pub_sourceid)){
  i <- i+1
  which_controls <- controls$pub_sourceid == pub
  which_cases <- cases$pub_sourceid == pub
  # within journals, remove any controls that are also cases
  controls.keep[which_controls] <- !(controls$auth_id[which_controls] %in% unique(cases$auth_id[which_cases]))
  # within journals, remove any duplicate cases
  cases.keep[which_cases] <- !duplicated(cases$auth_id[which_cases])
  #if((i %% 100)==0) print(i)
}
# put cases and controls back together
authors4 <- rbind(cases[cases.keep,],controls[controls.keep,])

cat("\n\n------- Flow chart Step 5 (Excluded)------- \n\n")
cat("Number of controls excluded = ",sum(1-authors2$case) - sum(1-authors4$case))
#cat("Number of unique control authors excluded = ",length(unique(authors2$auth_id[authors2$case==0]))-length(unique(authors4$auth_id[authors4$case==0])))

# Keep only top 10 controls
authors5 <- authors4 %>%
  group_by(pub_id) %>%
  mutate(match_score_rank = frank(-Match_Score))
authors5 <- authors5[order(authors5$pub_id,authors5$Match_Score,decreasing=T),]
# View(authors5[,c("pub_id","case","Match_Score","match_score_rank")])
authors5 <- authors5[authors5$match_score_rank <= 11, ] # keep 11 because the top match is always the case themself

cat("\n\n------- Flow chart Step 6 (Included)------- \n\n")
cat("Number of controls included after keeping only top 10 = ",sum(1-authors5$case))
#cat("Number of unique control authors included after keeping only top 10 = ",length(unique(authors5$auth_id[authors5$case==0])))

# Include only matched sets that have a case and at least one control
# i.e., remove invalid matched sets
checkFun <- function(status){
  if(length(status)==1) return(FALSE)
  if(sum(status) == 1 & length(status) > 1){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
authors5$pub_id <- factor(authors5$pub_id,levels=unique(authors5$pub_id))
included.pubs <- tapply(authors5$case,authors5$pub_id,FUN=checkFun,simplify = T)
included.pubs <- data.frame(pub_id=names(included.pubs),include=included.pubs)
authors6 <- merge(authors5,included.pubs,by="pub_id",all.x=T)
authors6 <- authors6[authors6$include,]
authors6$include <- NULL

cat("\n\n------- Flow chart Step 7 (Included)------- \n\n")
cat("Number of matched sets = ",sum(authors6$case))

# For remaining controls, categorize match score into deciles
authors6$match_quantile <- numeric(length=nrow(authors6)) + NA
# find quantiles of match score for controls
authors6$match_quantile[authors6$case==0] <- quantcut(authors6$Match_Score[authors6$case==0],q=10)
# set "quantile" to 11 for cases, so they will always be included
authors6$match_quantile[authors6$case==1] <- 11

# # Which controls are also focal article authors
# controls_focal <- read.csv(file="./data/Control_Equals_Cited.csv")
# controls_focal$focal_author <- 1
# controls_focal <- controls_focal[!duplicated(controls_focal[,c("Cited_Author","Publication")]),]
# controls_focal <- subset(controls_focal,Cited_Author %in% unique(authors$AuthorID))
# authors <- merge(authors,controls_focal[,c("Cited_Author","Publication","focal_author")],
#                  by.x=c("AuthorID","Publication"),by.y=c("Cited_Author","Publication"),
#                  all.x=T)
# authors$focal_author[is.na(authors$focal_author) & authors$case ==0] <- 0

### put this together in final dataset
# "icc" stands for "intra-citing commentary"
icc_df_all <- authors6

# Remove authors with no gender designation
icc_df <- icc_df_all[icc_df_all$Gender != "unknown",]
icc_df$Publication <- factor(icc_df$pub_id,levels=unique(icc_df$pub_id))
included.pubs <- tapply(icc_df$case,icc_df$pub_id,FUN=checkFun,simplify = T)
included.pubs <- data.frame(pub_id=names(included.pubs),include=included.pubs)
icc_df <- merge(icc_df,included.pubs,by="pub_id",all.x=T)
icc_df <- icc_df[icc_df$include,]
icc_df$include <- NULL

cat("\n\n------- Flow chart Step 7 (Excluded)------- \n\n")
cat("Matched sets excluded due to missing gender = ",sum(authors6$case) - sum(icc_df$case))

cat("\n\n------- Flow chart Step 8 (Included)------- \n\n")
cat("Matched sets with complete gender information = ",sum(icc_df$case))

cat("\n\n------- Flow chart Step 9 (Included)------- \n\n")
cat("Number of journals in one-stage meta-analysis = ",length(unique(icc_df$pub_sourceid)))

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

icc_df$Publication <- factor(icc_df$pub_id,levels=unique(icc_df$pub_id))
icc_df$pub_sourceid <- factor(icc_df$pub_sourceid,levels=unique(icc_df$pub_sourceid))
icc_df$Gender <- factor(icc_df$Gender,unique(icc_df$Gender))
included.journals <- tapply(1:nrow(icc_df),
                            icc_df$pub_sourceid,
                            check_matched_sets,
                            publications=icc_df$pub_id,
                            gender=icc_df$Gender,
                            case=icc_df$case,min_num_matched_sets=min_num_matched_sets)
included.journals <- data.frame(pub_sourceid=names(included.journals),include.journal=included.journals)
icc_df <- merge(icc_df,included.journals,by="pub_sourceid")

cat("\n\n------- Flow chart Step 9 (Excluded)------- \n\n")
cat("Journals with insufficient data for journal-specific estimate = ",length(included.journals$include.journal) - sum(included.journals$include.journal))

cat("\n\n------- Flow chart Step 10 (Included)------- \n\n")
cat("Journals included in two-stage meta-analysis = ",sum(included.journals$include.journal))
cat("\n\nMatched included in two-stage meta-analysis = ",sum(icc_df$case[icc_df$include.journal]))

# Create factor variables
icc_df$Gender <- factor(icc_df$Gender, levels=c("male","female"))
icc_df$pub_id <- factor(icc_df$pub_id, levels=unique(icc_df$pub_id))
icc_df$auth_id <- factor(icc_df$auth_id,levels=unique(icc_df$auth_id))

icc_df_all$Gender <- factor(icc_df_all$Gender, levels=c("male","female","unknown"))
icc_df_all$pub_id <- factor(icc_df_all$pub_id, levels=unique(icc_df_all$pub_id))
icc_df_all$auth_id <- factor(icc_df_all$auth_id,levels=unique(icc_df_all$auth_id))

# Measures of author seniority
icc_df$years_in_scopus <- 2019 - icc_df$First_Year_in_Scopus
icc_df_all$years_in_scopus <- 2019 - icc_df_all$First_Year_in_Scopus
icc_df$years_in_scopus_quintile <- as.numeric(quantcut(icc_df$years_in_scopus,q=5))
icc_df$h_index_quintile <- as.numeric(quantcut(icc_df$H_Index,q=5))
icc_df$n_pubs_quintile <- as.numeric(quantcut(icc_df$Total_Publications_In_Scopus,q=5))

#### Compute percentiles for seniority measures ####

# icc_df
YiS <- icc_df[,c("years_in_scopus","auth_id")]
YiS <- YiS[!duplicated(YiS),]
YiS$years_in_scopus_ptile <- percent_rank(YiS$years_in_scopus)/0.1
icc_df <- merge(icc_df,YiS[,c("auth_id","years_in_scopus_ptile")],by="auth_id",all.x=T,all.y=F)

h_index <- icc_df[,c("H_Index","auth_id")]
h_index <- h_index[!duplicated(h_index),]
h_index$h_index_ptile <- percent_rank(h_index$H_Index)/0.1
icc_df <- merge(icc_df,h_index[,c("auth_id","h_index_ptile")],by="auth_id",all.x=T,all.y=F)

n_pubs <- icc_df[,c("Total_Publications_In_Scopus","auth_id")]
n_pubs <- n_pubs[!duplicated(n_pubs),]
n_pubs$n_pubs_ptile <- percent_rank(n_pubs$Total_Publications_In_Scopus)/0.1
icc_df <- merge(icc_df,n_pubs[,c("auth_id","n_pubs_ptile")],by="auth_id",all.x=T,all.y=F)

# icc_df_all
YiS <- icc_df_all[,c("years_in_scopus","auth_id")]
YiS <- YiS[!duplicated(YiS),]
YiS$years_in_scopus_ptile <- percent_rank(YiS$years_in_scopus)/0.1
icc_df_all <- merge(icc_df_all,YiS[,c("auth_id","years_in_scopus_ptile")],by="auth_id",all.x=T,all.y=F)

h_index <- icc_df_all[,c("H_Index","auth_id")]
h_index <- h_index[!duplicated(h_index),]
h_index$h_index_ptile <- percent_rank(h_index$H_Index)/0.1
icc_df_all <- merge(icc_df_all,h_index[,c("auth_id","h_index_ptile")],by="auth_id",all.x=T,all.y=F)

n_pubs <- icc_df_all[,c("Total_Publications_In_Scopus","auth_id")]
n_pubs <- n_pubs[!duplicated(n_pubs),]
n_pubs$n_pubs_ptile <- percent_rank(n_pubs$Total_Publications_In_Scopus)/0.1
icc_df_all <- merge(icc_df_all,n_pubs[,c("auth_id","n_pubs_ptile")],by="auth_id",all.x=T,all.y=F)

#### Merge in author country data ####
country <- readRDS("./data/author_countries_raw.rds")
country$country <- as.character(country$country)
# UNSD country codes, downloaded April 23, 2019 from https://unstats.un.org/unsd/methodology/m49/overview/ 
country_codes <- read.csv("./data/UNSD_country_codes.csv")
# change three letter country code to lower case
country_codes$country <- tolower(country_codes$ISO.alpha3.Code)
country_codes <- country_codes[country_codes$country != "",c("country","Region.Name")]
# Merge in
country <- merge(country,country_codes,by="country",all.x=T,all.y=F)
# Do some manually
asian_countries_missed <- c("hkg","twn") # Hong Kong, Taiwan
country$Region.Name[country$country %in% asian_countries_missed] <- "Asia"
country$Asia <- country$Region.Name == "Asia"
# Which authors published in an Asian country at least once in their first year?
authors_asian <- tapply(country$Asia,as.character(country$auth_id),sum,na.rm=T)
authors_asian <- data.frame(auth_id=names(authors_asian),asia=as.integer(authors_asian>=1))
# merge into authors dataset
icc_df_all <- merge(icc_df_all,authors_asian,by="auth_id",all.x=T,all.y=F)

#### Merge in journal topics ####
# Read in journal topics
journal_topics <- read.csv(file="./data/All_Journals_ASJC.csv")
# Keep only medical or multidisciplinary topics (some journals have extra ASJC codes outside our range of interest)
journal_topics <- journal_topics[(journal_topics$AJSC_Codes < 2800 & journal_topics$AJSC_Codes >= 2700) |
                                   journal_topics$AJSC_Codes==1000,]
journal_names <- journal_topics[!duplicated(journal_topics$pub_sourceid),c("pub_sourceid","sourcetitle")]
topic_names <- read.csv(file="./data/ASJC Codes with levels.csv",
                        sep=";") # cloned from github.com/plreyes/Scopus.git
journal_topics <- merge(x=journal_topics,y=topic_names,
                        by.x="AJSC_Codes",by.y="Code",
                        all.x=T)

# Create dataframe of topics by journal
journal_topics_low <- dcast(journal_topics,pub_sourceid ~ Low,fun.aggregate = length, value.var="Low")
journal_topics_low <- merge(journal_topics_low,journal_names,by="pub_sourceid")
journal_topics_low <- journal_topics_low[,c(ncol(journal_topics_low),1:(ncol(journal_topics_low)-1))]

############## Save data #############
saveRDS(icc_df_all,file="./data/processed_data_all.rds")
saveRDS(icc_df,file="./data/processed_data_no_missing.rds")
saveRDS(journal_topics_low,file="./data/journal_topics.rds")

#######################################################################
sink()
#######################################################################

