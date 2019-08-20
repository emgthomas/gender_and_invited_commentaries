# Raw and processed datasets

This directory contains raw, de-identified data and processed datasets used for analysis.

## All_Journals_ASJC.csv
List of All Science Journal Classification (ASJC) codes for each journal included in our dataset.

## ASJC Codes with levels.csv
Names of topics corresponding to each ASJC code, for different levels of the ASJC hierarchy. This file was originally cloned from github.com/plreyes/Scopus.git, and some corrections were made.

## authors_asian.rds
Data frame indicating which authors have Asian country of origin. We do not provide full country of origin information to ensure authors cannot be identified.

## authors_raw.rds
List of authors in each matched set. Author ID variable (auth_id) was randomly generated.

## imputed_data.rds
Contains ten datasets with missing gender data randomly imputed, generated as part of multiple imputation analyses.

## journal_topics.rds
List of topics for each journal.

## processed_data_no_missing.rds
This data frame is in a format ready for input into the clogit function (survival package) to run conditional logistic regression analyses. Authors with missing gender have been excluded.

## processed_data_all.rds
Same dataset as above, but authors with missing gender have not been excluded. For use in multiple imputation analyses.

## publications_raw.rds
List of all intra-citing commentary articles included in our study, and journals in which they were published. Article titles have been removed. Publication ID variable (pub_id) was randomly generated.
