# ---------------------------------------------------------- #
# ----------------------- Master file ---------------------- #
# ---------------------------------------------------------- #

# Set the working directory to the ~~root directory of the git repository~~
setwd("/Users/emt380/Documents/PhD_Papers/Gender_bias/R_code/gender_and_invited_commentaries/github")

# Run files in appropriate order
source("./code/functions.R")
source("./code/data_cleaning.R",print.eval=T) 
source("./code/descriptive_analyses.R",print.eval=T)
source("./code/main_analyses.R",print.eval=T)
source("./code/sensitivity_analyses.R",print.eval=T)
source("./code/secondary_analyses.R",print.eval=T)
