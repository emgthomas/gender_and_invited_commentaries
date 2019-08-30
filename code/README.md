# Code for reproducing data processing and statistical analyses

This repository contains code to run all data processing steps and analyses, and to reproduce all figures and tables presented manuscript.

Package dependencies: data.table (1.10.4.3), dplyr (0.7.6), gtools (3.5.0),  reshape2 (1.4.3), gmodels (2.18.1), survival (2.41.3), splines (3.4.2), forestplot (1.7.2), plotly (4.7.1), RColorBrewer (1.1.2), lmtest (0.9.35), metafor (2.0.0), lme4 (1.1.14), mice (3.4.0)

The repository contains the following files:

## master.R
This file will source all other files, thereby reproducing all figures and tables. To run this file, first *edit the working directory in master.R to point to the root of this git repository*. Each .R file can also be run independently.

## data_cleaning.R
Produces datasets needed for analysis using raw data on authors (all case authors plus up to 50 controls per case) and intra-citing commentary articles. Produces a text file (../results/data_cleaning.txt) that reports numbers needed to reproduce Figure 1.

## descriptive_analyses.R
Produces a text file (../results/descriptive_analyses.txt) that reports numbers needed to reproduce Table 1, eTable 1, and other descriptive statistics provided in the Results section of the main manuscript.

## main_analyses.R
Produces a text file (../results/main_analyses.txt) that reports results of main analyses, including all one-stage meta analysis conditional logistic regression models (eTables 2 through 5). Produces Figure 2, eFigures 3 through 7. Runs journal-specific models and produces eTable 6 and data for interactive eFigure 8.

## sensitivity_analyses.R
Produces a text file (../results/sensitivity_analyses.txt) that reports results of sensitivity analyses, including eTables 8 through 9 and eFigure 10.

## secondary_analyses.R
Produces a text file (../results/secondary_analyses.txt) that reports results of secondary analyses, including subgroup analyses by journal topic and model allowing for effect modification by journal Cite Score. Produces Figure 3, Figure 4, and eFigure 9.
