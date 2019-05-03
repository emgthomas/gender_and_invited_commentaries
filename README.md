# Gender differences in authorship of invited commentary articles in medical journals: a matched case-control study

This repository contains code and and data to reproduce all analyses, figures, and tables presented in the above manuscript. Data are de-identified. We provide the matched dataset. Code for performing matching is not provided; matching of cases and controls was performed internally at Elsevier.

All analyses were performed in R version 3.4.2. Sourcing the file code/master.R will re-run all analyses and save all tables and figures in the results directory, but first *edit the working directory in master.R to point to the root of this git repository*.

The repository contains the following directories:

## code
All code for processing data, running analyses, and reproducing tables, plots, and statistics quoted in the text of the manuscript. See the README file in this directory, which provides (1) package dependencies; (2) for each R script in the directory, a list of which figures and tables in the main manuscript and supplementary material can be reproduced with that script.

## data
Raw and processed datasets used in analyses. We provide matched author-level data and data on all intra-citing commentary articles. Names of authors and publications are removed to ensure data are de-identified. However, actual journal names are included.

## results
Figures and text files containing tables and all statistics reported in the main manuscript and supplementary material.

## shiny_app
Code for reproducing eFigure 1, which shows odds ratios plotted against journal Cite Score for all journals.
