# Analysis of Pre-chemotherapy TB Studies

This directory contains the code and data to run analyses and produce results for TITLE by Rodriguez, CA and Leavitt SV, et al. This paper describes a meta-analysis of disease prognosis (including survival and self-cure rates) of pre-chemotherapy tuberculosis studies stratified by disease severity.

## Data

### pre_chemo_data.xlsx

This is an Excel spreadsheet with the cleaned, extracted life-table data for each of the studies. The tabs are labeled with the numeric ID given to each study.

### study_id.csv

This table details the concordance between the numeric study IDs and the papers they refer to (first author and year).

### mortality_data.csv

This is the individual-level mortality data created from pre_chemo_data.xlsx by data_prep.R.

### cure_data.csv

This is the individual-level self-cure data created from pre_chemo_data.xlsx by data_prep.R.


***

## Programs

### utils.R

This script contains many functions called by the other programs to format the data for analysis and results for presentation.

### data_prep.R

This script takes the Excel spreadsheet with the extracted life-table data (pre_chemo_data.xlsx) and formats it into the individual-level data for both mortality (mortality_data.csv) and self-cure (cure_data.csv) which will be used in analysis.

### mortality_analysis.R

This script runs all of the Bayesian mortality survival models: all-cause and TB-specific mortality both with and without a fixed effect for disease severity. It also runs the study type stratified (santorium/hospital verses non-sanatorium) sensitivity analysis.

### mortality_results.R

This script takes the results of the models run in mortality_analysis.R and formats the output to create the figures and tables included in the main text and supplement of the manuscript.

### cure_analysis.R

This script runs the Bayesian logistic regression model for self-cure and extracts the data from the results for the table used in the manuscript.


