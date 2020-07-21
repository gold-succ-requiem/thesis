# obtains data from ClinicalTrials.gov, hopefully including ADMET and ADE data
# input: ClinicalTrials.gov query
# output: CSV of study details

#=========

# load packages
library(dplyr)
library(jsonlite)
library(RCurl)

#=========

# FUNCTION: obtain JSON
## obtain study records: phase 3+, completed
## in JSON format: to circumvent security issues from https; and because CSV is not supported
fetch.json <- function(url) {
    getURL(url) %>%
        fromJSON()
}

#=========

# Pass URL to functions
## convert JSON to dataframe via as.data.frame()
## convert and write to CSV via write.csv()
## NOTE: will write to directory "../data"
u <- "https://clinicaltrials.gov/api/query/full_studies?expr=&recrs=e&recrs=l&phase=2&phase=3&min_rnk=1&max_rnk=100&fmt=JSON"
fetch.json(u) %>%
    as.matrix(as.data.frame()) %>%
    write.csv(.,"../data/clinicaltrials.csv")

#=========
