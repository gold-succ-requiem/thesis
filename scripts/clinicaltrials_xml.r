# obtains data from ClinicalTrials.gov, XML -> CSV, test
# input: ClinicalTrials.gov query
# output: CSV of study details

# load packages
library(dplyr)
library(tibble)
library(tidyselect)
library(XML)
library(xml2)

# Function to generate table of trial's reported adverse event counts
count.fun <- function(x){
    # Navigate to <serious_events/category> level
    records <- x %>%
        read_xml() %>%
        xml_find_all(., "////serious_events//category")
    
    # Obtain all group IDs at this level
    group.id <- records %>% 
        xml_find_all("//counts") %>% 
        xml_attr(., "group_id") %>% 
        as_tibble() %>% 
        rename("Group ID" = value)
    
    # Obtain all affected counts at this level
    affected <- records %>% 
        xml_find_all("//counts") %>% 
        xml_attr(., "subjects_affected") %>% 
        as_tibble() %>% 
        rename("Subjects affected" = value)
    
    # Obtain all at-risk counts at this level
    at.risk <- records %>% 
        xml_find_all("//counts") %>% 
        xml_attr(., "subjects_at_risk") %>% 
        as_tibble() %>% 
        rename("Subjects at risk" = value)
    
    # Obtain all adverse event subtitles at this level
    subtitle <- records %>% 
        xml_find_all("//sub_title") %>% 
        xml_text() %>% 
        as_tibble() %>% 
        rename("Subtitle" = value)
    
    # Merge tibbles
    #tab <- cbind(group.id, affected) %>%
    #    cbind(., at.risk)
    tab <- cbind(subtitle, group.id) %>%
        cbind(., affected) %>%
        cbind(., at.risk)
}

# Prepare vector of file names to load
## Test set - vector of two handpicked trials 
#pg <- c("../data/clinicaltrials2/NCT00000125.xml", "../data/clinicaltrials2/NCT00088699.xml")

## Real deal - vector of all trials in directory
trials <- list.files(path = "../data/clinicaltrials2", pattern = "*.xml", full.names = T, recursive = F)

# Apply function to generate counts table
## Create list of tables to append
tab.list <- list()

## Iteratively generate sub-tables
for (i in trials) {
    tab.out <- count.fun(i)
    tab.list[[i]] <- tab.out
}

## Bind sub-tables to final table
fin.tab <- do.call(rbind, tab.list)
