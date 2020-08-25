# obtains data from ClinicalTrials.gov, XML -> CSV, test
# input: ClinicalTrials.gov query
# output: CSV of study details

# load packages
library(dplyr)
library(tibble)
library(tidyselect)
library(XML)
library(xml2)

# Function to obtain trial's reported events and descriptions
event.fun <- function(x){
    # Navigate to <category> node
    records <- x %>%
        read_xml() %>%
        xml_find_all(., "////serious_events//category")
    
    # Obtain adverse event subtitles
    subtitle <- records %>% 
        xml_find_all("//sub_title") %>% 
        xml_text() %>% 
        as_tibble() %>% 
        rename("Subtitle" = value)
    
    # Obtain adverse event descriptions
    #desc <- records %>%
    #    xml_find_first("/event/description") %>%
    #    xml_text() %>%
    #    as_tibble()
}

# Prepare vector of file names to load
## Test set - vector of two handpicked trials 
#pg <- c("../data/clinicaltrials2/NCT00000125.xml", "../data/clinicaltrials2/NCT00088699.xml")

## Real deal - vector of all trials in directory
trials <- list.files(path = "../data/clinicaltrials2", pattern = "*.xml", full.names = T, recursive = F)

# Apply function to generate table of side effects only -- to be added to separate script?
## Iteratively fill table
### Create list of tables of append
tab2.list <- list()

### Iteratively generate sub-tables
fin.tab2 <- for (i in trials) {
    tab2.out <- event.fun(i)
    tab2.list[[i]] <- tab2.out
}

### Bind sub-tables to final table
fin.tab2 <- do.call(rbind, tab2.list)