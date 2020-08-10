# yo angelo!
library(dplyr)
library(rclinicaltrials)

# Sanity check
# clinicaltrials_count(query = c("phase=2", "phase=3", "type=Intr", "rslt=With"))

# Downloads list of lists of data
dl <- clinicaltrials_download(query = c("phase=2", "phase=3", "type=Intr", "rslt=With"), count = NULL, include_results = T)

# Merges arm and intervention dataframes into csv
dl$study_information %>% 
    c('arms', 'interventions') %>%
    merge(.,., by = "nct_id") %>%
    write.csv(., file = "clinicaltrials.csv")

# To-do list
# - filter rows by including "intervention_type == 'Drug'" only
# - edit two columns from study info df
# - run through katana
