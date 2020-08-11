# Load packages
library(dplyr)
library(rclinicaltrials, lib.loc = "/home/z5232927/packages")

# Sanity check
# clinicaltrials_count(query = c("phase=2", "phase=3", "type=Intr", "rslt=With"))

# Downloads list of lists of data
dl <- clinicaltrials_download(query = c("phase=2", "phase=3", "type=Intr", "rslt=With"), count = 10000, include_results = T)

# Merges arm and intervention dataframes into csv
merge(dl$study_information$arms, dl$study_information$interventions, by = "nct_id") %>%
    write.csv(., file = "clinicaltrials.csv")

# To-do list
# - filter rows by including "intervention_type == 'Drug'" only
# - edit two columns from study info df
# - run through katana
