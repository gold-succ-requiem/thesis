# TABLE OF SIDER EFFECTS AND FREQUENCY
# Load packages
library(data.table)
library(dplyr)
library(magrittr)
library(readr)
library(tibble)
library(tidyr)

# Load data
path <- "../data/"
meddra <- fread(paste(path, "meddra_freq.tsv", sep = "")) %>%
    as_tibble()
drugbank <- fread(paste(path, "drugBank_Sider_mapping_5_1_7.csv", sep = "")) %>%
    as_tibble() %>%
    select("DrugBank ID", STITCH_ID)

# Create SIDER table
meddra.sub <- meddra %>% 
    select(stitch.id.stereo, umls.concept.id.meddra, meddra.concept, is.placebo) %>%
    add_column(freq = rowMeans(select(meddra, starts_with("freq.")))) %>%
    subset(meddra.concept == "PT" & is.placebo == "") %>%
    select(-one_of(c('meddra.concept', 'is.placebo'))) %>%
    pivot_wider(names_from = umls.concept.id.meddra, values_from = freq, values_fn = median) %>%
    merge(., drugbank, by.x = "stitch.id.stereo", by.y = "STITCH_ID", all = T) %T>%
    saveRDS(., file = "meddra.sub.rds")
## How to add column of DrugBank IDs
## mapping to drugbank data frame

# Match stereo STITCH IDs to CIDs -- test here
