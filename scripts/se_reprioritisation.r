# SIDE EFFECT REPRIORITISATION
# Load packages
library(data.table)
library(dplyr)

# Load side effects
se <- readRDS("meddra.sub.rds") %>%
    mutate(., median = apply(select(., starts_with("C")), 1, median, na.rm = T)) %>%
    select(., `DrugBank ID`, median)

# Original table
# Load fused network -- target + pathway
# Pivot long
# Omit triangular
# Order by similarity ascending
fused <- readRDS("W.1.rds")[[1]] %>%
    inset(upper.tri(., diag = T), 0) %>%
    as.data.frame() %>%
    rownames_to_column(., var = "DrugBank_ID1") %>%
    as_tibble() %>%
    pivot_longer(., !DrugBank_ID1, names_to = "DrugBank_ID2", values_to = "score") %>%
    subset(., score > 0) %T>%
    saveRDS(., "fused.1.1.rds")

# Reprioritised table
# Based on fused table
# New column named median
# If DrugBank_ID2 = `DrugBank ID`, then median = median
# Scores as multiplier
# Reorder
fused.penal <- fused %>%
    full_join(., se, by = c("DrugBank_ID2" = "DrugBank ID"), keep = T) %>%
    select(-`DrugBank ID`) %>%
    drop_na() %>%
    rowwise() %>%
    mutate(score_penal = multiply_by(score, median)) %T>%
    saveRDS(., "fused.repr.1.1.rds")
