# SIMILARITY NETWORK FUSION
# Load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(SNFtool)
library(tibble)

# Vector of RepoDB unique drug IDs
repodb.ids <- fread("../data/repoDB_full.csv") %>%
    pull(., drug_id) %>%
    unique()

# Vector of CSV file names
files <- c("chem_similarity", "target_similarity", "path_similarity_v2", "GO_Sim_BP_combined")

# SNF variables
K <- 20 # K neighbours
alpha <- 0.5 # empirical hyperparameter
t <- 20 # t iterations

# Preprocessing function
preprocess <- function(x) {
    # Function takes a similarity matrix CSV and preprocesses it for SNF. Specifically sets matrix rownames, while also conforming to drug IDs in RepoDB, and removing NA values.
    
    # Set path variables
    #x <- files[2]
    path <- "../data/"
    csv <- ".csv"
    
    # Load similarity data
    # Filter rows based on matching drug IDs in repodb.ids
    # Convert first column to row names
    # Remove NA
    # Tranpose matrix to make colnames as rownames and repeat above steps
    fread(paste(path, x, sep = "")) %>%
        filter(V1 %in% repodb.ids) %>%
        column_to_rownames(var = "V1") %>%
        filter(Reduce(`+`, lapply(., is.na)) != ncol(.)) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = "V1") %>%
        filter(V1 %in% repodb.ids) %>%
        column_to_rownames(var = "V1") %>%
        filter(Reduce(`+`, lapply(., is.na)) != ncol(.)) %>%
        as.matrix() #%>%
        #View()
        #return()
    ## find way to set var as first index
}

# SNF function over list of matrices
ifelse.SNF <- function(i) {
    ifelse(length(W.all.combn[[i]]) == 1, W.all.combn[[i]],SNF(W.all.combn[[i]], K, t))
}

# Build list of preprocessed CSV files
# Convert to dist mat
# Cut off at threshold
W.all <- lapply(files, preprocess) %>%
    lapply(., `*`, -1) %>%
    lapply(., `+`, 1)

# Collect row names
common.rows <- Reduce(intersect, lapply(W.all, row.names))

# Subset dist mats by common drug IDs
# Convert to affinity matrix
W.all.repodb <- lapply(W.all, function(x) { x[row.names(x) %in% common.rows,] }) %>%
    lapply(., function(x) {x[ ,which(colnames(x) %in% common.rows)]}) %>%
    lapply(., affinityMatrix, K, alpha)

# List of all
W.all.combn <- do.call(c, lapply(1:length(W.all.repodb), function(x) {combn(W.all.repodb, x, simplify = F)}))

# Fuse affinity matrices
#W <- lapply(1:length(W.all.combn), ifelse.SNF)
W <- lapply(1:length(W.all.combn), function(i) {SNF(W.all.combn[[i]], K, t)})

# Output
saveRDS(W, "W.rds")
saveRDS(W.all.combn, "W.all.combn.rds")
