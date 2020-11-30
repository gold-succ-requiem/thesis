# SIMILARITY NETWORK FUSION
# Load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(rlist)
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
sigma <- 0.5 # local model variance
t <- 20 # t iterations

# Preprocessing function
preprocess <- function(x) {
    # Function takes a similarity matrix CSV and preprocesses it for SNF. Specifically sets matrix rownames, while also conforming to drug IDs in RepoDB, and removing NA values.
    # 1 = structure
    # 2 = sequence
    # 3 = path
    # 4 = BP
    
    # Set path variables
    #x <- files[2]
    path <- "../data/"
    csv <- ".csv"
    
    # Load similarity data
    # Filter rows based on matching drug IDs in repodb.ids
    # Convert first column to row names
    # Remove NA
    # Tranpose matrix to make colnames as rownames and repeat above steps
    fread(paste(path, x, csv, sep = "")) %>%
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

# Function for filtering by quantile
filter.quart <- function(sim.mat, quart) {
    # Adds sparsity to matrix based on threshold quartile
    # sim.mat = a similarity matrix
    # Quantile options: 1 = 0th; 2 = 25th; 3 = 50th; 4 = 75th; 5 = 100th
    
    # Test
    #sim.mat <- W.all.repodb[[1]]
    #quart <- 4
    
    # Collect quartiles -- test if input is integer in [1,5] first
    if (quart %in% seq(1, 5) && quart %% 1 == 0) {
        mt.quart <- quantile(sim.mat)[quart]
    } else {
        stop("argument 'quart' must be integer between 1 and 5 inclusive")
    }
    
    # Filter matrix
    sim.mat[sim.mat < mt.quart] <- 0
    
    # Return matrix as output
    return(sim.mat)
}

# Build list of preprocessed CSV files
W.all <- lapply(files, preprocess)

# Collect row names
common.rows <- Reduce(intersect, lapply(W.all, row.names))

# Subset dist mats by common drug IDs
W.all.repodb <- lapply(W.all, function(x) { x[row.names(x) %in% common.rows,] }) %>%
    lapply(., function(x) {x[ ,which(colnames(x) %in% common.rows)]})

# List of all
W.all.combn <- do.call(c, lapply(1:length(W.all.repodb), function(x) {combn(W.all.repodb, x, simplify = F)}))

# Fuse affinity matrices, filter by quartile threshold
W <- lapply(1:length(W.all.combn), function(i) {SNF(W.all.combn[[i]], K, t)}) %>%
    inset(., 1:4, W.all.combn[1:4]) %>%
    list.flatten() %>%
    rapply(., filter.quart, classes = "ANY", how = "list", quart = 1) %T>%
    saveRDS(., "W.1.rds")

# Output -- experiment
# "W.combn.x.rds" <- list of all
# "W.x.rds" <- fused
