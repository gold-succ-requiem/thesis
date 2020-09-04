# Testing similarity matrix and SNF
## Load packages
library(data.table)
library(dplyr)
library(SNFtool)

# Chemical structure
## Load structures (SDF)
struct.sim <- fread("$TMPDIR/chem_similarity.csv") %>%
    as.matrix()
rownames(struct.sim) <- struct.sim[,1]

# Target structure
## Load protein.fasta download as some matrix
seq.sim <- fread("$TMPDIR/protein_similarity.csv") %>%
    as.matrix()

# Similarity network fusion
## Set parameters
K <- 20 # K neighbours
alpha <- 0.5 # hyperparameter
t <- 20 # t iterations

## Construct similarity graphs
W.struct <- affinityMatrix(struct.sim, K, alpha)
W.seq <- affinityMatrix(seq.sim, K, alpha)

## Fuse similarity graphs
W <- list(W.struct, W.seq) %>%
    SNF(., K, t)

## Output as sanity check
fwrite(W, file = "$TMPDIR/fused_similarity.csv")
