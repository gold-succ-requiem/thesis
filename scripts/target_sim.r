# Compute similarity matrix -- target structure
# Load packages
library(data.table)
library(protr)

# Load protein.fasta download as some matrix
#seqs <- readFASTA("../data/protein.fasta")
seqs <- readFASTA("./scripts/protein.fasta")

# Use protr to generate similarity matrix -- needs list
seqs.sim <- parSeqSim(seqs, cores = 6)

# Export matrix to .csv
fwrite(seqs.sim, file = "seqs.sim.csv")
