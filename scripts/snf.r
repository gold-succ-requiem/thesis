# Testing similarity matrix and SNF
## Load packages
library(data.table)
library(dplyr)
library(grDevices)
library(SNFtool)

# Array of RepoDB unique drug IDs
repodb.ids <- fread("../data/repoDB_full.csv") %>%
    select(drug_id) %>%
    distinct(., drug_id, .keep_all = T)

# Chemical structure
## Load structure sims as matrix
struct.sim <- fread("../data/chem_similarity.csv") %>%
    as.matrix()

## Convert first column to row names
rownames(struct.sim) <- struct.sim[,1]
struct.sim <- struct.sim[,-1]

## Subset based on RepoDB
struct.sim.repodb <- struct.sim[ ,which(colnames(struct.sim) %in% repodb.ids$drug_id)]
struct.sim.repodb <- subset(struct.sim.repodb, rownames(struct.sim.repodb) %in% repodb.ids$drug_id)

## Make elements numeric -- so `affinityMatrix()` won't throw a hissy fit
struct.sim.repodb <- `dimnames<-`(`dim<-`(as.numeric(struct.sim.repodb), dim(struct.sim.repodb)), dimnames(struct.sim.repodb))

# Target structure
## Load target similarities as matrix
seq.sim <- fread("../data/protein_similarity.csv") %>%
    as.matrix()

## Divide matrix by 100 to make everything =< 1
seq.sim <- seq.sim * (1/100)

## Subset based on RepoDB
seq.sim.repodb <- seq.sim[ ,which(colnames(seq.sim) %in% repodb.ids$drug_id)]
seq.sim.repodb <- subset(seq.sim.repodb, rownames(seq.sim.repodb) %in% repodb.ids$drug_id)

# Pathway similarity
path.sim <- fread("../data/path_similarity_v2.csv") %>%
    as.matrix()

## Convert first column to row names
rownames(path.sim) <- path.sim[,1]
path.sim <- path.sim[,-1]

## Subset based on RepoDB
path.sim.repodb <- path.sim[ ,which(colnames(path.sim) %in% repodb.ids$drug_id)]
path.sim.repodb <- subset(path.sim.repodb, rownames(path.sim.repodb) %in% repodb.ids$drug_id)

## Make elements numeric -- so `affinityMatrix()` won't throw a hissy fit
#path.sim.100 <- apply(path.sim.100, 2, as.numeric)

# GO BP similarity
bp.sim <- fread("../data/GO_Sim_BP_combined.csv") %>%
    as.matrix()

## Convert first column to row names
rownames(bp.sim) <- bp.sim[,1]
bp.sim <- bp.sim[,-1]

## Subset based on RepoDB
bp.sim.repodb <- bp.sim[ ,which(colnames(bp.sim) %in% repodb.ids$drug_id)]
bp.sim.repodb <- subset(bp.sim.repodb, rownames(bp.sim.repodb) %in% repodb.ids$drug_id)

## Make elements numeric -- so `affinityMatrix()` won't throw a hissy fit
bp.sim.repodb <- `dimnames<-`(`dim<-`(as.numeric(bp.sim.repodb), dim(bp.sim.repodb)), dimnames(bp.sim.repodb))

# GO CC similarity
cc.sim <- fread("../data/GO_Sim_CC_combined.csv") %>%
    as.matrix()

## Convert first column to row names
rownames(cc.sim) <- cc.sim[,1]
cc.sim <- cc.sim[,-1]

## Omit rows all containing NA AND cols all containing NA, not actually needed
cc.sim <- cc.sim[rowSums(is.na(cc.sim)) != ncol(cc.sim), ]
cc.sim <- cc.sim[ ,colSums(is.na(cc.sim)) < nrow(cc.sim)]

## Subset first 100 rows and cols
cc.sim.100 <- cc.sim[1:100, 1:100]

## Make elements numeric -- so `affinityMatrix()` won't throw a hissy fit
cc.sim.100 <- apply(cc.sim.100, 2, as.numeric)

# GO MF similarity
mf.sim <- fread("../data/GO_Sim_MF_combined.csv") %>%
    as.matrix()

## Convert first column to row names
rownames(mf.sim) <- mf.sim[,1]
mf.sim <- mf.sim[,-1]

## Omit rows all containing NA AND cols all containing NA, not actually needed
mf.sim <- mf.sim[rowSums(is.na(mf.sim)) != ncol(mf.sim), ]
mf.sim <- mf.sim[ ,colSums(is.na(mf.sim)) < nrow(mf.sim)]

## Subset first 100 rows and cols
mf.sim.100 <- mf.sim[1:100, 1:100]

## Make elements numeric -- so `affinityMatrix()` won't throw a hissy fit
mf.sim.100 <- apply(mf.sim.100, 2, as.numeric)

# Similarity network fusion
## Set parameters
K <- 20 # K neighbours
alpha <- 0.5 # hyperparameter
t <- 20 # t iterations

## Construct similarity graphs
W.struct <- affinityMatrix(struct.sim.repodb, K, alpha)
W.bp <- affinityMatrix(bp.sim.repodb, K, alpha)

## Fuse similarity graphs
W <- list(W.struct, W.seq, W.path, W.bp, W.cc, W.mf) %>%
    SNF(., K, t)

## Output as sanity check
fwrite(W, file = "../data/fused_similarity.csv")

# Define parameters for clustering functions
C <- 2
nn <- 5

# Labels from `SNFtool`
labels.snftool <- spectralClustering(W, C)

## Display clusters
png("../data/snf_clusters.png")
displayClusters(W, labels.snftool)
dev.off()

# Spectral clustering and visualise eigenvalues using `kknn`
labels.kknn <- specClust(W, centers = C, nn = nn)

## Visualise eigenvalues
labels.kknn$centers %>%
    plot()

# Spectral clustering, PCA, and visualisation using `Spectrum`
#labels.spectrum <- Spectrum(W, showpca = T, fontsize = 8, dotsize = 2)
labels.spectrum <- Spectrum(list(W.struct, W.seq, W.path, W.bp, W.cc, W.mf), showpca = T, fontsize = 8, dotsize = 2)
pca(labels.spectrum$similarity_matrix, labels = labels.spectrum$assignments)

# Label comparison for quality assurance
labels.compare <- cbind(labels.snftool, labels.kknn$cluster)
labels.acc <- sum(labels.compare[,1] == labels.compare[,2]) / dim(labels.compare)[1]
