# Heatmap and graph visualisation of SNF
# Load packages
library(data.table)
library(dplyr)
library(gplots)
library(qgraph)

# Load SNF to table and normalise
snf <- fread("../data/new_net_info_V7_pval.csv") %>%
    column_to_rownames(var = "rn") %>%
    as.matrix()
snf.triang <- snf[upper.tri(snf, diag = T)] <- NA
#snf <- data.table(
#    row = rep(seq_len(nrow(snf)), ncol(snf)), 
#    col = rep(seq_len(ncol(snf)), each = nrow(snf)), 
#    value = c(snf)
#    ) %>%
#    subset(., !(value == 0))

# Load RepoDB based on common rows
repodb <- fread("../data/repoDB_drugPair_Transformation.csv") %>%
    `[`(, c(2, 1, 3)) %>%
    left_join(., snf.melt, by = c("DrugBank.ID1" = "Var2", "DrugBank.ID2" = "Var1"))

snf.melt.adj <- snf.melt %>%
    filter(Var2 %in% repodb$DrugBank.ID1) %>%
    filter(Var1 %in% repodb$DrugBank.ID2)

## Visualisation -- drug structure, original
#struct.sim.melt <- melt(struct.sim.repodb)
#ggplot(struct.sim.melt, aes(x = Var1, y = Var2, fill = value)) + 
#    geom_tile() + 
#    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
#    labs(x = "Drug ID", y = "Drug ID") +
#    guides(fill = guide_legend(title = "Similarity")) + 
#    ggsave("../data/struct_heatmap.png")

## Visualisation -- drug structure, affinity matrix
#struct.aff.melt <- melt(snf.all[[1]])
#ggplot(snf.melt, aes(x = Var1, y = Var2, fill = value)) + 
#    geom_tile() + 
#    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
#    labs(x = "Drug ID", y = "Drug ID") +
#    guides(fill = guide_legend(title = "Similarity")) + 
#    ggsave("../data/struct_aff_heatmap.png")

# Visualisation -- Fused network heatmap
png("../data/snf_heatmap.png")
heatmap.2(
    snf, 
    dendrogram = "none",
    symm = T,
    #na.color = par("bg"),
    scale="none", 
    Rowv=NA, 
    #Colv=NA, 
    #margins(5,5),
    #cexRow=0.5,cexCol=1.0,
    key=TRUE,keysize=1.5, 
    trace="none",
    labRow = NULL,
    labCol = NULL
    )
dev.off()

# Visualisation -- Fused network graph
png("../data/snf_graph.png")
qgraph(snf)
dev.off()
