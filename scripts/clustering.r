# CLUSTERING
# Load packages
library(data.table)
library(dplyr)
library(igraph)
library(magrittr)
library(plyr)

# Generate edgelist igraph object
# net <- readRDS("W.3.rds") %>%
#     extract2(., 1) %>%
#     inset(upper.tri(., diag = T), 0) %>%
#     as.data.frame() %>%
#     rownames_to_column(., var = "from") %>%
#     as_tibble() %>%
#     pivot_longer(., !from, names_to = "to", values_to = "weight") %>%
#     subset(., weight > 0) %>%
#     graph_from_data_frame(., directed = F) %>%
#     View()
net <- readRDS("W.1.rds") %>%
    extract2(., 1) %>%
    graph_from_adjacency_matrix(., mode = "undirected", weighted = T, diag = F)

# Mode
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

# Round
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# Function for network descriptives
get.net.des <- function(igraph.obj) {
    # Obtains network descriptives and stores in a neat list for later analysis
    # Requires igraph package
    
    # Test variables -- comment out when sourcing
    igraph.obj <- net
    
    # Generate list
    l <- list()
    
    # Density
    l[["density"]] <- edge_density(igraph.obj, loops = F)
    
    # Transitivity -- aka clustering coefficient
    # Weighted as per Barrat et al. 2004
    cc <- transitivity(igraph.obj, type = "barrat")
    l[["cc.local"]] <- cc
    l[["cc.global"]] <- mean(cc)
    #l[["cc.local"]] <- transitivity(igraph.obj, type = "local")
    
    # Node degree and distribution
    deg <- degree(igraph.obj, mode = "all")
    l[["deg"]] <- deg
    
    ## Histogram
    l[["deg.hist"]] <- hist(deg, breaks=1:vcount(igraph.obj)-1, main = "Degree distribution") %>%
        with(., recordPlot())
    
    ## Cum freq
    deg.dist <- degree_distribution(igraph.obj, cumulative = T, mode = "all")
    l[["deg.dist"]] <- plot(x = 0:max(deg), y = 1 - deg.dist, pch = 19, cex = 1.2, col = "orange", main = "Degree distribution", xlab = "Degree", ylab = "Cumulative Frequency") %>%
        with(., recordPlot())
    
    # Community detection -- edge betweenness
    # Runtime is too long ...?
    #ceb <- cluster_edge_betweenness(igraph.obj, weights = E(igraph.obj)$weight, directed = F, edge.betweenness = F)
    
    # Community detection -- leading eigenvector
    cle <- cluster_leading_eigen(igraph.obj, weights = E(igraph.obj)$weight)
    
    # Number of communities
    l[["cle.no"]] <- length(cle)
    
    # Community membership + plot
    mem <- membership(cle)
    l[["cle.mem"]] <- mem
    
    # Assign colours for visualisation
    V(igraph.obj)$color <- mem + 1
    
    # Cluster membership frequencies
    # Custom functions implemented for nicer visualisation
    # Calling list() conforms barplot object to len = 1, to be passed to recordPlot()
    l[["mem.plot"]] <- barplot(table(mem), col = 1:length(mem), xlab = "Clusters", ylab = "Drugs", ylim = c(0, roundUpNice(sum(mem == Mode(mem))))) %>%
        list() %>%
        with(., recordPlot())
    
    # Modularity
    l[["mod"]] <- modularity(igraph.obj, mem)
    
    # Community plot
    l[["com.plot"]] <- plot(cle, igraph.obj, vertex.label = NA) %>%
        with(., recordPlot())
    
    # Community plot, separate clusters
    l[["com.sep"]] <- delete_edges(igraph.obj, E(igraph.obj)[crossing(cle, igraph.obj)]) %>%
        plot(., vertex.label = NA) %>%
        with(., recordPlot())
    
    # Community dendogram
    l[["com.dend"]] <- dendPlot(cle, mode = "hclust") %>%
        with(., recordPlot())
    
    # Return list
    return(l)
}

# Obtain network descriptives for network
net.des <- list()
net.des <- get.net.des(net) %T>%
    saveRDS(., "net.des.1.1.rds")
