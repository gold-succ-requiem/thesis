# Assessing significance of clustering coefficient

# Compute p-value
# Takes vector of values
# Assumes normal distribution
p.val <- function(x) {
    a <- 5
    n <- length(x)
    s <- sd(x)
    xbar <- mean(x)
    z <- (xbar - a) / (s/sqrt(n))
    pval <- 2*pnorm(-abs(z))
    
    hist(x)
    return(c(z, pval))
}

# Execute p-value
# Threshold is 0.005
readRDS("net.des.3.8.rds") %>%
    extract2(., "cc.local") %>%
    p.val()
