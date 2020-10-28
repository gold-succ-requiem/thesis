# Load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(pbapply)
library(rje)
library(ROCit)
library(sets)
library(tibble)
library(tidyr)

#lapply(1:15, function(i) {paste("validation", i, "rds", sep = ".")})
# Set file variables
predFile <- "W.rds"
#predDat <- readRDS("W.rds") %>%
    extract2(5)
repoDBFile <- "../data/repoDB_full.csv"
li <- list()

# validate.fun for validating input
validate.fun <- function(predFile, labelFile) {
    # Peforms validation of aggregated matrices against RepoDB.
    
    # Empty list
    valid <- list()
    
    # Variables
    threshold <- 0.01
    
    # Load prediction data files
    #predDat = as.data.frame(fread(file = predDatFile, header = TRUE, sep = ",", quote = "\""))
    #predDat = fread(predDatFile)
    #predDat[1:2] = t(apply(predDat[1:2],1,sort))
    #predDat = predDat[!duplicated(predDat[1:2]),]
    #predDat = as.numeric(predDat)
    #predDat = fread(predDatFile) %>% 
    #  column_to_rownames(var = "rn") %>%
    #  as.matrix()
    predDat <- readRDS(predFile) %>%
        extract2(5)
    
    # Convert to table and eliminate duplicates
    predDat[upper.tri(predDat, diag = T)] = 0
    predDat <- as.data.frame.table(predDat) %>%
        filter(., Freq != 0)
    
    # Write to list
    valid[["pred.dat"]] <- predDat
    #fwrite(predDat, file="../data/repoDB_PredDat_V7.csv", sep=",", row.names = F, quote = T)
    
    # Load RepoDB for drug IDs
    repoDB <- fread(file = labelFile, header = TRUE, sep = ",", quote = "\"")
    
    # Filter predicted Drug-pairs for which drug data available in repoDB
    # I've commented this out as I've already done this
    #predDat = predDat[which((predDat$ID1 %in% repoDB$drug_id) & (predDat$ID2 %in% repoDB$drug_id)),]
    # predDat = predDat[which(predDat$p_value < 0.01),]
    
    # # Get TP and TN info from repoDB
    # generate Predicted labels
    # P_vector = ifelse(predDat$adjP_value < threshold, 1, 0)
    # OR, use "raw" p-values
    
    # P_vector = 1 - predDat$adjP_value
    #P_vector = predDat$rowMeans # this isn't appropriate, because, its just a combined score, doesn't indicate probability
    P_vector <- predDat$Freq
    
    getL.fun <- function(x){
        r1 = subset(repoDB, drug_id %in% x[[1]]) #repoDB[which(repoDB$drug_id %in% x[1]),]
        r2 = subset(repoDB, drug_id %in% x[[2]]) #repoDB[which(repoDB$drug_id %in% x[2]),]
        r1.OR.r2.aprv = (length(r1$stat[which(r1$stat == "Approved")]) > 0) | (length(r2$stat[which(r2$stat == "Approved")]) > 0)
        # r2.aprv = r2$stat[which(r2$stat == "Approved")]
        if(length(intersect(r1$ind_id, r2$ind_id))>0 & r1.OR.r2.aprv) # drug1 and drug2 have similar indication and drug2 is generally an approved drug
            return(1)
        else if (length(intersect(r1$ind_id, r2$ind_id))==0 & r1.OR.r2.aprv) # drug2 is approved, but drug1 and drug2 DON'T have similar indication
            return (NA)
        else
            return(0)
        # return(ifelse(length(intersect(r1$ind_id, r2$ind_id))>0 & length(r2.aprv) > 0, 1, 0))
    }
    
    getL <- pbapply(predDat, 1, getL.fun)
    
    # Create list of predicted values and true labels
    #return(list(pred=P_vector,class=getL))
    val.R <- list(pred = P_vector, class = getL)
    
    # for full Chemical similarity
    # val.R <- repoDB_validation(predDatFile = "../data/",repoDBFile = "../data/repoDB_full.csv", threshold = 0.01)
    # ROCit_obj = rocit(score=val.R$pred, class=val.R$class)
    # plot(ROCit_obj)
    
    
    # for full data
    #val.R <- repoDB_validation(predDatFile = "../data/new_net_info_V7_pval.csv",repoDBFile = "../data/repoDB_full.csv", threshold = 0.01)
    newPred = val.R$pred[-which(is.na(val.R$class))]
    newClass = val.R$class[-which(is.na(val.R$class))]
    
    # ------------ Simple ROC 
    ROCit_obj = rocit(score=newPred, class=newClass)
    valid[["auc"]] <- ROCit_obj$AUC
    
    # ------------ ROC + confidence interval
    score = newPred
    class = newClass
    rocit_emp <- rocit(score = score, class = class, method = "emp")
    rocit_bin <- rocit(score = score, class = class, method = "bin")
    # --------------------------
    ciROC_emp90 <- ciROC(rocit_emp, level = 0.9)
    set.seed(200)
    ciROC_bin90 <- ciROC(rocit_bin, level = 0.9, nboot = 200)
    
    #valid[["roc"]] <- plot(ciROC_emp90, col = 1, legend = FALSE)
    #lines(ciROC_bin90$TPR~ciROC_bin90$FPR, col = 2, lwd = 2)
    #lines(ciROC_bin90$LowerTPR~ciROC_bin90$FPR, col = 2, lty = 2)
    #lines(ciROC_bin90$UpperTPR~ciROC_bin90$FPR, col = 2, lty = 2)
    #legend("bottomright", c("Empirical ROC",
    #                        "Binormal ROC",
    #                        "90% CI (Empirical)", 
    #                        "90% CI (Binormal)"),
    #       lty = c(1,1,2,2), col = 
    #           c(1,2,1,2), lwd = c(2,2,1,1))
    
    # ----------------- KS plot
    # KS plot shows the cumulative density functions F(c) and G(c) in the positive 
    # and negative populations. If the positive population have higher value, then 
    # negative curve (F(c)) ramps up quickly. The KS statistic is the maximum difference 
    # of F(c) and G(c). (Source: https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html)
    #kplot <- ksplot(ROCit_obj)
    #valid[["kstat"]] <- kplot[["KS stat"]]
    #valid[["kstat.cutoff"]] <- kplot[["KS Cutoff"]]
    
    # ---------------- Gain table
    # For description: https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html)
    rocit_emp <- rocit(score = score, class = class, negref = "FP")
    gtable_custom <- gainstable(rocit_emp, breaks = seq(1,100,15))
    #valid[["gtable"]] <- plot(gtable_custom, type = 1)
    
    # ----------------- Precision-vs-Recall curve
    # ACC: Overall accuracy of classification.
    # MIS: Misclassification rate.
    # SENS: Sensitivity.
    # SPEC: Specificity.
    # PREC: Precision.
    # REC: Recall. Same as sensitivity.
    # PPV: Positive predictive value.
    # NPV: Positive predictive value.
    # TPR: True positive rate.
    # FPR: False positive rate.
    # TNR: True negative rate.
    # FNR: False negative rate.
    # pDLR: Positive diagnostic likelihood ratio.
    # nDLR: Negative diagnostic likelihood ratio.
    # FSCR: F-score,.
    
    measure <- measureit(score = score, class = class, measure = c("PREC", "REC", "FSCR", "ACC"))
    valid[["measure.names"]] <- names(measure)
    
    # Plot PREC
    valid[["prec"]] <- plot(measure$PREC~measure$REC, type = "l")
    
    # Plot ACC
    valid[["fscr"]] <- measure$FSCR
    valid[["acc"]] <- plot(measure$ACC~measure$Cutoff, type = "l")
    
    # Check to see if entire function works
    valid[["sanity"]] <- c(1:3)
    
    # Return list
    return(valid)
}

# plot.fun
plot.fun <- function(x) {
    # Plots objects contained in list generated by `validate.fun`.
    
    data <- "../data/"
    ext <- ".png"
    
    png(paste(data, x, ext, sep = "")) %>%
        with(plot(x, type = "l"))
    dev.off()
}

# Run function
li <- validate.fun(predFile, repoDBFile)

#li <- lapply(1:length(predDat), function(i) {
#    l <- list()
#    l <- validate.fun(predDat[[i]], repoDBFile)
#})

# TEST ZONE
# Define function
foo.fun <- function(m, n) {
    # Takes matrix and scalar inputs, and saves outcomes of matrix operations in a list.
    # m is matrix
    # n is scalar
    
    # Empty list
    l <- list()
    
    # Adds 1 to all matrix elements
    l[["add"]] <- m + 1
    
    # Multiplies matrix elements by n
    l[["times"]] <- m * n
    
    # Transposes matrix
    l[["trans"]] <- t(m)
    
    # Return output
    return(l)
}

# Define variables
n <- 2
mt <- matrix(data = 1:9, nrow = 3, ncol = 3)

# Run function
li2 <- foo.fun(mt, n)


