# Load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(magrittr)
library(pbapply)
library(rje)
library(rlist)
library(ROCit)
library(sets)
library(tibble)
library(tidyr)

#lapply(1:15, function(i) {paste("validation", i, "rds", sep = ".")})
# Set file variables
#predFile <- "W.rds"
#predFile <- readRDS("W.1.rds") #%>%
    #extract2(5)
#repoDBFile <- "../data/repoDB_full.csv"
#predDat.single <- readRDS("W.combn.1.rds")
predDat <- readRDS("W.4.rds")
li <- list()

# Prepare predDat
#predDat <- lapply(1:4, function(i) {predFile[[i]] <- predDat.single[[i]]})
# predDat <- inset(predDat.single, 5:15, NULL) %>%
#     list.append(., predFile[5:15]) %>%
#     list.flatten() #%>%
    #View()

# Function for validating against RepoDB
repodb.validate <- function(predDat) {
    # Peforms validation of aggregated matrices against RepoDB.
    
    # Empty list
    valid <- list()
    
    # Variables
    threshold <- 0.01
    #predDat <- predDat[[1]]
    
    # Load prediction data files
    #predDat = as.data.frame(fread(file = predDatFile, header = TRUE, sep = ",", quote = "\""))
    #predDat = fread(predDatFile)
    #predDat[1:2] = t(apply(predDat[1:2],1,sort))
    #predDat = predDat[!duplicated(predDat[1:2]),]
    #predDat = as.numeric(predDat)
    #predDat = fread(predDatFile) %>% 
    #  column_to_rownames(var = "rn") %>%
    #  as.matrix()
    #predDat <- readRDS(predFile) %>%
    #    extract2(5) %>%
    #    View()
    
    # Convert to table and eliminate duplicates
    predDat[upper.tri(predDat, diag = T)] = 0
    predDat <- as.data.frame.table(predDat) %>%
        filter(., Freq != 0)
    
    # Write to list
    valid[["pred.dat"]] <- predDat
    #fwrite(predDat, file="../data/repoDB_PredDat_V7.csv", sep=",", row.names = F, quote = T)
    
    # Load RepoDB for drug IDs
    repoDB <- fread(file = "../data/repoDB_full.csv", header = TRUE, sep = ",", quote = "\"")
    
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
    score = val.R$pred[-which(is.na(val.R$class))]
    class = val.R$class[-which(is.na(val.R$class))]
    
    # ------------ Simple ROC 
    ROCit_obj = rocit(score = score, class = class)
    valid[["auc"]] <- ROCit_obj$AUC
    
    # ------------ ROC + confidence interval
    #score = newPred
    #class = newClass
    rocit_emp <- rocit(score = score, class = class, method = "emp")
    rocit_bin <- rocit(score = score, class = class, method = "bin")
    # --------------------------
    ciROC_emp90 <- ciROC(rocit_emp, level = 0.9)
    set.seed(200)
    ciROC_bin90 <- ciROC(rocit_bin, level = 0.9, nboot = 200)

    # valid[["roc"]] <- plot(ciROC_emp90, col = 1, legend = FALSE) %>%
    #     with(., recordPlot())
    
    # Plot ROC curve with lines, record to list
    valid[["roc"]] <- {
        plot(ciROC_emp90, col = 1, legend = FALSE)
        lines(ciROC_bin90$TPR~ciROC_bin90$FPR, col = 2, lwd = 2)
        lines(ciROC_bin90$LowerTPR~ciROC_bin90$FPR, col = 2, lty = 2)
        lines(ciROC_bin90$UpperTPR~ciROC_bin90$FPR, col = 2, lty = 2)
        legend("bottomright", c("Empirical ROC",
                                "Binormal ROC",
                                "90% CI (Empirical)", 
                                "90% CI (Binormal)"),
            lty = c(1,1,2,2), col = 
                c(1,2,1,2), lwd = c(2,2,1,1))
        } %>%
        with(., recordPlot())
    # lines(ciROC_bin90$TPR~ciROC_bin90$FPR, col = 2, lwd = 2)
    # lines(ciROC_bin90$LowerTPR~ciROC_bin90$FPR, col = 2, lty = 2)
    # lines(ciROC_bin90$UpperTPR~ciROC_bin90$FPR, col = 2, lty = 2)
    # legend("bottomright", c("Empirical ROC",
    #                         "Binormal ROC",
    #                         "90% CI (Empirical)", 
    #                         "90% CI (Binormal)"),
    #        lty = c(1,1,2,2), col = 
    #            c(1,2,1,2), lwd = c(2,2,1,1))
    
    # ----------------- KS plot
    # KS plot shows the cumulative density functions F(c) and G(c) in the positive 
    # and negative populations. If the positive population have higher value, then 
    # negative curve (F(c)) ramps up quickly. The KS statistic is the maximum difference 
    # of F(c) and G(c). (Source: https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html)
    kplot <- ksplot(ROCit_obj)
    valid[["kstat"]] <- kplot[["KS stat"]] %>%
        with(., recordPlot())
    valid[["kstat.cutoff"]] <- kplot[["KS Cutoff"]] %>%
        with(., recordPlot())
    
    # ---------------- Gain table
    # For description: https://cran.r-project.org/web/packages/ROCit/vignettes/my-vignette.html)
    #rocit_emp <- rocit(score = score, class = class, negref = "FP")
    #gtable_custom <- gainstable(rocit_emp, breaks = seq(1,100,15))
    #valid[["gtable"]] <- plot(gtable_custom, type = 1) %>%
    #    with(., recordPlot())
    
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
    valid[["prec"]] <- plot(measure$PREC~measure$REC, type = "l") %>%
        with(., recordPlot())
    
    # Plot ACC
    valid[["fscr"]] <- measure$FSCR
    valid[["acc"]] <- plot(measure$ACC~measure$Cutoff, type = "l") %>%
        with(., recordPlot())
    
    # Return list
    return(valid)
}

# Run function
li <- lapply(c(1:2, 4:length(predDat)), function(i) {
    l <- list()
    l <- repodb.validate(predDat[[i]]) %>%
        list.append()
    return(l)
})
saveRDS(li, "li.4.rds")

# Experiments
# li.x.rds

# v <- c(1:5)
# v %>%
#     magrittr::extract(c(1:2, 4:5)) %>%
#     View()
