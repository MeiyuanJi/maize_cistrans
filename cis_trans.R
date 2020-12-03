
suppressMessages(library("ggplot2"))
suppressMessages(library("argparse"))
suppressMessages(library("RColorBrewer"))

parse_arg <- function() {
    parser <- ArgumentParser(description = "classify gene regulatory patterns based on log2foldchange between parents and in their reciprocal F1 hybrid. ")
    parser$add_argument("-par", "--parents", help = "differential expression results for parents. ")
    parser$add_argument("-F1", "--F1", help = "differential expression result for F1 hybrid. ")
    parser$add_argument("-p", "--pvalue", help = "pcutoff for significance. ", default = 0.05)
    parser$add_argument("-O", "--output", help = "prefix of output file (for the merged-results). ")
    parser$parse_args()
}

argv <- parse_arg()
par <- argv$parents
F1 <- argv$F1
pre <- argv$output
p <- as.numeric(argv$pvalue)
### load and read DESeq2 results
par_df <- read.table(par, header = T, sep = "\t", row.names = 1)
F1_df <- read.table(F1, header = T, sep = "\t", row.names = 1)
A <- par_df[, c("log2FoldChange", "lfcSE", "padj")]
names(A) <- c("A_lfc", "A_SE", "A_p")
B <- F1_df[, c("log2FoldChange", "lfcSE", "padj")]
names(B) <- c("B_lfc", "B_SE", "B_p")
AB <- merge(A, B, by = 0, all = T) # 0 specifices the row.names
colnames(AB) <- c("gene", colnames(AB)[2:length(colnames(AB))])
AB2 <- AB[, -1]
rownames(AB2) <- AB[, 1]
AB <- AB2
############################################################
#################### Student's t test ######################
############################################################
### Student's t test ###
### Compare if A-B = 0
### Given variables: log2FoldChange, Standard Error of lfc (lfcSE)
### SE = S/sqrt(n)
# T-statistic = (X1^2 - X2^2 - HypothesizedDifference)/SE(x1-x2)
# When equal variance:
# SE(x1-x2)^2 = ((n1-1)*S1^2 + (n2-1)*S2^2)/(n1+n2-2) * (1/n1 + 1/n2)
# df = n1 + n2 -2
# When unequal variance:
# SE(x1-x2)^2 = S1^2/n1 + S2^2/n2
# df = (S1^2/n1 + S2^2/n2)^2/(S1^2/n1)^2/(n1-1) + (S2^2/n2)^2/(n2-1)
### Given df and t-statistic, we can get the probability of the null hypothesis is true.
t.test2 <- function(lfc1, lfc2, se1, se2, n1, n2, equal.variance = F) {
    if (equal.variance == F) {
        # Aspin-Welch Unequal-Variance t-test
        se <- sqrt((se1^2) + (se2^2))
        df <- (se1^2 + se2^2)^2/((se1^2)^2/(n1-1) + (se2^2)^2/(n2-1))
    } else {
        se <- sqrt(((n1-1)*se1^2*n + (n2-1)*se2^2*n)/((n1+n2-2)*(1/n1 + 1/n2)))
        df <- n1 + n2 - 2
    }
    t <- (lfc1 - lfc2)/se
    res <- c(lfc1-lfc2, se, t, 2*pt(-abs(t), df))
    names(res) <- c("Difference of means", "Std Error", "t", "p-value")
    return(res)
}

criteria <- as.data.frame(rbind(c("A!=0;B!=0;A=B", "Pure Cis"), c("A!=0;B!=0;A!=B", "Cis+Trans"), c("A!=0;B=0;A!=B", "Pure Trans"), c("A=0;B!=0;A!=B", "Compensatory"), c("A=0;B=0;A=B", "Conserved")))
names(criteria) <- c("class", "category")

A.n = 4
B.n = 4
pcutoff = p
##
AB <- AB[complete.cases(AB), ] # drop rows with NA
AB$A_B <- AB$A_lfc - AB$B_lfc
AB$AB_p <- apply(AB, 1, function(x) t.test2(x["A_lfc"], x["B_lfc"], x["A_SE"], x["B_SE"], A.n, B.n)["p-value"])
#AB$AB_padj <- p.adjust(AB$AB_p, method = "BH")
AB$A <- ifelse(AB$A_p < pcutoff, "A!=0", "A=0")
AB$B <- ifelse(AB$B_p < pcutoff, "B!=0", "B=0")
AB$AB <- ifelse(AB$AB_p < pcutoff, "A!=B", "A=B")
AB$class <- paste(AB$A, AB$B, AB$AB, sep = ";")
AB$category <- as.character(criteria$category[match(AB$class, criteria$class)])
AB$category[is.na(AB$category)] <- "Ambiguous"
AB$category[AB$category == "Cis+Trans" & (AB$B_lfc)*(AB$A_B) > 0] <- "Cis+Trans:enhancing"
AB$category[AB$category == "Cis+Trans" & (AB$B_lfc)*(AB$A_B) < 0] <- "Cis+Trans:compensating"
write.table(AB, file = paste0(pre, "_raw.txt"), sep = "\t", quote = F, row.names = T)

AB_clean <- AB[, c("A_lfc", "B_lfc", "AB_p", "class", "category")]
write.table(AB_clean, file = paste0(pre, ".txt"), sep = "\t", quote = F, row.names = T)

