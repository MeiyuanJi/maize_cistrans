suppressMessages(library("DESeq2"))
suppressMessages(library("argparse"))

parse_arg <- function() {
    # accept arguments from terminal
    parser <- ArgumentParser(description = "Usage: Rscript F1vspar_inheritance.R -cnt <count> -sample <sample> -par <par> -F1 <F1> -condition <condition> -O <output>")
    parser$add_argument("-cnt", "--cnt", help = "read count data. ")
    parser$add_argument("-sample", "--sample", help = "sample table matched to count data. ")
    parser$add_argument("-par", "--parents", nargs = "+", help = "parents genotype (in P1vsP2). ")
    parser$add_argument("-F1", "--F1", help = "F1 genotype. ")
    parser$add_argument("-condition", "--condition", help = "condition for differential expression comparison. ")
    parser$add_argument("-O", "--output", help = "the prefix for the output name.")
    parser$parse_args()
}

factorsize <- function(cnt) {
    # get factor size for all samples based on library size
    htseq_mat <- as.matrix(cnt)
    col <- data.frame(samples = colnames(cnt))
    dds <- DESeqDataSetFromMatrix(htseq_mat, colData = col, design = ~1) # ~1 means no design
    dds <- estimateSizeFactors(dds)
    factor_size <- sizeFactors(dds)
    factor_size <- as.data.frame(factor_size)
    #print(factor_size)
}

getMid <- function(count, coldata, index.P1, index.P2) {
    # get the average expression for P1 and P2
    meanLibSize = mean(coldata$factor_size[c(index.P1, index.P2)])
    g <- expand.grid(index.P1, index.P2)
    i <- 1
    mid <- (count[,g$Var1[i]]/coldata$factor_size[g$Var1[i]] + count[,g$Var2[i]]/coldata$factor_size[g$Var2[i]])*meanLibSize/2
    for (i in 2:nrow(g)) {
        mid <- cbind(mid,(count[,g$Var1[i]]/coldata$factor_size[g$Var1[i]] + count[,g$Var2[i]]/coldata$factor_size[g$Var2[i]])*meanLibSize/2)
    }
    mid <- data.frame(mid)
    names(mid) <- paste0("mid.", paste0(coldata$genotype[g$Var1], coldata$condition[g$Var1], coldata$replicate[g$Var1], "-", coldata$genotype[g$Var2], coldata$condition[g$Var2], coldata$replicate[g$Var2]))
    mid <- round(mid, 0)
    return(mid)
    #print(head(mid))
}

F1vsmid <- function(count, index.F1, mid) {
    cnt <- cbind(count[, index.F1], mid)
    info <- data.frame(sample = names(cnt), genome = c(rep("F1", 4), rep("mid", 16)))
    dds <- DESeqDataSetFromMatrix(countData = cnt, colData = info, design = ~genome)
    res <- results(DESeq(dds), alpha = 0.05, contrast = c("genome", "F1", "mid"))
    res <- as.data.frame(res)
    #print(head(res))
    #write.table(res, file = paste0(out, ".txt"), sep = "\t", quote = F)
}

F1vspar <- function(count, index.F1, index.P, out) {
    cnt <- cbind(count[, index.F1], count[, index.P])
    info <- data.frame(sample = names(cnt), genome = c(rep("F1", 4), rep("Par", 4)))
    dds <- DESeqDataSetFromMatrix(countData = cnt, colData = info, design = ~genome)
    res <- results(DESeq(dds), alpha = 0.05, contrast = c("genome", "F1", "Par"))
    res <- as.data.frame(res)
    #print(head(res))
    #write.table(res, file = paste0(out, ".txt"), sep = "\t", quote = F)
}

P1vsP2 <- function(count, index.P1, index.P2, out) {
    cnt <- cbind(count[, index.P1], count[, index.P2])
    info <- data.frame(sample = names(cnt), genome = c(rep("P1", 4), rep("P2", 4)))
    dds <- DESeqDataSetFromMatrix(countData = cnt, colData = info, design = ~genome)
    res <- results(DESeq(dds), alpha = 0.05, contrast = c("genome", "P1", "P2"))
    res <- as.data.frame(res)
    #print(head(res))
    #write.table(res, file = paste0(out, ".txt"), sep = "\t", quote = F)
}

classDominance<-function(F1vsMid, F1vsP1, F1vsP2, P1vsP2, Pname) {
    # Assign inheritance mode function based on differential expression
    log2fc.threshold <- 0
    # Hybrid vs Mid parental val
    F1vsMid <- data.frame(F1vsMid[,c("log2FoldChange", "padj")])
    names(F1vsMid) <- c("F1vsMid.lfc", "F1vsMid.padj")
    # Hybrid vs Parent 1
    F1vsP1 <- data.frame(F1vsP1[,c("log2FoldChange", "padj")])
    names(F1vsP1) <- c("F1vsP1.lfc", "F1vsP1.padj")
    # Hybrid vs Parent 2
    F1vsP2 <- data.frame(F1vsP2[,c("log2FoldChange", "padj")])
    names(F1vsP2) <- c("F1vsP2.lfc", "F1vsP2.padj")
    # Parent 1 vs Parent 2
    P1vsP2 <- data.frame(P1vsP2[,c("log2FoldChange", "padj")])
    names(P1vsP2) <- c("P1vsP2.lfc", "P1vsP2.padj")
    ### combine DE results
    tbl <- cbind(F1vsMid, F1vsP1, F1vsP2, P1vsP2)
    ### 
    tbl$mid <- ifelse(abs(tbl$F1vsMid.lfc) > log2fc.threshold & tbl$F1vsMid.padj < 0.05, "F1!=Mid", "F1=Mid")
    tbl$P1 <- ifelse(abs(tbl$F1vsP1.lfc) > log2fc.threshold & tbl$F1vsP1.padj < 0.05, "F1!=P1", "F1=P1")
    tbl$P2 <- ifelse(abs(tbl$F1vsP2.lfc) > log2fc.threshold & tbl$F1vsP2.padj < 0.05, "F1!=P2", "F1=P2")
    tbl$P1P2 <- ifelse(abs(tbl$P1vsP2.lfc) > log2fc.threshold & tbl$P1vsP2.padj < 0.05, ifelse(tbl$P1vsP2.lfc > log2fc.threshold & tbl$P1vsP2.padj < 0.05, "P1>P2","P1<P2"), "P1=P2")
    # together
    tbl$class <- paste(tbl$P1P2, tbl$mid, tbl$P1, tbl$P2, sep=";")
    # assign category
    tbl$category = "7.Other"
    tbl$category[grep("F1=Mid;F1!=P1;F1!=P2",tbl$class)] = "1.Additivity"
    tbl$category[grep("P1>P2;F1!=Mid;F1=P1;F1!=P2",tbl$class)] = paste0("2.",Pname[1], "-dominant,higher")
    tbl$category[grep("P1<P2;F1!=Mid;F1=P1;F1!=P2",tbl$class)] = paste0("2.",Pname[1], "-dominant,lower")
    tbl$category[grep("P1>P2;F1!=Mid;F1!=P1;F1=P2",tbl$class)] = paste0("3.",Pname[2], "-dominant,lower")
    tbl$category[grep("P1<P2;F1!=Mid;F1!=P1;F1=P2",tbl$class)] = paste0("3.",Pname[2], "-dominant,higher")
    tbl$category[grepl("P1=P2;F1!=Mid;F1!=P1;F1!=P2",tbl$class) & tbl$F1vsP1.lfc > 0 & tbl$F1vsP2.lfc > 0] = paste0("4.Transgressive Up:", Pname[1],"=", Pname[2])
    tbl$category[grepl("P1=P2;F1!=Mid;F1!=P1;F1!=P2",tbl$class) & tbl$F1vsP1.lfc < 0 & tbl$F1vsP2.lfc < 0] = paste0("4.Transgressive Down:",Pname[1],"=", Pname[2])
    #tbl$category[grep("P1=P2;F1!=Mid;F1!=P1;F1!=P2", tbl$class) & tbl$F1vsP1.lfc > 0 & tbl$F1vsP2.lfc > 0] = paste0("4.Transgressive Up:", Pname[1], "=", Pname[2])
    #tbl$category[grep("P1=P2;F1!=Mid;F1!=P1;F1!=P2", tbl$class) & tbl$F1vsP1.lfc < 0 & tbl$F1vsP2.lfc < 0] = paste0("4.Transgressive Down:", Pname[1], "=", Pname[2])
    tbl$category[grepl("P1>P2;F1!=Mid;F1!=P1;F1!=P2",tbl$class) & tbl$F1vsP1.lfc > 0 & tbl$F1vsP2.lfc > 0] = paste0("4.Transgressive Up:",Pname[1]," higher")
    tbl$category[grepl("P1<P2;F1!=Mid;F1!=P1;F1!=P2",tbl$class) & tbl$F1vsP1.lfc > 0 & tbl$F1vsP2.lfc > 0] = paste0("4.Transgressive Up:",Pname[2]," higher")
    tbl$category[grepl("P1>P2;F1!=Mid;F1!=P1;F1!=P2",tbl$class) & tbl$F1vsP1.lfc <0 & tbl$F1vsP2.lfc <0] = paste0("5.Transgressive Down:",Pname[1]," higher")
    tbl$category[grepl("P1<P2;F1!=Mid;F1!=P1;F1!=P2",tbl$class) & tbl$F1vsP1.lfc <0 & tbl$F1vsP2.lfc <0] = paste0("5.Transgressive Down: ",Pname[2]," higher")
    #tbl <- tbl[, c("class", "category")]
    return(tbl)
}

### arguments
args <- parse_arg()
cnt <- args$cnt
sample <- args$sample
output <- args$output
par1 <- args$parents[1]
par2 <- args$parents[2]
F1 <- args$F1
condition <- args$condition

countdf <- read.table(cnt, sep = "\t", header = T, row.names = 1)
coldf <- read.table(sample, sep = "\t", header = T, row.names = 1)

### factor size for libraries
fz <- factorsize(countdf)
coldf <- cbind(coldf, fz)
coldf$sample <- rownames(coldf)
rownames(coldf) <- 1:nrow(coldf)

### the average expression of P1 and P2
par1 <- which(coldf$genotype == par1 & coldf$condition == condition)
par2 <- which(coldf$genotype == par2 & coldf$condition == condition)

mid <- getMid(countdf, coldf, par1, par2)

### F1 expression vs. the average expression of P1&P2
F1 <- which(coldf$genotype == F1 & coldf$condition == condition)
print("Compare between F1 and mid (P1&P2)")
resF1vsmid <- F1vsmid(countdf, F1, mid)

### F1 expression vs. its parent
print("Compare between F1 and its parent (P1/P2)")
resF1vsP1 <- F1vspar(countdf, F1, par1)
resF1vsP2 <- F1vspar(countdf, F1, par2)

### parental lines comparison
print("Compare between P1 and P2 ")
resP1vsP2 <- P1vsP2(countdf, par1, par2)

### heritability
heritability <- classDominance(resF1vsmid, resF1vsP1, resF1vsP2, resP1vsP2, c(args$parents[1], args$parents[2]))
write.table(heritability, file = paste0(output, ".txt"), sep = "\t", quote = F)
